#pragma once
#include "field.hpp"
#include "poly.hpp"
#include <string>
#include <vector>
#include <cassert>

/*
 * Reed-Solomon codec, "original view" (Reed & Solomon 1960):
 *
 *   - A length-k message is treated as k field symbols.
 *   - We interpolate the unique polynomial P of degree < k through the points
 *     (x_0, m_0), ..., (x_{k-1}, m_{k-1}), where x_i = i+1 are fixed, distinct
 *     evaluation points.
 *   - The codeword is P evaluated at all n = k + parity points. Because the
 *     first k evaluations reproduce the message symbols, the code is
 *     systematic in its first k positions.
 *
 * Any k of the n evaluations determine P, so the code corrects up to
 * t = floor(parity / 2) symbol errors.
 *
 * Decoding uses the Gao decoder: it interpolates the received values, then runs
 * a partial extended-Euclidean algorithm against prod (x - x_i), stopping once
 * the remainder drops below degree (n+k)/2. Dividing that remainder by the
 * tracked Bezout coefficient recovers P exactly when at most t errors occurred.
 */
class ReedSolomon {
public:
	int k;      // message length in symbols
	int parity; // number of parity symbols
	int n;      // codeword length = k + parity
	std::vector<Fp> xs; // n distinct evaluation points

	ReedSolomon(int k_, int parity_) : k(k_), parity(parity_), n(k_ + parity_) {
		assert(k > 0 && parity >= 0);
		assert(n < (int)Fp::MOD && "not enough distinct points in the field");
		xs.resize(n);
		for (int i = 0; i < n; ++i) xs[i] = Fp(i + 1); // points 1, 2, ..., n
	}

	// Number of symbol errors this configuration can correct.
	int correctable() const { return parity / 2; }

	// Encode k message symbols into an n-symbol codeword.
	std::vector<Fp> encode_symbols(const std::vector<Fp>& msg) const {
		assert((int)msg.size() == k);
		std::vector<Fp> xk(xs.begin(), xs.begin() + k);
		std::vector<Fp> P = poly::interpolate(xk, msg);
		std::vector<Fp> code(n);
		for (int i = 0; i < n; ++i) code[i] = poly::eval(P, xs[i]);
		return code;
	}

	// Convenience: encode a string, one symbol per byte.
	std::vector<Fp> encode(const std::string& m) const {
		std::vector<Fp> msg(m.size());
		for (size_t i = 0; i < m.size(); ++i)
			msg[i] = Fp((long long)(unsigned char)m[i]);
		return encode_symbols(msg);
	}

	struct DecodeResult {
		bool ok;                    // false if too many errors to correct
		std::vector<Fp> message;    // recovered k message symbols (if ok)
	};

	// Decode an n-symbol received word, correcting up to correctable() errors.
	DecodeResult decode_symbols(const std::vector<Fp>& received) const {
		assert((int)received.size() == n);

		// g0 = prod_i (x - x_i), the vanishing polynomial of the point set.
		std::vector<Fp> g0 = { Fp::one() };
		for (int i = 0; i < n; ++i) {
			std::vector<Fp> factor = { -xs[i], Fp::one() };
			g0 = poly::mul(g0, factor);
		}

		// g1 = interpolation of the received values (= P plus error).
		std::vector<Fp> g1 = poly::interpolate(xs, received);

		// Partial extended-Euclidean algorithm on (g0, g1). We track only the
		// Bezout coefficient v that multiplies g1, since the recovered message
		// polynomial is remainder / v.
		std::vector<Fp> r_prev = g0, r_cur = g1;
		std::vector<Fp> v_prev = {}, v_cur = { Fp::one() };
		const int stop = n + k; // stop when 2*deg(r_cur) < n + k
		while (2 * poly::degree(r_cur) >= stop) {
			std::vector<Fp> q, rem;
			poly::divmod(r_prev, r_cur, q, rem);
			r_prev = r_cur;
			r_cur = rem;
			std::vector<Fp> v_next = poly::sub(v_prev, poly::mul(q, v_cur));
			v_prev = v_cur;
			v_cur = v_next;
			if (poly::degree(r_cur) < 0) break; // exact remainder, nothing left
		}

		// Recover P = r_cur / v_cur. If the division is inexact, or the result
		// has too high a degree, the error weight exceeded what we can correct.
		if (poly::degree(v_cur) < 0) return { false, {} };
		std::vector<Fp> P, rem;
		poly::divmod(r_cur, v_cur, P, rem);
		if (poly::degree(rem) >= 0) return { false, {} };
		if (poly::degree(P) > k - 1) return { false, {} };

		std::vector<Fp> msg(k);
		for (int i = 0; i < k; ++i) msg[i] = poly::eval(P, xs[i]);
		return { true, msg };
	}

	struct DecodeString {
		bool ok;
		std::string message;
	};

	// Convenience: decode back to a string.
	DecodeString decode(const std::vector<Fp>& received) const {
		DecodeResult r = decode_symbols(received);
		if (!r.ok) return { false, "" };
		std::string s(k, '\0');
		for (int i = 0; i < k; ++i) s[i] = (char)(unsigned char)(r.message[i].v);
		return { true, s };
	}
};
