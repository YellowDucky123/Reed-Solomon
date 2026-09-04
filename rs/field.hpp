#pragma once
#include <cstdint>

/*
 * Fp — arithmetic in the prime field GF(p), p = 524287 = 2^19 - 1 (a Mersenne
 * prime). Kept self-contained on purpose: every value is a plain integer in
 * [0, p), and since p < 2^20 any product fits comfortably in a uint64_t, so we
 * never need Montgomery form, GMP, or a runtime init step. This makes the field
 * exact and trivial to test — the "Fp as source of truth" model.
 */
struct Fp {
	static const uint64_t MOD = 524287ULL; // 2^19 - 1, prime

	uint64_t v;

	Fp() : v(0) {}
	Fp(long long x) {
		long long m = x % (long long)MOD;
		if (m < 0) m += (long long)MOD;
		v = (uint64_t)m;
	}

	static Fp zero() { return Fp(0); }
	static Fp one()  { return Fp(1); }
	static Fp raw(uint64_t val) { Fp f; f.v = val; return f; }

	Fp operator+(const Fp& o) const { return raw((v + o.v) % MOD); }
	Fp operator-(const Fp& o) const { return raw((v + MOD - o.v) % MOD); }
	Fp operator*(const Fp& o) const { return raw((v * o.v) % MOD); }
	Fp operator-() const { return raw((MOD - v) % MOD); }

	bool operator==(const Fp& o) const { return v == o.v; }
	bool operator!=(const Fp& o) const { return v != o.v; }

	// Fast exponentiation.
	Fp pow(uint64_t e) const {
		Fp base = *this, result = one();
		while (e) {
			if (e & 1) result = result * base;
			base = base * base;
			e >>= 1;
		}
		return result;
	}

	// Multiplicative inverse via Fermat's little theorem (MOD is prime, so
	// a^(p-2) == a^{-1}). Inverting zero is undefined and returns zero.
	Fp inverse() const { return pow(MOD - 2); }
};
