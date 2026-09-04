#pragma once
#include "field.hpp"
#include <vector>
#include <algorithm>
#include <cassert>

/*
 * Polynomials over Fp, stored little-endian: coeff[i] is the coefficient of
 * x^i, so the degree is (size - 1) once trailing zeros are trimmed. The zero
 * polynomial is the empty vector, with degree -1.
 */
namespace poly {

using std::vector;

// Drop trailing zero coefficients so size-1 is the true degree.
inline void trim(vector<Fp>& a) {
	while (!a.empty() && a.back() == Fp::zero()) a.pop_back();
}

// Degree of the polynomial; -1 for the zero polynomial.
inline int degree(const vector<Fp>& a) {
	for (int i = (int)a.size() - 1; i >= 0; --i)
		if (a[i] != Fp::zero()) return i;
	return -1;
}

// Evaluate p(x) by Horner's method.
inline Fp eval(const vector<Fp>& a, Fp x) {
	Fp r = Fp::zero();
	for (int i = (int)a.size() - 1; i >= 0; --i) r = r * x + a[i];
	return r;
}

inline vector<Fp> add(const vector<Fp>& a, const vector<Fp>& b) {
	vector<Fp> r(std::max(a.size(), b.size()), Fp::zero());
	for (size_t i = 0; i < a.size(); ++i) r[i] = r[i] + a[i];
	for (size_t i = 0; i < b.size(); ++i) r[i] = r[i] + b[i];
	trim(r);
	return r;
}

inline vector<Fp> sub(const vector<Fp>& a, const vector<Fp>& b) {
	vector<Fp> r(std::max(a.size(), b.size()), Fp::zero());
	for (size_t i = 0; i < a.size(); ++i) r[i] = r[i] + a[i];
	for (size_t i = 0; i < b.size(); ++i) r[i] = r[i] - b[i];
	trim(r);
	return r;
}

// Schoolbook O(n*m) multiplication — fine for the sizes an RS codeword uses.
inline vector<Fp> mul(const vector<Fp>& a, const vector<Fp>& b) {
	if (a.empty() || b.empty()) return {};
	vector<Fp> r(a.size() + b.size() - 1, Fp::zero());
	for (size_t i = 0; i < a.size(); ++i)
		for (size_t j = 0; j < b.size(); ++j)
			r[i + j] = r[i + j] + a[i] * b[j];
	trim(r);
	return r;
}

inline vector<Fp> scalar_mul(const vector<Fp>& a, Fp s) {
	vector<Fp> r(a.size());
	for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] * s;
	trim(r);
	return r;
}

// Polynomial long division: computes q, r such that a = q*b + r with
// deg(r) < deg(b). b must be nonzero.
inline void divmod(const vector<Fp>& a, const vector<Fp>& b,
                   vector<Fp>& q, vector<Fp>& r) {
	int db = degree(b);
	assert(db >= 0 && "division by zero polynomial");
	r = a;
	trim(r);
	Fp binv = b[db].inverse();
	int qsize = degree(a) - db + 1;
	q.assign(qsize > 0 ? qsize : 0, Fp::zero());
	while (degree(r) >= db) {
		int dr = degree(r);
		Fp coeff = r[dr] * binv;
		int shift = dr - db;
		q[shift] = coeff;
		for (int i = 0; i <= db; ++i)
			r[shift + i] = r[shift + i] - coeff * b[i];
		trim(r);
	}
	trim(q);
}

// Lagrange interpolation: the unique polynomial of degree < n passing through
// the n points (xs[i], ys[i]). The xs must be distinct.
inline vector<Fp> interpolate(const vector<Fp>& xs, const vector<Fp>& ys) {
	assert(xs.size() == ys.size());
	int n = (int)xs.size();
	vector<Fp> result; // zero polynomial
	for (int i = 0; i < n; ++i) {
		// Basis polynomial L_i(x) = prod_{j != i} (x - xs[j]) / (xs[i] - xs[j]).
		vector<Fp> num = { Fp::one() };
		Fp denom = Fp::one();
		for (int j = 0; j < n; ++j) {
			if (j == i) continue;
			vector<Fp> factor = { -xs[j], Fp::one() }; // (x - xs[j])
			num = mul(num, factor);
			denom = denom * (xs[i] - xs[j]);
		}
		result = add(result, scalar_mul(num, ys[i] * denom.inverse()));
	}
	trim(result);
	return result;
}

} // namespace poly
