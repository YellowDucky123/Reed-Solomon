#include <bits/stdc++.h>
#include <cmath>
#include <libff/algebra/field_utils/bigint.hpp>
#include <libff/algebra/fields/prime_base/fp.hpp>
//#include <gmpxx.h> // compile with -lgmpxx -lgmp

#pragma once

using namespace std;

typedef complex<double> cd;

// The modulus must live at namespace scope with external linkage: Fp_model
// takes it as a `const bigint<n>&` non-type template parameter, so it needs
// linkage (hence `extern`, defined once in field.cpp) and cannot be a class
// member.
const mp_size_t RS_LIMBS = 4;
extern libff::bigint<RS_LIMBS> rs_modulus;
typedef libff::Fp_model<RS_LIMBS, rs_modulus> Fp;

class Poly {
public:
	Poly() {};

	/*
	 * Encode a message into a polynomial with lagrange interpolation
	 */
	vector<Fp> encode_message(string m, vector<Fp>& agreed_x) {
		int k = m.length();

		// Using each character in the message as m_i
		// m_i is just a symbol, could be whatever is defined
		vector<Fp> values(k);
		for (int i = 0; i < k; ++i) {
			values[i] = Fp(static_cast<long>(m[i]));
		}

		vector<Fp> polynomial = interpolate(values, k, agreed_x);
		return polynomial;
	}


	// Only interpolate on the positions of the message not on parity
	vector<Fp> interpolate(vector<Fp> y, int k, vector<Fp>& agreed_x) {
		return lagrange_interpolation(y, agreed_x);
	}

	/*
	 * Evaluate a polynomial at a point
	 */
	Fp evaluate_polynomial(vector<Fp> p, Fp x) {
		if (p.empty()) {
			return Fp::zero();
		}

		Fp value = p[0];
		for (size_t i = 1; i < p.size(); ++i) {
			value = (value * x) + p[i];
		}

		return value;
	}

	/*
	 * Gets the polynomial where its 1 at point x and 0 on all others
	 * takes the index of the agreed_x, not the value of the x
	 */
	vector<Fp> polyAtPoint(int idx, vector<Fp>& agreed_x) {
		vector<vector<Fp>> polynomial;
		for (size_t i = 0; i < agreed_x.size(); ++i) {
			if ((int)i == idx) continue;

			vector<Fp> zero_p = {Fp::one(), -agreed_x[i]};
			polynomial.push_back(zero_p);
		}

		if (polynomial.empty()) return {Fp::one()};

		vector<Fp> p = polynomial[0];
		for (size_t i = 1; i < polynomial.size(); ++i) {
			p = multiply_polys(p, polynomial[i]);
		}

		Fp value_at_point = evaluate_polynomial(p, agreed_x[idx]);
		Fp inv = value_at_point.inverse();
		for (Fp& coeff : p) {
			coeff = coeff * inv;
		}

		return p;
	}

	/*
	 * Interpolate the polynomial with lagrange interpolation
	 */
	vector<Fp> lagrange_interpolation(vector<Fp> y, vector<Fp>& agreed_x) {
		int n = agreed_x.size();
		vector<vector<Fp>> polynomials(n);

		for (int i = 0; i < n; ++i) {
			polynomials[i] = polyAtPoint(i, agreed_x);
		}

		// Each polynomial for each point gets scalared by their y values in those points
		for (int i = 0; i < n; ++i) {
			vector<Fp>& p_i = polynomials[i];
			Fp y_val = y[i];
			for (Fp& p_ij : p_i) {
				p_ij = p_ij * y_val;
			}
		}

		size_t interpolated_degree = polynomials[0].size();
		vector<Fp> interpolated(interpolated_degree, Fp::zero());
		for (int i = 0; i < n; ++i) {
			for (size_t j = 0; j < interpolated_degree; ++j) {
				interpolated[j] = interpolated[j] + polynomials[i][j];
			}
		}

		return interpolated;
	}

	/*
	 * Helper to multiply two Fp polynomials directly for Lagrange
	 */
	vector<Fp> multiply_polys(const vector<Fp>& p1, const vector<Fp>& p2) {
		vector<Fp> result(p1.size() + p2.size() - 1, Fp::zero());
		for (size_t i = 0; i < p1.size(); ++i) {
			for (size_t j = 0; j < p2.size(); ++j) {
				result[i + j] = result[i + j] + (p1[i] * p2[j]);
			}
		}
		return result;
	}

	/*
	 * Fast Fourier Transform
	 */
	vector<cd> fft(vector<cd> a, bool invert = false) {
		int n = a.size();

		if (n == 1) {
			return a;
		}

		vector<cd> w(n);
		for (int i = 0; i < n; ++i) {
			// Sign is negative for forward FFT, positive for inverse IFFT
			double alpha = 2 * M_PI * i / n * (invert ? 1 : -1);
			w[i] = cd(cos(alpha), sin(alpha));
		}

		vector<cd> A0(n / 2), A1(n / 2);

		for (int i = 0; i < n / 2; ++i) {
			A0[i] = a[i * 2];
			A1[i] = a[i * 2 + 1];
		}

		vector<cd> y0 = fft(A0, invert);
		vector<cd> y1 = fft(A1, invert);

		vector<cd> y(n);
		for (int k = 0; k < n / 2; ++k) {
			y[k] = y0[k] + w[k] * y1[k];
			y[k + n / 2] = y0[k] - w[k] * y1[k];
		}

		return y;
	}

	vector<double> Polynomial_multiplication(vector<int>& poly1, vector<int>& poly2) {
		int n = 1;
		while (n < (int)(poly1.size() + poly2.size())) {
			n <<= 1;
		}

		vector<cd> a(n, 0), b(n, 0);
		for (size_t i = 0; i < poly1.size(); i++) a[i] = poly1[i];
		for (size_t i = 0; i < poly2.size(); i++) b[i] = poly2[i];

		vector<cd> fa = fft(a, false);
		vector<cd> fb = fft(b, false);

		vector<cd> fc(n);
		for (int i = 0; i < n; i++) {
			fc[i] = fa[i] * fb[i];
		}

		vector<cd> c = fft(fc, true);

		vector<double> result(n);
		for (int i = 0; i < n; i++) {
			result[i] = round(c[i].real() / n);
		}

		while (result.size() > 1 && result.back() == 0) {
			result.pop_back();
		}

		return result;
	}
};
