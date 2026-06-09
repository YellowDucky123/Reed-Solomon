#include <bits/stdc++.h>
#include <cmath>
#include <libff/algebra/field_utils/bigint.hpp>
#include <libff/algebra/fields/prime_base/fp.hpp>
//#include <gmpxx.h> // compile with -lgmpxx -lgmp

#pragma once

using namespace std;

typedef complex<double> cd;

class Poly {
	const mp_size_t limbs = 4;
	const libff:bigint<limbs> my_prime(524287);
	typedef libff::Fp_model<limbs, my_prime> Fp; 

public:
	Poly() {};

	/*
	 * Encode a message into a polynomial with lagrange interpolation
	 */
	vector<Fp> encode_message(string m, vector<Fp>& agreed_x) {
		int k = m.length();
		Fp::num_bits = my_prime.num_bits();
	
		// Using each character in the message as m_i
		// m_i is just a symbol, could be whatever is defined
		vector<Fp> values(k);
		for (int i = 0; i < k; ++i) {
			values[i] = Fp y(static_cast<int>(m[i]));
		}

		Vector<Fp> polynomial = interpolate(values, k, agreed_x);
		return polynomial;
	}


	// Only interpolate on the positions of the message not on parity
	vector<Fp> interpolate(vector<Fp> y, int& k, vector<Fp>& agreed_x) {
		vector<tuple<int, int>> points;
		for (size_t i = 0; i < k; ++i) {
			points.push_back(make_tuple(agreed_x[i], y[i]));
		}

		return lagrange_interpolation(points, agreed_x);
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
	vector<double> polyAtPoint(int idx, vector<Fp>& agreed_x) {
		vector<vector<double>> polynomial;
		for (int i = 0; i < agreed_x.size(); ++i) {
			if (i == idx) continue;

			vector<double> zero_p = {1.0, -agreed_x[i]}
			polynomial.push_back(zero_p);
		}
		
		if (polynomial.empty()) return {1.0};

		vector<double> p = polynomial[0];
		for (size_t i = 1; i < polynomial.size(); ++i) {
			p = multiply_double_polys(p, polynomials[i]);
		}

		double value_at_point = evaluate_polynomial(p, agreed_x[idx]);
		for (double& coeff : p) {
			coeff /= value_at_point;
		}

		return p;
	}

	/* 
	 * Interpolate the polynomial with lagrange interpolation
	 */
	vector<double> lagrange_interpolation(vector<tuple<int, int>> points, vector<Fp>& agreed_x) {
		int n = agreed_x.size();
		vector<vector<double>> polynomials(n);

		for (int i = 0; i < n; ++i) {
			polynomials[i] = polyAtPoint(i, agreed_x);
		}

		// Each polynomial for each point gets scalared by their y values in those points
		for (int i = 0; i < n; ++i) {
			vector<double>& p_i = polynomials[i]; 
			double y_val = get<1>(points[i]);
			for (double& p_ij : p_i) {
				p_ij *= y_val;
			}
		}

		int interpolated_degree = polynomials[0].size();
		vector<double> interpolated(interpolated_degree, 0.0);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < interpolated_degree; ++j) {
				interpolated[j] += polynomials[i][j];
			}
		}

		return interpolated;
	}

	/*
	 * Helper to multiply two double polynomials directly for Lagrange
	 */ 
    vector<double> multiply_double_polys(const vector<double>& p1, const vector<double>& p2) {
        vector<double> result(p1.size() + p2.size() - 1, 0.0);
        for (size_t i = 0; i < p1.size(); ++i) {
            for (size_t j = 0; j < p2.size(); ++j) {
                result[i + j] += p1[i] * p2[j];
            }
        }
        return result;
    }

	/* 
	 * Fast Fourier Transform
	 */
	vector<cd> fft(vector<cd> a) {
		int n = a.size();

		if (n == 1) {
			return a;
		}

		vector<cd> w(n);

		for (int i = 0; i < n; ++i) {
			double alpha = -2 * M_PI * i / n;
			w[i] = cd(cos(alpha), sin(alpha));
		}

		vector<cd> A0(n / 2), A1(n / 2);

		for (int i = 0; i < n; ++i) {
			A0[i] = a[i * 2];
			A1[i] = a[i * 2 + 1];
		}

		vector<cd> y0 = fft(A0);

		vector<cd> y1 = fft(A1);

		vector<cd> y(n);

		for (int k = 0; k < n / 2; ++k) {
			y[k] = y0[k] + w[k] * y1[k];
			y[k + n / 2] = y0[k] - w[k] * y1[k];
		}

		return y;
	}

	vector<double> Polynomial_multiplication(vector<int>& poly1, vector<int>& poly2) {
	    int n = 1;
		while (n < poly1.size() + poly2.size()) {
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
	}
