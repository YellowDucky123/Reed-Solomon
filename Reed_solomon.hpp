#pragma once
#include <libff/algebra/field_utils/bigint.hpp>
#include <libff/algebra/fields/prime_base/fp.hpp>
#include "Poly.hpp"

// Fp and rs_modulus are defined in Poly.hpp (at namespace scope, as the
// Fp_model modulus template parameter requires).

class Reed_solomon {

public:
static vector<Fp> encode(string m, int parity) {
	int k = m.length();
	int n = k + parity;

	Poly poly = Poly();
	vector<Fp> output(n);

	// Agreed evaluation points: 1, 2, 3, ... n in the field.
	vector<Fp> agreed_x(n);
	Fp element = Fp::one();
	Fp step = Fp::one();
	for (int i = 0; i < n; ++i) {
		agreed_x[i] = element;
		element = element + step;
	}

	vector<Fp> encoded = poly.encode_message(m, agreed_x);

	// Get the values associated at each x to be sent out as the codeword
	for (int i = 0; i < n; ++i) {
		output[i] = poly.evaluate_polynomial(encoded, agreed_x[i]);
	}

	// output is values of the polynomial at each x_agreed positions
	return output;
}

// polynomial is codeword because codeword is coefficients of the polynomial
static string decode(vector<Fp> polynomial) {
	// Berlekamp-massey, syndrome decode
	return "";
}
};
