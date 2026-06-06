#pragma once
#include <libff/algebra/field_utils/bigint.hpp>
#include <libff/algebra/fields/prime_base/fp.hpp>
#include "Poly.hpp"

class Reed_solomon {
public:
static vector<double> encode(string m, int parity) {
	Poly poly = Poly();
	vector<Fp> encoded = poly.encode_message(m, parity);


	Fp value = poly.evaluate_polynomial(encoded, );
}
}
