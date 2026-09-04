#include "Poly.hpp"

// Single definition of the field modulus declared `extern` in Poly.hpp.
libff::bigint<RS_LIMBS> rs_modulus = libff::bigint<RS_LIMBS>("524287");
