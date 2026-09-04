#include <gtest/gtest.h>
#include "poly.hpp"
#include <random>

using poly::vector;

static vector<Fp> randpoly(std::mt19937_64& rng, int deg) {
	vector<Fp> p(deg + 1);
	for (auto& c : p) c = Fp((long long)(rng() % Fp::MOD));
	poly::trim(p);
	if (p.empty()) p = { Fp::one() };
	return p;
}

TEST(Poly, DegreeAndTrim) {
	vector<Fp> p = { Fp(1), Fp(2), Fp(0), Fp(0) };
	poly::trim(p);
	EXPECT_EQ(p.size(), 2u);
	EXPECT_EQ(poly::degree(p), 1);
	vector<Fp> z = {};
	EXPECT_EQ(poly::degree(z), -1);
}

TEST(Poly, EvalHorner) {
	// p(x) = 3 + 2x + x^2, little-endian {3, 2, 1}
	vector<Fp> p = { Fp(3), Fp(2), Fp(1) };
	EXPECT_EQ(poly::eval(p, Fp(0)).v, 3u);
	EXPECT_EQ(poly::eval(p, Fp(1)).v, 6u);
	EXPECT_EQ(poly::eval(p, Fp(2)).v, 11u);
}

TEST(Poly, MultiplyKnown) {
	// (1 + x)^2 = 1 + 2x + x^2
	vector<Fp> a = { Fp(1), Fp(1) };
	vector<Fp> r = poly::mul(a, a);
	ASSERT_EQ(r.size(), 3u);
	EXPECT_EQ(r[0].v, 1u);
	EXPECT_EQ(r[1].v, 2u);
	EXPECT_EQ(r[2].v, 1u);
}

TEST(Poly, DivModIdentity) {
	std::mt19937_64 rng(12345);
	for (int t = 0; t < 500; ++t) {
		auto a = randpoly(rng, rng() % 9);
		auto b = randpoly(rng, rng() % 6);
		vector<Fp> q, r;
		poly::divmod(a, b, q, r);
		// a == q*b + r and deg(r) < deg(b)
		auto recon = poly::add(poly::mul(q, b), r);
		auto aa = a;
		poly::trim(aa);
		EXPECT_EQ(recon, aa);
		EXPECT_LT(poly::degree(r), poly::degree(b));
	}
}

TEST(Poly, InterpolationRecoversPoints) {
	std::mt19937_64 rng(999);
	int n = 6;
	vector<Fp> xs(n), ys(n);
	for (int i = 0; i < n; ++i) {
		xs[i] = Fp(i + 1);
		ys[i] = Fp((long long)(rng() % Fp::MOD));
	}
	auto p = poly::interpolate(xs, ys);
	EXPECT_LE(poly::degree(p), n - 1);
	for (int i = 0; i < n; ++i)
		EXPECT_EQ(poly::eval(p, xs[i]).v, ys[i].v) << "i=" << i;
}
