#include <gtest/gtest.h>
#include "field.hpp"

TEST(Field, AddSubMul) {
	Fp a(10), b(20);
	EXPECT_EQ((a + b).v, 30u);
	EXPECT_EQ((b - a).v, 10u);
	EXPECT_EQ((a * b).v, 200u);
}

TEST(Field, Wraparound) {
	Fp a((long long)Fp::MOD - 1);
	EXPECT_EQ((a + Fp(1)).v, 0u);          // (p-1) + 1 == 0
	EXPECT_EQ((Fp(0) - Fp(1)).v, Fp::MOD - 1); // 0 - 1 == p-1
}

TEST(Field, NegativeConstruction) {
	EXPECT_EQ(Fp(-1).v, Fp::MOD - 1);
	EXPECT_EQ(Fp(-(long long)Fp::MOD).v, 0u);
	EXPECT_EQ(Fp((long long)Fp::MOD + 5).v, 5u);
}

TEST(Field, InverseIsReciprocal) {
	for (long long x = 1; x < 2000; ++x) {
		Fp a(x);
		EXPECT_EQ((a * a.inverse()).v, 1u) << "x=" << x;
	}
}

TEST(Field, PowMatchesRepeatedMul) {
	Fp a(7);
	Fp acc = Fp::one();
	for (int e = 0; e < 50; ++e) {
		EXPECT_EQ(a.pow(e).v, acc.v) << "e=" << e;
		acc = acc * a;
	}
}

TEST(Field, NoOverflowOnLargeProduct) {
	Fp a((long long)Fp::MOD - 1);
	// (p-1)^2 mod p == 1; verifies the uint64 product path doesn't overflow.
	EXPECT_EQ((a * a).v, 1u);
}
