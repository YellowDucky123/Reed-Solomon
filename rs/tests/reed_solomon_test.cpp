#include <gtest/gtest.h>
#include "reed_solomon.hpp"
#include <random>
#include <algorithm>

TEST(RS, EncodeLengthAndSystematic) {
	ReedSolomon rs(5, 4);
	auto code = rs.encode("hello");
	ASSERT_EQ(code.size(), 9u);
	std::string h = "hello";
	for (int i = 0; i < 5; ++i)
		EXPECT_EQ(code[i].v, (uint64_t)(unsigned char)h[i]) << "i=" << i;
}

TEST(RS, RoundTripNoErrors) {
	ReedSolomon rs(5, 4);
	auto code = rs.encode("hello");
	auto d = rs.decode(code);
	ASSERT_TRUE(d.ok);
	EXPECT_EQ(d.message, "hello");
}

TEST(RS, ExhaustiveSingleError) {
	ReedSolomon rs(4, 2); // t = 1
	std::string msg = "Wxyz";
	auto code = rs.encode(msg);
	for (size_t p = 0; p < code.size(); ++p) {
		for (int delta = 1; delta < 40; ++delta) {
			auto recv = code;
			recv[p] = recv[p] + Fp(delta);
			auto d = rs.decode(recv);
			ASSERT_TRUE(d.ok) << "pos " << p << " delta " << delta;
			EXPECT_EQ(d.message, msg) << "pos " << p << " delta " << delta;
		}
	}
}

TEST(RS, CorrectsUpToTErrors) {
	ReedSolomon rs(8, 6); // t = 3
	std::string msg = "ABCDEFGH";
	auto code = rs.encode(msg);
	std::mt19937_64 rng(2024);
	int t = rs.correctable();
	ASSERT_EQ(t, 3);
	for (int trial = 0; trial < 200; ++trial) {
		auto recv = code;
		std::vector<int> pos;
		while ((int)pos.size() < t) {
			int p = (int)(rng() % recv.size());
			if (std::find(pos.begin(), pos.end(), p) == pos.end()) pos.push_back(p);
		}
		for (int p : pos) // add a nonzero field element => guaranteed corruption
			recv[p] = recv[p] + Fp((long long)(1 + rng() % (Fp::MOD - 1)));
		auto d = rs.decode(recv);
		ASSERT_TRUE(d.ok) << "trial " << trial;
		EXPECT_EQ(d.message, msg) << "trial " << trial;
	}
}

TEST(RS, DetectsTooManyErrors) {
	ReedSolomon rs(4, 2); // t = 1, so 2 errors exceed capacity
	std::string msg = "test";
	auto code = rs.encode(msg);
	std::mt19937_64 rng(7);
	int failures = 0, silent = 0;
	for (int trial = 0; trial < 200; ++trial) {
		auto recv = code;
		std::vector<int> pos;
		while ((int)pos.size() < 2) { // 2 errors, one over capacity
			int p = (int)(rng() % recv.size());
			if (std::find(pos.begin(), pos.end(), p) == pos.end()) pos.push_back(p);
		}
		for (int p : pos)
			recv[p] = recv[p] + Fp((long long)(1 + rng() % (Fp::MOD - 1)));
		auto d = rs.decode(recv);
		// The decoder must never silently return a WRONG message: it either
		// reports failure, or (rarely) the corrupted word coincides with the
		// original codeword — never a confidently-wrong answer.
		if (!d.ok) failures++;
		else { EXPECT_EQ(d.message, msg) << "silent miscorrection, trial " << trial; silent++; }
	}
	EXPECT_GT(failures, 0); // at least some are detected as uncorrectable
}
