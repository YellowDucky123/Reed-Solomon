# Reed–Solomon codec (`rs/`)

A complete, self-contained Reed–Solomon encoder/decoder over the prime field
`GF(524287)` (`2^19 - 1`), with a GoogleTest suite. Header-only.

## Files

| File | Purpose |
|------|---------|
| `field.hpp` | `Fp` — exact arithmetic in `GF(524287)`. Plain integer mod-p; products fit in `uint64_t`, so no Montgomery form / GMP / runtime init. |
| `poly.hpp` | Polynomials over `Fp` (little-endian): eval, add/sub, multiply, long division, Lagrange interpolation. |
| `reed_solomon.hpp` | `ReedSolomon` — encode + decode. |
| `demo.cpp` | End-to-end example: encode → corrupt → recover. |
| `tests/` | GoogleTest unit tests for the field, polynomials, and the codec. |

## The scheme (original / evaluation view)

- A length-`k` message is `k` field symbols (one per byte for the string API).
- Encoding interpolates the unique polynomial `P` of degree `< k` through
  `(x_i, m_i)` for `i < k`, where `x_i = i+1` are fixed distinct points, then
  evaluates `P` at all `n = k + parity` points. The first `k` evaluations equal
  the message, so the code is **systematic** in its first `k` positions.
- Corrects up to `t = floor(parity / 2)` symbol errors.

Decoding uses the **Gao decoder**: interpolate the received values into `g1`,
then run a partial extended-Euclidean algorithm against
`g0 = prod (x - x_i)`, stopping when the remainder degree drops below
`(n + k) / 2`. Dividing that remainder by the tracked Bézout coefficient
recovers `P` exactly when at most `t` errors occurred; an inexact division
signals "too many errors". This works for arbitrary distinct points (no
primitive element / consecutive powers required).

## Build & run

GoogleTest must be installed (`libgtest-dev`).

```sh
cd rs
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure   # or: ./build/rs_tests
./build/rs_demo
```

> Note: on this machine an unrelated conda `libstdc++` shadows the system one at
> runtime. If you hit a `GLIBCXX_3.4.32 not found` error, prefix commands with
> `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`.

## Status

- Field arithmetic, polynomial ops, encode, and decode: **done and tested** (16
  tests passing — including exhaustive single-error correction, random
  up-to-`t`-error correction, and too-many-errors detection).
- Not implemented: erasure decoding, byte-packing of >8-bit symbols, and a
  performance pass (interpolation/multiply are schoolbook `O(n^2)`).
