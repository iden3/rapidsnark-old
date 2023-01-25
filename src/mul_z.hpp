#ifndef MUL_Z_HPP
#define MUL_Z_HPP

#include "assert.h"
#include <sstream>
#include <gmp.h>
#include "fft2.hpp"

template<typename Engine>
class MulZ {
    using FrElement = typename Engine::FrElement;

    Engine &E;

    FrElement Z1[4];
    FrElement Z2[4];
    FrElement Z3[4];
public:
    MulZ(Engine &E, FFT<typename Engine::Fr> *fft);

    tuple <FrElement, FrElement> mul2(const FrElement &a, const FrElement &b,
                                      const FrElement &ap, const FrElement &bp,
                                      int64_t p);

    tuple <FrElement, FrElement> mul4(const FrElement &a, const FrElement &b,
                                      const FrElement &c, const FrElement &d,
                                      const FrElement &ap, const FrElement &bp,
                                      const FrElement &cp, const FrElement &dp,
                                      int64_t p);
};

#include "mul_z.cpp"

#endif
