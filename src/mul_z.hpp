#ifndef MUL_Z_HPP
#define MUL_Z_HPP

#include "assert.h"
#include <sstream>
#include <gmp.h>
#include "fft.hpp"

template<typename Engine>
class MulZ {
    using FrElement = typename Engine::FrElement;

    Engine &E;

    FrElement Z1[4];
    FrElement Z2[4];
    FrElement Z3[4];
public:
    MulZ(Engine &E, FFT<typename Engine::Fr> *fft);

    tuple <FrElement, FrElement> mul2( FrElement &a,  FrElement &b,
                                       FrElement &ap,  FrElement &bp,
                                      int64_t p);

    tuple <FrElement, FrElement> mul4( FrElement &a,  FrElement &b,
                                       FrElement &c,  FrElement &d,
                                       FrElement &ap,  FrElement &bp,
                                       FrElement &cp,  FrElement &dp,
                                      int64_t p);
};

#include "mul_z.cpp"

#endif
