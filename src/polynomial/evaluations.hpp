#ifndef EVALUATIONS_HPP
#define EVALUATIONS_HPP

#include "assert.h"
#include <sstream>
#include <gmp.h>
#include "fft.hpp"
#include "polynomial.hpp"

template<typename Engine>
class Evaluations {
    using FrElement = typename Engine::FrElement;

    u_int64_t length;

    Engine &E;
    FFT<typename Engine::Fr> *fft;

    void initialize(u_int64_t length);

public:
    FrElement *eval;

    Evaluations(u_int64_t length);

    Evaluations(FrElement evaluations[]);

    Evaluations(Polynomial<FrElement> *polynomial);

    ~Evaluations();

    FrElement getEvaluation(u_int64_t index) const;

    u_int64_t getLength() const;
};

#endif