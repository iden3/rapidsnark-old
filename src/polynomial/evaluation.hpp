#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "assert.h"
#include <sstream>
#include <gmp.h>
#include "fft.hpp"

template<typename Engine>
class Evaluation {
    using FrElement = typename Engine::FrElement;

    Engine &E;
    FFT<typename Engine::Fr> *fft;

    void initialize(u_int64_t length);
public:
    FrElement *eval;

    u_int64_t length;

    Evaluation(u_int64_t length);

    Evaluation(FrElement evaluations[], u_int64_t length);

    ~Evaluation();

    //static Evaluation<FrElement> fromPolynomial(Polynomial<FrElement> polynomial);

    //FrElement getEvaluation(unsigned long index);

    //unisgned long length();
};


#endif //EVALUATION_HPP