#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "assert.h"
#include <sstream>
#include <gmp.h>
#include "fft.hpp"
#include "polynomial.hpp"

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

    Evaluation(FrElement evaluations[]);

    Evaluation(Polynomial<FrElement> *polynomial);

    ~Evaluation();

    FrElement getEvaluation(u_int64_t index) const;

    //unisgned long length();
};


#endif //EVALUATION_HPP