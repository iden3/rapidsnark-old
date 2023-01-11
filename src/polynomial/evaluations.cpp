#include "evaluations.hpp"


template<typename Engine>
void Evaluations<Engine>::initialize(u_int64_t length) {
    fft = new FFT<typename Engine::Fr>(length);

    eval = new FrElement[length];
    memset(eval, 0, sizeof(eval));
    this->length = length;
}

template<typename Engine>
Evaluations<Engine>::Evaluations(Engine &_E, u_int64_t length) : E(_E) {
    this->initialize(length);
}

template<typename Engine>
Evaluations<Engine>::Evaluations(Engine &_E, FrElement *evaluations) : E(_E) {
    u_int64_t len = sizeof(*evaluations) / sizeof(Engine::FrElement);
    initialize(len);

    memcpy(eval, evaluations, sizeof(eval));
}

template<typename Engine>
Evaluations<Engine>::Evaluations(Engine &_E, Polynomial<Engine> &polynomial) : E(_E) {
    u_int64_t extendedLength = polynomial.getLength() * 4;

    //Extend polynomial
    initialize(extendedLength);

#pragma omp parallel for
    for (int64_t index = 0; index < polynomial.getDegree() + 1; ++index) {
        eval[index] = polynomial.getCoef(index);
    }

    //Coefficients to evaluations
    fft->fft(eval, extendedLength);
}

template<typename Engine>
Evaluations<Engine>::~Evaluations() {
    delete[] this->eval;
}

template<typename Engine>
typename Engine::FrElement Evaluations<Engine>::getEvaluation(u_int64_t index) const {
    if (index > length - 1) {
        throw std::runtime_error("Evaluations::getEvaluation: invalid index");
    }
    return eval[index];
}

template<typename Engine>
u_int64_t Evaluations<Engine>::getLength() const {
    return length;
}