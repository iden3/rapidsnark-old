#include "evaluations.hpp"


template<typename Engine>
void Evaluations<Engine>::initialize(u_int64_t length) {
    E = Engine::engine;
    fft = new FFT<typename Engine::Fr>(length);

    eval = new FrElement[length];
    memset(eval, 0, sizeof(eval));
    this->length = length;
}

template<typename Engine>
Evaluations<Engine>::Evaluations(u_int64_t length) {
    this->initialize(length);
}

template<typename Engine>
Evaluations<Engine>::Evaluations(FrElement evaluations[]) {
    u_int64_t len = sizeof(*evaluations) / sizeof(Engine::FrElement);
    initialize(len);

    memcpy(eval, evaluations, sizeof(eval));
}

template<typename Engine>
Evaluations<Engine>::Evaluations(Polynomial<FrElement> *polynomial) {
    u_int64_t extendedLength = polynomial->length * 4;

    //Extend polynomial
    initialize(extendedLength);

#pragma omp parallel for
    for (int64_t index = 0; index < polynomial->degree + 1; ++index) {
        eval[index] = polynomial[index];
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