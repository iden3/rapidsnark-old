#include "evaluations.hpp"
#include "thread_utils.hpp"


template<typename Engine>
void Evaluations<Engine>::initialize(u_int64_t length) {
    eval = new FrElement[length];
    memset(eval, 0, length * sizeof(FrElement));
    this->length = length;
}

template<typename Engine>
Evaluations<Engine>::Evaluations(Engine &_E, u_int64_t length) : E(_E) {
    this->initialize(length);
}

template<typename Engine>
Evaluations<Engine>::Evaluations(Engine &_E, FrElement *evaluations, u_int64_t length) : E(_E) {
    initialize(length);

    int nThreads = omp_get_max_threads() / 2;
    ThreadUtils<Engine>::parcpy(eval,
                                evaluations,
                                length * sizeof(FrElement), nThreads);

}

template<typename Engine>
Evaluations<Engine>::Evaluations(Engine &_E, FFT<typename Engine::Fr> *fft, Polynomial<Engine> &polynomial, u_int32_t extensionLength) : E(_E) {
    //Extend polynomial
    initialize(extensionLength);

    //TODO improve performance
    #pragma omp parallel for
    for (int64_t index = 0; index < polynomial.getDegree() + 1; ++index) {
        eval[index] = polynomial.getCoef(index);
    }

    //Coefficients to evaluations
    fft->fft(eval, extensionLength);
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