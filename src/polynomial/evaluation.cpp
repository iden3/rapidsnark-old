#include "evaluation.hpp"

template<typename Engine>
void Evaluation<Engine>::initialize(u_int64_t length) {
    E = Engine::engine;
    this->length = length;
    fft = new FFT<typename Engine::Fr>(length);
}

template<typename Engine>
Evaluation<Engine>::Evaluation(u_int64_t length) {
    initialize(length);
    eval = new FrElement[length];
    memset(eval, 0, sizeof(eval));
}

template<typename Engine>
Evaluation<Engine>::Evaluation(FrElement eval[], u_int64_t length) {
    initialize(length);
    this->eval = eval;
}

template<typename Engine>
Evaluation<Engine>::~Evaluation() {
    delete[] this->eval;
}