#include "polynomial.hpp"

template<typename Engine>
void Polynomial<Engine>::initialize(u_int64_t length) {
    E = Engine::engine;
    this->length = length;
    fft = new FFT<typename Engine::Fr>(length);
}

template<typename Engine>
Polynomial<Engine>::Polynomial(u_int64_t length) {
    initialize(length);
    coef = new FrElement[length];
    memset(coef, 0, sizeof(coef));
    degree = 0;
}

template<typename Engine>
Polynomial<Engine>::Polynomial(FrElement coef[], u_int64_t length) {
    initialize(length);
    this->coef = coef;
    fixDegree();
}

template<typename Engine>
Polynomial<Engine>::~Polynomial() {
    delete[] this->coef;
}

template<typename Engine>
void Polynomial<Engine>::fixDegree() {
    u_int64_t degree;
    for (degree = length - 1; degree != 0 && this->E.fr.isZero(coef[degree]); degree--);
    this->degree = degree;
}

//template <typename Engine>
//static Polynomial<Engine::FrElement> Polynomial<Engine>::fromEvaluations(FrElement (*buffer)[]) {
//    //fft->fft(polynomial.data(), polynomial.size());
//}
/*
template <typename Engine>
static Polynomial<FrElement> fromCoefficientsArray(FrElement (*buffer)[]) {
}

template <typename Engine>
bool isEqual(const Polynomial& other) const {
    if (degree != other.degree) {
        return false;
    }
    for (int i = 0; i <= degree; i++) {
        if (coefficients[i] != other.coefficients[i]) {
            return false;
        }
    }
    return true;
}

//TODO     blindCoefficients(blindingFactors)

template <typename Engine>
int getCoef(int index) const {
    if (index < 0 || index > degree) {
        std::cerr << "Error: invalid index" << std::endl;
        exit(1);
    }
    return coefficients[index];
}

template <typename Engine>
void setCoef(int index, int value) {
    if (index < 0 || index > degree) {
        std::cerr << "Error: invalid index" << std::endl;
        exit(1);
    }
    coefficients[index] = value;
}

//TODO     static async to4T(buffer, domainSize, blindingFactors, Fr) {

template <typename Engine>
int length() const {
    return degree + 1;
}

template <typename Engine>
int degree() const {
    return degree;
}

template <typename Engine>
int evaluate(int x) const {
    int result = 0;
    for (int i = 0; i <= degree; i++) {
        result = result * x + coefficients[i];
    }
    return result;
}

template <typename Engine>
Polynomial add(const Polynomial& other) const {
    int newDegree = std::max(degree, other.degree);
    int* newCoefficients = new int[newDegree + 1];
    for (int i = 0; i <= newDegree; i++) {
        newCoefficients[i] = getCoef(i) + other.getCoef(i);
    }
    Polynomial result(newCoefficients, newDegree);
    delete[] newCoefficients;
    return result;
}

template <typename Engine>
Polynomial sub(const Polynomial& other) const {
    int newDegree = std::max(degree, other.degree);
    int* newCoefficients = new int[newDegree + 1];
    for (int i = 0; i <= newDegree; i++) {
        newCoefficients[i] = getCoef(i) - other.getCoef(i);
    }
    Polynomial result(newCoefficients, newDegree);
    delete[] newCoefficients;
    return result;
}

template <typename Engine>
void mulScalar(int scalar) {
    for (int i = 0; i <= degree; i++) {
        coefficients[i] *= scalar;
    }
}

//TODO template <typename Engine>
void addScalar(int scalar) {
}

//TODO template <typename Engine>
void subScalar(int scalar) {
}

//TODO byXSubValue(value) {

//TODO     divBy(polynomial) {

//TODO divZh

//TODO byX

//TODO static lagrangePolynomialInterpolation

//TODO static zerofierPolynomial(xArr, Fr)

//TODO print()


template <typename Engine>
int degree() const {
    for (int i = degree; i >= 0; i--) {
        if (coefficients[i] != 0) {
            return i;
        }
    }
    return 0;
}


*/