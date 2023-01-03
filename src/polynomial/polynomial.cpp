#include "polynomial.hpp"


template<typename Engine>
void Polynomial<Engine>::initialize(u_int64_t length) {
    E = Engine::engine;
    fft = new FFT<typename Engine::Fr>(length);

    coef = new FrElement[length];
    memset(coef, 0, sizeof(coef));
    this->length = length;
    degree = 0;
}

template<typename Engine>
Polynomial<Engine>::Polynomial(u_int64_t length) {
    this->initialize(length);
}

template<typename Engine>
Polynomial<Engine>::Polynomial(FrElement elements[]) {
    //TODO checks
    u_int64_t len = sizeof(*elements) / sizeof(Engine::FrElement);
    initialize(len);
    memcpy(coef, elements, sizeof(elements));
    fixDegree();
}

template<typename Engine>
Polynomial<Engine>::Polynomial(FrElement elements[], u_int64_t domainSize, u_int64_t nBlindCoefficients) {
    //TODO checks
    initialize(domainSize + nBlindCoefficients);
    memcpy(coef, elements, sizeof(elements));

    fft->ifft(coef, domainSize);

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


template<typename Engine>
bool Polynomial<Engine>::isEqual(const Polynomial<Engine> &other) const {
    if (degree != other.degree) {
        return false;
    }

    for (int i = 0; i <= degree; i++) {
        if (E.fr.noeq(coef[i], other.coef[i])) {
            return false;
        }
    }
    return true;
}

template<typename Engine>
void Polynomial<Engine>::blindCoefficients(FrElement blindingFactors[]){
    u_int64_t lenBlindingFactors = sizeof(*blindingFactors) / sizeof(Engine::FrElement);

    for (int i = 0; i < lenBlindingFactors; i++) {
        coef[length - lenBlindingFactors + i] = E.fr.add(coef[length - lenBlindingFactors + i], blindingFactors[i]);
        coef[i] = E.fr.sub(coef[i], blindingFactors[i]);
    }
}

template <typename Engine>
typename Engine::FrElement Polynomial<Engine>::getCoef(u_int64_t index) const {
    if (index > degree) {
        throw std::runtime_error("Polynomial::getCoef: invalid index");
    }
    return coef[index];
}

template <typename Engine>
void Polynomial<Engine>::setCoef(u_int64_t index, FrElement value) {
    if (index > degree) {
        throw std::runtime_error("Polynomial::getCoef: invalid index");
    }
    coef[index] = value;
}


//TODO     static async to4T(buffer, domainSize, blindingFactors, Fr) {
/*
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
*/
template <typename Engine>
void  Polynomial<Engine>::add(Polynomial<FrElement> &other, FrElement &blindingValue) {
/*    int newDegree = std::max(degree, other.degree);
    int* newCoefficients = new int[newDegree + 1];
    for (int i = 0; i <= newDegree; i++) {
        newCoefficients[i] = getCoef(i) + other.getCoef(i);
    }
    Polynomial result(newCoefficients, newDegree);
    delete[] newCoefficients;
    return result;*/
}
/*
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
*/

template <typename Engine>
void  Polynomial<Engine>::divZh() {
}

/*
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