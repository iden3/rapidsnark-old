#include "polynomial.hpp"

#include "logger.hpp"

using namespace CPlusPlusLogging;

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
void Polynomial<Engine>::blindCoefficients(FrElement blindingFactors[]) {
    u_int64_t lenBlindingFactors = sizeof(*blindingFactors) / sizeof(Engine::FrElement);

    for (int i = 0; i < lenBlindingFactors; i++) {
        coef[length - lenBlindingFactors + i] = E.fr.add(coef[length - lenBlindingFactors + i], blindingFactors[i]);
        coef[i] = E.fr.sub(coef[i], blindingFactors[i]);
    }
}

template<typename Engine>
typename Engine::FrElement Polynomial<Engine>::getCoef(u_int64_t index) const {
    if (index > degree) {
        throw std::runtime_error("Polynomial::getCoef: invalid index");
    }
    return coef[index];
}

template<typename Engine>
void Polynomial<Engine>::setCoef(u_int64_t index, FrElement value) {
    if (index > degree) {
        throw std::runtime_error("Polynomial::getCoef: invalid index");
    }
    coef[index] = value;
}

//TODO     static async to4T(buffer, domainSize, blindingFactors, Fr) {

template<typename Engine>
u_int64_t Polynomial<Engine>::getLength() const {
    return length;
}

template<typename Engine>
u_int64_t Polynomial<Engine>::getDegree() const {
    return degree;
}

template<typename Engine>
typename Engine::FrElement Polynomial<Engine>::evaluate(FrElement point) const {
    FrElement result = E.fr.zero;
    for (u_int64_t i = 0; i <= degree; i++) {
        result = E.fr.add(E.fr.mul(result, point), coef[i]);
    }
    return result;
}


template<typename Engine>
void Polynomial<Engine>::add(Polynomial<Engine> &polynomial) {
    FrElement *newCoef;
    bool resize = polynomial.length > this->length;

    if (resize) {
        newCoef = new FrElement[polynomial.length];
    }

    u_int64_t thisLength = this->length;
    u_int64_t polyLength = polynomial.length;
    // TODO parallelize
    for (u_int64_t i = 0; i < std::max(thisLength, polyLength); i++) {
        FrElement a = i < thisLength ? this->coef[i] : this.Fr.zero;
        FrElement b = i < polyLength ? polynomial->coef[i] : this.Fr.zero;
        FrElement sum = E.fr.add(a, b);

        if (resize) {
            newCoef[i] = sum;
        } else {
            this->coef[i] = sum;
        }
    }

    if (resize) {
        delete this->coef;
        this->coef = newCoef;
    }

    fixDegree();
}

template<typename Engine>
void Polynomial<Engine>::sub(Polynomial<Engine> &polynomial) {
    FrElement *newCoef;
    bool resize = polynomial.length > this->length;

    if (resize) {
        newCoef = new FrElement[polynomial.length];
    }

    u_int64_t thisLength = this->length;
    u_int64_t polyLength = polynomial.length;
    // TODO parallelize
    for (u_int64_t i = 0; i < std::max(thisLength, polyLength); i++) {
        FrElement a = i < thisLength ? this->coef[i] : this.Fr.zero;
        FrElement b = i < polyLength ? polynomial->coef[i] : this.Fr.zero;
        FrElement sub = E.fr.sub(a, b);

        if (resize) {
            newCoef[i] = sub;
        } else {
            this->coef[i] = sub;
        }
    }

    if (resize) {
        delete this->coef;
        this->coef = newCoef;
    }

    fixDegree();
}

template<typename Engine>
void Polynomial<Engine>::mulScalar(FrElement &value) {
    for (u_int64_t i = 0; i < this->length; i++) {
        this.coef[i] = E.fr.mul(this->coef[i], value);
    }
}

template<typename Engine>
void Polynomial<Engine>::addScalar(FrElement &value) {
    FrElement currentValue = (0 == this->length) ? E.fr.zero : this->coef[0];
    this->coef[0] = E.fr.add(currentValue, value);
}

template<typename Engine>
void Polynomial<Engine>::subScalar(FrElement &value) {
    FrElement currentValue = (0 == this->length) ? E.fr.zero : this->coef[0];
    this->coef[0] = E.fr.sub(currentValue, value);
}

// Multiply current polynomial by the polynomial (X - value)
template<typename Engine>
void Polynomial<Engine>::byXSubValue(FrElement &value) {
    bool resize = E.fr.neq(E.fr.zero, this->coef[this->length - 1]);

    u_int64_t length = resize ? this->length + 1 : this->length;
    Polynomial<Engine> *pol = new Polynomial<Engine>(length);

    // Step 0: Set current coefficients to the new buffer shifted one position
    memcpy(pol->coef[1], this->coef, sizeof(this->coef));

    // Step 1: multiply each coefficient by (-value)
    this->mulScalar(E.fr.neg(value));

    // Step 2: Add current polynomial to destination polynomial
    pol->add(this);

    // Swap buffers
    delete this->coef;
    this->coef = pol->coef;

    fixDegree();
}

// Euclidean division
template<typename Engine>
Polynomial<Engine> Polynomial<Engine>::divBy(Polynomial<Engine> &polynomial) {
    u_int64_t degreeA = this->degree;
    u_int64_t degreeB = polynomial->degree;

    FrElement *buffPolR = this->coef;

    this.coef = new FrElement[this->length];

    for (u_int64_t i = degreeA - degreeB; i >= 0; i--) {
        this->coef[i] = E.fr.div(buffPolR[i + degreeB], polynomial[degreeB]);

        for (u_int64_t j = 0; j <= degreeB; j++) {
            buffPolR[i + j] = E.fr.sub(buffPolR[i + j], E.fr.mul(this->coef[i], polynomial->coef[j]));
        }
    }

    Polynomial<Engine> *polR = new Polynomial<Engine>(0);
    polR->coef = buffPolR;

    fixDegree();
    polR->fixDegree();

    return polR;
}

template<typename Engine>
void Polynomial<Engine>::divZh(u_int64_t domainSize) {
    for (u_int64_t i = 0; i < domainSize; i++) {
        this->coef[i] = E.fr.neg(this->coef.slice[i]);
    }

    for (u_int64_t i = domainSize; i < this->length; i++) {
        coef[i] = E.fr.sub(coef[i - domainSize], coef[i]);

        if (i > (domainSize * 3 - 4)) {
            if (!E.fr.isZero(coef[i])) {
                throw std::runtime_error("Polynomial is not divisible");
            }
        }
    }
    fixDegree();
}

template<typename Engine>
void Polynomial<Engine>::byX() {
    bool resize = E.fr.neq(E.fr.zero, this->coef[this->length - 1]);

    if (resize) {
        FrElement *newCoef = new FrElement[this->length + 1];
        memcpy(newCoef[1], coef[0], sizeof(coef));
        coef = newCoef;
    } else {
        memcpy(coef[1], coef[0], sizeof(coef));
    }

    coef[0] = E.fr.zero;
    fixDegree();
}

template<typename Engine>
Polynomial<Engine>*
Polynomial<Engine>::lagrangePolynomialInterpolation(FrElement xArr[], FrElement yArr[]) {
    Polynomial<Engine> *polynomial = computeLagrangePolynomial(0);
    u_int64_t len = sizeof(xArr);

    for (u_int64_t i = 1; i < len; i++) {
        polynomial->add(computeLagrangePolynomial(i));
    }

    return polynomial;
}

template<typename Engine>
Polynomial<Engine> Polynomial<Engine>::computeLagrangePolynomial(u_int64_t i, Polynomial<Engine> &pol, FrElement xArr[],
                                                                 FrElement yArr[]) {
    Engine &E = Engine::engine;
    Polynomial<Engine> *polynomial = nullptr;
    u_int64_t len = sizeof(xArr);

    for (u_int64_t j = 0; j < len; j++) {
        if (j == i) continue;

        if (nullptr == polynomial) {
            polynomial = new Polynomial<Engine>(len);
            polynomial[0] = E.fr.neg(xArr[j]);
            polynomial[1] = E.fr.one;
        } else {
            polynomial->byXSubValue(xArr[j]);
        }
    }

    FrElement denominator = polynomial->evaluate(xArr[i]);
    denominator = E.fr.inv(denominator);
    FrElement mulFactor = E.fr.mul(yArr[i], denominator);

    polynomial->mulScalar(mulFactor);

    return polynomial;
}

template<typename Engine>
Polynomial<Engine>* Polynomial<Engine>::zerofierPolynomial(FrElement xArr[]) {
    Engine &E = Engine::engine;
    u_int64_t len = sizeof(xArr);
    Polynomial<Engine> polynomial = new Polynomial<Engine>(len);

    // Build a zerofier polynomial with the following form:
    // zerofier(X) = (X-xArr[0])(X-xArr[1])...(X-xArr[n])
    polynomial[0] = E.fr.neg(xArr[0]);
    polynomial[1] = E.fr.one;

    for (u_int64_t i = 1; i < len; i++) {
        polynomial->byXSubValue(xArr[i]);
    }

    return polynomial;
}

template<typename Engine>
void Polynomial<Engine>::print() {
    std::ostringstream res;

    for (u_int64_t i = this->degree; i >= 0; i--) {
        FrElement coef = coef[i];
        if (E.fr.neq(E.fr.zero, coef)) {
            if (E.fr.isNegative(coef)) {
                res << " - ";
            } else if (i != this->degree) {
                res << " + ";
            }
            res << E.fr.toString(coef);
            if (i > 0) {
                if (i > 1) {
                    res << "x^" << i;
                } else {
                    res << "x";
                }
            }
        }
    }
    LOG_TRACE(res);
}
