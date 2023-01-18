#include "polynomial.hpp"
#include "thread_utils.hpp"

#include "logger.hpp"

using namespace CPlusPlusLogging;

template<typename Engine>
void Polynomial<Engine>::initialize(u_int64_t length, u_int64_t blindLength) {
    u_int64_t totalLength = length + blindLength;
    coef = new FrElement[totalLength];
    memset(coef, 0, totalLength * sizeof(FrElement));
    this->length = totalLength;
    degree = 0;
}

template<typename Engine>
Polynomial<Engine>::Polynomial(Engine &_E, u_int64_t length, u_int64_t blindLength) : E(_E) {
    this->initialize(length, blindLength);
}

template<typename Engine>
Polynomial<Engine>* Polynomial<Engine>::fromCoefficients(Engine &_E, FrElement* coefficients, u_int64_t length, u_int64_t blindLength) {
    Polynomial<Engine>* pol = new Polynomial<Engine>(_E, length, blindLength);

    //memcpy(pol->coef, coefficients, length * sizeof(FrElement));
    int nThreads = omp_get_max_threads() / 2;
    ThreadUtils<Engine>::parcpy(pol->coef, coefficients, length * sizeof(FrElement), nThreads);
    pol->fixDegree();

    return pol;
}

template<typename Engine>
Polynomial<Engine>* Polynomial<Engine>::fromEvaluations(Engine &_E, FFT<typename Engine::Fr> *fft, FrElement* evaluations, u_int64_t length, u_int64_t blindLength) {
    Polynomial<Engine>* pol = new Polynomial<Engine>(_E, length, blindLength);

    int nThreads = omp_get_max_threads() / 2;
    ThreadUtils<Engine>::parcpy(pol->coef, evaluations, length * sizeof(FrElement), nThreads);

    fft->ifft(pol->coef, length);

    pol->fixDegree();

    return pol;
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
void Polynomial<Engine>::blindCoefficients(FrElement blindingFactors[], u_int32_t length) {
    const u_int32_t polLength = this->length;

    for (int i = 0; i < length; i++) {
        coef[polLength - length + i] = E.fr.add(coef[polLength - length + i], blindingFactors[i]);
        coef[i] = E.fr.sub(coef[i], blindingFactors[i]);
    }
    fixDegree();
}

template<typename Engine>
typename Engine::FrElement Polynomial<Engine>::getCoef(u_int64_t index) const {
    if (index > length) {
        return E.fr.zero();
        //throw std::runtime_error("Polynomial::getCoef: invalid index");
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
    FrElement result = E.fr.zero();

    for (u_int64_t i = degree + 1; i > 0; i--) {
        result = E.fr.add(coef[i - 1], E.fr.mul(result, point));
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
        FrElement a = i < thisLength ? this->coef[i] : E.fr.zero();
        FrElement b = i < polyLength ? polynomial.coef[i] : E.fr.zero();
        FrElement sum;
        E.fr.add(sum, a, b);

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
        E.fr.mul(this->coef[i], this->coef[i], value);
    }
}

template<typename Engine>
void Polynomial<Engine>::addScalar(FrElement &value) {
    FrElement currentValue = (0 == this->length) ? E.fr.zero : this->coef[0];
    E.fr.add(this->coef[0], currentValue, value);
}

template<typename Engine>
void Polynomial<Engine>::subScalar(FrElement &value) {
    FrElement currentValue = (0 == this->length) ? E.fr.zero : this->coef[0];
    E.fr.sub(this->coef[0], currentValue, value);
}

// Multiply current polynomial by the polynomial (X - value)
template<typename Engine>
void Polynomial<Engine>::byXSubValue(FrElement &value) {
    bool resize = !E.fr.eq(E.fr.zero(), this->coef[this->length - 1]);

    u_int64_t length = resize ? this->length + 1 : this->length;
    Polynomial<Engine> *pol = new Polynomial<Engine>(E, length);

    // Step 0: Set current coefficients to the new buffer shifted one position
    int nThreads = omp_get_max_threads() / 2;
    ThreadUtils<Engine>::parcpy(&pol->coef[1], this->coef, (resize ? this->length : this->length -1)  * sizeof(FrElement), nThreads);
    pol->fixDegree();

    // Step 1: multiply each coefficient by (-value)
    FrElement negValue = E.fr.neg(value);
    this->mulScalar(negValue);

    // Step 2: Add current polynomial to destination polynomial
    pol->add(*this);

    // Swap buffers
    delete[] this->coef;
    this->coef = pol->coef;

    fixDegree();
}

// Euclidean division
template<typename Engine>
Polynomial<Engine>* Polynomial<Engine>::divBy(Polynomial<Engine> &polynomial) {
    u_int64_t degreeA = this->degree;
    u_int64_t degreeB = polynomial.degree;

    Polynomial<Engine> *polR = fromCoefficients(this->E, this->coef, this->length);

    delete[] this->coef;
    initialize(this->length);

    for (int64_t i = degreeA - degreeB; i >= 0; i--) {
        E.fr.div(this->coef[i], polR->getCoef(i + degreeB), polynomial.getCoef(degreeB));

        for (u_int64_t j = 0; j <= degreeB; j++) {
            polR->coef[i + j] = E.fr.sub(polR->getCoef(i + j), E.fr.mul(this->coef[i], polynomial.coef[j]));
        }
    }

    fixDegree();
    polR->fixDegree();

    return polR;
}

template<typename Engine>
void Polynomial<Engine>::divZh(u_int64_t domainSize) {
    for (u_int64_t i = 0; i < domainSize; i++) {
        E.fr.neg(this->coef[i], this->coef[i]);
    }

    for (u_int64_t i = domainSize; i < this->length; i++) {
        E.fr.sub(coef[i], coef[i - domainSize], coef[i]);

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
    int nThreads = omp_get_max_threads() / 2;

    if (resize) {
        FrElement *newCoef = new FrElement[this->length + 1];
        ThreadUtils<Engine>::parcpy(newCoef[1], coef[0], sizeof(coef), nThreads);
//        memcpy(newCoef[1], coef[0], sizeof(coef));
        coef = newCoef;
    } else {
        ThreadUtils<Engine>::parcpy(coef[1], coef[0], sizeof(coef), nThreads);
//        memcpy(coef[1], coef[0], sizeof(coef));
    }

    coef[0] = E.fr.zero;
    fixDegree();
}

template<typename Engine>
Polynomial<Engine>*
Polynomial<Engine>::lagrangePolynomialInterpolation(FrElement xArr[], FrElement yArr[], u_int32_t length) {
    Polynomial<Engine> *polynomial = computeLagrangePolynomial(0, xArr, yArr, length);

    for (u_int64_t i = 1; i < length; i++) {
        Polynomial<Engine> *polynomialI = computeLagrangePolynomial(i, xArr, yArr, length);
        polynomial->add(*polynomialI);
    }

    return polynomial;
}

template<typename Engine>
Polynomial<Engine>* Polynomial<Engine>::computeLagrangePolynomial(u_int64_t i, FrElement xArr[], FrElement yArr[], u_int32_t length) {
    Engine &E = Engine::engine;
    Polynomial<Engine> *polynomial = NULL;

    for (u_int64_t j = 0; j < length; j++) {
        if (j == i) continue;

        if (NULL == polynomial) {
            polynomial = new Polynomial<Engine>(E, length + 1);
            polynomial->coef[0] = E.fr.neg(xArr[j]);
            polynomial->coef[1] = E.fr.one();
            polynomial->fixDegree();

        } else {
            polynomial->byXSubValue(xArr[j]);
        }
    }

    FrElement denominator = polynomial->evaluate(xArr[i]);
    E.fr.inv(denominator, denominator);
    FrElement mulFactor = E.fr.mul(yArr[i], denominator);

    polynomial->mulScalar(mulFactor);

    return polynomial;
}

template<typename Engine>
Polynomial<Engine>* Polynomial<Engine>::zerofierPolynomial(FrElement xArr[], u_int32_t length) {
    Engine &E = Engine::engine;
    Polynomial<Engine> *polynomial = new Polynomial<Engine>(E, length + 1);

    // Build a zerofier polynomial with the following form:
    // zerofier(X) = (X-xArr[0])(X-xArr[1])...(X-xArr[n])
    E.fr.neg(polynomial->coef[0], xArr[0]);
    polynomial->coef[1] = E.fr.one();

    polynomial->fixDegree();

    for (u_int64_t i = 1; i < length; i++) {
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