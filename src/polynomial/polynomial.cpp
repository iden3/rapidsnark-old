#include "polynomial.hpp"
#include "thread_utils.hpp"
#include "evaluations.hpp"
#include <math.h>
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
Polynomial<Engine> *
Polynomial<Engine>::fromCoefficients(Engine &_E, FrElement *coefficients, u_int64_t length, u_int64_t blindLength) {
    Polynomial<Engine> *pol = new Polynomial<Engine>(_E, length, blindLength);

    int nThreads = omp_get_max_threads() / 2;
    ThreadUtils::parcpy(pol->coef, coefficients, length * sizeof(FrElement), nThreads);
    pol->fixDegree();

    return pol;
}

template<typename Engine>
Polynomial<Engine> *
Polynomial<Engine>::fromEvaluations(Engine &_E, FFT<typename Engine::Fr> *fft, FrElement *evaluations, u_int64_t length,
                                    u_int64_t blindLength) {
    Polynomial<Engine> *pol = new Polynomial<Engine>(_E, length, blindLength);

    int nThreads = omp_get_max_threads() / 2;
    ThreadUtils::parcpy(pol->coef, evaluations, length * sizeof(FrElement), nThreads);

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
typename Engine::FrElement Polynomial<Engine>::fastEvaluate(FrElement point) const {
    int nThreads = omp_get_max_threads();

    uint64_t nCoefs = this->degree + 1;
    uint64_t coefsThread = nCoefs / nThreads;
    uint64_t residualCoefs = nCoefs - coefsThread * nThreads;

    FrElement res[nThreads];
    FrElement xN[nThreads];

    xN[0] = E.fr.one();

#pragma omp parallel for
    for (int i = 0; i < nThreads; i++) {
        res[i] = E.fr.zero();

        uint64_t nCoefs = i == (nThreads - 1) ? coefsThread + residualCoefs : coefsThread;
        for (u_int64_t j = nCoefs; j > 0; j--) {
            res[i] = E.fr.add(coef[(i * coefsThread) + j - 1], E.fr.mul(res[i], point));

            if (i == 0) xN[0] = E.fr.mul(xN[0], point);
        }
    }

    for (uint i = 1; i < nThreads; i++) {
        res[0] = E.fr.add(res[0], E.fr.mul(xN[i - 1], res[i]));
        xN[i] = E.fr.mul(xN[i - 1], xN[0]);
    }

    return res[0];
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

#pragma omp parallel for
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

#pragma omp parallel for
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
#pragma omp parallel for
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
    ThreadUtils::parcpy(&pol->coef[1], this->coef, (resize ? this->length : this->length - 1) * sizeof(FrElement),
                        nThreads);
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
Polynomial<Engine> *Polynomial<Engine>::divBy(Polynomial<Engine> &polynomial) {
    u_int64_t degreeA = this->degree;
    u_int64_t degreeB = polynomial.degree;

    Polynomial<Engine> *polR = new Polynomial<Engine>(this->E, this->length);
    FrElement *ptr = this->coef;
    this->coef = polR->coef;
    polR->coef = ptr;

    FrElement val = polynomial.coef[degreeB];
    for (int64_t i = degreeA - degreeB; i >= 0; i--) {
        E.fr.div(this->coef[i], polR->coef[i + degreeB], val);

        for (u_int64_t j = 0; j <= degreeB; j++) {
            polR->coef[i + j] = E.fr.sub(polR->coef[i + j], E.fr.mul(this->coef[i], polynomial.coef[j]));
        }
    }

    fixDegree();
    polR->fixDegree();

    return polR;
}

template<typename Engine>
void Polynomial<Engine>::divByMonic(uint32_t m, FrElement beta) {
    u_int64_t d = this->degree;

    Polynomial<Engine> *polResult = new Polynomial<Engine>(this->E, this->length);
    FrElement *bArr = new FrElement[m];

    for (uint32_t i = 0; i < m; i++) {
        polResult->coef[(d - i) - m] = this->coef[d - i];
        bArr[i] = this->coef[d - i];
    }

    uint32_t nThreads = m;
    #pragma omp parallel for
    for(int k = 0; k<nThreads;k++) {
        for (int i = d - 2 * m - k; i >= 0; i = i - nThreads) {
            if(i<0) break;
            uint32_t idx = k;
            bArr[idx] = E.fr.add(this->coef[i + m], E.fr.mul(bArr[idx], beta));
            polResult->coef[i] = bArr[idx];
        }
    }

    // Swap buffers
    delete[] this->coef;
    this->coef = polResult->coef;

    fixDegree();
}

template<typename Engine>
Polynomial<Engine> *Polynomial<Engine>::divByVanishing(uint32_t m, FrElement beta) {
    if(this->degree < m) {
        throw std::runtime_error("divByVanishing polynomial divisor must be of degree lower than the dividend polynomial");
    }

    Polynomial<Engine> *polR = new Polynomial<Engine>(this->E, this->length);
    FrElement *ptr = this->coef;
    this->coef = polR->coef;
    polR->coef = ptr;

    #pragma omp parallel for
    for (int k = 0; k < m; k++) {
        for (int32_t i = this->length - 1 - k; i >= m; i = i - m) {
            FrElement leadingCoef = polR->coef[i];
            if (E.fr.eq(E.fr.zero(), leadingCoef)) continue;

            polR->coef[i] = E.fr.zero();
            polR->coef[i - m] = E.fr.add(polR->coef[i - m], E.fr.mul(beta, leadingCoef));
            this->coef[i - m] = E.fr.add(this->coef[i - m], leadingCoef);
        }
    }

    fixDegree();
    polR->fixDegree();

    return polR;
}

template<typename Engine>
void Polynomial<Engine>::divZh(u_int64_t domainSize) {
#pragma omp parallel for
    for (u_int64_t i = 0; i < domainSize; i++) {
        E.fr.neg(this->coef[i], this->coef[i]);
    }

    int nThreads = pow(2, log2(omp_get_max_threads()));
    int nElementsThread = domainSize / nThreads;
    int nChunks = this->length / domainSize;

    for (uint i = 0; i < nChunks - 1; i++) {
#pragma omp parallel for
        for (uint k = 0; k < nThreads; k++) {
            for (uint j = 0; j < nElementsThread; j++) {
                int id = k;
                u_int64_t idxBase = id * nElementsThread + j;
                u_int64_t idx0 = idxBase + i * domainSize;
                u_int64_t idx1 = idxBase + (i + 1) * domainSize;
                E.fr.sub(coef[idx1], coef[idx0], coef[idx1]);

                if (i > (domainSize * 3 - 4)) {
                    if (!E.fr.isZero(coef[idx1])) {
                        throw std::runtime_error("Polynomial is not divisible");
                    }
                }
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
        ThreadUtils::parcpy(newCoef[1], coef[0], sizeof(coef), nThreads);
        coef = newCoef;
    } else {
        ThreadUtils::parcpy(coef[1], coef[0], sizeof(coef), nThreads);
    }

    coef[0] = E.fr.zero;
    fixDegree();
}

template<typename Engine>
Polynomial<Engine> *
Polynomial<Engine>::lagrangePolynomialInterpolation(FrElement xArr[], FrElement yArr[], u_int32_t length) {
    Polynomial<Engine> *polynomial = computeLagrangePolynomial(0, xArr, yArr, length);

    for (u_int64_t i = 1; i < length; i++) {
        Polynomial<Engine> *polynomialI = computeLagrangePolynomial(i, xArr, yArr, length);
        polynomial->add(*polynomialI);
    }

    return polynomial;
}

template<typename Engine>
Polynomial<Engine> *
Polynomial<Engine>::computeLagrangePolynomial(u_int64_t i, FrElement xArr[], FrElement yArr[], u_int32_t length) {
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

    FrElement denominator = polynomial->fastEvaluate(xArr[i]);
    E.fr.inv(denominator, denominator);
    FrElement mulFactor = E.fr.mul(yArr[i], denominator);

    polynomial->mulScalar(mulFactor);

    return polynomial;
}

template<typename Engine>
Polynomial<Engine> *Polynomial<Engine>::zerofierPolynomial(FrElement xArr[], u_int32_t length) {
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