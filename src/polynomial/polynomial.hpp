#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include "assert.h"
#include <sstream>
#include <gmp.h>
#include "fft.hpp"


template<typename Engine>
class Polynomial {
    using FrElement = typename Engine::FrElement;
    using G1Point = typename Engine::G1Point;

    u_int64_t length;
    u_int64_t degree;

    Engine &E;
    FFT<typename Engine::Fr> *fft;

    void initialize(u_int64_t length, u_int64_t blindLength = 0);

    static Polynomial<Engine>* computeLagrangePolynomial(u_int64_t i, FrElement xArr[], FrElement yArr[], u_int32_t length);
public:
    FrElement *coef;

    Polynomial(Engine &_E, u_int64_t length, u_int64_t blindLength = 0);

    // From coefficients
    static Polynomial<Engine>* fromCoefficients(Engine &_E, FrElement *coefficients, u_int64_t length, u_int64_t blindLength = 0);

    // From evaluations
    static Polynomial<Engine>* fromEvaluations(Engine &_E, FrElement *evaluations, u_int64_t length, u_int64_t blindLength = 0);

    ~Polynomial();

    void fixDegree();

    bool isEqual(const Polynomial<Engine> &other) const;

    void blindCoefficients(FrElement blindingFactors[], u_int32_t length);

    typename Engine::FrElement getCoef(u_int64_t index) const;

    void setCoef(u_int64_t index, FrElement value);

    u_int64_t getLength() const;

    u_int64_t getDegree() const;

    typename Engine::FrElement evaluate(FrElement point) const;

    void add(Polynomial<Engine> &polynomial);

    void sub(Polynomial<Engine> &polynomial);

    void mulScalar(FrElement &value);

    void addScalar(FrElement &value);

    void subScalar(FrElement &value);

    // Multiply current polynomial by the polynomial (X - value)
    void byXSubValue(FrElement &value);

    // Euclidean division, returns reminder polygon
    Polynomial<Engine>* divBy(Polynomial<Engine> &polynomial);

    void divZh(u_int64_t domainSize);

    void byX();

    static Polynomial<Engine>* lagrangePolynomialInterpolation(FrElement xArr[], FrElement yArr[], u_int32_t length);

    static Polynomial<Engine>* zerofierPolynomial(FrElement xArr[], u_int32_t length);

    void print();
};

#include "polynomial.cpp"

#endif
