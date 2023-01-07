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

    void initialize(u_int64_t length);

    static Polynomial<Engine> computeLagrangePolynomial(u_int64_t i, Polynomial<Engine> &pol, FrElement xArr[], FrElement yArr[]);
public:
    FrElement *coef;

    Polynomial(u_int64_t length);

    // From coefficients
    Polynomial(FrElement elements[]);

    // From evaluations
    Polynomial(FrElement elements[], u_int64_t domainSize, u_int64_t nBlindCoefficients = 0);

    ~Polynomial();

    void fixDegree();

    bool isEqual(const Polynomial<Engine> &other) const;

    void blindCoefficients(FrElement blindingFactors[]);

    typename Engine::FrElement getCoef(u_int64_t index) const;

    void setCoef(u_int64_t index, FrElement value);

//    static to4T(buffer, domainSize, blindingFactors, Fr); //TODO return dues coses!!!!!!

    u_int64_t getLength() const;

    u_int64_t getDegree() const;

    typename Engine::FrElement evaluate(FrElement point) const;

    void add(Polynomial<FrElement> &polynomial);

    void sub(Polynomial<FrElement> &polynomial);

    void mulScalar(FrElement &value);

    void addScalar(FrElement &value);

    void subScalar(FrElement &value);

    // Multiply current polynomial by the polynomial (X - value)
    void byXSubValue(FrElement &value);

    // Euclidean division
    Polynomial<Engine> divBy(Polynomial<Engine> &polynomial);

    void divZh(u_int64_t domainSize);

    void byX();

    static Polynomial<Engine> lagrangePolynomialInterpolation(Polynomial<Engine> &pol, FrElement xArr[], FrElement yArr[]);

    static Polynomial<Engine> zerofierPolynomial(FrElement xArr[]);

    void print();
};

#endif
