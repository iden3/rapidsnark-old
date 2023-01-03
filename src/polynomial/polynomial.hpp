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

    Engine &E;
    FFT<typename Engine::Fr> *fft;

    void initialize(u_int64_t length);

public:
    FrElement *coef;

    u_int64_t length;
    u_int64_t degree;

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
//
//    unsigned long length();
//
//    unsigned long degree();
//
//    FrElement evaluate(FrElement x);

    void add(Polynomial<FrElement> &other, FrElement &blindingValue);

//    void sub(Polynomial<FrElement>, FrElement blindingValue);
//
//    void mulScalar(FrElement value);
//
//    void addScalar(FrElement value);
//
//    void subScalar(FrElement value);
//
//    // Multiply current polynomial by the polynomial (X - value)
//    void byXSubValue(FrElement value);
//
//    // Euclidean division
//    void divBy(Polynomial<FrElement> polynomial);
//
//    // Divide polynomial by X - value
//    void divByXSubValue(FrElement value);

    void divZh();

//    void byX();
//
//    // Compute a new polynomial f(x^n) from f(x)
//    // f(x)   = a_0 + a_1·x + a_2·x^2 + ... + a_j·x^j
//    // f(x^n) = a_0 + a_1·x^n + a_2·x^2n + ... + a_j·x^jn
//    static expX(Polynomial<FrElement> polynomial, unsigned long n, bool truncate = false);
//
//    std::vector<Polynomial<FrElement>> split(unsigned int numPols, unsigned long degPols, std::vector<FrElement> blindingFactors);
//
//    void truncate();
//
//    static lagrangePolynomialInterpolation(std::vector<FrElement> xArr, std::vector<FrElement> yArr);
//
//    static zerofierPolynomial(std::vector<FrElement> xArr);
//
//    void print();
};

#endif //POLYNOMIAL_HPP
