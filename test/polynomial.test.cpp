#include <gtest/gtest.h>
#include <iostream>

#include <alt_bn128.hpp>
#include "../src/polynomial/polynomial.hpp"
#include <immintrin.h>

extern AltBn128::Engine E;

Polynomial<AltBn128::Engine>* getPolynomial(u_int32_t n) {
    Polynomial<AltBn128::Engine>* pol = new Polynomial<AltBn128::Engine>(E, n);
    for (u_int32_t i = 0; i < n; i++)
    {
        pol->setCoef(i, E.fr.set(i + 1));
    }

    return pol;
}

// POLYNOMIAL BASIC OPERATIONS
TEST(POLYNOMIAL, PolynomialConstructor) {
    Polynomial<AltBn128::Engine>* polA = new Polynomial<AltBn128::Engine>(E, 8);

    ASSERT_EQ(polA->getLength(), 8);
    ASSERT_EQ(polA->getDegree(), 0);

    Polynomial<AltBn128::Engine>* polB = new Polynomial<AltBn128::Engine>(E, 8, 2);

    ASSERT_EQ(polB->getLength(), 10);
    ASSERT_EQ(polB->getDegree(), 0);

    AltBn128::FrElement bufferC[8];
    Polynomial<AltBn128::Engine>* polC = new Polynomial<AltBn128::Engine>(E, bufferC, 8);

    ASSERT_EQ(polC->getLength(), 8);
    ASSERT_EQ(polC->getDegree(), 0);

    AltBn128::FrElement bufferD[10];
    Polynomial<AltBn128::Engine>* polD = new Polynomial<AltBn128::Engine>(E, bufferD, 8, 2);

    ASSERT_EQ(polD->getLength(), 10);
    ASSERT_EQ(polD->getDegree(), 0);
}

TEST(POLYNOMIAL, fromPolynomial) {
    auto polA = getPolynomial(8);
    auto polB = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);

    ASSERT_EQ(polA->getLength(), polB->getLength());
    ASSERT_EQ(polA->getDegree(), polB->getDegree());
    ASSERT_TRUE(polA->isEqual(*polB));

    polA->setCoef(0, E.fr.set(33));
    ASSERT_FALSE(polA->isEqual(*polB));

    polA = getPolynomial(8);
    auto polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA, 2);

    ASSERT_EQ(polA->getLength() + 2, polC->getLength());
    ASSERT_EQ(polA->getDegree(), polC->getDegree());
    ASSERT_TRUE(polA->isEqual(*polC));

    polC->setCoef(9, E.fr.set(33));
    ASSERT_FALSE(polA->isEqual(*polC));

    polA = getPolynomial(8);
    AltBn128::FrElement bufferD[8];
    auto polD = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA, bufferD);

    ASSERT_EQ(polA->getLength(), polD->getLength());
    ASSERT_EQ(polA->getDegree(), polD->getDegree());
    ASSERT_TRUE(polA->isEqual(*polD));

    polA->setCoef(0, E.fr.set(33));
    ASSERT_FALSE(polA->isEqual(*polD));

    polA = getPolynomial(8);
    AltBn128::FrElement bufferE[10];
    auto polE = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA, bufferE, 2);

    ASSERT_EQ(polA->getLength() + 2, polE->getLength());
    ASSERT_EQ(polA->getDegree(), polE->getDegree());
    ASSERT_TRUE(polA->isEqual(*polE));

    polE->setCoef(9, E.fr.set(33));
    ASSERT_FALSE(polA->isEqual(*polE));
}

TEST(POLYNOMIAL, PolynomialFromReservedBuffer) {
    AltBn128::FrElement buffer[8];
    Polynomial pol = Polynomial(E, &buffer[0], 8);

    ASSERT_EQ(pol.getLength(), 8);
    ASSERT_EQ(pol.getDegree(), 0);
}

TEST(POLYNOMIAL, isEqual) {
    AltBn128::FrElement buffer[8] = {
        E.fr.set(1),
        E.fr.set(2),
        E.fr.set(3),
        E.fr.set(4),
        E.fr.set(5),
        E.fr.set(6),
        E.fr.zero(), E.fr.zero()
    };

    auto polA = Polynomial<AltBn128::Engine>::fromCoefficients(E, buffer, 8);
    auto polB = Polynomial<AltBn128::Engine>::fromCoefficients(E, buffer, 8);

    ASSERT_TRUE(polA->isEqual(*polB));
    ASSERT_TRUE(polB->isEqual(*polA));

    polB->setCoef(7, E.fr.set(33));

    ASSERT_FALSE(polA->isEqual(*polB));
    ASSERT_FALSE(polB->isEqual(*polA));
}

TEST(POLYNOMIAL, getCoef) {
    auto pol = getPolynomial(8);

    ASSERT_TRUE(E.fr.eq(E.fr.zero(), pol->getCoef(8)));
    ASSERT_TRUE(E.fr.eq(pol->coef[6], pol->getCoef(6)));
}

TEST(POLYNOMIAL, setCoef) {
    auto polA = getPolynomial(8);

    EXPECT_ANY_THROW(polA->setCoef(8, E.fr.one()));
    EXPECT_NO_THROW(polA->setCoef(7, E.fr.one()));
    EXPECT_NO_THROW(polA->setCoef(0, E.fr.one()));

    AltBn128::FrElement buffer[8] = {
    E.fr.set(1),
    E.fr.set(2),
    E.fr.set(3),
    E.fr.set(4),
    E.fr.set(5),
    E.fr.set(6),
    E.fr.zero(), E.fr.zero()
};

    auto polB = Polynomial<AltBn128::Engine>::fromCoefficients(E, buffer, 8);

    ASSERT_EQ(polB->getDegree(), 5);

    polB->setCoef(7, E.fr.one());

    ASSERT_EQ(polB->getDegree(), 7);

    polB->setCoef(7, E.fr.zero());

    ASSERT_EQ(polB->getDegree(), 5);
}

TEST(POLYNOMIAL, fixDegree) {
    AltBn128::FrElement buffer[8] = {
        E.fr.set(1),
        E.fr.set(2),
        E.fr.set(3),
        E.fr.set(4),
        E.fr.set(5),
        E.fr.set(6),
        E.fr.zero(), E.fr.zero()
    };

    auto polA = Polynomial<AltBn128::Engine>::fromCoefficients(E, buffer, 8);

    ASSERT_EQ(polA->getDegree(), 5);

    polA->fixDegreeFrom(7);

    ASSERT_EQ(polA->getDegree(), 5);

    polA->fixDegreeFrom(5);

    ASSERT_EQ(polA->getDegree(), 5);

    polA->setCoef(7, E.fr.set(33));

    ASSERT_EQ(polA->getDegree(), 7);

    for(int i = 0; i < 8; i++) {
        polA->setCoef(i, E.fr.zero());
    }

    ASSERT_EQ(polA->getDegree(), 0);
}

TEST(POLYNOMIAL, AddWhenPowerOfTwoDegree) {
    // Check correctness of a+b polynomial operation where a degree is less than b degree
    auto polA = getPolynomial(8);
    auto polB = getPolynomial(16);
    auto polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->add(*polB);

    for(int i = 0; i < 16; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.add(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 15);

    delete polA;
    delete polB;
    delete polC;

    // Check correctness of a+b polynomial operation where b degree is less than a degree
    polA = getPolynomial(16);
    polB = getPolynomial(8);
    polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->add(*polB);


    for(int i = 0; i < 16; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.add(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 15);

    delete polA;
    delete polB;
    delete polC;

    // Check correctness of a+b polynomial operation where b degree is equal than a degree
    polA = getPolynomial(16);
    polB = getPolynomial(16);
    polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->add(*polB);

    for(int i = 0; i < 16; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.add(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 15);

    delete polA;
    delete polB;
    delete polC;
}

TEST(POLYNOMIAL, AddWhenNOTPowerOfTwoDegree) {
    // Check correctness of a+b polynomial operation where a degree is less than b degree and is not a power of two
    auto polA = getPolynomial(7);
    auto polB = getPolynomial(8);
    auto polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->add(*polB);

    for(int i = 0; i < 8; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.add(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 7);

    delete polA;
    delete polB;
    delete polC;

    // Check correctness of a+b polynomial operation where b degree is less than a degree and is not a power of two
    polA = getPolynomial(8);
    polB = getPolynomial(7);
    polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->add(*polB);


    for(int i = 0; i < 8; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.add(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 7);

    delete polA;
    delete polB;
    delete polC;

    // Check correctness of a+b polynomial operation where b degree is equal than a degree and is not a power of two
    polA = getPolynomial(7);
    polB = getPolynomial(7);
    polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->add(*polB);

    for(int i = 0; i < 7; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.add(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 6);

    delete polA;
    delete polB;
    delete polC;
}

TEST(POLYNOMIAL, SubtractWhenPowerOfTwoDegree) {
    // Check correctness of a-b polynomial operation where a degree is less than b degree
    auto polA = getPolynomial(8);
    auto polB = getPolynomial(16);
    auto polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->sub(*polB);

    for(int i = 0; i < 16; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.sub(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 15);

    delete polA;
    delete polB;
    delete polC;

    // Check correctness of a-b polynomial operation where b degree is less than a degree
    polA = getPolynomial(16);
    polB = getPolynomial(8);
    polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->sub(*polB);

    for(int i = 0; i < 16; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.sub(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 15);

    delete polA;
    delete polB;
    delete polC;

    // Check correctness of a-b polynomial operation where b degree is equal than a degree
    polA = getPolynomial(16);
    polB = getPolynomial(16);
    polC = Polynomial<AltBn128::Engine>::fromPolynomial(E, *polA);
    polC->sub(*polB);

    for(int i = 0; i < 16; i++) {
        ASSERT_TRUE(E.fr.eq(polC->getCoef(i), E.fr.sub(polA->getCoef(i), polB->getCoef(i))));
    }

    ASSERT_EQ(polC->getDegree(), 0);

    delete polA;
    delete polB;
    delete polC;
}

TEST(POLYNOMIAL, mulScalar) {
    auto pol = getPolynomial(8);
    auto scalar = E.fr.set(7);
    auto polRes = Polynomial<AltBn128::Engine>::fromPolynomial(E, *pol);
    
    polRes->mulScalar(scalar);

    for(int i = 0; i < 8; i++) {
        ASSERT_TRUE(E.fr.eq(polRes->getCoef(i), E.fr.mul(scalar, pol->getCoef(i))));
    }

    ASSERT_EQ(pol->getDegree(), 7);
}

TEST(POLYNOMIAL, addScalar) {
    auto pol = getPolynomial(8);
    auto scalar = E.fr.set(7);
    auto polRes = Polynomial<AltBn128::Engine>::fromPolynomial(E, *pol);
    
    polRes->addScalar(scalar);

    for(int i = 0; i < 8; i++) {
        if(i == 0) {
            ASSERT_TRUE(E.fr.eq(polRes->getCoef(i), E.fr.add(pol->getCoef(i), scalar)));
        } else {
            ASSERT_TRUE(E.fr.eq(polRes->getCoef(i), pol->getCoef(i)));
        }
    }

        ASSERT_EQ(pol->getDegree(), 7);
}

TEST(POLYNOMIAL, subScalar) {
    auto pol = getPolynomial(8);
    auto scalar = E.fr.set(7);
    auto polRes = Polynomial<AltBn128::Engine>::fromPolynomial(E, *pol);
    
    polRes->subScalar(scalar);

    for(int i = 0; i < 8; i++) {
        if(i == 0) {
            ASSERT_TRUE(E.fr.eq(polRes->getCoef(i), E.fr.sub(pol->getCoef(i), scalar)));
        } else {
            ASSERT_TRUE(E.fr.eq(polRes->getCoef(i), pol->getCoef(i)));
        }
    }

    ASSERT_EQ(pol->getDegree(), 7);
}

TEST(POLYNOMIAL, byXSubValue) {
    AltBn128::FrElement bufferA[8] = {
        E.fr.set(1),
        E.fr.set(2),
        E.fr.set(3),
        E.fr.set(4),
        E.fr.set(5),
        E.fr.zero(), E.fr.zero(), E.fr.zero()
    };

    AltBn128::FrElement bufferB[8] = {
        E.fr.set(-7),
        E.fr.set(-13),
        E.fr.set(-19),
        E.fr.set(-25),
        E.fr.set(-31),
        E.fr.set(5),
        E.fr.zero(), E.fr.zero()
    };

    auto polA = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferA, 8);
    auto polB = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferB, 8);

    AltBn128::FrElement element = E.fr.set(7);
    polA->byXSubValue(element);

    ASSERT_TRUE(polA->isEqual(*polB));
}

TEST(POLYNOMIAL, byXNSubValue) {
    AltBn128::FrElement bufferA[8] = {
        E.fr.set(1),
        E.fr.set(2),
        E.fr.set(3),
        E.fr.set(4),
        E.fr.set(5),
        E.fr.zero(), E.fr.zero(), E.fr.zero()
    };

    AltBn128::FrElement bufferB[8] = {
        E.fr.set(-7),
        E.fr.set(-13),
        E.fr.set(-19),
        E.fr.set(-25),
        E.fr.set(-31),
        E.fr.set(5),
        E.fr.zero(), E.fr.zero()
    };

    auto polA = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferA, 8);
    auto polB = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferB, 8);

    AltBn128::FrElement element = E.fr.set(7);
    polA->byXNSubValue(1, element);

    ASSERT_TRUE(polA->isEqual(*polB));

    AltBn128::FrElement bufferC[8] = {
        E.fr.set(-7),
        E.fr.set(-14),
        E.fr.set(-20),
        E.fr.set(-26),
        E.fr.set(-32),
        E.fr.set(4),
        E.fr.set(5), E.fr.zero()
    };

    delete polA;
    polA = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferA, 8);
    auto polC = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferC, 8);

    polA->byXNSubValue(2, element);

    ASSERT_TRUE(polA->isEqual(*polC));
}

TEST(POLYNOMIAL, divByXSubValue) {
    AltBn128::FrElement bufferA[8] = {
        E.fr.set(1),
        E.fr.set(2),
        E.fr.set(3),
        E.fr.set(4),
        E.fr.set(5),
        E.fr.zero(), E.fr.zero(), E.fr.zero()
    };

    AltBn128::FrElement bufferB[8] = {
        E.fr.set(-7),
        E.fr.set(-13),
        E.fr.set(-19),
        E.fr.set(-25),
        E.fr.set(-31),
        E.fr.set(5),
        E.fr.zero(), E.fr.zero()
    };

    auto polA = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferA, 8);
    auto polB = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferB, 8);

    AltBn128::FrElement element = E.fr.set(7);
    polB->divByXSubValue(element);

    ASSERT_TRUE(polA->isEqual(*polB));
}

TEST(POLYNOMIAL, divZh) {
    AltBn128::FrElement bufferA[8] = {
        E.fr.set(-1),
        E.fr.set(-2),
        E.fr.set(-2),
        E.fr.set(-2),
        E.fr.set(3),
        E.fr.set(4),
        E.fr.zero(), E.fr.zero()
    };

    AltBn128::FrElement bufferB[4] = {
        E.fr.set(1),
        E.fr.set(2),
        E.fr.set(3),
        E.fr.set(4)
    };

    auto polA = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferA, 8);
    auto polB = Polynomial<AltBn128::Engine>::fromCoefficients(E, bufferB, 4);

    polA->divZh(2, 1);

    ASSERT_TRUE(polA->isEqual(*polB));

    auto polC = getPolynomial(64);
    auto polD = getPolynomial(64);

    auto one = E.fr.one();
    polC->byXNSubValue(16, one);

    polC->divZh(16, 1);

    ASSERT_TRUE(polC->isEqual(*polD));
}

TEST(POLYNOMIAL, divByZerofier) {
    auto polA = getPolynomial(64);
    auto polB = getPolynomial(64);

    auto seven = E.fr.set(7);
    polA->byXNSubValue(16, seven);

    polA->divByZerofier(16, seven);

    ASSERT_TRUE(polA->isEqual(*polB));
}

TEST(POLYNOMIAL, byX) {
    auto polA = getPolynomial(8);
    auto polB = getPolynomial(8);

    polA->byX();

    ASSERT_TRUE(E.fr.isZero(polA->getCoef(0)));

    for(int i = 1; i < 8; i++) {
        ASSERT_TRUE(E.fr.eq(polA->getCoef(i), polB->getCoef(i-1)));
    }

    polA = getPolynomial(8);
    polA->setCoef(7, E.fr.zero());

    polA->byX();

    ASSERT_TRUE(E.fr.isZero(polA->getCoef(0)));

    for(int i = 1; i < 8; i++) {
        ASSERT_TRUE(E.fr.eq(polA->getCoef(i), polB->getCoef(i-1)));
    }

}
