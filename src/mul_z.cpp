#include "mul_z.hpp"


template<typename Engine>
MulZ<Engine>::MulZ(Engine &_E, FFT<typename Engine::Fr> *fft) : E(_E) {
    FrElement w2 = fft->root(2, 1);

    FrElement negTwo, negEight, two, four, tmp;
    E.fr.fromString(negTwo, "-2");
    E.fr.fromString(negEight, "-8");
    E.fr.fromString(two, "2");
    E.fr.fromString(four, "4");

    Z1[0] = E.fr.zero();
    E.fr.add(Z1[1], E.fr.negOne(), w2);
    Z1[2] = negTwo;
    E.fr.sub(Z1[3], E.fr.negOne(), w2);

    Z2[0] = E.fr.zero();
    E.fr.mul(tmp, negTwo, w2);
    E.fr.add(Z2[1], E.fr.zero(), tmp);
    Z2[2] = four;
    E.fr.mul(tmp, negTwo, w2);
    E.fr.sub(Z2[3], E.fr.zero(), tmp);

    Z3[0] = E.fr.zero();
    E.fr.mul(tmp, two, w2);
    E.fr.add(Z3[1], two, tmp);
    Z3[2] = negEight;
    E.fr.mul(tmp, two, w2);
    E.fr.sub(Z3[3], two, tmp);
}

template<typename Engine>
tuple<typename Engine::FrElement, typename Engine::FrElement> MulZ<Engine>::mul2(FrElement &a,
                                                                                 FrElement &b,
                                                                                 FrElement &ap,
                                                                                 FrElement &bp,
                                                                                 int64_t p) {
    FrElement a_b, a_bp, ap_b, ap_bp;

    E.fr.mul(a_b, a, b);
    E.fr.mul(a_bp, a, bp);
    E.fr.mul(ap_b, ap, b);
    E.fr.mul(ap_bp, ap, bp);

    FrElement r = a_b;

    FrElement a0, a1;
    E.fr.add(a0, a_bp, ap_b);
    a1 = ap_bp;

    FrElement rz = a0;
    if (p >= 0) {
        FrElement tmp;
        E.fr.mul(tmp, Z1[p], a1);
        E.fr.add(rz, rz, tmp);
    }

    return make_tuple(r, rz);
}

template<typename Engine>
tuple<typename Engine::FrElement, typename Engine::FrElement> MulZ<Engine>::mul4( FrElement &a,
                                                                                  FrElement &b,
                                                                                  FrElement &c,
                                                                                  FrElement &d,
                                                                                  FrElement &ap,
                                                                                  FrElement &bp,
                                                                                  FrElement &cp,
                                                                                  FrElement &dp,
                                                                                 int64_t p) {
    FrElement a_b, a_bp, ap_b, ap_bp;
    FrElement c_d, c_dp, cp_d, cp_dp;

    E.fr.mul(a_b, a, b);
    E.fr.mul(a_bp, a, bp);
    E.fr.mul(ap_b, ap, b);
    E.fr.mul(ap_bp, ap, bp);

    E.fr.mul(c_d, c, d);
    E.fr.mul(c_dp, c, dp);
    E.fr.mul(cp_d, cp, d);
    E.fr.mul(cp_dp, cp, dp);

    FrElement r;
    E.fr.mul(r, a_b, c_d);

    FrElement a0;
    FrElement tmp;

    E.fr.mul(a0, ap_b, c_d);
    E.fr.mul(tmp, a_bp, c_d);
    E.fr.add(a0, a0, tmp);
    E.fr.mul(tmp, a_b, cp_d);
    E.fr.add(a0, a0, tmp);
    E.fr.mul(tmp, a_b, c_dp);
    E.fr.add(a0, a0, tmp);

    FrElement a1;
    E.fr.mul(a1, ap_bp, c_d);
    E.fr.mul(tmp, ap_b, cp_d);
    E.fr.add(a1, a1, tmp);
    E.fr.mul(tmp, ap_b, c_dp);
    E.fr.add(a1, a1, tmp);
    E.fr.mul(tmp, a_bp, cp_d);
    E.fr.add(a1, a1, tmp);
    E.fr.mul(tmp, a_bp, c_dp);
    E.fr.add(a1, a1, tmp);
    E.fr.mul(tmp, a_b, cp_dp);
    E.fr.add(a1, a1, tmp);

    FrElement a2;
    E.fr.mul(a2, a_bp, cp_dp);
    E.fr.mul(tmp, ap_b, cp_dp);
    E.fr.add(a2, a2, tmp);
    E.fr.mul(tmp, ap_bp, c_dp);
    E.fr.add(a2, a2, tmp);
    E.fr.mul(tmp, ap_bp, cp_d);
    E.fr.add(a2, a2, tmp);

    FrElement a3;
    E.fr.mul(a3, ap_bp, cp_dp);

    FrElement rz = a0;
    if (p >= 0) {
        E.fr.mul(tmp, Z1[p], a1);
        E.fr.add(rz, rz, tmp);
        E.fr.mul(tmp, Z2[p], a2);
        E.fr.add(rz, rz, tmp);
        E.fr.mul(tmp, Z3[p], a3);
        E.fr.add(rz, rz, tmp);
    }

    return make_tuple(r, rz);
}