namespace Groth16 {

template <typename Engine>
std::unique_ptr<Prover<Engine>> makeProver(
    u_int32_t nVars, 
    u_int32_t nPublic, 
    u_int32_t domainSize, 
    u_int64_t nCoeffs, 
    void *vk_alpha1,
    void *vk_beta_1,
    void *vk_beta_2,
    void *vk_delta_1,
    void *vk_delta_2,
    void *coefs, 
    void *pointsA, 
    void *pointsB1, 
    void *pointsB2, 
    void *pointsC, 
    void *pointsH
) {
    Prover<Engine> *p = new Prover<Engine>(
        Engine::engine, 
        nVars, 
        nPublic, 
        domainSize, 
        nCoeffs, 
        *(typename Engine::G1PointAffine *)vk_alpha1,
        *(typename Engine::G1PointAffine *)vk_beta_1,
        *(typename Engine::G2PointAffine *)vk_beta_2,
        *(typename Engine::G1PointAffine *)vk_delta_1,
        *(typename Engine::G2PointAffine *)vk_delta_2,
        (Coef<Engine> *)((uint64_t)coefs + 4), 
        (typename Engine::G1PointAffine *)pointsA,
        (typename Engine::G1PointAffine *)pointsB1,
        (typename Engine::G2PointAffine *)pointsB2,
        (typename Engine::G1PointAffine *)pointsC,
        (typename Engine::G1PointAffine *)pointsH
    );
    return std::unique_ptr< Prover<Engine> >(p);
}

template <typename Engine>
std::unique_ptr<Proof<Engine>> Prover<Engine>::prove(typename Engine::FrElement *wtns) {

    std::cout << "Start Initializing a b c A\n";
    auto a = new typename Engine::FrElement[domainSize];
    auto b = new typename Engine::FrElement[domainSize];
    auto c = new typename Engine::FrElement[domainSize];

    #pragma omp parallel for
    for (u_int32_t i=0; i<domainSize; i++) {
        E.fr.copy(a[i], E.fr.zero());
        E.fr.copy(b[i], E.fr.zero());
    }

    std::cout << "Processing coefs\n";
    #define NLOCKS 1024
    omp_lock_t locks[NLOCKS];
    for (int i=0; i<NLOCKS; i++) omp_init_lock(&locks[i]);
    #pragma omp parallel for 
    for (u_int64_t i=0; i<nCoefs; i++) {
        typename Engine::FrElement *ab = (coefs[i].m == 0) ? a : b;
        typename Engine::FrElement aux;

        E.fr.mul(
            aux,
            wtns[coefs[i].s],
            coefs[i].coef
        );

        omp_set_lock(&locks[coefs[i].c % NLOCKS]);
        E.fr.add(
            ab[coefs[i].c],
            ab[coefs[i].c],
            aux
        );
        omp_unset_lock(&locks[coefs[i].c % NLOCKS]);
    }
    for (int i=0; i<NLOCKS; i++) omp_destroy_lock(&locks[i]);


    std::cout << "Calculating c\n";
    #pragma omp parallel for
    for (u_int32_t i=0; i<domainSize; i++) {
        E.fr.mul(
            c[i],
            a[i],
            b[i]
        );
    }

    std::cout << "Initializing fft\n";
    u_int32_t domainPower = fft->log2(domainSize);

    std::cout << "Start iFFT A\n";
    fft->ifft(a, domainSize);
    cout << "a After ifft:" << endl;
    cout << E.fr.toString(a[0]) << endl;
    cout << E.fr.toString(a[1]) << endl;
    std::cout << "Start Shift A\n";
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(a[i], a[i], fft->root(domainPower+1, i));
    }
    cout << "a After shift:" << endl;
    cout << E.fr.toString(a[0]) << endl;
    cout << E.fr.toString(a[1]) << endl;
    std::cout << "Start FFT A\n";
    fft->fft(a, domainSize);
    cout << "a After fft:" << endl;
    cout << E.fr.toString(a[0]) << endl;
    cout << E.fr.toString(a[1]) << endl;

    std::cout << "Start iFFT B\n";
    fft->ifft(b, domainSize);
    cout << "b After ifft:" << endl;
    cout << E.fr.toString(b[0]) << endl;
    cout << E.fr.toString(b[1]) << endl;
    std::cout << "Start Shift B\n";
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(b[i], b[i], fft->root(domainPower+1, i));
    }
    cout << "b After shift:" << endl;
    cout << E.fr.toString(b[0]) << endl;
    cout << E.fr.toString(b[1]) << endl;
    std::cout << "Start FFT B\n";
    fft->fft(b, domainSize);
    cout << "b After fft:" << endl;
    cout << E.fr.toString(b[0]) << endl;
    cout << E.fr.toString(b[1]) << endl;

    std::cout << "Start iFFT C\n";
    fft->ifft(c, domainSize);
    cout << "c After ifft:" << endl;
    cout << E.fr.toString(c[0]) << endl;
    cout << E.fr.toString(c[1]) << endl;
    std::cout << "Start Shift C\n";
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(c[i], c[i], fft->root(domainPower+1, i));
    }
    cout << "c After shift:" << endl;
    cout << E.fr.toString(c[0]) << endl;
    cout << E.fr.toString(c[1]) << endl;
    std::cout << "Start FFT C\n";
    fft->fft(c, domainSize);
    cout << "c After fft:" << endl;
    cout << E.fr.toString(c[0]) << endl;
    cout << E.fr.toString(c[1]) << endl;

    std::cout << "Start ABC\n";
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(a[i], a[i], b[i]);
        E.fr.sub(a[i], a[i], c[i]);
        E.fr.fromMontgomery(a[i], a[i]);
    }
    cout << "abc:" << endl;
    cout << E.fr.toString(a[0]) << endl;
    cout << E.fr.toString(a[1]) << endl;

    delete b;
    delete c;

    std::cout << "Start Multiexp H\n";
    typename Engine::G1Point pih;
    E.g1.multiMulByScalar(pih, pointsH, (uint8_t *)a, sizeof(a[0]), domainSize);
    cout << "pih: " << E.g1.toString(pih) << "\n";

    delete a;

    std::cout << "Start Multiexp A\n";
    uint32_t sW = sizeof(wtns[0]);
    typename Engine::G1Point pi_a;
    E.g1.multiMulByScalar(pi_a, pointsA, (uint8_t *)wtns, sW, nVars);
    cout << "pi_a: " << E.g1.toString(pi_a) << "\n";

    std::cout << "Start Multiexp B1\n";
    typename Engine::G1Point pib1;
    E.g1.multiMulByScalar(pib1, pointsB1, (uint8_t *)wtns, sW, nVars);
    cout << "pib1: " << E.g1.toString(pib1) << "\n";

    std::cout << "Start Multiexp B2\n";
    typename Engine::G2Point pi_b;
    E.g2.multiMulByScalar(pi_b, pointsB2, (uint8_t *)wtns, sW, nVars);
    cout << "pi_b: " << E.g2.toString(pi_b) << "\n";

    std::cout << "Start Multiexp C\n";
    typename Engine::G1Point pi_c;
    E.g1.multiMulByScalar(pi_c, pointsC, (uint8_t *)((uint64_t)wtns + (nPublic +1)*sW), sW, nVars-nPublic-1);
    cout << "pi_c: " << E.g1.toString(pi_c) << "\n";

    typename Engine::FrElement r;
    typename Engine::FrElement s;
    typename Engine::FrElement rs;

    // TODO set to random
    E.fr.copy(r, E.fr.zero());
    E.fr.copy(s, E.fr.zero());

    typename Engine::G1Point p1;
    typename Engine::G2Point p2;

    E.g1.add(pi_a, pi_a, vk_alpha1);
    E.g1.mulByScalar(p1, vk_delta1, (uint8_t *)&r, sizeof(r));
    E.g1.add(pi_a, pi_a, p1);

    E.g2.add(pi_b, pi_b, vk_beta2);
    E.g2.mulByScalar(p2, vk_delta2, (uint8_t *)&s, sizeof(s));
    E.g2.add(pi_b, pi_b, p2);

    E.g1.add(pib1, pib1, vk_beta1);
    E.g1.mulByScalar(p1, vk_delta1, (uint8_t *)&s, sizeof(s));
    E.g1.add(pib1, pib1, p1);

    E.g1.add(pi_c, pi_c, pih);

    E.g1.mulByScalar(p1, pi_a, (uint8_t *)&s, sizeof(s));
    E.g1.add(pi_c, pi_c, p1);

    E.g1.mulByScalar(p1, pib1, (uint8_t *)&r, sizeof(r));
    E.g1.add(pi_c, pi_c, p1);

    E.fr.mul(rs, r, s);
    E.fr.toMontgomery(rs, rs);

    E.g1.mulByScalar(p1, vk_delta1, (uint8_t *)&rs, sizeof(rs));
    E.g1.add(pi_c, pi_c, p1);

    Proof<Engine> *p = new Proof<Engine>(Engine::engine);
    E.g1.copy(p->A, pi_a);
    E.g2.copy(p->B, pi_b);
    E.g1.copy(p->C, pi_c);

    return std::unique_ptr<Proof<Engine>>(p);
}

template <typename Engine>
std::string Proof<Engine>::toJsonStr() {

    std::ostringstream ss;
    ss << "{ \"pi_a\":[\"" << E.f1.toString(A.x) << "\",\"" << E.f1.toString(A.y) << "\",\"1\"], ";
    ss << " \"pi_b\": [[\"" << E.f1.toString(B.x.a) << "\",\"" << E.f1.toString(B.x.b) << "\"],[\"" << E.f1.toString(B.y.a) << "\",\"" << E.f1.toString(B.y.b) << "\"], [\"1\",\"0\"]], ";
    ss << " \"pi_c\": [\"" << E.f1.toString(C.x) << "\",\"" << E.f1.toString(C.y) << "\",\"1\"], ";
    ss << " \"protocol\":\"groth16\" }";
        
    return ss.str();
}

template <typename Engine>
json Proof<Engine>::toJson() {

    json p;

    p["pi_a"] = {};
    p["pi_a"].push_back(E.f1.toString(A.x) );
    p["pi_a"].push_back(E.f1.toString(A.y) );
    p["pi_a"].push_back("1" );


    json x2;
    x2.push_back(E.f1.toString(B.x.a));
    x2.push_back(E.f1.toString(B.x.b));
    json y2;
    y2.push_back(E.f1.toString(B.y.a));
    y2.push_back(E.f1.toString(B.y.b));
    json z2;
    z2.push_back("1");
    z2.push_back("0");
    p["pi_b"] = {};
    p["pi_b"].push_back(x2);
    p["pi_b"].push_back(y2);
    p["pi_b"].push_back(z2);

    p["pi_c"] = {};
    p["pi_c"].push_back(E.f1.toString(C.x) );
    p["pi_c"].push_back(E.f1.toString(C.y) );
    p["pi_c"].push_back("1" );

    p["protocol"] = "groth16";
            
    return p;
}

} // namespace