//#include "fflonk_prover.hpp"

#include "logger.hpp"
#include "curve_utils.hpp"
#include "zkey.hpp"
#include "zkey_fflonk.hpp"
#include "wtns_utils.hpp"
#include <sodium.h>
#include "mul_z.hpp"
#include "keccak_256_transcript.hpp"

using namespace CPlusPlusLogging;

namespace Fflonk {

    template<typename Engine>
    FflonkProver<Engine>::FflonkProver() {
        E = &Engine::engine;
        curveName = CurveUtils::getCurveNameByEngine(this->E);
    }

    template<typename Engine>
    FflonkProver<Engine>::~FflonkProver() {
        delete fft;
    }

    template<typename Engine>
    std::tuple<json, json> FflonkProver<Engine>::prove(BinFileUtils::BinFile *fdZkey, BinFileUtils::BinFile *fdWtns) {
        try {
            LOG_TRACE("FFLONK PROVER STARTED");

            LOG_TRACE("> Reading witness file");

            this->fdZkey = fdZkey;
            this->fdWtns = fdWtns;

            auto wtns = WtnsUtils::loadHeader(fdWtns);

            LOG_TRACE("> Reading zkey file");
            zkey = Zkey::FflonkZkeyHeader::loadFflonkZkeyHeader(fdZkey);

            if (zkey->protocolId != Zkey::FFLONK_PROTOCOL_ID) {
                throw std::invalid_argument("zkey file is not fflonk");
            }

            fft = new FFT<typename Engine::Fr>(zkey->domainSize * 9 + 18);
            zkeyPower = fft->log2(zkey->domainSize);

            mulZ = new MulZ<Engine>(fft);

            LOG_TRACE("> Computing omegaBuffer");
            omegaBuffer = new FrElement[zkey->domainSize * 16];
            FrElement omega = fft->root(zkeyPower + 4, 1);
            omegaBuffer[0] = E.fr.one;
            for (int64_t i = 1; i < zkey->domainSize * 16; ++i) {
                omegaBuffer[i] = E.fr.mul(omegaBuffer[i - 1], omega);
            }

            //TODO compare zkey field with wtns field
//            if (!Scalar.eq(zkey.r, wtns.q)) {
//                throw std::invalid_argument("Curve of the witness does not match the curve of the proving key");
//            }

            //TODO compare zkey field with current field
//            if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
//                throw std::invalid_argument("different wtns curve");
//            }

            //TODO
//            if (wtns.nWitness != zkey.nVars - zkey->nAdditions) {
//                std::ostringstream ss;
//                ss << "Invalid witness length. Circuit: " << zkey->nVars << ", witness: " << wtns.nWitness << ", " << zkey->nAdditions;
//                throw std::invalid_argument(ss);
//            }

            sDomain = zkey->domainSize * sizeof(FrElement);

            std::ostringstream ss;
            ss << "----------------------------\n"
               << "  FFLONK PROVE SETTINGS\n"
               << "  Curve:         " << curveName << "\n"
               << "  Circuit power: " << zkey->power << "\n"
               << "  Domain size:   " << zkey->domainSize << "\n"
               << "  Vars:          " << zkey->nVars << "\n"
               << "  Public vars:   " << zkey->nPublic << "\n"
               << "  Constraints:   " << zkey->nConstraints << "\n"
               << "  Additions:     " << zkey->nAdditions << "\n"
               << "----------------------------\n";
            LOG_TRACE(ss);

            //Read witness data
            LOG_TRACE("> Reading witness file data");
            FrElement *buffWitness = (FrElement *) fdWtns->getSectionData(2);

            // First element in plonk is not used and can be any value. (But always the same).
            // We set it to zero to go faster in the exponentiations.
            buffWitness[0] = E.fr.zero();

            buffInternalWitness = new FrElement[zkey->nAdditions];

            buffers = new std::map<string, FrElement[]>;
            polynomials = new std::map<std::string, Polynomial<Engine>>;
            evaluations = new std::map<std::string, Evaluations<Engine>>;
            toInverse = new std::map<std::string, FrElement>;
            challenges = new std::map<std::string, FrElement>;
            roots = new std::map<std::string, FrElement[]>;

            // To divide prime fields the Extended Euclidean Algorithm for computing modular inverses is needed.
            // NOTE: This is the equivalent of compute 1/denominator and then multiply it by the numerator.
            // The Extended Euclidean Algorithm is expensive in terms of computation.
            // For the special case where we need to do many modular inverses, there's a simple mathematical trick
            // that allows us to compute many inverses, called Montgomery batch inversion.
            // More info: https://vitalik.ca/general/2018/07/21/starks_part_3.html
            // Montgomery batch inversion reduces the n inverse computations to a single one
            // To save this (single) inverse computation on-chain, will compute it in proving time and send it to the verifier.
            // The verifier will have to check:
            // 1) the denominator is correct multiplying by himself non-inverted -> a * 1/a == 1
            // 2) compute the rest of the denominators using the Montgomery batch inversion
            // The inversions are:
            //   · denominator needed in step 8 and 9 of the verifier to multiply by 1/Z_H(xi)
            //   · denominator needed in step 10 and 11 of the verifier
            //   · denominator needed in the verifier when computing L_i^{S1}(X) and L_i^{S2}(X)
            //   · L_i i=1 to num public inputs, needed in step 6 and 7 of the verifier to compute L_1(xi) and PI(xi)
            //toInverse property is the variable to store the values to be inverted

            proof = new SnarkProof<Engine>("fflonk");

            ss << "> Reading Section " << Zkey::ZKEY_FF_ADDITIONS_SECTION << ". Additions\n";
            LOG_TRACE(ss);

            calculateAdditions(fdZkey);

            ss << "> Reading Section " << Zkey::ZKEY_FF_SIGMA1_SECTION << "," << Zkey::ZKEY_FF_SIGMA2_SECTION
               << "," << Zkey::ZKEY_FF_SIGMA3_SECTION << ". Sigma1, Sigma2 & Sigma 3\n";
            LOG_TRACE(ss);

            LOG_TRACE("··· Reading Sigma polynomials ");
            polynomials["Sigma1"] = new Polynomial<Engine>(zkey->domainSize);
            polynomials["Sigma2"] = new Polynomial<Engine>(zkey->domainSize);
            polynomials["Sigma3"] = new Polynomial<Engine>(zkey->domainSize);

            memcpy(polynomials["Sigma1"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION), sDomain);
            memcpy(polynomials["Sigma2"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION), sDomain);
            memcpy(polynomials["Sigma3"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION), sDomain);

            LOG_TRACE("··· Reading Sigma evaluations ");
            evaluations["Sigma1"] = new Evaluations<Engine>(zkey->domainSize);
            evaluations["Sigma2"] = new Evaluations<Engine>(zkey->domainSize);
            evaluations["Sigma3"] = new Evaluations<Engine>(zkey->domainSize);

            memcpy(evaluations["Sigma1"].eval, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION) + sDomain, sDomain * 4);
            memcpy(evaluations["Sigma2"].eval, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION) + sDomain, sDomain * 4);
            memcpy(evaluations["Sigma3"].eval, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION) + sDomain, sDomain * 4);

            ss << "> Reading Section " << Zkey::ZKEY_FF_PTAU_SECTION << ". Powers of Tau\n";
            LOG_TRACE(ss);
            PTau = new G1Point[zkey->domainSize * 16];
            memset(PTau, 0, sizeof(PTau));

            // domainSize * 9 + 18 = SRS length in the zkey saved in setup process.
            // it corresponds to the maximum SRS length needed, specifically to commit C2
            // notice that the reserved buffers size is zkey->domainSize * 16 * sG1 because a power of two buffer size is needed
            // the remaining buffer not filled from SRS are set to 0
            memcpy(PTau, (G1Point *) fdZkey->getSectionData(Zkey::ZKEY_FF_PTAU_SECTION),
                   (zkey->domainSize * 9 + 18) * sizeof(G1Point));

            // START FFLONK PROVER PROTOCOL

            // ROUND 1. Compute C1(X) polynomial
            LOG_TRACE("> ROUND 1");
            round1();

            delete polynomials["T0"];
            delete evaluations["QL"];
            delete evaluations["QR"];
            delete evaluations["QM"];
            delete evaluations["QO"];
            delete evaluations["QC"];

            // ROUND 2. Compute C2(X) polynomial
            LOG_TRACE("> ROUND 2");
            round2();

            delete buffers["A"];
            delete buffers["B"];
            delete buffers["C"];
            delete evaluations["A"];
            delete evaluations["B"];
            delete evaluations["C"];
            delete evaluations["Sigma1"];
            delete evaluations["Sigma2"];
            delete evaluations["Sigma3"];
            delete evaluations["lagrange1"];
            delete evaluations["Z"];

            // ROUND 3. Compute opening evaluations
            LOG_TRACE("> ROUND 3");
            round3();

            delete fdZkey;
            delete polynomials["A"];
            delete polynomials["B"];
            delete polynomials["C"];
            delete polynomials["Z"];
            delete polynomials["T1"];
            delete polynomials["T2"];
            delete polynomials["Sigma1"];
            delete polynomials["Sigma2"];
            delete polynomials["Sigma3"];
            delete polynomials["QL"];
            delete polynomials["QR"];
            delete polynomials["QM"];
            delete polynomials["QC"];
            delete polynomials["QO"];

            // ROUND 4. Compute W(X) polynomial
            LOG_TRACE("> ROUND 4");
            round4();

            // ROUND 5. Compute W'(X) polynomial
            LOG_TRACE("> ROUND 5");
            round5();

            delete polynomials["C1"];
            delete polynomials["C2"];
            delete polynomials["R1"];
            delete polynomials["R2"];
            delete polynomials["F"];
            delete polynomials["L"];
            delete polynomials["ZT"];
            delete polynomials["ZTS2"];

            proof.addEvaluation("inv", getMontgomeryBatchedInverse());

            // Prepare public inputs
            json publicSignals;
            for (u_int32_t i = 1; i <= zkey->nPublic; i++) {
                publicSignals.push_back(E.Fr.toString(E.Fr.toMontgomery(buffWitness[i])));
            }

            delete fdWtns;

            LOG_TRACE("FFLONK PROVE FINISHED");

            return { proof->toJson(), publicSignals};
        }
        catch (const std::exception &e) {
            std::cerr << "EXCEPTION: " << e.what() << "\n";
            exit(EXIT_FAILURE);
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::calculateAdditions() {
        Zkey::Addition<Engine> *additionsBuff = fdZkey->getSectionData(Zkey::ZKEY_FF_ADDITIONS_SECTION);

        for (u_int32_t i = 0; i < zkey->nAdditions; i++) {
            // Get witness value
            FrElement witness1 = getWitness(additionsBuff[i]->signalId1);
            FrElement witness2 = getWitness(additionsBuff[i]->signalId2);

            //Calculate final result
            buffInternalWitness[i] = E.fr.add(E.fr.mul(additionsBuff[i]->factor1, witness1),
                                              E.fr.mul(additionsBuff[i]->factor2, witness2));
        }
    }

    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getWitness(u_int64_t idx) {
        u_int32_t diff = zkey->nVars - zkey->nAdditions;

        if (idx < diff) {
            return buffWitness[index];
        } else if (idx < zkey->nVars) {
            return buffInternalWitness[idx - diff];
        }

        return E.fr.zero();
    }

    // ROUND 1
    template<typename Engine>
    void FflonkProver<Engine>::round1() {
        // STEP 1.1 - Generate random blinding scalars (b_1, ..., b9) ∈ F
        blindingFactors = new FrElement[10];

        //0 index not used, set to zero
        //TODO check it fill all bytes with random values!!!!
        randombytes_buf((void *) &blindingFactors[0], sizeof(blindingFactors));

        // STEP 1.2 - Compute wire polynomials a(X), b(X) and c(X)
        computeWirePolynomials();

        // STEP 1.3 - Compute the quotient polynomial T0(X)
        computeT0();

        // STEP 1.4 - Compute the FFT-style combination polynomial C1(X)
        computeC1();

        // The first output of the prover is ([C1]_1)
        proof.addPolynomial("C1", multiExponentiation(polynomials["C1"], "C1"));
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeWirePolynomials() {
        // Build A, B and C evaluations buffer from zkey and witness files
        buffers["A"] = new FrElement[zkey->domainSize];
        buffers["B"] = new FrElement[zkey->domainSize];
        buffers["C"] = new FrElement[zkey->domainSize];

        memset(buffers["A"], 0, sizeof(buffers["A"]));
        memset(buffers["B"], 0, sizeof(buffers["B"]));
        memset(buffers["C"], 0, sizeof(buffers["C"]));

        // Read zkey sections and fill the buffers
        memcpy(buffers["A"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_A_MAP_SECTION), zkey->domainSize);
        memcpy(buffers["B"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_B_MAP_SECTION), zkey->domainSize);
        memcpy(buffers["C"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_C_MAP_SECTION), zkey->domainSize);

        #pragma omp parallel
        {
            #pragma omp sections
            {
                #pragma omp section
                {
                    computeWirePolynomial("A", {blindingFactors[2], blindingFactors[1]});
                }

                #pragma omp section
                {
                    computeWirePolynomial("B", {blindingFactors[4], blindingFactors[3]});
                }

                #pragma omp section
                {
                    computeWirePolynomial("C", {blindingFactors[6], blindingFactors[5]});
                }
            }
        }

        // Check degrees
        if (polynomials["A"].degree() >= zkey->domainSize + 2) {
            throw std::runtime_error("A Polynomial is not well calculated");
        }
        if (polynomials["B"].degree() >= zkey->domainSize + 2) {
            throw std::runtime_error("B Polynomial is not well calculated");
        }
        if (polynomials["C"].degree() >= zkey->domainSize + 2) {
            throw std::runtime_error("C Polynomial is not well calculated");
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeWirePolynomial(std::string polName, FrElement blindingFactors[]) {

        // Compute all witness from signal ids and set them to the polynomial buffers
        #pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->nConstraints; ++i) {
            FrElement aux;
            FrElement witness = getWitness(buffers[polName][i]);
            buffers[polName][i] = E.fr.toMontgomery(witness);
        }

        // Create the polynomial
        // and compute the coefficients of the wire polynomials from evaluations
        u_int32_t bFactorsLen = sizeof(*blindingFactors) / sizeof(Engine::FrElement);
        polynomials[polName] = new Polynomial<Engine>(buffers[polName], zkey->domainSize, bFactorsLen);

        // Compute the extended evaluations of the wire polynomials
        evaluations[polName] = new Evaluations<Engine>(polynomials[polName]);

        // Blind polynomial coefficients with blinding scalars blindingFactors
        polynomials[polName]->blindCoefficients(blindingFactors);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeT0() {
        std::ostringstream ss;
        ss << "> Reading sections"
           << Zkey::ZKEY_FF_QL_SECTION << "," << Zkey::ZKEY_FF_QR_SECTION << ","
           << Zkey::ZKEY_FF_QM_SECTION << "," << Zkey::ZKEY_FF_QO_SECTION << "," << Zkey::ZKEY_FF_QC_SECTION
           << ". Q selectors\n";
        LOG_TRACE(ss);

        // Reserve memory for Q's evaluations
        evaluations["QL"] = new Evaluations<Engine>(zkey->domainSize * 4);
        evaluations["QR"] = new Evaluations<Engine>(zkey->domainSize * 4);
        evaluations["QM"] = new Evaluations<Engine>(zkey->domainSize * 4);
        evaluations["QO"] = new Evaluations<Engine>(zkey->domainSize * 4);
        evaluations["QC"] = new Evaluations<Engine>(zkey->domainSize * 4);

        // Read Q's evaluations from zkey file
        memcpy(evaluations["QL"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QR"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QR_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QM"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QM_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QO"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QO_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QC"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QC_SECTION) + sDomain, sDomain * 4);

        // Read Lagrange polynomials & evaluations from zkey file
        evaluations["lagrangePolynomials"] = new Evaluations<Engine>(zkey->domainSize * 5);
        memcpy(evaluations["lagrangePolynomials"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_LAGRANGE_SECTION), sDomain * 5);

        // Reserve memory for buffers T0 and T0z
        buffers["T0"] =  new FrElement[zkey->domainSize * 4];
        buffers["T0z"] =  new FrElement[zkey->domainSize * 4];

        LOG_TRACE("> Computing T0");
        #pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->domainSize * 4; i++) {
            if ((0 != i) && (i % 5000 == 0)) {
                ss.str("");
                ss << "> Computing t0 evaluation" << i / (zkey->domainSize * 4);
                LOG_TRACE(ss);
            }

            FrElement omega =  omegaBuffer[i * 4];

            // Get related evaluations to compute current T0 evaluation
            FrElement a = evaluations["A"]->getEvaluation(i);
            FrElement b = evaluations["B"]->getEvaluation(i);
            FrElement c = evaluations["C"]->getEvaluation(i);

            FrElement ql = evaluations["QL"]->getEvaluation(i);
            FrElement qr = evaluations["QR"]->getEvaluation(i);
            FrElement qm = evaluations["QM"]->getEvaluation(i);
            FrElement qo = evaluations["QO"]->getEvaluation(i);
            FrElement qc = evaluations["QC"]->getEvaluation(i);

            // Compute blinding factors
            FrElement az = E.fr.add(E.fr.mul(blindingFactors[1], omega), blindingFactors[2]);
            FrElement bz = E.fr.add(E.fr.mul(blindingFactors[3], omega), blindingFactors[4]);
            FrElement cz = E.fr.add(E.fr.mul(blindingFactors[5], omega), blindingFactors[6]);

            // Compute current public input
            FrElement pi = E.Fr.zero;
            for (u_int32_t j = 0; j < zkey->nPublic; j++) {
                u_int32_t offset = (j * 5 * zkey->domainSize) + zkey->domainSize + i;

                FrElement lPol = evaluations["lagrange1"]->getEvaluation(offset);
                FrElement aVal = buffers["A"][j];

                pi = E.Fr.sub(pi, E.Fr.mul(lPol, aVal));
            }

            // T0(X) = [q_L(X)·a(X) + q_R(X)·b(X) + q_M(X)·a(X)·b(X) + q_O(X)·c(X) + q_C(X) + PI(X)] · 1/Z_H(X)
            // Compute first T0(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            // expression 1 -> q_L(X)·a(X)
            FrElement e1 = E.Fr.mul(a, ql);
            FrElement e1z = E.Fr.mul(az, ql);

            // expression 2 -> q_R(X)·b(X)
            FrElement e2 = E.Fr.mul(b, qr);
            FrElement e2z = E.Fr.mul(bz, qr);

            // expression 3 -> q_M(X)·a(X)·b(X)
            auto [e3, e3z] = mulZ.mul2(a, b, az, bz, i % 4);
            e3 = E.Fr.mul(e3, qm);
            e3z = E.Fr.mul(e3z, qm);

            // expression 4 -> q_O(X)·c(X)
            FrElement e4 = E.Fr.mul(c, qo);
            FrElement e4z = E.Fr.mul(cz, qo);

            // t0 = expressions 1 + expression 2 + expression 3 + expression 4 + qc + pi
            FrElement t0 = E.Fr.add(e1, E.Fr.add(e2, E.Fr.add(e3, E.Fr.add(e4, E.Fr.add(qc, pi)))));
            FrElement t0z = E.Fr.add(e1z, E.Fr.add(e2z, E.Fr.add(e3z, e4z)));

            buffers["T0"] = t0;
            buffers["T0z"] = t0z;
        }

        // Compute the coefficients of the polynomial T0(X) from buffers.T0
        LOG_TRACE("··· Computing T0 ifft");
        polynomials["T0"] = new Polynomial<Engine>(buffers["T0"], zkey->domainSize * 4);

        // Divide the polynomial T0 by Z_H(X)
        polynomials["T0"]->divZh();

        // Compute the coefficients of the polynomial T0z(X) from buffers.T0z
        LOG_TRACE("··· Computing T0z ifft");
        polynomials["T0z"] = new Polynomial<Engine>(buffers["T0z"], zkey->domainSize * 4);

        // Add the polynomial T0z to T0 to get the final polynomial T0
        polynomials["T0"]->add(polynomials["T0z"]);

        // Check degree
        if (polynomials["T0"]->degree() >= 2 * zkey->domainSize + 2) {
            throw std::runtime_error("T0 Polynomial is not well calculated");
        }

        delete buffers["T0"];
        delete buffers["T0z"];
        delete polynomials["T0z"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeC1() {
        LOG_TRACE("··· Computing C1");

        // C1(X) := a(X^4) + X · b(X^4) + X^2 · c(X^4) + X^3 · T0(X^4)
        // Get X^n · f(X) by shifting the f(x) coefficients n positions,
        // the resulting polynomial will be degree deg(f(X)) + n
        u_int64_t lengthA = polynomials["A"].length;
        u_int64_t lengthB = polynomials["B"].length;
        u_int64_t lengthC = polynomials["C"].length;
        u_int64_t lengthT0 = polynomials["T0"].length;

        // Compute degree of the new polynomial C1 to reserve the buffer memory size
        // Will be the next power of two to bound the maximum(deg(A_4), deg(B_4)+1, deg(C_4)+2, deg(T0_4)+3)
        u_int64_t degreeA = polynomials["A"]->degree;
        u_int64_t degreeB = polynomials["B"]->degree;
        u_int64_t degreeC = polynomials["C"]->degree;
        u_int64_t degreeT0 = polynomials["T0"]->degree;

        u_int64_t maxLength = std::max(lengthA, std::max(lengthB, std::max(lengthC, lengthT0)));
        u_int64_t maxDegree = std::max(degreeA * 4 + 1, std::max(degreeB * 4 + 2, std::max(degreeC * 4 + 3, degreeT0 * 4 + 3)));

        u_int64_t lengthBuffer = 2 ** (fft->log2(maxDegree - 1) + 1);

        polynomials["C1"] = new Polynomial<Engine>(lengthBuffer);

        for (u_int64_t i = 0; i < maxLength; i++) {
            polynomials["C1"].coef[i * 4] = polynomials["A"]->getCoef(i);
            polynomials["C1"].coef[i * 4 + 1] = polynomials["B"]->getCoef(i);
            polynomials["C1"].coef[i * 4 + 2] = polynomials["C"]->getCoef(i);
            polynomials["C1"].coef[i * 4 + 3] = polynomials["T0"]->getCoef(i);
        }

        // Check degree
        if (polynomials["C1"]->degree >= 8 * zkey->domainSize + 8) {
            throw std::runtime_error("C1 Polynomial is not well calculated");
        }
    }


    // ROUND 2
    template<typename Engine>
    void FflonkProver<Engine>::round2() {
        // STEP 2.1 - Compute permutation challenge beta and gamma ∈ F
        // Compute permutation challenge beta
        Keccak256Transcript transcript = new Keccak256Transcript<Engine>();
        for (u_int32_t i = 0; i < zkey->nPublic; i++) {
            transcript->addScalar(buffers["A"][i]);
        }
        transcript->addPolCommitment(proof.getPolynomial("C1"));

        challenges["beta"] = transcript->getChallenge();
        std::ostringstream ss;
        ss << "challenges.beta: " << E.Fr.toString(challenges["beta"]);
        LOG_TRACE(ss);

        // Compute permutation challenge gamma
        transcript->reset();
        transcript->addScalar(challenges["beta"]);
        challenges["gamma"] = transcript->getChallenge();
        ss.str("");
        ss << "challenges.gamma: " << E.Fr.toString(challenges["gamma"]);

        // STEP 2.2 - Compute permutation polynomial z(X)
        computeZ();

        // STEP 2.3 - Compute quotient polynomial T1(X) and T2(X)
        computeT1();
        computeT2();

        // STEP 2.4 - Compute the FFT-style combination polynomial C2(X)
        computeC2();

        // The second output of the prover is ([C2]_1)
        proof.addPolynomial("C2", multiExponentiation(polynomials["C2"], "C2"));
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZ() {
        FrElement numArr = new FrElement[zkey->domainSize];
        FrElement denArr = new FrElement[zkey->domainSize];

        // Set the first values to 1
        numArr[0] = E.Fr.one;
        denArr[0] = E.Fr.one;

        LOG_TRACE("> Computing Z");

        std::ostringstream ss;
        for (u_int64_t i = 0; i < zkey->domainSize; i++) {
            if ((0 != i) && (i % 5000 == 0)) {
                ss.str("");
                ss << "> Computing Z evaluation" << i / zkey->domainSize;
                LOG_TRACE(ss);
            }

            FrElement omega =  omegaBuffer[i * 16];

            // Z(X) := numArr / denArr
            // numArr := (a + beta·ω + gamma)(b + beta·ω·k1 + gamma)(c + beta·ω·k2 + gamma)
            FrElement betaw = E.Fr.mul(challenges["beta"], omega);

            FrElement num1 = buffers["A"][i];
            num1 = E.Fr.add(num1, betaw);
            num1 = E.Fr.add(num1, challenges["gamma"]);

            FrElement num2 = buffers["B"][i];
            num2 = E.Fr.add(num2, E.Fr.mul(zkey->k1, betaw));
            num2 = E.Fr.add(num2, challenges["gamma"]);

            FrElement num3 = buffers["C"][i];
            num3 = E.Fr.add(num3, E.Fr.mul(zkey->k2, betaw));
            num3 = E.Fr.add(num3, challenges["gamma"]);

            FrElement num = E.Fr.mul(num1, E.Fr.mul(num2, num3));

            // denArr := (a + beta·sigma1 + gamma)(b + beta·sigma2 + gamma)(c + beta·sigma3 + gamma)
            FrElement den1 = buffers["A"][i];
            den1 = E.Fr.add(den1, E.Fr.mul(challenges["beta"], evaluations["Sigma1"]->getEvaluation(i * 4)));
            den1 = E.Fr.add(den1, challenges["gamma"]);

            FrElement den2 = buffers["B"][i];
            den2 = E.Fr.add(den2, E.Fr.mul(challenges["beta"], evaluations["Sigma2"]->getEvaluation(i * 4)));
            den2 = E.Fr.add(den2, challenges["gamma"]);

            FrElement den3 = buffers["C"][i];
            den3 = E.Fr.add(den3, E.Fr.mul(challenges["beta"], evaluations["Sigma3"]->getEvaluation(i * 4)));
            den3 = E.Fr.add(den3, challenges["gamma"]);

            FrElement den = E.Fr.mul(den1, E.Fr.mul(den2, den3));

            // Multiply current num value with the previous one saved in numArr
            num = E.Fr.mul(numArr[i], num);
            numArr[(i + 1) % zkey->domainSize] = num;

            // Multiply current den value with the previous one saved in denArr
            den = E.Fr.mul(denArr[i], den);
            denArr[(i + 1) % zkey->domainSize] = den;
        }
        // Compute the inverse of denArr to compute in the next command the
        // division numArr/denArr by multiplying num · 1/denArr
        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //denArr = await Fr.batchInverse(denArr);

        // Multiply numArr · denArr where denArr was inversed in the previous command
        #pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->domainSize; i++) {
            buffers["Z"][i] = E.Fr.mul(numArr[i], denArr[i]);
        }

        if (E.Fr.neq(numArr[0], E.Fr.one)) {
            throw std::runtime_error("Copy constraints does not match");
        }

        // Compute polynomial coefficients z(X) from buffers.Z
        LOG_TRACE("··· Computing Z ifft");
        polynomials["Z"] = new Polynomial<Engine>(buffers["Z"], zkey->domainSize);

        // Compute extended evaluations of z(X) polynomial
        evaluations["Z"] = new Evaluations<Engine>(polynomials["Z"]);

        // Blind z(X) polynomial coefficients with blinding scalars b
        polynomials["Z"]->blindCoefficients({blindingFactors[9], blindingFactors[8], blindingFactors[7]});

        // Check degree
        if (polynomials["Z"]->degree >= zkey->domainSize + 3) {
            throw std::runtime_error("Z Polynomial is not well calculated");
        }

        delete buffers["Z"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeT1() {
        buffers["T1"] = new FrElement[zkey->domainSize * 4];
        buffers["T1z"] = new FrElement[zkey->domainSize * 4];

        LOG_TRACE("> Computing T1");

        std::ostringstream ss;
        for (u_int64_t i = 0; i < zkey->domainSize * 4; i++) {
            if ((0 != i) && (i % 5000 == 0)) {
                ss.str("");
                ss << "> Computing t1 evaluation" << i / (zkey->domainSize * 4);
                LOG_TRACE(ss);
            }

            FrElement omega =  omegaBuffer[i];
            FrElement omega2 =  E.Fr.square(omega);

            FrElement z = evaluations["Z"]->getEvaluation(i);
            FrElement zp = E.Fr.add(E.Fr.add(
                    E.Fr.mul(blindingFactors[7], omega2), E.Fr.mul(blindingFactors[8], omega)), blindingFactors[9]);

            // T1(X) := (z(X) - 1) · L_1(X)
            // Compute first T1(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            FrElement lagrange1 = evaluations["lagrange1"]->getEvaluation(zkey->domainSize + i);
            FrElement t1 = E.Fr.mul(E.Fr.sub(z, E.Fr.one), lagrange1);
            FrElement t1z = E.Fr.mul(zp, lagrange1);

            buffers["T1"][i] = t1;
            buffers["T1z"][i] = t1z;
        }

        // Compute the coefficients of the polynomial T1(X) from buffers.T1
        LOG_TRACE("··· Computing T1 ifft");
        polynomials["T1"] = new Polynomial<Engine>(buffers["T1"], zkey->domainSize * 4);

        // Divide the polynomial T1 by Z_H(X)
        polynomials["T1"]->divZh();

        // Compute the coefficients of the polynomial T1z(X) from buffers.T1z
        LOG_TRACE("··· Computing T1z ifft");
        polynomials["T1z"] = new Polynomial<Engine>(buffers["T1z"], zkey->domainSize * 4);

        // Add the polynomial T0z to T0 to get the final polynomial T0
        polynomials["T1"]->add(polynomials["T1z"]);

        // Check degree
        if (polynomials["T1"]->degree >= zkey->domainSize + 2) {
            throw std::runtime_error("T1 Polynomial is not well calculated");
        }

        delete buffers.T1;
        delete buffers.T1z;
        delete polynomials.T1z;
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeT2() {
        buffers["T2"] = new FrElement[zkey->domainSize * 4];
        buffers["T2z"] = new FrElement[zkey->domainSize * 4];

        LOG_TRACE("> Computing T2");

        std::ostringstream ss;
        for (u_int64_t i = 0; i < zkey->domainSize * 4; i++) {
            if ((0 != i) && (i % 5000 == 0)) {
                ss.str("");
                ss << "> Computing t2 evaluation" << i / (zkey->domainSize * 4);
                LOG_TRACE(ss);
            }

            FrElement omega =  omegaBuffer[i];
            FrElement omega2 =  E.Fr.square(omega);
            FrElement omegaW = E.Fr.mul(omega, fft->root(zkeyPower, 1));
            FrElement omegaW2 = E.Fr.square(omegaW);

            FrElement a = evaluations["A"]->getEvaluation(i);
            FrElement b = evaluations["B"]->getEvaluation(i);
            FrElement c = evaluations["C"]->getEvaluation(i);
            FrElement z = evaluations["Z"]->getEvaluation(i);
            FrElement zW = evaluations["Z"]->getEvaluation((zkey->domainSize * 4 + 4 + i) % (zkey->domainSize * 4));

            FrElement ap = E.Fr.add(E.Fr.mul(blindingFactors[1], omega), blindingFactors[2]);
            FrElement bp = E.Fr.add(E.Fr.mul(blindingFactors[3], omega), blindingFactors[4]);
            FrElement cp = E.Fr.add(E.Fr.mul(blindingFactors[5], omega), blindingFactors[6]);
            FrElement zp = E.Fr.add(E.Fr.add(E.Fr.mul(blindingFactors[7], omega2), E.Fr.mul(blindingFactors[8], omega)), blindingFactors[9]);
            FrElement zWp = E.Fr.add(E.Fr.add(E.Fr.mul(blindingFactors[7], omegaW2), E.Fr.mul(blindingFactors[8], omegaW)), blindingFactors[9]);

            FrElement sigma1 = evaluations["Sigma1"]->getEvaluation(i);
            FrElement sigma2 = evaluations["Sigma2"]->getEvaluation(i);
            FrElement sigma3 = evaluations["Sigma3"]->getEvaluation(i);

            // T2(X) := [ (a(X) + beta·X + gamma)(b(X) + beta·k1·X + gamma)(c(X) + beta·k2·X + gamma)z(X)
            //           -(a(X) + beta·sigma1(X) + gamma)(b(X) + beta·sigma2(X) + gamma)(c(X) + beta·sigma3(X) + gamma)z(Xω)] · 1/Z_H(X)
            // Compute first T2(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)

            // expression 1 -> (a(X) + beta·X + gamma)(b(X) + beta·k1·X + gamma)(c(X) + beta·k2·X + gamma)z(X)
            FrElement betaX = E.Fr.mul(challenges["beta"], omega);

            FrElement e11 = E.Fr.add(a, betaX);
            e11 = E.Fr.add(e11, challenges["gamma"]);

            FrElement e12 = E.Fr.add(b, E.Fr.mul(betaX, zkey->k1));
            e12 = E.Fr.add(e12, challenges["gamma"]);

            FrElement e13 = E.Fr.add(c, E.Fr.mul(betaX, zkey->k2));
            e13 = E.Fr.add(e13, challenges["gamma"]);

            auto [e1, e1z] = mulZ.mul4(e11, e12, e13, z, ap, bp, cp, zp, i % 4);

            // expression 2 -> (a(X) + beta·sigma1(X) + gamma)(b(X) + beta·sigma2(X) + gamma)(c(X) + beta·sigma3(X) + gamma)z(Xω)
            FrElement e21 = E.Fr.add(a, E.Fr.mul(challenges["beta"], sigma1));
            e21 = E.Fr.add(e21, challenges["gamma"]);

            FrElement e22 = E.Fr.add(b, E.Fr.mul(challenges["beta"], sigma2));
            e22 = E.Fr.add(e22, challenges["gamma"]);

            FrElement e23 = E.Fr.add(c, E.Fr.mul(challenges["beta"], sigma3));
            e23 = E.Fr.add(e23, challenges["gamma"]);

            auto [e2, e2z] = mulZ.mul4(e21, e22, e23, zW, ap, bp, cp, zWp, i % 4);

            FrElement t2 = E.Fr.sub(e1, e2);
            FrElement t2z = E.Fr.sub(e1z, e2z);

            buffers["T2"][i] = t2;
            buffers["T2z"][i] = t2z;
        }

        // Compute the coefficients of the polynomial T2(X) from buffers.T2
        LOG_TRACE("··· Computing T2 ifft");
        polynomials["T2"] = new Polynomial<Engine>(buffers["T2"], zkey->domainSize * 4);

        // Divide the polynomial T2 by Z_H(X)
        polynomials["T2"]->divZh();

        // Compute the coefficients of the polynomial T2z(X) from buffers.T2z
        LOG_TRACE("··· Computing T2z ifft");
        polynomials["T2z"] = new Polynomial<Engine>(buffers["T2z"], zkey->domainSize * 4);

        // Add the polynomial T2z to T2 to get the final polynomial T2
        polynomials["T2"]->add(polynomials["T2z"]);

        // Check degree
        if (polynomials["T2"]->degree >= 3 * zkey->domainSize + 6) {
            throw std::runtime_error("T2 Polynomial is not well calculated");
        }

        delete buffers["T2"];
        delete buffers["T2z"];
        delete polynomials["T2z"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeC2() {
        LOG_TRACE("··· Computing C2");

        // C2(X) := z(X^3) + X · T1(X^3) + X^2 · T2(X^3)
        // Get X^n · f(X) by shifting the f(x) coefficients n positions,
        // the resulting polynomial will be degree deg(f(X)) + n
        u_int64_t lengthZ = polynomials["Z"]->length;
        u_int64_t lengthT1 = polynomials["T1"]->length;
        u_int64_t lengthT2 = polynomials["T2"]->length;
        // Compute degree of the new polynomial C2(X) to reserve the buffer memory size
        // Will be the maximum(deg(Z_3), deg(T1_3)+1, deg(T2_3)+2)
        u_int64_t degreeZ = polynomials["Z"]->degree;
        u_int64_t degreeT1 = polynomials["T1"]->degree;
        u_int64_t degreeT2 = polynomials["T2"]->degree;

        u_int64_t maxLength = std::max(lengthZ, std::max(lengthT1, lengthT2));
        u_int64_t maxDegree = std::max(degreeZ * 3 + 1, std::max(degreeT1 * 3 + 2, degreeT2 * 3 + 3));

        u_int64_t lengthBuffer = 2 ** (fft->log2(maxDegree - 1) + 1);

        polynomials["C2"] = new Polynomial<Engine>(lengthBuffer);

        for (u_int64_t i = 0; i < maxLength; i++) {
            polynomials["C2"].coef[i * 4] = polynomials["Z"]->getCoef(i);
            polynomials["C2"].coef[i * 4 + 1] = polynomials["T1"]->getCoef(i);
            polynomials["C2"].coef[i * 4 + 2] = polynomials["T2"]->getCoef(i);
        }

        // Check degree
        if (polynomials["C2"]->degree >= 9 * zkey->domainSize + 18) {
            throw std::runtime_error("C2 Polynomial is not well calculated");
        }
    }

    // ROUND 3
    template<typename Engine>
    void FflonkProver<Engine>::round3() {
        // STEP 3.1 - Compute evaluation challenge xi ∈ S
        Keccak256Transcript transcript = new Keccak256Transcript<Engine>();
        transcript->addPolCommitment(proof.getPolynomial("C2"));

        // Obtain a xi_seeder from the transcript
        // To force h1^4 = xi, h2^3 = xi and h_3^2 = xiω
        // we compute xi = xi_seeder^12, h1 = xi_seeder^3, h2 = xi_seeder^4 and h3 = xi_seeder^6
        FrElement xiSeed = transcript->getChallenge();
        FrElement xiSeed2 = E.Fr.square(xiSeed);

        // Compute omega3 and omega4
        roots["w4"] = new FrElement[4];
        roots["w4"][0] = E.Fr.one;
        roots["w4"][1] = zkey->w4;
        roots["w4"][2] = E.Fr.square(zkey->w4);
        roots["w4"][3] = E.Fr.mul(roots["w4"][2], zkey->w4);

        roots["w3"] = new FrElement[3];
        roots["w3"][0] = E.Fr.one;
        roots["w3"][1] = zkey->w3;
        roots["w3"][2] = E.Fr.square(zkey->w3);

        // Compute h1 = xi_seeder^3
        roots["S1h1"] = new FrElement[4];
        roots["S1h1"][0] = E.Fr.mul(xiSeed2, xiSeed);
        roots["S1h1"][1] = E.Fr.mul(roots["S1h1"][0], roots["w4"][1]);
        roots["S1h1"][2] = E.Fr.mul(roots["S1h1"][0], roots["w4"][2]);
        roots["S1h1"][3] = E.Fr.mul(roots["S1h1"][0], roots["w4"][3]);

        roots["S2h2"] = new FrElement[3];
        roots["S2h2"][0] = E.Fr.square(xiSeed2);
        roots["S2h2"][1] = E.Fr.mul(roots["S2h2"][0], roots["w3"][1]);
        roots["S2h2"][2] = E.Fr.mul(roots["S2h2"][0], roots["w3"][2]);

        roots["S2h3"] = new FrElement[3];
        // Multiply h3 by third-root-omega to obtain h_3^3 = xiω
        roots["S2h3"][0] = E.Fr.mul(roots["S2h2"][0], zkey->wr);
        roots["S2h3"][1] = E.Fr.mul(roots["S2h3"][0], roots["w3"][1]);
        roots["S2h3"][2] = E.Fr.mul(roots["S2h3"][0], roots["w3"][2]);

        // Compute xi = xi_seeder^12
        challenges.xi = E.Fr.mul(E.Fr.square(roots["S2h2"][0]), roots["S2h2"][0]);

        std::ostringstream ss;
        ss << "challenges.xi: " << E.Fr.toString(challenges["xi"]);
        LOG_TRACE(ss);

        // Reserve memory for Q's polynomials
        polynomials["QL"] = new Polynomial<Engine>(zkey->domainSize);
        polynomials["QR"] = new Polynomial<Engine>(zkey->domainSize);
        polynomials["QM"] = new Polynomial<Engine>(zkey->domainSize);
        polynomials["QO"] = new Polynomial<Engine>(zkey->domainSize);
        polynomials["QC"] = new Polynomial<Engine>(zkey->domainSize);

        // Read Q's evaluations from zkey file
        memcpy(polynomials["QL"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QR"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QM"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QO"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QC"].coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);

        // STEP 3.2 - Compute opening evaluations and add them to the proof (third output of the prover)
        proof->addEvaluationCommitment("ql", polynomials["QL"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qr", polynomials["QR"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qm", polynomials["QM"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qo", polynomials["QO"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qc", polynomials["QC"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("s1", polynomials["Sigma1"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("s2", polynomials["Sigma2"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("s3", polynomials["Sigma3"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("a", polynomials["A"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("b", polynomials["B"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("c", polynomials["C"]->evaluate(challenges["xi"]));
        proof->addEvaluationCommitment("z", polynomials["Z"]->evaluate(challenges["xi"]));

        FrElement xiw = E.Fr.mul(challenges["xi"], fft->root(zkeyPower + 4, 1));
        proof->addEvaluationCommitment("zw", polynomials["Z"]->evaluate(challenges["xiw"]));
        proof->addEvaluationCommitment("t1w", polynomials["T1"]->evaluate(challenges["xiw"]));
        proof->addEvaluationCommitment("t2w", polynomials["T2"]->evaluate(challenges["xiw"]));
    }

    // ROUND 4
    template<typename Engine>
    void FflonkProver<Engine>::round4() {
        // STEP 4.1 - Compute challenge alpha ∈ F
        Keccak256Transcript<Engine> transcript = new Keccak256Transcript<Engine>();
        transcript->addScalar(proof->addEvaluationCommitment("ql"));
        transcript->addScalar(proof->addEvaluationCommitment("qr"));
        transcript->addScalar(proof->addEvaluationCommitment("qm"));
        transcript->addScalar(proof->addEvaluationCommitment("qo"));
        transcript->addScalar(proof->addEvaluationCommitment("qc"));
        transcript->addScalar(proof->addEvaluationCommitment("s1"));
        transcript->addScalar(proof->addEvaluationCommitment("s2"));
        transcript->addScalar(proof->addEvaluationCommitment("s3"));
        transcript->addScalar(proof->addEvaluationCommitment("a"));
        transcript->addScalar(proof->addEvaluationCommitment("b"));
        transcript->addScalar(proof->addEvaluationCommitment("c"));
        transcript->addScalar(proof->addEvaluationCommitment("z"));
        transcript->addScalar(proof->addEvaluationCommitment("zw"));
        transcript->addScalar(proof->addEvaluationCommitment("t1w"));
        transcript->addScalar(proof->addEvaluationCommitment("t2w"));
        challenges["alpha"] = transcript->getChallenge();

        std::ostringstream ss;
        ss << "challenges.alpha: " << E.Fr.toString(challenges["alpha"]);
        LOG_TRACE(ss);

        // STEP 4.2 - Compute F(X)
        computeR1();
        computeR2();

        computeF();
        computeZT();

        Polynomial<Engine> polRemainder = polynomials["F"]->divBy(polynomials["ZT"]);

        // Check degrees
        if (polRemainder->degree > 0) {
            std::ostringstream ss;
            ss << "Degree of f(X)/ZT(X) remainder is " << polRemainder->degree << " and should be 0";
            throw std::runtime_error(ss);
        }
        if (polynomials["F"]->degree >= 9 * zkey->domainSize + 12) {
            throw std::runtime_error("Degree of f(X)/ZT(X) is not correct");
        }

        // The fourth output of the prover is ([W1]_1), where W1:=(f/Z_t)(x)
        proof.addPolynomial("W1", multiExponentiation(polynomials["W1"], "W1"));
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeR1() {
        // COMPUTE R1
        // Compute the coefficients of R1(X) from 4 evaluations using lagrange interpolation. R1(X) ∈ F_{<4}[X]
        // We decide to use Lagrange interpolations because the R1 degree is very small (deg(R1)===3),
        // and we were not able to compute it using current ifft implementation because the omega are different
        LOG_TRACE("> Computing R1");

        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        polynomials["R1"] = Polynomial<Engine>::lagrangePolynomialInterpolation(
                {roots["S1h1"][0], roots["S1h1"][1], roots["S1h1"][2], roots["S1h1"][3]},
                {polynomials["C1"]->evaluate(roots["S1h1"][0]), polynomials["C1"]->evaluate(roots["S1h1"][1]),
                polynomials["C1"]->evaluate(roots["S1h1"][2]), polynomials["C1"]->evaluate(roots["S1h1"][3])});

        // Check the degree of r1(X) < 4
        if (polynomials["R1"]->degree > 3) {
            throw std::runtime_error("R1 Polynomial is not well calculated");
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeR2() {
        // COMPUTE R2
        // Compute the coefficients of r2(X) from 6 evaluations using lagrange interpolation. r2(X) ∈ F_{<6}[X]
        // We decide to use Lagrange interpolations because the R2.degree is very small (deg(R2)===5),
        // and we were not able to compute it using current ifft implementation because the omega are different
        LOG_TRACE("> Computing R2");
        polynomials["R2"] = Polynomial<Engine>::lagrangePolynomialInterpolation(
                {roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                 roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]},
                {polynomials.C2.evaluate(roots["S2h2"][0]), polynomials.C2.evaluate(roots["S2h2"][1]),
                 polynomials.C2.evaluate(roots["S2h2"][2]), polynomials.C2.evaluate(roots["S2h3"][0]),
                 polynomials.C2.evaluate(roots["S2h3"][1]), polynomials.C2.evaluate(roots["S2h3"][2])});

        // Check the degree of r2(X) < 6
        if (polynomials["R2"]->degree > 5) {
            throw std::runtime_error("R2 Polynomial is not well calculated");
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeF() {
        buffers["F"] = new FrElement(zkey->domainSize * 16);

        LOG_TRACE("> Computing C1 & C2 fft");

        evaluations["C1"] = new Evaluations<Engine>(polynomials["C1"]);
        evaluations["C2"] = new Evaluations<Engine>(polynomials["C2"]);

        LOG_TRACE("> Computing F");
        // COMPUTE F(X)
        std::ostringstream ss;
        for (u_int64_t i = 0; i < zkey->domainSize * 16; i++) {
            if ((0 != i) && (i % 5000 == 0)) {
                ss.str("");
                ss << "> Computing F evaluation" << i / (zkey->domainSize * 16);
                LOG_TRACE(ss);
            }

            FrElement omega =  omegaBuffer[i * 16];

            FrElement c1 = evaluations["C1"]->getEvaluation(i * 4);
            FrElement c2 = evaluations["C2"]->getEvaluation(i * 4);
            FrElement r1 = polynomials["R1"]->evaluate(omega);
            FrElement r2 = polynomials["R2"]->evaluate(omega);

            // f1 = (X - h2) (X - h2w3) (X - h2w3_2) (X - h3) (X - h3w3) (X - h3w3_2) (C1(X) - R1(X))
            FrElement f1 = E.Fr.sub(omega, roots["S2h2"][0]);
            f1 = E.Fr.mul(f1, E.Fr.sub(omega, roots["S2h2"][1]));
            f1 = E.Fr.mul(f1, E.Fr.sub(omega, roots["S2h2"][2]));
            f1 = E.Fr.mul(f1, E.Fr.sub(omega, roots["S2h3"][0]));
            f1 = E.Fr.mul(f1, E.Fr.sub(omega, roots["S2h3"][1]));
            f1 = E.Fr.mul(f1, E.Fr.sub(omega, roots["S2h3"][2]));
            f1 = E.Fr.mul(f1, E.Fr.sub(c1, r1));

            // f2 = alpha (X - h1) (X - h1w4) (X - h1w4_2) (X - h1w4_3) (C2(X) - R2(X))
            FrElement f2 = E.Fr.mul(challenges["alpha"], E.Fr.sub(omega, roots["S1h1"][0]));
            f2 = E.Fr.mul(f2, E.Fr.sub(omega, roots["S1h1"][1]));
            f2 = E.Fr.mul(f2, E.Fr.sub(omega, roots["S1h1"][2]));
            f2 = E.Fr.mul(f2, E.Fr.sub(omega, roots["S1h1"][3]));
            f2 = E.Fr.mul(f2, E.Fr.sub(c2, r2));

            FrElement f = E.Fr.add(f1, f2);

            buffers["F"][i] = f;
        }

        LOG_TRACE("··· Computing F ifft");
        polynomials["F"] = new Polynomial<Engine>(buffers["F"], zkey->domainSize * 16);

        // Check degree
        if (polynomials["F"]->degree >= 9 * zkey->domainSize + 22) {
            throw std::runtime_error("F Polynomial is not well calculated");
        }

        delete buffers["F"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZT() {
        polynomials["ZT"] = Polynomial<Engine>::zerofierPolynomial(
                {roots["S1h1"][0], roots["S1h1"][1], roots["S1h1"][2], roots["S1h1"][3],
                 roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                 roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]});
    }

    // ROUND 5
    template<typename Engine>
    void FflonkProver<Engine>::round5() {
        // STEP 5.1 - Compute random evaluation point y ∈ F
        // STEP 4.1 - Compute challenge alpha ∈ F
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>();
        transcript->addPolCommitment(proof->getPolynomialCommitment("W1"));

        challenges["y"] = transcript->getChallenge();

        std::ostringstream ss;
        ss << "challenges.y: " << E.Fr.toString(challenges["y"]);
        LOG_TRACE(ss);

        // STEP 5.2 - Compute L(X)
        computeL();
        computeZTS2();

        FrElement ZTS2Y = polynomials["ZTS2"]->evaluate(challenges["y"]);
        ZTS2Y = E.Fr.inv(ZTS2Y);
        polynomials["L"].mulScalar(ZTS2Y);

        Polynomial<Engine> *polDividend = new Polynomial<Engine>({E.Fr.neg(challenges["y"]), E.Fr.one});
        Polynomial<Engine> *polRemainder = polynomials["L"]->divBy(polDividend);

        // Check degrees
        if (polRemainder->degree > 0) {
            ss.str("");
            ss << "Degree of L(X)/(ZTS2(y)(X-y)) remainder is " << polRemainder->degree << " and should be 0";
            throw std::runtime_error(ss);
        }
        if (polynomials["L"]->degree >= 9 * zkey->domainSize + 17) {
            throw std::runtime_error("Degree of L(X)/(ZTS2(y)(X-y)) is not correct");
        }

        // The fifth output of the prover is ([W2]_1), where W2:=(f/Z_t)(x)
        proof.addPolynomial("W2", multiExponentiation(polynomials["L"], "W2"));
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeL() {
        buffers["L"] = new FrElement[zkey->domainSize * 16];

        FrElement evalR1Y = polynomials["R1"]->evaluate(challenges["y"]);
        FrElement evalR2Y = polynomials["R2"]->evaluate(challenges["y"]);
        FrElement evalZTY = polynomials["ZT"]->evaluate(challenges["y"]);

        FrElement preL1 = E.Fr.sub(challenges["y"], roots["S2h2"][0]);
        preL1 = E.Fr.mul(preL1, E.Fr.sub(challenges["y"], roots["S2h2"][1]));
        preL1 = E.Fr.mul(preL1, E.Fr.sub(challenges["y"], roots["S2h2"][2]));
        preL1 = E.Fr.mul(preL1, E.Fr.sub(challenges["y"], roots["S2h3"][0]));
        preL1 = E.Fr.mul(preL1, E.Fr.sub(challenges["y"], roots["S2h3"][1]));
        preL1 = E.Fr.mul(preL1, E.Fr.sub(challenges["y"], roots["S2h3"][2]));
        toInverse["yBatch"] = preL1;

        FrElement preL2 = E.Fr.mul(challenges["alpha"], E.Fr.sub(challenges["y"], roots["S1h1"][0]));
        preL2 = E.Fr.mul(preL2, E.Fr.sub(challenges["y"], roots["S1h1"][1]));
        preL2 = E.Fr.mul(preL2, E.Fr.sub(challenges["y"], roots["S1h1"][2]));
        preL2 = E.Fr.mul(preL2, E.Fr.sub(challenges["y"], roots["S1h1"][3]));

        LOG_TRACE("> Computing F fft");

        evaluations["F"] = new Evaluations<Engine>(polynomials["F"]);

        LOG_TRACE("> Computing L");

        // Set initial omega
        std::ostringstream ss;
        for (u_int64_t i = 0; i < zkey->domainSize * 16; i++) {
            if ((0 != i) && (i % 5000 == 0)) {
                ss.str("");
                ss << "> Computing L evaluation" << i / (zkey->domainSize * 16);
                LOG_TRACE(ss);
            }

            FrElement omega =  omegaBuffer[i * 16];

            FrElement c1 = evaluations["C1"]->getEvaluation(i * 4);
            FrElement c2 = evaluations["C2"]->getEvaluation(i * 4);
            FrElement f = evaluations["F"]->getEvaluation(i * 4);

            // l1 = (y - h2) (y - h2w3) (y - h2w3_2) (y - h3) (y - h3w3) (y - h3w3_2) (C1(X) - R1(y))
            FrElement l1 = E.Fr.mul(preL1, E.Fr.sub(c1, evalR1Y));

            // l2 = alpha (y - h1) (y - h1w4) (y - h1w4_2) (y - h1w4_3) (C2(X) - R2(y))
            FrElement l2 = E.Fr.mul(preL2, E.Fr.sub(c2, evalR2Y));

            // l3 = ZT(y) (f(X)/ZT(X))
            // Recall f is already a f(X)/ZT(X)
            FrElement l3 = E.Fr.mul(evalZTY, f);

            FrElement l = E.Fr.sub(E.Fr.add(l1, l2), l3);

            buffers["L"][i] = l;
        }

        LOG_TRACE("··· Computing L ifft");
        polynomials["L"] = new Polynomial<Engine>(buffers["L"], zkey->domainSize * 16);

        // Check degree
        if (polynomials["L"]->degree >= 9 * zkey->domainSize + 18) {
            throw std::runtime_error("L Polynomial is not well calculated");
        }

        delete buffers["L"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZTS2() {
        polynomials["ZTS2"] = Polynomial<Engine>::zerofierPolynomial(
                {roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                 roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]});
    }

    //TODO !!!!!!!!!!
    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getMontgomeryBatchedInverse() {
        return E.fr.zero;
    }

    template<typename Engine>
    typename Engine::G1Point FflonkProver<Engine>::expTau(const FrElement *polynomial, int64_t from, int64_t count) {
//        G1P value;
//        FrElements pol = polynomialFromMontgomery(polynomial, from, count);
//
//        E.g1.multiMulByScalar(value, ptau.data(), (uint8_t *) pol.data(), sizeof(pol[0]), pol.size());
//        return value;
    }
}