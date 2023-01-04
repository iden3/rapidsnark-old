#include "fflonk_prover.hpp"

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
        E = Engine::engine;
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
            evaluations = new std::map<std::string, Evaluation<Engine>>;
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
            evaluations["Sigma1"] = new Evaluation<Engine>(zkey->domainSize);
            evaluations["Sigma2"] = new Evaluation<Engine>(zkey->domainSize);
            evaluations["Sigma3"] = new Evaluation<Engine>(zkey->domainSize);

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

        for (u_int64_t i = 0; i < zkey->nAdditions; i++) {
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
        u_int64_t diff = zkey->nVars - zkey->nAdditions;

        if (idx < diff) {
            return buffWitness[index];
        } else if (idx < zkey->nVars) {
            return buffInternalWitness[idx - diff];
        }

        return E.fr.zero();
    }

    //// ROUND 1
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
//        // Build A, B and C evaluations buffer from zkey and witness files
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
        for (u_int64_t i = 0; i < zkey->nConstraints; ++i) {
            FrElement aux;
            FrElement witness = getWitness(buffers[polName][i]);
            buffers[polName][i] = E.fr.toMontgomery(witness);
        }

        // Create the polynomial
        // and compute the coefficients of the wire polynomials from evaluations
        u_int64_t bFactorsLen = sizeof(*blindingFactors) / sizeof(Engine::FrElement);
        polynomials[polName] = new Polynomial<Engine>(buffers[polName], zkey->domainSize, bFactorsLen);

        // Compute the extended evaluations of the wire polynomials
        evaluations[polName] = new Evaluation<Engine>(polynomials[polName]);

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
        evaluations["QL"] = new Evaluation<Engine>(zkey->domainSize * 4);
        evaluations["QR"] = new Evaluation<Engine>(zkey->domainSize * 4);
        evaluations["QM"] = new Evaluation<Engine>(zkey->domainSize * 4);
        evaluations["QO"] = new Evaluation<Engine>(zkey->domainSize * 4);
        evaluations["QC"] = new Evaluation<Engine>(zkey->domainSize * 4);

        // Read Q's evaluations from zkey file
        memcpy(evaluations["QL"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QR"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QR_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QM"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QM_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QO"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QO_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QC"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QC_SECTION) + sDomain, sDomain * 4);

        // Read Lagrange polynomials & evaluations from zkey file
        evaluations["lagrangePolynomials"] = new Evaluation<Engine>(zkey->domainSize * 5);
        memcpy(evaluations["lagrangePolynomials"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_LAGRANGE_SECTION), sDomain * 5);

        // Reserve memory for buffers T0 and T0z
        buffers["T0"] =  new FrElement[zkey->domainSize * 4];
        buffers["T0z"] =  new FrElement[zkey->domainSize * 4];

        LOG_TRACE("> Computing T0");
        #pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize * 4; i++) {
            if ((0 != i) && (i % 5000 == 0)) {
                ss << "> Computing t0 evaluation" << i / (zkey->domainSize * 4);
                LOG_TRACE(ss);
            }

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
            FrElement az = E.fr.add(E.fr.mul(blindingFactors[1], omegaBuffer[i] * 4), blindingFactors[2]);
            FrElement bz = E.fr.add(E.fr.mul(blindingFactors[3], omegaBuffer[i] * 4), blindingFactors[4]);
            FrElement cz = E.fr.add(E.fr.mul(blindingFactors[5], omegaBuffer[i] * 4), blindingFactors[6]);

            // Compute current public input
            FrElement pi = E.Fr.zero;
            for (u_int64_t j = 0; j < zkey->nPublic; j++) {
                u_int32_t offset = (j * 5 * zkey->domainSize) + zkey->domainSize + i;

                FrElement lPol = evaluations["lagrange1"]->getEvaluation(offset);
                FrElement aVal = buffers["A"][j];

                pi = E.Fr.sub(pi, E.Fr.mul(lPol, aVal));
            }

            //T0(X) = [q_L(X)·a(X) + q_R(X)·b(X) + q_M(X)·a(X)·b(X) + q_O(X)·c(X) + q_C(X) + PI(X)] · 1/Z_H(X)
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
        u_int64_t lengthA = polynomials["A"]->length;
        u_int64_t lengthB = polynomials["B"]->length;
        u_int64_t lengthC = polynomials["C"]->length;
        u_int64_t lengthT0 = polynomials["T0"]->length;

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


    //// ROUND 2
    template<typename Engine>
    void FflonkProver<Engine>::round2() {
        // STEP 2.1 - Compute permutation challenge beta and gamma ∈ F
        // Compute permutation challenge beta
        Keccak256Transcript transcript = new Keccak256Transcript<Engine>();
        for (u_int32_t i = 0; i < zkey->nPublic; i++) {
            transcript.addScalar(buffers["A"][i]);
        }
        transcript.addPolCommitment(proof.getPolynomial("C1"));

        challenges["beta"] = transcript.getChallenge();
        std::ostringstream ss;
        ss << "challenges.beta: " << E.Fr.toString(challenges["beta"]);
        LOG_TRACE(ss);

        // Compute permutation challenge gamma
        transcript.reset();
        transcript.addScalar(challenges["beta"]);
        challenges["gamma"] = transcript.getChallenge();
        std::ostringstream ss;
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

    }

    template<typename Engine>
    void FflonkProver<Engine>::computeT1() {

    }

    template<typename Engine>
    void FflonkProver<Engine>::computeT2() {

    }

    template<typename Engine>
    void FflonkProver<Engine>::computeC2() {

    }

    //// ROUND 3
    template<typename Engine>
    void FflonkProver<Engine>::round3() {
    }

    //// ROUND 4
    template<typename Engine>
    void FflonkProver<Engine>::round4() {
    }

    //// ROUND 5
    template<typename Engine>
    void FflonkProver<Engine>::round5() {
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