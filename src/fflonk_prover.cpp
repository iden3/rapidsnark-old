#include "fflonk_prover.hpp"

#include "curve_utils.hpp"
#include "zkey.hpp"
#include "zkey_fflonk.hpp"
#include "wtns_utils.hpp"
#include <sodium.h>
#include "mul_z.hpp"
#include "keccak_256_transcript.hpp"
#include "logger.hpp"

using namespace CPlusPlusLogging;

namespace Fflonk {

    template<typename Engine>
    FflonkProver<Engine>::FflonkProver(Engine &_E) : E(_E) {
        curveName = CurveUtils::getCurveNameByEngine();
    }

    template<typename Engine>
    FflonkProver<Engine>::~FflonkProver() {
        delete fft;
    }

    template<typename Engine>
    std::tuple <json, json> FflonkProver<Engine>::prove(BinFileUtils::BinFile *fdZkey, BinFileUtils::BinFile *fdWtns) {
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

            mulZ = new MulZ<Engine>(E, fft);

            LOG_TRACE("> Computing omegaBuffer");
            omegaBuffer = new FrElement[zkey->domainSize * 16];
            FrElement omega = fft->root(zkeyPower + 4, 1);
            omegaBuffer[0] = E.fr.one();
            for (int64_t i = 1; i < zkey->domainSize * 16; ++i) {
                E.fr.mul(omegaBuffer[i], omegaBuffer[i - 1], omega);
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

//            buffers = new std::map<string, FrElement[]>;
//            polynomials = new std::map<std::string, Polynomial<Engine>>;
//            evaluations = new std::map<std::string, Evaluations<Engine>>;
//            toInverse = new std::map<std::string, FrElement>;
//            challenges = new std::map<std::string, FrElement>;
//            roots = new std::map<std::string, FrElement[]>;

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

            proof = new SnarkProof<Engine>(E, "fflonk");

            ss << "> Reading Section " << Zkey::ZKEY_FF_ADDITIONS_SECTION << ". Additions\n";
            LOG_TRACE(ss);

            calculateAdditions(fdZkey);

            ss << "> Reading Section " << Zkey::ZKEY_FF_SIGMA1_SECTION << "," << Zkey::ZKEY_FF_SIGMA2_SECTION
               << "," << Zkey::ZKEY_FF_SIGMA3_SECTION << ". Sigma1, Sigma2 & Sigma 3\n";
            LOG_TRACE(ss);

            LOG_TRACE("··· Reading Sigma polynomials ");
            polynomials["Sigma1"] = new Polynomial<Engine>(E, zkey->domainSize);
            polynomials["Sigma2"] = new Polynomial<Engine>(E, zkey->domainSize);
            polynomials["Sigma3"] = new Polynomial<Engine>(E, zkey->domainSize);

            memcpy(polynomials["Sigma1"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION), sDomain);
            memcpy(polynomials["Sigma2"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION), sDomain);
            memcpy(polynomials["Sigma3"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION), sDomain);

            LOG_TRACE("··· Reading Sigma evaluations ");
            evaluations["Sigma1"] = new Evaluations<Engine>(E, zkey->domainSize);
            evaluations["Sigma2"] = new Evaluations<Engine>(E, zkey->domainSize);
            evaluations["Sigma3"] = new Evaluations<Engine>(E, zkey->domainSize);

            memcpy(evaluations["Sigma1"]->eval, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION) + sDomain, sDomain * 4);
            memcpy(evaluations["Sigma2"]->eval, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION) + sDomain, sDomain * 4);
            memcpy(evaluations["Sigma3"]->eval, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION) + sDomain, sDomain * 4);

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

//            proof->addEvaluationCommitment("inv", getMontgomeryBatchedInverse());

            // Prepare public inputs
            json publicSignals;
            for (u_int32_t i = 1; i <= zkey->nPublic; i++) {
                E.fr.toMontgomery(buffWitness[i], buffWitness[i]);
                publicSignals.push_back(E.fr.toString(buffWitness[i]).c_str());
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
    void FflonkProver<Engine>::calculateAdditions(BinFileUtils::BinFile *fdZkey) {
        Zkey::Addition<Engine> *additionsBuff = (Zkey::Addition<Engine> *)fdZkey->getSectionData(Zkey::ZKEY_FF_ADDITIONS_SECTION);

        for (u_int32_t i = 0; i < zkey->nAdditions; i++) {
            // Get witness value
            FrElement witness1 = getWitness(additionsBuff[i].signalId1);
            FrElement witness2 = getWitness(additionsBuff[i].signalId2);

            //Calculate final result
            E.fr.mul(witness1, additionsBuff[i].factor1, witness1);
            E.fr.mul(witness2, additionsBuff[i].factor2, witness2);
            E.fr.add(buffInternalWitness[i], witness1, witness2);
        }
    }

    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getWitness(u_int64_t idx) {
        u_int32_t diff = zkey->nVars - zkey->nAdditions;

        if (idx < diff) {
            return buffWitness[idx];
        } else if (idx < zkey->nVars) {
            return buffInternalWitness[idx - diff];
        }

        return E.fr.zero();
    }

    // ROUND 1
    template<typename Engine>
    void FflonkProver<Engine>::round1() {
        // STEP 1.1 - Generate random blinding scalars (b_1, ..., b9) ∈ F
//        blindingFactors = new FrElement[10];

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
        G1Point c1 = multiExponentiation(polynomials["C1"], "C1");
        proof->addPolynomialCommitment("C1", c1);
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
                    FrElement bFactors[2] = {blindingFactors[2], blindingFactors[1]};
                    computeWirePolynomial("A", bFactors);
                }

#pragma omp section
                {
                    FrElement bFactors[2] = {blindingFactors[4], blindingFactors[3]};
                    computeWirePolynomial("B", bFactors);
                }

#pragma omp section
                {
                    FrElement bFactors[2] = {blindingFactors[6], blindingFactors[5]};
                    computeWirePolynomial("C", bFactors);
                }
            }
        }

        // Check degrees
        if (polynomials["A"]->getDegree() >= zkey->domainSize + 2) {
            throw std::runtime_error("A Polynomial is not well calculated");
        }
        if (polynomials["B"]->getDegree() >= zkey->domainSize + 2) {
            throw std::runtime_error("B Polynomial is not well calculated");
        }
        if (polynomials["C"]->getDegree() >= zkey->domainSize + 2) {
            throw std::runtime_error("C Polynomial is not well calculated");
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeWirePolynomial(std::string polName, FrElement blindingFactors[]) {

        // Compute all witness from signal ids and set them to the polynomial buffers
#pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->nConstraints; ++i) {
            FrElement aux;
//TODO remove            FrElement witness = getWitness(buffers[polName][i]);
            FrElement witness = E.fr.one();
            E.fr.toMontgomery(buffers[polName][i], witness);
        }

        // Create the polynomial
        // and compute the coefficients of the wire polynomials from evaluations
//TODO remove        u_int32_t bFactorsLen = sizeof(*blindingFactors) / sizeof(Engine::FrElement);
        u_int32_t bFactorsLen = 1;
        polynomials[polName] = new Polynomial<Engine>(E, buffers[polName], zkey->domainSize, bFactorsLen);

        // Compute the extended evaluations of the wire polynomials
        evaluations[polName] = new Evaluations<Engine>(E, *polynomials[polName]);

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
        evaluations["QL"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QR"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QM"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QO"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QC"] = new Evaluations<Engine>(E, zkey->domainSize * 4);

        // Read Q's evaluations from zkey file
        memcpy(evaluations["QL"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QR"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QR_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QM"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QM_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QO"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QO_SECTION) + sDomain, sDomain * 4);
        memcpy(evaluations["QC"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QC_SECTION) + sDomain, sDomain * 4);

        // Read Lagrange polynomials & evaluations from zkey file
        evaluations["lagrangePolynomials"] = new Evaluations<Engine>(E, zkey->domainSize * 5);
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
            FrElement az, bz, cz;
            E.fr.mul(az, blindingFactors[1], omega);
            E.fr.add(az, az, blindingFactors[2]);
            E.fr.mul(bz, blindingFactors[3], omega);
            E.fr.add(bz, bz, blindingFactors[4]);
            E.fr.mul(cz, blindingFactors[5], omega);
            E.fr.add(cz, cz, blindingFactors[6]);

            // Compute current public input
            FrElement pi = E.fr.zero();
            for (u_int32_t j = 0; j < zkey->nPublic; j++) {
                u_int32_t offset = (j * 5 * zkey->domainSize) + zkey->domainSize + i;

                FrElement lPol = evaluations["lagrange1"]->getEvaluation(offset);
                FrElement aVal = buffers["A"][j];
                E.fr.mul(lPol, lPol, aVal);
                E.fr.sub(pi, pi, lPol);
            }

            // T0(X) = [q_L(X)·a(X) + q_R(X)·b(X) + q_M(X)·a(X)·b(X) + q_O(X)·c(X) + q_C(X) + PI(X)] · 1/Z_H(X)
            // Compute first T0(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            // expression 1 -> q_L(X)·a(X)
            FrElement e1, e1z;
            E.fr.mul(e1, a, ql);
            E.fr.mul(e1z, az, ql);

            // expression 2 -> q_R(X)·b(X)
            FrElement e2, e2z;
            E.fr.mul(e2, b, qr);
            E.fr.mul(e2z, bz, qr);

            // expression 3 -> q_M(X)·a(X)·b(X)
            auto [e3, e3z] = mulZ->mul2(a, b, az, bz, i % 4);
            E.fr.mul(e2, e3, qm);
            E.fr.mul(e3z, e3z, qm);

            // expression 4 -> q_O(X)·c(X)
            FrElement e4, e4z;
            E.fr.mul(e4, c, qo);
            E.fr.mul(e4z, cz, qo);

            // t0 = expressions 1 + expression 2 + expression 3 + expression 4 + qc + pi
            FrElement t0, t0z;
            E.fr.add(t0, qc, pi);
            E.fr.add(t0, t0, e4);
            E.fr.add(t0, t0, e3);
            E.fr.add(t0, t0, e2);
            E.fr.add(t0, t0, e1);

            E.fr.add(t0z, e3z, e4z);
            E.fr.add(t0z, t0z, e2z);
            E.fr.add(t0z, t0z, e1z);

            buffers["T0"][i] = t0;
            buffers["T0z"][i] = t0z;
        }

        // Compute the coefficients of the polynomial T0(X) from buffers.T0
        LOG_TRACE("··· Computing T0 ifft");
        polynomials["T0"] = new Polynomial<Engine>(E, buffers["T0"], zkey->domainSize * 4);

        // Divide the polynomial T0 by Z_H(X)
        polynomials["T0"]->divZh(zkey->domainSize);

        // Compute the coefficients of the polynomial T0z(X) from buffers.T0z
        LOG_TRACE("··· Computing T0z ifft");
        polynomials["T0z"] = new Polynomial<Engine>(E, buffers["T0z"], zkey->domainSize * 4);

        // Add the polynomial T0z to T0 to get the final polynomial T0
        polynomials["T0"]->add(*polynomials["T0z"]);

        // Check degree
        if (polynomials["T0"]->getDegree() >= 2 * zkey->domainSize + 2) {
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
        u_int64_t lengthA = polynomials["A"]->getLength();
        u_int64_t lengthB = polynomials["B"]->getLength();
        u_int64_t lengthC = polynomials["C"]->getLength();
        u_int64_t lengthT0 = polynomials["T0"]->getLength();

        // Compute degree of the new polynomial C1 to reserve the buffer memory size
        // Will be the next power of two to bound the maximum(deg(A_4), deg(B_4)+1, deg(C_4)+2, deg(T0_4)+3)
        u_int64_t degreeA = polynomials["A"]->getDegree();
        u_int64_t degreeB = polynomials["B"]->getDegree();
        u_int64_t degreeC = polynomials["C"]->getDegree();
        u_int64_t degreeT0 = polynomials["T0"]->getDegree();

        u_int64_t maxLength = std::max(lengthA, std::max(lengthB, std::max(lengthC, lengthT0)));
        u_int64_t maxDegree = std::max(degreeA * 4 + 1, std::max(degreeB * 4 + 2, std::max(degreeC * 4 + 3, degreeT0 * 4 + 3)));

        u_int64_t lengthBuffer = std::pow(2, fft->log2(maxDegree - 1) + 1);

        polynomials["C1"] = new Polynomial<Engine>(E, lengthBuffer);

        for (u_int64_t i = 0; i < maxLength; i++) {
            polynomials["C1"]->coef[i * 4] = polynomials["A"]->getCoef(i);
            polynomials["C1"]->coef[i * 4 + 1] = polynomials["B"]->getCoef(i);
            polynomials["C1"]->coef[i * 4 + 2] = polynomials["C"]->getCoef(i);
            polynomials["C1"]->coef[i * 4 + 3] = polynomials["T0"]->getCoef(i);
        }

        // Check degree
        if (polynomials["C1"]->getDegree() >= 8 * zkey->domainSize + 8) {
            throw std::runtime_error("C1 Polynomial is not well calculated");
        }
    }


    // ROUND 2
    template<typename Engine>
    void FflonkProver<Engine>::round2() {
        // STEP 2.1 - Compute permutation challenge beta and gamma ∈ F
        // Compute permutation challenge beta
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>(E);
        for (u_int32_t i = 0; i < zkey->nPublic; i++) {
            transcript->addScalar(buffers["A"][i]);
        }
        transcript->addPolCommitment(proof->getPolynomialCommitment("C1"));

        challenges["beta"] = transcript->getChallenge();
        std::ostringstream ss;
        ss << "challenges.beta: " << E.fr.toString(challenges["beta"]);
        LOG_TRACE(ss);

        // Compute permutation challenge gamma
        transcript->reset();
        transcript->addScalar(challenges["beta"]);
        challenges["gamma"] = transcript->getChallenge();
        ss.str("");
        ss << "challenges.gamma: " << E.fr.toString(challenges["gamma"]);

        // STEP 2.2 - Compute permutation polynomial z(X)
        computeZ();

        // STEP 2.3 - Compute quotient polynomial T1(X) and T2(X)
        computeT1();
        computeT2();

        // STEP 2.4 - Compute the FFT-style combination polynomial C2(X)
        computeC2();

        // The second output of the prover is ([C2]_1)
        G1Point c2 = multiExponentiation(polynomials["C2"], "C2");
        proof->addPolynomialCommitment("C2", c2);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZ() {
        FrElement *numArr = new FrElement[zkey->domainSize];
        FrElement *denArr = new FrElement[zkey->domainSize];

        // Set the first values to 1
        numArr[0] = E.fr.one();
        denArr[0] = E.fr.one();

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
            FrElement betaw;
            E.fr.mul(betaw, challenges["beta"], omega);

            FrElement num1;
            E.fr.add(num1, buffers["A"][i], betaw);
            E.fr.add(num1, num1, challenges["gamma"]);

            FrElement num2;
            E.fr.mul(num2, *((FrElement *)zkey->k1), betaw);
            E.fr.add(num2, num2,buffers["B"][i]);
            E.fr.add(num2, num2, challenges["gamma"]);

            FrElement num3;
            E.fr.mul(num3, *((FrElement *)zkey->k2), betaw);
            E.fr.add(num3, num3, buffers["C"][i]);
            E.fr.add(num3, num3, challenges["gamma"]);

            FrElement num;
            E.fr.mul(num, num2, num3);
            E.fr.mul(num, num, num1);

            // denArr := (a + beta·sigma1 + gamma)(b + beta·sigma2 + gamma)(c + beta·sigma3 + gamma)
            FrElement den1 = evaluations["Sigma1"]->getEvaluation(i * 4);
            E.fr.mul(den1, den1, challenges["beta"]);
            E.fr.add(den1, den1, buffers["A"][i]);
            E.fr.add(den1, den1, challenges["gamma"]);

            FrElement den2 = evaluations["Sigma2"]->getEvaluation(i * 4);
            E.fr.mul(den2, den2, challenges["beta"]);
            E.fr.add(den2, den2, buffers["B"][i]);
            E.fr.add(den2, den2, challenges["gamma"]);

            FrElement den3 = evaluations["Sigma3"]->getEvaluation(i * 4);
            E.fr.mul(den3, den3, challenges["beta"]);
            E.fr.add(den3, den3, buffers["C"][i]);
            E.fr.add(den3, den3, challenges["gamma"]);

            FrElement den;
            E.fr.mul(den, den2, den3);
            E.fr.mul(den, den, den1);

            // Multiply current num value with the previous one saved in numArr
            E.fr.mul(num, num, numArr[i]);
            numArr[(i + 1) % zkey->domainSize] = num;

            // Multiply current den value with the previous one saved in denArr
            E.fr.mul(den, denArr[i], den);
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
            E.fr.mul(buffers["Z"][i], numArr[i], denArr[i]);
        }

        if (!E.fr.eq(numArr[0], E.fr.one())) {
            throw std::runtime_error("Copy constraints does not match");
        }

        // Compute polynomial coefficients z(X) from buffers.Z
        LOG_TRACE("··· Computing Z ifft");
        polynomials["Z"] = new Polynomial<Engine>(E, buffers["Z"], zkey->domainSize);

        // Compute extended evaluations of z(X) polynomial
        evaluations["Z"] = new Evaluations<Engine>(E, *polynomials["Z"]);

        // Blind z(X) polynomial coefficients with blinding scalars b
        FrElement bFactors[3] = {blindingFactors[9], blindingFactors[8], blindingFactors[7]};
        polynomials["Z"]->blindCoefficients(bFactors);

        // Check degree
        if (polynomials["Z"]->getDegree() >= zkey->domainSize + 3) {
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

            FrElement omega = omegaBuffer[i];
            FrElement omega2;
            E.fr.square(omega2, omega);

            FrElement z = evaluations["Z"]->getEvaluation(i);
            FrElement zp, zptmp;
            E.fr.mul(zp, blindingFactors[7], omega2);
            E.fr.mul(zptmp, blindingFactors[8], omega);
            E.fr.add(zp, zp, zptmp);
            E.fr.add(zp, zp, blindingFactors[9]);

            // T1(X) := (z(X) - 1) · L_1(X)
            // Compute first T1(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            FrElement lagrange1 = evaluations["lagrange1"]->getEvaluation(zkey->domainSize + i);
            FrElement t1, t1z;
            E.fr.sub(t1, z, E.fr.one());
            E.fr.mul(t1, t1, lagrange1);
            E.fr.mul(t1z, zp, lagrange1);

            buffers["T1"][i] = t1;
            buffers["T1z"][i] = t1z;
        }

        // Compute the coefficients of the polynomial T1(X) from buffers.T1
        LOG_TRACE("··· Computing T1 ifft");
        polynomials["T1"] = new Polynomial<Engine>(E, buffers["T1"], zkey->domainSize * 4);

        // Divide the polynomial T1 by Z_H(X)
        polynomials["T1"]->divZh(zkey->domainSize);

        // Compute the coefficients of the polynomial T1z(X) from buffers.T1z
        LOG_TRACE("··· Computing T1z ifft");
        polynomials["T1z"] = new Polynomial<Engine>(E, buffers["T1z"], zkey->domainSize * 4);

        // Add the polynomial T0z to T0 to get the final polynomial T0
        polynomials["T1"]->add(*polynomials["T1z"]);

        // Check degree
        if (polynomials["T1"]->getDegree() >= zkey->domainSize + 2) {
            throw std::runtime_error("T1 Polynomial is not well calculated");
        }

        delete buffers["T1"];
        delete buffers["T1z"];
        delete polynomials["T1z"];
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
            FrElement omega2;
            E.fr.square(omega2, omega);
            FrElement omegaW;
            E.fr.mul(omegaW, omega, fft->root(zkeyPower, 1));
            FrElement omegaW2;
            E.fr.square(omegaW2, omegaW);

            FrElement a = evaluations["A"]->getEvaluation(i);
            FrElement b = evaluations["B"]->getEvaluation(i);
            FrElement c = evaluations["C"]->getEvaluation(i);
            FrElement z = evaluations["Z"]->getEvaluation(i);
            FrElement zW = evaluations["Z"]->getEvaluation((zkey->domainSize * 4 + 4 + i) % (zkey->domainSize * 4));

            FrElement ap;
            E.fr.mul(ap, blindingFactors[1], omega);
            E.fr.add(ap, ap, blindingFactors[2]);
            FrElement bp;
            E.fr.mul(bp, blindingFactors[3], omega);
            E.fr.add(bp, bp, blindingFactors[4]);
            FrElement cp;
            E.fr.mul(cp, blindingFactors[5], omega);
            E.fr.add(cp, cp, blindingFactors[6]);
            FrElement zp, zptmp;
            E.fr.mul(zp, blindingFactors[7], omega2);
            E.fr.mul(zptmp, blindingFactors[8], omega);
            E.fr.add(zptmp, zptmp, blindingFactors[9]);
            E.fr.add(zp, zp, zptmp);
            FrElement zWp, zWptmp;
            E.fr.mul(zWp, blindingFactors[7], omegaW2);
            E.fr.mul(zWptmp, blindingFactors[8], omegaW);
            E.fr.add(zWptmp, zWptmp, blindingFactors[9]);
            E.fr.add(zWp, zWp, zWptmp);

            FrElement sigma1 = evaluations["Sigma1"]->getEvaluation(i);
            FrElement sigma2 = evaluations["Sigma2"]->getEvaluation(i);
            FrElement sigma3 = evaluations["Sigma3"]->getEvaluation(i);

            // T2(X) := [ (a(X) + beta·X + gamma)(b(X) + beta·k1·X + gamma)(c(X) + beta·k2·X + gamma)z(X)
            //           -(a(X) + beta·sigma1(X) + gamma)(b(X) + beta·sigma2(X) + gamma)(c(X) + beta·sigma3(X) + gamma)z(Xω)] · 1/Z_H(X)
            // Compute first T2(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)

            // expression 1 -> (a(X) + beta·X + gamma)(b(X) + beta·k1·X + gamma)(c(X) + beta·k2·X + gamma)z(X)
            FrElement betaX;
            E.fr.mul(betaX, challenges["beta"], omega);

            FrElement e11;
            E.fr.add(e11, a, betaX);
            E.fr.add(e11, e11, challenges["gamma"]);

            FrElement e12;
            E.fr.mul(e12, betaX, *((FrElement *)zkey->k1));
            E.fr.add(e12, e12, b);
            E.fr.add(e12, e12, challenges["gamma"]);

            FrElement e13;
            E.fr.mul(e13, betaX, *((FrElement *)zkey->k2));
            E.fr.add(e13, e13, c);
            E.fr.add(e13, e13, challenges["gamma"]);

            auto [e1, e1z] = mulZ->mul4(e11, e12, e13, z, ap, bp, cp, zp, i % 4);

            // expression 2 -> (a(X) + beta·sigma1(X) + gamma)(b(X) + beta·sigma2(X) + gamma)(c(X) + beta·sigma3(X) + gamma)z(Xω)
            FrElement e21;
            E.fr.mul(e21, challenges["beta"], sigma1);
            E.fr.add(e21, e21, a);
            E.fr.add(e21, e21, challenges["gamma"]);

            FrElement e22;
            E.fr.mul(e22, challenges["beta"], sigma2);
            E.fr.add(e22, e22, b);
            E.fr.add(e22, e22, challenges["gamma"]);

            FrElement e23;
            E.fr.mul(e23, challenges["beta"], sigma3);
            E.fr.add(e23, e23, c);
            E.fr.add(e23, e23, challenges["gamma"]);

            auto [e2, e2z] = mulZ->mul4(e21, e22, e23, zW, ap, bp, cp, zWp, i % 4);

            FrElement t2, t2z;
            E.fr.sub(t2, e1, e2);
            E.fr.sub(t2z, e1z, e2z);

            buffers["T2"][i] = t2;
            buffers["T2z"][i] = t2z;
        }

        // Compute the coefficients of the polynomial T2(X) from buffers.T2
        LOG_TRACE("··· Computing T2 ifft");
        polynomials["T2"] = new Polynomial<Engine>(E, buffers["T2"], zkey->domainSize * 4);

        // Divide the polynomial T2 by Z_H(X)
        polynomials["T2"]->divZh(zkey->domainSize);

        // Compute the coefficients of the polynomial T2z(X) from buffers.T2z
        LOG_TRACE("··· Computing T2z ifft");
        polynomials["T2z"] = new Polynomial<Engine>(E, buffers["T2z"], zkey->domainSize * 4);

        // Add the polynomial T2z to T2 to get the final polynomial T2
        polynomials["T2"]->add(*polynomials["T2z"]);

        // Check degree
        if (polynomials["T2"]->getDegree() >= 3 * zkey->domainSize + 6) {
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
        u_int64_t lengthZ = polynomials["Z"]->getLength();
        u_int64_t lengthT1 = polynomials["T1"]->getLength();
        u_int64_t lengthT2 = polynomials["T2"]->getLength();
        // Compute degree of the new polynomial C2(X) to reserve the buffer memory size
        // Will be the maximum(deg(Z_3), deg(T1_3)+1, deg(T2_3)+2)
        u_int64_t degreeZ = polynomials["Z"]->getDegree();
        u_int64_t degreeT1 = polynomials["T1"]->getDegree();
        u_int64_t degreeT2 = polynomials["T2"]->getDegree();

        u_int64_t maxLength = std::max(lengthZ, std::max(lengthT1, lengthT2));
        u_int64_t maxDegree = std::max(degreeZ * 3 + 1, std::max(degreeT1 * 3 + 2, degreeT2 * 3 + 3));

        u_int64_t lengthBuffer = std::pow(2, fft->log2(maxDegree - 1) + 1);

        polynomials["C2"] = new Polynomial<Engine>(E, lengthBuffer);

        for (u_int64_t i = 0; i < maxLength; i++) {
            polynomials["C2"]->coef[i * 4] = polynomials["Z"]->getCoef(i);
            polynomials["C2"]->coef[i * 4 + 1] = polynomials["T1"]->getCoef(i);
            polynomials["C2"]->coef[i * 4 + 2] = polynomials["T2"]->getCoef(i);
        }

        // Check degree
        if (polynomials["C2"]->getDegree() >= 9 * zkey->domainSize + 18) {
            throw std::runtime_error("C2 Polynomial is not well calculated");
        }
    }

    // ROUND 3
    template<typename Engine>
    void FflonkProver<Engine>::round3() {
        // STEP 3.1 - Compute evaluation challenge xi ∈ S
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>(E);
        transcript->addPolCommitment(proof->getPolynomialCommitment("C2"));

        // Obtain a xi_seeder from the transcript
        // To force h1^4 = xi, h2^3 = xi and h_3^2 = xiω
        // we compute xi = xi_seeder^12, h1 = xi_seeder^3, h2 = xi_seeder^4 and h3 = xi_seeder^6
        FrElement xiSeed = transcript->getChallenge();
        FrElement xiSeed2;
        E.fr.square(xiSeed2, xiSeed);

        // Compute omega3 and omega4
        roots["w4"] = new FrElement[4];
        roots["w4"][0] = E.fr.one();
        roots["w4"][1] = *((FrElement *)zkey->w4);
        E.fr.square(roots["w4"][2], roots["w4"][1]);
        E.fr.mul(roots["w4"][3], roots["w4"][2], roots["w4"][1]);

        roots["w3"] = new FrElement[3];
        roots["w3"][0] = E.fr.one();
        roots["w3"][1] = *((FrElement *)zkey->w3);
        E.fr.square(roots["w3"][2], roots["w3"][1]);

        // Compute h1 = xi_seeder^3
        roots["S1h1"] = new FrElement[4];
        E.fr.mul(roots["S1h1"][0], xiSeed2, xiSeed);
        E.fr.mul(roots["S1h1"][1], roots["S1h1"][0], roots["w4"][1]);
        E.fr.mul(roots["S1h1"][2], roots["S1h1"][0], roots["w4"][2]);
        E.fr.mul(roots["S1h1"][3], roots["S1h1"][0], roots["w4"][3]);

        roots["S2h2"] = new FrElement[3];
        E.fr.square(roots["S2h2"][0], xiSeed2);
        E.fr.mul(roots["S2h2"][1], roots["S2h2"][0], roots["w3"][1]);
        E.fr.mul(roots["S2h2"][2], roots["S2h2"][0], roots["w3"][2]);

        roots["S2h3"] = new FrElement[3];
        // Multiply h3 by third-root-omega to obtain h_3^3 = xiω
        E.fr.mul(roots["S2h3"][0], roots["S2h2"][0], *((FrElement *)zkey->wr));
        E.fr.mul(roots["S2h3"][1], roots["S2h3"][0], roots["w3"][1]);
        E.fr.mul(roots["S2h3"][2], roots["S2h3"][0], roots["w3"][2]);

        // Compute xi = xi_seeder^12
        E.fr.square(challenges["xi"], roots["S2h2"][0]);
        E.fr.mul(challenges["xi"], challenges["xi"], roots["S2h2"][0]);

        std::ostringstream ss;
        ss << "challenges.xi: " << E.fr.toString(challenges["xi"]);
        LOG_TRACE(ss);

        // Reserve memory for Q's polynomials
        polynomials["QL"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QR"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QM"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QO"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QC"] = new Polynomial<Engine>(E, zkey->domainSize);

        // Read Q's evaluations from zkey file
        memcpy(polynomials["QL"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QR"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QM"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QO"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);
        memcpy(polynomials["QC"]->coef, (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION), sDomain);

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

        FrElement xiw;
        E.fr.mul(xiw, challenges["xi"], fft->root(zkeyPower + 4, 1));
        proof->addEvaluationCommitment("zw", polynomials["Z"]->evaluate(challenges["xiw"]));
        proof->addEvaluationCommitment("t1w", polynomials["T1"]->evaluate(challenges["xiw"]));
        proof->addEvaluationCommitment("t2w", polynomials["T2"]->evaluate(challenges["xiw"]));
    }

    // ROUND 4
    template<typename Engine>
    void FflonkProver<Engine>::round4() {
        // STEP 4.1 - Compute challenge alpha ∈ F
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>(E);
        transcript->addScalar(proof->getEvaluationCommitment("ql"));
        transcript->addScalar(proof->getEvaluationCommitment("qr"));
        transcript->addScalar(proof->getEvaluationCommitment("qm"));
        transcript->addScalar(proof->getEvaluationCommitment("qo"));
        transcript->addScalar(proof->getEvaluationCommitment("qc"));
        transcript->addScalar(proof->getEvaluationCommitment("s1"));
        transcript->addScalar(proof->getEvaluationCommitment("s2"));
        transcript->addScalar(proof->getEvaluationCommitment("s3"));
        transcript->addScalar(proof->getEvaluationCommitment("a"));
        transcript->addScalar(proof->getEvaluationCommitment("b"));
        transcript->addScalar(proof->getEvaluationCommitment("c"));
        transcript->addScalar(proof->getEvaluationCommitment("z"));
        transcript->addScalar(proof->getEvaluationCommitment("zw"));
        transcript->addScalar(proof->getEvaluationCommitment("t1w"));
        transcript->addScalar(proof->getEvaluationCommitment("t2w"));
        challenges["alpha"] = transcript->getChallenge();

        std::ostringstream ss;
        ss << "challenges.alpha: " << E.fr.toString(challenges["alpha"]);
        LOG_TRACE(ss);

        // STEP 4.2 - Compute F(X)
        computeR1();
        computeR2();

        computeF();
        computeZT();

        Polynomial<Engine>* polRemainder = polynomials["F"]->divBy(*polynomials["ZT"]);

        // Check degrees
        if (polRemainder->getDegree() > 0) {
            std::ostringstream ss;
            ss << "Degree of f(X)/ZT(X) remainder is " << polRemainder->getDegree() << " and should be 0";
            throw std::runtime_error(ss.str());
        }
        if (polynomials["F"]->getDegree() >= 9 * zkey->domainSize + 12) {
            throw std::runtime_error("Degree of f(X)/ZT(X) is not correct");
        }

        // The fourth output of the prover is ([W1]_1), where W1:=(f/Z_t)(x)
        G1Point w1 = multiExponentiation(polynomials["W1"], "W1");
        proof->addPolynomialCommitment("W1", w1);
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
        FrElement xArr[4] = {roots["S1h1"][0], roots["S1h1"][1], roots["S1h1"][2], roots["S1h1"][3]};
        FrElement yArr[4] = {polynomials["C1"]->evaluate(roots["S1h1"][0]), polynomials["C1"]->evaluate(roots["S1h1"][1]),
                             polynomials["C1"]->evaluate(roots["S1h1"][2]), polynomials["C1"]->evaluate(roots["S1h1"][3])};
        polynomials["R1"] = Polynomial<Engine>::lagrangePolynomialInterpolation(xArr, yArr);

        // Check the degree of r1(X) < 4
        if (polynomials["R1"]->getDegree() > 3) {
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
        FrElement xArr[6] = {roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                             roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]};
        FrElement yArr[6] = {polynomials["C2"]->evaluate(roots["S2h2"][0]), polynomials["C2"]->evaluate(roots["S2h2"][1]),
                             polynomials["C2"]->evaluate(roots["S2h2"][2]), polynomials["C2"]->evaluate(roots["S2h3"][0]),
                             polynomials["C2"]->evaluate(roots["S2h3"][1]), polynomials["C2"]->evaluate(roots["S2h3"][2])};

        polynomials["R2"] = Polynomial<Engine>::lagrangePolynomialInterpolation(xArr, yArr);

        // Check the degree of r2(X) < 6
        if (polynomials["R2"]->getDegree() > 5) {
            throw std::runtime_error("R2 Polynomial is not well calculated");
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeF() {
        buffers["F"] = new FrElement[zkey->domainSize * 16];

        LOG_TRACE("> Computing C1 & C2 fft");

        evaluations["C1"] = new Evaluations<Engine>(E, *polynomials["C1"]);
        evaluations["C2"] = new Evaluations<Engine>(E, *polynomials["C2"]);

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
            FrElement f1, ftmp;
            E.fr.sub(f1, omega, roots["S2h2"][0]);
            E.fr.sub(ftmp, omega, roots["S2h2"][1]);
            E.fr.mul(f1, f1, ftmp);
            E.fr.sub(ftmp, omega, roots["S2h2"][2]);
            E.fr.mul(f1, f1, ftmp);
            E.fr.sub(ftmp, omega, roots["S2h3"][0]);
            E.fr.mul(f1, f1, ftmp);
            E.fr.mul(f1, f1, ftmp);
            E.fr.sub(ftmp, omega, roots["S2h3"][1]);
            E.fr.mul(f1, f1, ftmp);
            E.fr.sub(ftmp, omega, roots["S2h3"][2]);
            E.fr.mul(f1, f1, ftmp);
            E.fr.sub(ftmp, c1, r1);
            E.fr.mul(f1, f1, ftmp);

            // f2 = alpha (X - h1) (X - h1w4) (X - h1w4_2) (X - h1w4_3) (C2(X) - R2(X))
            FrElement f2;
            E.fr.sub(ftmp, omega, roots["S1h1"][0]);
            E.fr.mul(f2, challenges["alpha"], ftmp);
            E.fr.sub(ftmp, omega, roots["S1h1"][1]);
            E.fr.mul(f2, f2, ftmp);
            E.fr.sub(ftmp, omega, roots["S1h1"][2]);
            E.fr.mul(f2, f2, ftmp);
            E.fr.sub(ftmp, omega, roots["S1h1"][3]);
            E.fr.mul(f2, f2, ftmp);
            E.fr.sub(ftmp, c2, r2);
            E.fr.mul(f2, f2, ftmp);

            FrElement f;
            E.fr.add(f, f1, f2);

            buffers["F"][i] = f;
        }

        LOG_TRACE("··· Computing F ifft");
        polynomials["F"] = new Polynomial<Engine>(E, buffers["F"], zkey->domainSize * 16);

        // Check degree
        if (polynomials["F"]->getDegree() >= 9 * zkey->domainSize + 22) {
            throw std::runtime_error("F Polynomial is not well calculated");
        }

        delete buffers["F"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZT() {
        FrElement arr[10] = {roots["S1h1"][0], roots["S1h1"][1], roots["S1h1"][2], roots["S1h1"][3],
                             roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                             roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]};

        polynomials["ZT"] = Polynomial<Engine>::zerofierPolynomial(arr);
    }

    // ROUND 5
    template<typename Engine>
    void FflonkProver<Engine>::round5() {
        // STEP 5.1 - Compute random evaluation point y ∈ F
        // STEP 4.1 - Compute challenge alpha ∈ F
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>(E);
        transcript->addPolCommitment(proof->getPolynomialCommitment("W1"));

        challenges["y"] = transcript->getChallenge();

        std::ostringstream ss;
        ss << "challenges.y: " << E.fr.toString(challenges["y"]);
        LOG_TRACE(ss);

        // STEP 5.2 - Compute L(X)
        computeL();
        computeZTS2();

        FrElement ZTS2Y = polynomials["ZTS2"]->evaluate(challenges["y"]);
        E.fr.inv(ZTS2Y, ZTS2Y);
        polynomials["L"]->mulScalar(ZTS2Y);

        FrElement dividendArr[2];
        E.fr.neg(dividendArr[0], challenges["y"]);
        dividendArr[1] = E.fr.one();
        Polynomial<Engine> *polDividend = new Polynomial<Engine>(E, dividendArr);
        Polynomial<Engine>* polRemainder = polynomials["L"]->divBy(*polDividend);

        // Check degrees
        if (polRemainder->getDegree() > 0) {
            ss.str("");
            ss << "Degree of L(X)/(ZTS2(y)(X-y)) remainder is " << polRemainder->getDegree() << " and should be 0";
            throw std::runtime_error(ss.str());
        }
        if (polynomials["L"]->getDegree() >= 9 * zkey->domainSize + 17) {
            throw std::runtime_error("Degree of L(X)/(ZTS2(y)(X-y)) is not correct");
        }

        // The fifth output of the prover is ([W2]_1), where W2:=(f/Z_t)(x)
        G1Point w2 = multiExponentiation(polynomials["L"], "W2");
        proof->addPolynomialCommitment("W2", w2);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeL() {
        buffers["L"] = new FrElement[zkey->domainSize * 16];

        FrElement evalR1Y = polynomials["R1"]->evaluate(challenges["y"]);
        FrElement evalR2Y = polynomials["R2"]->evaluate(challenges["y"]);
        FrElement evalZTY = polynomials["ZT"]->evaluate(challenges["y"]);

        FrElement preL1;
        FrElement preTmp;
        E.fr.sub(preL1, challenges["y"], roots["S2h2"][0]);
        E.fr.sub(preTmp, challenges["y"], roots["S2h2"][1]);
        E.fr.mul(preL1, preL1, preTmp);
        E.fr.sub(preTmp, challenges["y"], roots["S2h2"][2]);
        E.fr.mul(preL1, preL1, preTmp);
        E.fr.sub(preTmp, challenges["y"], roots["S2h3"][0]);
        E.fr.mul(preL1, preL1, preTmp);
        E.fr.sub(preTmp, challenges["y"], roots["S2h3"][1]);
        E.fr.mul(preL1, preL1, preTmp);
        E.fr.sub(preTmp, challenges["y"], roots["S2h3"][2]);
        E.fr.mul(preL1, preL1, preTmp);
        toInverse["yBatch"] = preL1;

        FrElement preL2;
        E.fr.sub(preL2, challenges["y"], roots["S1h1"][0]);
        E.fr.mul(preL2, preL2, challenges["alpha"]);
        E.fr.sub(preTmp, challenges["y"], roots["S1h1"][1]);
        E.fr.mul(preL2, preL2, preTmp);
        E.fr.sub(preTmp, challenges["y"], roots["S1h1"][2]);
        E.fr.mul(preL2, preL2, preTmp);
        E.fr.sub(preTmp, challenges["y"], roots["S1h1"][3]);
        E.fr.mul(preL2, preL2, preTmp);

        LOG_TRACE("> Computing F fft");

        evaluations["F"] = new Evaluations<Engine>(E, *polynomials["F"]);

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
            FrElement l1;
            E.fr.sub(l1, c1, evalR1Y);
            E.fr.mul(l1, preL1, l1);

            // l2 = alpha (y - h1) (y - h1w4) (y - h1w4_2) (y - h1w4_3) (C2(X) - R2(y))
            FrElement l2;
            E.fr.sub(l2, c2, evalR2Y);
            E.fr.mul(l2, preL2, l2);

            // l3 = ZT(y) (f(X)/ZT(X))
            // Recall f is already a f(X)/ZT(X)
            FrElement l3;
            E.fr.mul(l3, evalZTY, f);

            FrElement l;
            E.fr.add(l1, l1, l2);
            E.fr.sub(l, l1, l3);

            buffers["L"][i] = l;
        }

        LOG_TRACE("··· Computing L ifft");
        polynomials["L"] = new Polynomial<Engine>(E, buffers["L"], zkey->domainSize * 16);

        // Check degree
        if (polynomials["L"]->getDegree() >= 9 * zkey->domainSize + 18) {
            throw std::runtime_error("L Polynomial is not well calculated");
        }

        delete buffers["L"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZTS2() {
        FrElement arr[6] = {roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                            roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]};
        polynomials["ZTS2"] = Polynomial<Engine>::zerofierPolynomial(arr);
    }

    //TODO !!!!!!!!!!
    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getMontgomeryBatchedInverse() {
        return E.fr.zero();
    }

    template<typename Engine>
    typename Engine::G1Point FflonkProver<Engine>::multiExponentiation(const Polynomial<Engine> *polynomial, const std::string name) {
//        G1P value;
//        FrElements pol = polynomialFromMontgomery(polynomial, from, count);
//
//        E.g1.multiMulByScalar(value, ptau.data(), (uint8_t *) pol.data(), sizeof(pol[0]), pol.size());
//        return value;
    }
}
