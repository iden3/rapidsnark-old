#include "fflonk_prover.hpp"

#include "curve_utils.hpp"
#include "zkey.hpp"
#include "zkey_fflonk.hpp"
#include "wtns_utils.hpp"
#include <sodium.h>
#include "mul_z.hpp"
#include "keccak_256_transcript.hpp"
#include "logger.hpp"
#include "thread_utils.hpp"

using namespace CPlusPlusLogging;

namespace Fflonk {

    template<typename Engine>
    FflonkProver<Engine>::FflonkProver(Engine &_E) : E(_E) {
        curveName = CurveUtils::getCurveNameByEngine();
    }

    template<typename Engine>
    FflonkProver<Engine>::~FflonkProver() {
        if (NULL != fft) {
            delete fft;
        }
    }

    template<typename Engine>
    std::tuple <json, json> FflonkProver<Engine>::prove(BinFileUtils::BinFile *fdZkey, BinFileUtils::BinFile *fdWtns) {
        try {
            LOG_TRACE("FFLONK PROVER STARTED");

            LOG_TRACE("> Reading witness file");

            dump = new Dump::Dump<Engine>(E);

            this->fdZkey = fdZkey;
            this->fdWtns = fdWtns;
            auto beg = high_resolution_clock::now();
            auto wtns = WtnsUtils::loadHeader(fdWtns);
            processingTime.push_back(ProcessingTime("WtnsUtils load header", high_resolution_clock::now()));
            LOG_TRACE("> Reading zkey file");
            zkey = Zkey::FflonkZkeyHeader::loadFflonkZkeyHeader(fdZkey);
            processingTime.push_back(ProcessingTime("WtnloadFflonkZkeyHeader", high_resolution_clock::now()));

            if (zkey->protocolId != Zkey::FFLONK_PROTOCOL_ID) {
                throw std::invalid_argument("zkey file is not fflonk");
            }

            fft = new FFT<typename Engine::Fr>(zkey->domainSize * 16);
            zkeyPower = fft->log2(zkey->domainSize);
            processingTime.push_back(ProcessingTime("fft init", high_resolution_clock::now()));

            mulZ = new MulZ<Engine>(E, fft);

            if (mpz_cmp(zkey->rPrime, wtns->prime) != 0) {
                throw std::invalid_argument("Curve of the witness does not match the curve of the proving key");
            }

            //TODO compare zkey field with current field
//            if (mpz_cmp(zkeyHeader->rPrime, altBbn128r) != 0) {
//                throw std::invalid_argument( "zkey curve not supported" );
//            }

            if (wtns->nVars != zkey->nVars - zkey->nAdditions) {
                std::ostringstream ss;
                ss << "Invalid witness length. Circuit: " << zkey->nVars << ", witness: " << wtns->nVars << ", "
                   << zkey->nAdditions;
                throw std::invalid_argument(ss.str());
            }

            sDomain = zkey->domainSize * sizeof(FrElement);

            std::ostringstream ss;
            LOG_TRACE("----------------------------");
            LOG_TRACE("  FFLONK PROVE SETTINGS");
            ss.str("");
            ss << "  Curve:         " << curveName;
            LOG_TRACE(ss);
            ss.str("");
            ss << "  Circuit power: " << zkeyPower;
            LOG_TRACE(ss);
            ss.str("");
            ss << "  Domain size:   " << zkey->domainSize;
            LOG_TRACE(ss);
            ss.str("");
            ss << "  Vars:          " << zkey->nVars;
            LOG_TRACE(ss);
            ss.str("");
            ss << "  Public vars:   " << zkey->nPublic;
            LOG_TRACE(ss);
            ss.str("");
            ss << "  Constraints:   " << zkey->nConstraints;
            LOG_TRACE(ss);
            ss.str("");
            ss << "  Additions:     " << zkey->nAdditions;
            LOG_TRACE(ss);
            LOG_TRACE("----------------------------");

            //Read witness data
            LOG_TRACE("> Reading witness file data");
            buffWitness = (FrElement *) fdWtns->getSectionData(2);
            processingTime.push_back(ProcessingTime("Reading witness file data", high_resolution_clock::now()));

            // First element in plonk is not used and can be any value. (But always the same).
            // We set it to zero to go faster in the exponentiations.
            buffWitness[0] = E.fr.zero();

            buffInternalWitness = new FrElement[zkey->nAdditions];

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

            ss.str("");
            ss << "> Reading Section " << Zkey::ZKEY_FF_ADDITIONS_SECTION << ". Additions";
            LOG_TRACE(ss);
            processingTime.push_back(ProcessingTime("Read section Additions", high_resolution_clock::now()));

            calculateAdditions(fdZkey);
            processingTime.push_back(ProcessingTime("calculateAdditions", high_resolution_clock::now()));

            ss.str("");
            ss << "> Reading Section " << Zkey::ZKEY_FF_SIGMA1_SECTION << "," << Zkey::ZKEY_FF_SIGMA2_SECTION
               << "," << Zkey::ZKEY_FF_SIGMA3_SECTION << ". Sigma1, Sigma2 & Sigma 3";
            LOG_TRACE(ss);

            LOG_TRACE("··· Reading Sigma polynomials ");
            polynomials["Sigma1"] = new Polynomial<Engine>(E, zkey->domainSize);
            polynomials["Sigma2"] = new Polynomial<Engine>(E, zkey->domainSize);
            polynomials["Sigma3"] = new Polynomial<Engine>(E, zkey->domainSize);

            int nThreads = omp_get_max_threads() / 2;

            ThreadUtils::parcpy(polynomials["Sigma1"]->coef,
                                (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["Sigma2"]->coef,
                                (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["Sigma3"]->coef,
                                (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION),
                                sDomain, nThreads);

            processingTime.push_back(ProcessingTime("read Sigma pol", high_resolution_clock::now()));

            polynomials["Sigma1"]->fixDegree();
            polynomials["Sigma2"]->fixDegree();
            polynomials["Sigma3"]->fixDegree();

            LOG_TRACE("··· Reading Sigma evaluations ");
            evaluations["Sigma1"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
            evaluations["Sigma2"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
            evaluations["Sigma3"] = new Evaluations<Engine>(E, zkey->domainSize * 4);

            ThreadUtils::parcpy(evaluations["Sigma1"]->eval,
                                (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["Sigma2"]->eval,
                                (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["Sigma3"]->eval,
                                (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            processingTime.push_back(ProcessingTime("read Sigma eval", high_resolution_clock::now()));

            ss.str("");
            ss << "> Reading Section " << Zkey::ZKEY_FF_PTAU_SECTION << ". Powers of Tau";
            LOG_TRACE(ss);
            PTau = new G1PointAffine[zkey->domainSize * 16];
            memset(PTau, 0, sizeof(PTau));

            // domainSize * 9 + 18 = SRS length in the zkey saved in setup process.
            // it corresponds to the maximum SRS length needed, specifically to commit C2
            // notice that the reserved buffers size is zkey->domainSize * 16 * sG1 because a power of two buffer size is needed
            // the remaining buffer not filled from SRS are set to 0
            ThreadUtils::parcpy(this->PTau,
                                (G1PointAffine *) fdZkey->getSectionData(Zkey::ZKEY_FF_PTAU_SECTION),
                                (zkey->domainSize * 9 + 18) * sizeof(G1PointAffine), nThreads);
            processingTime.push_back(ProcessingTime("read PTau", high_resolution_clock::now()));

            // START FFLONK PROVER PROTOCOL

            // ROUND 1. Compute C1(X) polynomial
            LOG_TRACE("");
            LOG_TRACE("> ROUND 1");
            round1();
            processingTime.push_back(ProcessingTime("round1", high_resolution_clock::now()));

            delete polynomials["T0"];
            delete evaluations["QL"];
            delete evaluations["QR"];
            delete evaluations["QM"];
            delete evaluations["QO"];
            delete evaluations["QC"];

            // ROUND 2. Compute C2(X) polynomial
            LOG_TRACE("> ROUND 2");
            round2();
            processingTime.push_back(ProcessingTime("round2", high_resolution_clock::now()));

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
            processingTime.push_back(ProcessingTime("round3", high_resolution_clock::now()));

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
            processingTime.push_back(ProcessingTime("round4", high_resolution_clock::now()));

            delete polynomials["C0"];

            // ROUND 5. Compute W'(X) polynomial
            LOG_TRACE("> ROUND 5");
            round5();
            processingTime.push_back(ProcessingTime("round5", high_resolution_clock::now()));

            delete polynomials["C1"];
            delete polynomials["C2"];
            delete polynomials["R1"];
            delete polynomials["R2"];
            delete polynomials["F"];
            delete polynomials["L"];
            delete polynomials["ZT"];
            delete polynomials["ZTS2"];

            proof->addEvaluationCommitment("inv", getMontgomeryBatchedInverse());

            ss.str("\n");
            for (auto i = 0; i != processingTime.size(); i++) {
                long duration = i == 0 ? duration_cast<microseconds>(processingTime[i].duration - beg).count()
                                       : duration_cast<milliseconds>(
                                processingTime[i].duration - processingTime[i - 1].duration).count();
                ss << processingTime[i].label << ": " << duration << " ms\n";
            }
            LOG_TRACE(ss);

            // Prepare public inputs
            json publicSignals;
            for (u_int32_t i = 1; i <= zkey->nPublic; i++) {
                E.fr.toMontgomery(buffWitness[i], buffWitness[i]);
                publicSignals.push_back(E.fr.toString(buffWitness[i]).c_str());
            }

            LOG_TRACE("FFLONK PROVER FINISHED");

            return {proof->toJson(), publicSignals};
        }
        catch (const std::exception &e) {
            std::cerr << "EXCEPTION: " << e.what() << "\n";
            exit(EXIT_FAILURE);
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::calculateAdditions(BinFileUtils::BinFile *fdZkey) {
        LOG_TRACE("··· Computing additions");
        Zkey::Addition<Engine> *additionsBuff = (Zkey::Addition<Engine> *) fdZkey->getSectionData(
                Zkey::ZKEY_FF_ADDITIONS_SECTION);

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

        //0 index not used, set to zero
        //TODO check it fill all bytes with random values!!!!
        //randombytes_buf((void *) &blindingFactors[0], sizeof(blindingFactors));
        for (u_int32_t i = 0; i < zkey->nAdditions; i++) {
            blindingFactors[i] = E.fr.one();
        }
        processingTime.push_back(ProcessingTime("computeWirePolynomials inici", high_resolution_clock::now()));

        // STEP 1.2 - Compute wire polynomials a(X), b(X) and c(X)
        LOG_TRACE("> Computing A, B, C wire polynomials");
        computeWirePolynomials();
        processingTime.push_back(ProcessingTime("computeWirePolynomials", high_resolution_clock::now()));

        // STEP 1.3 - Compute the quotient polynomial T0(X)
        LOG_TRACE("> Computing T0 polynomial");
        computeT0();
        processingTime.push_back(ProcessingTime("computeT0", high_resolution_clock::now()));

        // STEP 1.4 - Compute the FFT-style combination polynomial C1(X)
        LOG_TRACE("> Computing C1 polynomial");
        computeC1();
        processingTime.push_back(ProcessingTime("compute C1", high_resolution_clock::now()));

        // The first output of the prover is ([C1]_1)
        LOG_TRACE("> Computing C1 multi exponentiation");
        G1Point C1 = multiExponentiation(polynomials["C1"]);
        proof->addPolynomialCommitment("C1", C1);
        dump->dump("[C1]_1", C1);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeWirePolynomials() {

        // Build A, B and C evaluations buffer from zkey and witness files
        buffers["A"] = new FrElement[zkey->domainSize];
        buffers["B"] = new FrElement[zkey->domainSize];
        buffers["C"] = new FrElement[zkey->domainSize];

        u_int32_t byteLength = sizeof(u_int32_t) * zkey->nConstraints;
        mapBuffers["A"] = new u_int32_t[zkey->nConstraints];
        mapBuffers["B"] = new u_int32_t[zkey->nConstraints];
        mapBuffers["C"] = new u_int32_t[zkey->nConstraints];

        memset(mapBuffers["A"], 0, byteLength);
        memset(mapBuffers["B"], 0, byteLength);
        memset(mapBuffers["C"], 0, byteLength);

        // Read zkey sections and fill the buffers
        int nThreads = omp_get_max_threads() / 2;
        ThreadUtils::parcpy(mapBuffers["A"],
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_A_MAP_SECTION),
                            byteLength, nThreads);
        ThreadUtils::parcpy(mapBuffers["B"],
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_B_MAP_SECTION),
                            byteLength, nThreads);
        ThreadUtils::parcpy(mapBuffers["C"],
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_C_MAP_SECTION),
                            byteLength, nThreads);

        FrElement bFactorsA[2] = {blindingFactors[2], blindingFactors[1]};
        FrElement bFactorsB[2] = {blindingFactors[4], blindingFactors[3]};
        FrElement bFactorsC[2] = {blindingFactors[6], blindingFactors[5]};

        computeWirePolynomial("A", bFactorsA);
        computeWirePolynomial("B", bFactorsB);
        computeWirePolynomial("C", bFactorsC);

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
            FrElement witness = getWitness(mapBuffers[polName][i]);
            E.fr.toMontgomery(buffers[polName][i], witness);
        }


        // Create the polynomial
        // and compute the coefficients of the wire polynomials from evaluations
        u_int32_t bFactorsLen = 2;
        std::ostringstream ss;
        ss << "··· Computing " << polName << " ifft";
        LOG_TRACE(ss);
        polynomials[polName] = Polynomial<Engine>::fromEvaluations(E, fft, buffers[polName], zkey->domainSize,
                                                                   bFactorsLen);

        // Compute the extended evaluations of the wire polynomials
        ss.str("");
        ss << "··· Computing " << polName << " fft";
        LOG_TRACE(ss);
        evaluations[polName] = new Evaluations<Engine>(E, fft, *polynomials[polName], zkey->domainSize * 4);

        // Blind polynomial coefficients with blinding scalars blindingFactors
        polynomials[polName]->blindCoefficients(blindingFactors, bFactorsLen);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeT0() {
        std::ostringstream ss;
        ss << "> Reading sections "
           << Zkey::ZKEY_FF_QL_SECTION << "," << Zkey::ZKEY_FF_QR_SECTION << ","
           << Zkey::ZKEY_FF_QM_SECTION << "," << Zkey::ZKEY_FF_QO_SECTION << "," << Zkey::ZKEY_FF_QC_SECTION
           << ". Q selectors";
        LOG_TRACE(ss);

        // Reserve memory for Q's evaluations
        evaluations["QL"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QR"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QM"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QO"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        evaluations["QC"] = new Evaluations<Engine>(E, zkey->domainSize * 4);

        // Read Q's evaluations from zkey file
        int nThreads = omp_get_max_threads() / 2;
        ThreadUtils::parcpy(evaluations["QL"]->eval,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION) + zkey->domainSize,
                            sDomain * 4, nThreads);
        ThreadUtils::parcpy(evaluations["QR"]->eval,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QR_SECTION) + zkey->domainSize,
                            sDomain * 4, nThreads);
        ThreadUtils::parcpy(evaluations["QM"]->eval,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QM_SECTION) + zkey->domainSize,
                            sDomain * 4, nThreads);
        ThreadUtils::parcpy(evaluations["QO"]->eval,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QO_SECTION) + zkey->domainSize,
                            sDomain * 4, nThreads);
        ThreadUtils::parcpy(evaluations["QC"]->eval,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QC_SECTION) + zkey->domainSize,
                            sDomain * 4, nThreads);

        // Read Lagrange polynomials & evaluations from zkey file
        evaluations["lagrange1"] = new Evaluations<Engine>(E, zkey->domainSize * 4);
        ThreadUtils::parcpy(evaluations["lagrange1"]->eval,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_LAGRANGE_SECTION) + zkey->domainSize,
                            sDomain * 4, nThreads);

        // Reserve memory for buffers T0 and T0z
        buffers["T0"] = new FrElement[zkey->domainSize * 4];
        buffers["T0z"] = new FrElement[zkey->domainSize * 4];

#pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->domainSize * 4; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "      t0 evaluation " << i << "/" << zkey->domainSize * 4;
//                LOG_TRACE(ss);
            }

            FrElement omega = fft->root(zkeyPower + 2, i);

            // Get related evaluations to compute current T0 evaluation
            FrElement a = evaluations["A"]->eval[i];
            FrElement b = evaluations["B"]->eval[i];
            FrElement c = evaluations["C"]->eval[i];

            FrElement ql = evaluations["QL"]->eval[i];
            FrElement qr = evaluations["QR"]->eval[i];
            FrElement qm = evaluations["QM"]->eval[i];
            FrElement qo = evaluations["QO"]->eval[i];
            FrElement qc = evaluations["QC"]->eval[i];

            // Compute blinding factors
            FrElement az = E.fr.add(E.fr.mul(blindingFactors[1], omega), blindingFactors[2]);
            FrElement bz = E.fr.add(E.fr.mul(blindingFactors[3], omega), blindingFactors[4]);
            FrElement cz = E.fr.add(E.fr.mul(blindingFactors[5], omega), blindingFactors[6]);

            // Compute current public input
            FrElement pi = E.fr.zero();
            for (u_int32_t j = 0; j < zkey->nPublic; j++) {
                u_int32_t offset = (j * 4 * zkey->domainSize) + i;
                FrElement lPol = evaluations["lagrange1"]->eval[offset];
                FrElement aVal = buffers["A"][j];

                pi = E.fr.sub(pi, E.fr.mul(lPol, aVal));
            }

            // T0(X) = [q_L(X)·a(X) + q_R(X)·b(X) + q_M(X)·a(X)·b(X) + q_O(X)·c(X) + q_C(X) + PI(X)] · 1/Z_H(X)
            // Compute first T0(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            // expression 1 -> q_L(X)·a(X)
            FrElement e1 = E.fr.mul(a, ql);
            FrElement e1z = E.fr.mul(az, ql);

            // expression 2 -> q_R(X)·b(X)
            FrElement e2 = E.fr.mul(b, qr);
            FrElement e2z = E.fr.mul(bz, qr);

            // expression 3 -> q_M(X)·a(X)·b(X)
            auto [e3, e3z] = mulZ->mul2(a, b, az, bz, i % 4);
            e3 = E.fr.mul(e3, qm);
            e3z = E.fr.mul(e3z, qm);

            // expression 4 -> q_O(X)·c(X)
            FrElement e4 = E.fr.mul(c, qo);
            FrElement e4z = E.fr.mul(cz, qo);

            // t0 = expressions 1 + expression 2 + expression 3 + expression 4 + qc + pi
            FrElement t0 = E.fr.add(e1, E.fr.add(e2, E.fr.add(e3, E.fr.add(e4, E.fr.add(qc, pi)))));
            FrElement t0z = E.fr.add(e1z, E.fr.add(e2z, E.fr.add(e3z, e4z)));

            buffers["T0"][i] = t0;
            buffers["T0z"][i] = t0z;
        }

        // Compute the coefficients of the polynomial T0(X) from buffers.T0
        LOG_TRACE("··· Computing T0 ifft");
        polynomials["T0"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T0"], zkey->domainSize * 4);

        // Divide the polynomial T0 by Z_H(X)
        processingTime.push_back(ProcessingTime("T0 divZh ini", high_resolution_clock::now()));
        LOG_TRACE("··· Computing T0 / ZH");
        polynomials["T0"]->divZh(zkey->domainSize);
        processingTime.push_back(ProcessingTime("T0 divZh fi", high_resolution_clock::now()));

        // Compute the coefficients of the polynomial T0z(X) from buffers.T0z
        LOG_TRACE("··· Computing T0z ifft");
        polynomials["T0z"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T0z"], zkey->domainSize * 4);

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
        u_int64_t maxDegree = std::max(degreeA * 4 + 1,
                                       std::max(degreeB * 4 + 2, std::max(degreeC * 4 + 3, degreeT0 * 4 + 3)));

        u_int64_t lengthBuffer = std::pow(2, fft->log2(maxDegree - 1) + 1);

        polynomials["C1"] = new Polynomial<Engine>(E, lengthBuffer);

        std::ostringstream ss;
#pragma omp parallel for
        for (u_int64_t i = 0; i < maxLength; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "      c1 coefficients " << i << "/" << maxLength;
                //LOG_TRACE(ss);
            }

            if (i <= degreeA) polynomials["C1"]->coef[i * 4] = polynomials["A"]->coef[i];
            if (i <= degreeB) polynomials["C1"]->coef[i * 4 + 1] = polynomials["B"]->coef[i];
            if (i <= degreeC) polynomials["C1"]->coef[i * 4 + 2] = polynomials["C"]->coef[i];
            if (i <= degreeT0) polynomials["C1"]->coef[i * 4 + 3] = polynomials["T0"]->coef[i];
        }

        polynomials["C1"]->fixDegree();

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
        LOG_TRACE("> Computing challenges beta and gamma");
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>(E);
        for (u_int32_t i = 0; i < zkey->nPublic; i++) {
            transcript->addScalar(buffers["A"][i]);
        }
        transcript->addPolCommitment(proof->getPolynomialCommitment("C1"));
        challenges["beta"] = transcript->getChallenge();
        std::ostringstream ss;
        ss << "··· challenges.beta: " << E.fr.toString(challenges["beta"]);
        LOG_TRACE(ss);

        // Compute permutation challenge gamma
        transcript->reset();
        transcript->addScalar(challenges["beta"]);
        challenges["gamma"] = transcript->getChallenge();
        ss.str("");
        ss << "··· challenges.gamma: " << E.fr.toString(challenges["gamma"]);
        LOG_TRACE(ss);

        // STEP 2.2 - Compute permutation polynomial z(X)
        LOG_TRACE("> Computing Z polynomial");
        computeZ();

        // STEP 2.3 - Compute quotient polynomial T1(X) and T2(X)
        LOG_TRACE("> Computing T1 polynomial");
        computeT1();
        LOG_TRACE("> Computing T2 polynomial");
        computeT2();

        // STEP 2.4 - Compute the FFT-style combination polynomial C2(X)
        LOG_TRACE("> Computing C2 polynomial");
        computeC2();

        // The second output of the prover is ([C2]_1)
        LOG_TRACE("> Computing C2 multi exponentiation");
        G1Point C2 = multiExponentiation(polynomials["C2"]);
        proof->addPolynomialCommitment("C2", C2);
        dump->dump("[C2]_1", C2);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZ() {
        FrElement *num = new FrElement[zkey->domainSize];
        FrElement *den = new FrElement[zkey->domainSize];
        FrElement *numArr = new FrElement[zkey->domainSize];
        FrElement *denArr = new FrElement[zkey->domainSize];
        buffers["Z"] = new FrElement[zkey->domainSize];

        LOG_TRACE("··· Computing Z evaluations");

        std::ostringstream ss;
#pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "> Computing Z evaluation " << i << "/" << zkey->domainSize;
                //LOG_TRACE(ss);
            }

            FrElement omega = fft->root(zkeyPower, i);

            // Z(X) := numArr / denArr
            // numArr := (a + beta·ω + gamma)(b + beta·ω·k1 + gamma)(c + beta·ω·k2 + gamma)
            FrElement betaw = E.fr.mul(challenges["beta"], omega);

            FrElement num1 = buffers["A"][i];
            num1 = E.fr.add(num1, betaw);
            num1 = E.fr.add(num1, challenges["gamma"]);

            FrElement num2 = buffers["B"][i];
            num2 = E.fr.add(num2, E.fr.mul(*((FrElement *) zkey->k1), betaw));
            num2 = E.fr.add(num2, challenges["gamma"]);

            FrElement num3 = buffers["C"][i];
            num3 = E.fr.add(num3, E.fr.mul(*((FrElement *) zkey->k2), betaw));
            num3 = E.fr.add(num3, challenges["gamma"]);

            num[i] = E.fr.mul(num1, E.fr.mul(num2, num3));

            // denArr := (a + beta·sigma1 + gamma)(b + beta·sigma2 + gamma)(c + beta·sigma3 + gamma)
            FrElement den1 = buffers["A"][i];
            den1 = E.fr.add(den1, E.fr.mul(challenges["beta"], evaluations["Sigma1"]->eval[i * 4]));
            den1 = E.fr.add(den1, challenges["gamma"]);

            FrElement den2 = buffers["B"][i];
            den2 = E.fr.add(den2, E.fr.mul(challenges["beta"], evaluations["Sigma2"]->eval[i * 4]));
            den2 = E.fr.add(den2, challenges["gamma"]);

            FrElement den3 = buffers["C"][i];
            den3 = E.fr.add(den3, E.fr.mul(challenges["beta"], evaluations["Sigma3"]->eval[i * 4]));
            den3 = E.fr.add(den3, challenges["gamma"]);

            den[i] = E.fr.mul(den1, E.fr.mul(den2, den3));
        }

        // Set the first values to 1
        numArr[0] = E.fr.one();
        denArr[0] = E.fr.one();

        for (u_int64_t i = 0; i < zkey->domainSize; i++) {
            // Multiply current num value with the previous one saved in numArr
            numArr[(i + 1) % zkey->domainSize] = E.fr.mul(numArr[i], num[i]);

            // Multiply current den value with the previous one saved in denArr
            denArr[(i + 1) % zkey->domainSize] = E.fr.mul(denArr[i], den[i]);
        }

        // Compute the inverse of denArr to compute in the next command the
        // division numArr/denArr by multiplying num · 1/denArr
        batchInverse(denArr, zkey->domainSize);

        // Multiply numArr · denArr where denArr was inverted in the previous command
#pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->domainSize; i++) {
            buffers["Z"][i] = E.fr.mul(numArr[i], denArr[i]);
        }

        if (!E.fr.eq(buffers["Z"][0], E.fr.one())) {
            throw std::runtime_error("Copy constraints does not match");
        }

        // Compute polynomial coefficients z(X) from buffers.Z
        LOG_TRACE("··· Computing Z ifft");
        polynomials["Z"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["Z"], zkey->domainSize, 3);

        // Compute extended evaluations of z(X) polynomial
        LOG_TRACE("··· Computing Z fft");
        evaluations["Z"] = new Evaluations<Engine>(E, fft, *polynomials["Z"], zkey->domainSize * 4);

        // Blind z(X) polynomial coefficients with blinding scalars b
        FrElement bFactors[3] = {blindingFactors[9], blindingFactors[8], blindingFactors[7]};
        polynomials["Z"]->blindCoefficients(bFactors, 3);

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

        LOG_TRACE("··· Computing T1 evaluations");

        std::ostringstream ss;
#pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize * 4; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "    T1 evaluation " << i << "/" << zkey->domainSize * 4;
                //LOG_TRACE(ss);
            }

            FrElement omega = fft->root(zkeyPower + 2, i);
            FrElement omega2 = E.fr.square(omega);

            FrElement z = evaluations["Z"]->eval[i];
            FrElement zp = E.fr.add(E.fr.add(
                    E.fr.mul(blindingFactors[7], omega2), E.fr.mul(blindingFactors[8], omega)), blindingFactors[9]);

            // T1(X) := (z(X) - 1) · L_1(X)
            // Compute first T1(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            FrElement lagrange1 = evaluations["lagrange1"]->eval[i];
            FrElement t1 = E.fr.mul(E.fr.sub(z, E.fr.one()), lagrange1);
            FrElement t1z = E.fr.mul(zp, lagrange1);

            buffers["T1"][i] = t1;
            buffers["T1z"][i] = t1z;

        }

        // Compute the coefficients of the polynomial T1(X) from buffers.T1
        LOG_TRACE("··· Computing T1 ifft");
        polynomials["T1"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T1"], zkey->domainSize * 4);

        // Divide the polynomial T1 by Z_H(X)
        processingTime.push_back(ProcessingTime("T1 divZh ini", high_resolution_clock::now()));
        polynomials["T1"]->divZh(zkey->domainSize);
        processingTime.push_back(ProcessingTime("T1 divZh fi", high_resolution_clock::now()));

        // Compute the coefficients of the polynomial T1z(X) from buffers.T1z
        LOG_TRACE("··· Computing T1z ifft");
        polynomials["T1z"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T1z"], zkey->domainSize * 4);

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

        LOG_TRACE("··· Computing T2 evaluations");

        std::ostringstream ss;
#pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize * 4; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "    T2 evaluation " << i << "/" << zkey->domainSize * 4;
                //LOG_TRACE(ss);
            }

            FrElement omega = fft->root(zkeyPower + 2, i);
            FrElement omega2 = E.fr.square(omega);
            FrElement omegaW = E.fr.mul(omega, fft->root(zkeyPower, 1));
            FrElement omegaW2 = E.fr.square(omegaW);

            FrElement a = evaluations["A"]->eval[i];
            FrElement b = evaluations["B"]->eval[i];
            FrElement c = evaluations["C"]->eval[i];
            FrElement z = evaluations["Z"]->eval[i];
            FrElement zW = evaluations["Z"]->eval[(zkey->domainSize * 4 + 4 + i) % (zkey->domainSize * 4)];

            FrElement ap = E.fr.add(E.fr.mul(blindingFactors[1], omega), blindingFactors[2]);
            FrElement bp = E.fr.add(E.fr.mul(blindingFactors[3], omega), blindingFactors[4]);
            FrElement cp = E.fr.add(E.fr.mul(blindingFactors[5], omega), blindingFactors[6]);
            FrElement zp = E.fr.add(E.fr.add(E.fr.mul(blindingFactors[7], omega2), E.fr.mul(blindingFactors[8], omega)),
                                    blindingFactors[9]);
            FrElement zWp = E.fr.add(
                    E.fr.add(E.fr.mul(blindingFactors[7], omegaW2), E.fr.mul(blindingFactors[8], omegaW)),
                    blindingFactors[9]);

            FrElement sigma1 = evaluations["Sigma1"]->eval[i];
            FrElement sigma2 = evaluations["Sigma2"]->eval[i];
            FrElement sigma3 = evaluations["Sigma3"]->eval[i];

            // T2(X) := [ (a(X) + beta·X + gamma)(b(X) + beta·k1·X + gamma)(c(X) + beta·k2·X + gamma)z(X)
            //           -(a(X) + beta·sigma1(X) + gamma)(b(X) + beta·sigma2(X) + gamma)(c(X) + beta·sigma3(X) + gamma)z(Xω)] · 1/Z_H(X)
            // Compute first T2(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)

            // expression 1 -> (a(X) + beta·X + gamma)(b(X) + beta·k1·X + gamma)(c(X) + beta·k2·X + gamma)z(X)
            FrElement betaX = E.fr.mul(challenges["beta"], omega);

            FrElement e11 = E.fr.add(a, betaX);
            e11 = E.fr.add(e11, challenges["gamma"]);

            FrElement e12 = E.fr.add(b, E.fr.mul(betaX, *((FrElement *) zkey->k1)));
            e12 = E.fr.add(e12, challenges["gamma"]);

            FrElement e13 = E.fr.add(c, E.fr.mul(betaX, *((FrElement *) zkey->k2)));
            e13 = E.fr.add(e13, challenges["gamma"]);

            auto [e1, e1z] = mulZ->mul4(e11, e12, e13, z, ap, bp, cp, zp, i % 4);

            // expression 2 -> (a(X) + beta·sigma1(X) + gamma)(b(X) + beta·sigma2(X) + gamma)(c(X) + beta·sigma3(X) + gamma)z(Xω)
            FrElement e21 = E.fr.add(a, E.fr.mul(challenges["beta"], sigma1));
            e21 = E.fr.add(e21, challenges["gamma"]);

            FrElement e22 = E.fr.add(b, E.fr.mul(challenges["beta"], sigma2));
            e22 = E.fr.add(e22, challenges["gamma"]);

            FrElement e23 = E.fr.add(c, E.fr.mul(challenges["beta"], sigma3));
            e23 = E.fr.add(e23, challenges["gamma"]);

            auto [e2, e2z] = mulZ->mul4(e21, e22, e23, zW, ap, bp, cp, zWp, i % 4);

            FrElement t2 = E.fr.sub(e1, e2);
            FrElement t2z = E.fr.sub(e1z, e2z);

            buffers["T2"][i] = t2;
            buffers["T2z"][i] = t2z;
        }

        // Compute the coefficients of the polynomial T2(X) from buffers.T2
        LOG_TRACE("··· Computing T2 ifft");
        polynomials["T2"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T2"], zkey->domainSize * 4);

        // Divide the polynomial T2 by Z_H(X)
        processingTime.push_back(ProcessingTime("T2 divZh ini", high_resolution_clock::now()));
        polynomials["T2"]->divZh(zkey->domainSize);
        processingTime.push_back(ProcessingTime("T2 divZh fi", high_resolution_clock::now()));

        // Compute the coefficients of the polynomial T2z(X) from buffers.T2z
        LOG_TRACE("··· Computing T2z ifft");
        polynomials["T2z"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T2z"], zkey->domainSize * 4);

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

        std::ostringstream ss;
#pragma omp parallel for
        for (u_int64_t i = 0; i < maxLength; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "    C2 coefficient " << i << "/" << maxLength;
                //LOG_TRACE(ss);
            }

            if (i <= degreeZ) polynomials["C2"]->coef[i * 3] = polynomials["Z"]->coef[i];
            if (i <= degreeT1) polynomials["C2"]->coef[i * 3 + 1] = polynomials["T1"]->coef[i];
            if (i <= degreeT2) polynomials["C2"]->coef[i * 3 + 2] = polynomials["T2"]->coef[i];
        }

        polynomials["C2"]->fixDegree();

        // Check degree
        if (polynomials["C2"]->getDegree() >= 9 * zkey->domainSize + 18) {
            throw std::runtime_error("C2 Polynomial is not well calculated");
        }
    }

    // ROUND 3
    template<typename Engine>
    void FflonkProver<Engine>::round3() {
        LOG_TRACE("> Computing challenge xi");
        // STEP 3.1 - Compute evaluation challenge xi ∈ S
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>(E);
        transcript->addPolCommitment(proof->getPolynomialCommitment("C2"));

        // Obtain a xi_seeder from the transcript
        // To force h1^4 = xi, h2^3 = xi and h_3^2 = xiω
        // we compute xi = xi_seeder^12, h1 = xi_seeder^3, h2 = xi_seeder^4 and h3 = xi_seeder^6
        FrElement xiSeed = transcript->getChallenge();
        FrElement xiSeed2;
        E.fr.square(xiSeed2, xiSeed);

        // Compute omega8, omega4 and omega3
        roots["w8"] = new FrElement[8];
        roots["w8"][0] = E.fr.one();
        for (uint i = 1; i < 8; i++) {
            roots["w8"][i] = E.fr.mul(roots["w8"][i - 1], *((FrElement *) zkey->w8));
        }

        // Compute omega3 and omega4
        roots["w4"] = new FrElement[4];
        roots["w4"][0] = E.fr.one();
        for (uint i = 1; i < 4; i++) {
            roots["w4"][i] = E.fr.mul(roots["w4"][i - 1], *((FrElement *) zkey->w4));
        }

        roots["w3"] = new FrElement[3];
        roots["w3"][0] = E.fr.one();
        roots["w3"][1] = *((FrElement *) zkey->w3);
        E.fr.square(roots["w3"][2], roots["w3"][1]);

        // Compute h0 = xiSeeder^3
        roots["S0h0"] = new FrElement[8];
        roots["S0h0"][0] = E.fr.mul(xiSeed2, xiSeed);
        for (uint i = 1; i < 8; i++) {
            roots["S0h0"][i] = E.fr.mul(roots["S0h0"][0], roots["w8"][i]);
        }

        // Compute h1 = xi_seeder^6
        roots["S1h1"] = new FrElement[4];
        roots["S1h1"][0] = E.fr.square(roots["S0h0"][0]);
        for (uint i = 1; i < 4; i++) {
            roots["S1h1"][i] = E.fr.mul(roots["S1h1"][0], roots["w4"][i]);
        }

        // Compute h2 = xi_seeder^8
        roots["S2h2"] = new FrElement[3];
        roots["S2h2"][0] = E.fr.mul(roots["S1h1"][0], xiSeed2);
        roots["S2h2"][1] = E.fr.mul(roots["S2h2"][0], roots["w3"][1]);
        roots["S2h2"][2] = E.fr.mul(roots["S2h2"][0], roots["w3"][2]);

        roots["S2h3"] = new FrElement[3];
        // Multiply h3 by third-root-omega to obtain h_3^3 = xiω
        // So, h3 = xi_seeder^8 ω^{1/3}
        roots["S2h3"][0] = E.fr.mul(roots["S2h2"][0], *((FrElement *) zkey->wr));
        roots["S2h3"][1] = E.fr.mul(roots["S2h3"][0], roots["w3"][1]);
        roots["S2h3"][2] = E.fr.mul(roots["S2h3"][0], roots["w3"][2]);

        // Compute xi = xi_seeder^24
        challenges["xi"] = E.fr.mul(E.fr.square(roots["S2h2"][0]), roots["S2h2"][0]);

        std::ostringstream ss;
        ss << "··· challenges.xi: " << E.fr.toString(challenges["xi"]);
        LOG_TRACE(ss);

        // Reserve memory for Q's polynomials
        polynomials["QL"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QR"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QM"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QO"] = new Polynomial<Engine>(E, zkey->domainSize);
        polynomials["QC"] = new Polynomial<Engine>(E, zkey->domainSize);

        // Read Q's evaluations from zkey file
        int nThreads = omp_get_max_threads() / 2;
        ThreadUtils::parcpy(polynomials["QL"]->coef,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION),
                            sDomain, nThreads);
        ThreadUtils::parcpy(polynomials["QR"]->coef,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QR_SECTION),
                            sDomain, nThreads);
        ThreadUtils::parcpy(polynomials["QM"]->coef,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QM_SECTION),
                            sDomain, nThreads);
        ThreadUtils::parcpy(polynomials["QO"]->coef,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QO_SECTION),
                            sDomain, nThreads);
        ThreadUtils::parcpy(polynomials["QC"]->coef,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_QC_SECTION),
                            sDomain, nThreads);

        polynomials["QL"]->fixDegree();
        polynomials["QR"]->fixDegree();
        polynomials["QM"]->fixDegree();
        polynomials["QO"]->fixDegree();
        polynomials["QC"]->fixDegree();

        // STEP 3.2 - Compute opening evaluations and add them to the proof (third output of the prover)
        LOG_TRACE("··· Computing evaluations");
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

        FrElement xiw = E.fr.mul(challenges["xi"], fft->root(zkeyPower, 1));
        proof->addEvaluationCommitment("zw", polynomials["Z"]->evaluate(xiw));
        proof->addEvaluationCommitment("t1w", polynomials["T1"]->evaluate(xiw));
        proof->addEvaluationCommitment("t2w", polynomials["T2"]->evaluate(xiw));
    }

    // ROUND 4
    template<typename Engine>
    void FflonkProver<Engine>::round4() {
        LOG_TRACE("> Computing challenge alpha");
        processingTime.push_back(ProcessingTime("inici round4()", high_resolution_clock::now()));

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
        processingTime.push_back(ProcessingTime("evaluations", high_resolution_clock::now()));

        std::ostringstream ss;
        ss << "··· challenges.alpha: " << E.fr.toString(challenges["alpha"]);
        LOG_TRACE(ss);

        // STEP 4.2 - Compute F(X)
        LOG_TRACE("> Reading C0 polynomial");
        polynomials["C0"] = new Polynomial<Engine>(E, zkey->domainSize * 8);
        int nThreads = omp_get_max_threads() / 2;
        ThreadUtils::parcpy(polynomials["C0"]->coef,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_C0_SECTION),
                            sDomain * 8, nThreads);

        LOG_TRACE("> Computing R0 polynomial");
        computeR0();
        processingTime.push_back(ProcessingTime("computeR0", high_resolution_clock::now()));

        LOG_TRACE("> Computing R1 polynomial");
        computeR1();
        processingTime.push_back(ProcessingTime("computeR1", high_resolution_clock::now()));

        LOG_TRACE("> Computing R2 polynomial");
        computeR2();
        processingTime.push_back(ProcessingTime("computeR2", high_resolution_clock::now()));

        LOG_TRACE("> Computing F polynomial");
        computeF();
        processingTime.push_back(ProcessingTime("computeF", high_resolution_clock::now()));
        LOG_TRACE("> Computing ZT polynomial");
        computeZT();
        processingTime.push_back(ProcessingTime("computeZT", high_resolution_clock::now()));

        LOG_TRACE("> Computing W = F / ZT polynomial");
        polynomials["F"]->divBy(*polynomials["ZT"]);
        processingTime.push_back(ProcessingTime("divBy ZT", high_resolution_clock::now()));

        // Check degrees
//        if (polRemainder->getDegree() > 0) {
//            std::ostringstream ss;
//            ss << "Degree of f(X)/ZT(X) remainder is " << polRemainder->getDegree() << " and should be 0";
//            throw std::runtime_error(ss.str());
//        }
        if (polynomials["F"]->getDegree() >= 9 * zkey->domainSize + 12) {
            throw std::runtime_error("Degree of f(X)/ZT(X) is not correct");
        }

        // The fourth output of the prover is ([W1]_1), where W1:=(f/Z_t)(x)
        LOG_TRACE("> Computing W1 multi exponentiation");
        G1Point W1 = multiExponentiation(polynomials["F"]);

        proof->addPolynomialCommitment("W1", W1);
        dump->dump("[W1]_1", W1);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeR0() {
        // COMPUTE R0
        // Compute the coefficients of R0(X) from 8 evaluations using lagrange interpolation. R0(X) ∈ F_{<8}[X]
        // We decide to use Lagrange interpolations because the R0 degree is very small (deg(R0)===7),
        // and we were not able to compute it using current ifft implementation because the omega are different
        FrElement xArr[8] = {roots["S0h0"][0], roots["S0h0"][1], roots["S0h0"][2], roots["S0h0"][3],
                             roots["S0h0"][4], roots["S0h0"][5], roots["S0h0"][6], roots["S0h0"][7]};
        FrElement yArr[8] = {polynomials["C0"]->evaluate(roots["S0h0"][0]),
                             polynomials["C0"]->evaluate(roots["S0h0"][1]),
                             polynomials["C0"]->evaluate(roots["S0h0"][2]),
                             polynomials["C0"]->evaluate(roots["S0h0"][3]),
                             polynomials["C0"]->evaluate(roots["S0h0"][4]),
                             polynomials["C0"]->evaluate(roots["S0h0"][5]),
                             polynomials["C0"]->evaluate(roots["S0h0"][6]),
                             polynomials["C0"]->evaluate(roots["S0h0"][7])};

        polynomials["R0"] = Polynomial<Engine>::lagrangePolynomialInterpolation(xArr, yArr, 8);

        // Check the degree of R0(X) < 8
        if (polynomials["R0"]->getDegree() > 7) {
            throw std::runtime_error("R0 Polynomial is not well calculated");
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeR1() {
        // COMPUTE R1
        // Compute the coefficients of R1(X) from 4 evaluations using lagrange interpolation. R1(X) ∈ F_{<4}[X]
        // We decide to use Lagrange interpolations because the R1 degree is very small (deg(R1)===3),
        // and we were not able to compute it using current ifft implementation because the omega are different
        FrElement xArr[4] = {roots["S1h1"][0], roots["S1h1"][1], roots["S1h1"][2], roots["S1h1"][3]};
        FrElement yArr[4] = {polynomials["C1"]->evaluate(roots["S1h1"][0]),
                             polynomials["C1"]->evaluate(roots["S1h1"][1]),
                             polynomials["C1"]->evaluate(roots["S1h1"][2]),
                             polynomials["C1"]->evaluate(roots["S1h1"][3])};

        polynomials["R1"] = Polynomial<Engine>::lagrangePolynomialInterpolation(xArr, yArr, 4);

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
        FrElement xArr[6] = {roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                             roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]};
        FrElement yArr[6] = {polynomials["C2"]->evaluate(roots["S2h2"][0]),
                             polynomials["C2"]->evaluate(roots["S2h2"][1]),
                             polynomials["C2"]->evaluate(roots["S2h2"][2]),
                             polynomials["C2"]->evaluate(roots["S2h3"][0]),
                             polynomials["C2"]->evaluate(roots["S2h3"][1]),
                             polynomials["C2"]->evaluate(roots["S2h3"][2])};

        polynomials["R2"] = Polynomial<Engine>::lagrangePolynomialInterpolation(xArr, yArr, 6);

        // Check the degree of r2(X) < 6
        if (polynomials["R2"]->getDegree() > 5) {
            throw std::runtime_error("R2 Polynomial is not well calculated");
        }
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeF() {
        buffers["F"] = new FrElement[zkey->domainSize * 16];

        LOG_TRACE("··· Reading C0 evaluations");
        evaluations["C0"] = new Evaluations<Engine>(E, zkey->domainSize * 16);
        int nThreads = omp_get_max_threads() / 2;
        ThreadUtils::parcpy(evaluations["C0"]->eval,
                            (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_C0_SECTION) + zkey->domainSize * 8,
                            sDomain * 16, nThreads);

        LOG_TRACE("··· Computing C1 fft");
        evaluations["C1"] = new Evaluations<Engine>(E, fft, *polynomials["C1"], zkey->domainSize * 16);
        LOG_TRACE("··· Computing C2 fft");
        evaluations["C2"] = new Evaluations<Engine>(E, fft, *polynomials["C2"], zkey->domainSize * 16);

        LOG_TRACE("··· Computing F evaluations");

        // COMPUTE F(X)
        std::ostringstream ss;
#pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize * 16; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "    F evaluation " << i << "/" << zkey->domainSize * 16;
                //LOG_TRACE(ss);
            }

            FrElement omega = fft->root(zkeyPower + 4, i);

            FrElement c0 = evaluations["C0"]->eval[i];
            FrElement c1 = evaluations["C1"]->eval[i];
            FrElement c2 = evaluations["C2"]->eval[i];
            FrElement r0 = polynomials["R0"]->evaluate(omega);
            FrElement r1 = polynomials["R1"]->evaluate(omega);
            FrElement r2 = polynomials["R2"]->evaluate(omega);

            // Compute the multiplier of all h0w8 terms
            FrElement preF0 = E.fr.sub(omega, roots["S0h0"][0]);
            for (uint i = 1; i < 8; i++) {
                preF0 = E.fr.mul(preF0, E.fr.sub(omega, roots["S0h0"][i]));
            }

            // Compute the multiplier of all h1w4 terms
            FrElement preF1 = E.fr.sub(omega, roots["S1h1"][0]);
            for (uint i = 1; i < 4; i++) {
                preF1 = E.fr.mul(preF1, E.fr.sub(omega, roots["S1h1"][i]));
            }

            // Compute the multiplier of all h2w3 and h2w3 terms
            FrElement preF2 = E.fr.sub(omega, roots["S2h2"][0]);
            for (uint i = 1; i < 3; i++) {
                preF2 = E.fr.mul(preF2, E.fr.sub(omega, roots["S2h2"][i]));
            }
            for (uint i = 0; i < 3; i++) {
                preF2 = E.fr.mul(preF2, E.fr.sub(omega, roots["S2h3"][i]));
            }

            // f0 = (X-h1) (X-h1w4) (X-h1w4_2) (X-h1w4_3) (X-h2) (X-h2w3) (X-h2w3_2) (X-h3) (X-h3w3) (X-h3w3_2) (C0(X) - R0(X))
            FrElement f0 = E.fr.mul(E.fr.mul(preF1, preF2), E.fr.sub(c0, r0));

            // f1 = alpha (X-h0) (X-h0w8) (X-h0w8_2) (X-h0w8_3) (X-h0w8_4) (X-h0w8_5) (X-h0w8_6) (X-h0w8_7) (X-h0w8_8)
            //            (X-h2) (X-h2w3) (X-h2w3_2) (X-h3) (X-h3w3) (X-h3w3_2) (C1(X) - R1(X))
            FrElement f1 = E.fr.mul(challenges["alpha"], E.fr.mul(preF0, preF2));
            f1 = E.fr.mul(f1, E.fr.sub(c1, r1));

            // f2 = alpha (X - h1) (X - h1w4) (X - h1w4_2) (X - h1w4_3) (C2(X) - R2(X))
            // f2 = alpha^2 (X-h0) (X-h0w8) (X-h0w8_2) (X-h0w8_3) (X-h0w8_4) (X-h0w8_5) (X-h0w8_6) (X-h0w8_7) (X-h0w8_8)
            //            (X-h1) (X-h1w4) (X-h1w4_2) (X-h1w4_3) (C2(X) - R2(X))
            FrElement f2 = E.fr.mul(E.fr.square(challenges["alpha"]), E.fr.mul(preF0, preF1));
            f2 = E.fr.mul(f2, E.fr.sub(c2, r2));

            FrElement f = E.fr.add(E.fr.add(f0, f1), f2);

            buffers["F"][i] = f;
        }

        LOG_TRACE("··· Computing F ifft");
        polynomials["F"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["F"], zkey->domainSize * 16);

        // Check degree
        if (polynomials["F"]->getDegree() >= 9 * zkey->domainSize + 30) {
            throw std::runtime_error("F Polynomial is not well calculated");
        }

        delete[] buffers["F"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZT() {
        FrElement arr[18] = {roots["S0h0"][0], roots["S0h0"][1], roots["S0h0"][2], roots["S0h0"][3],
                             roots["S0h0"][4], roots["S0h0"][5], roots["S0h0"][6], roots["S0h0"][7],
                             roots["S1h1"][0], roots["S1h1"][1], roots["S1h1"][2], roots["S1h1"][3],
                             roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                             roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]};

        polynomials["ZT"] = Polynomial<Engine>::zerofierPolynomial(arr, 18);
    }

    // ROUND 5
    template<typename Engine>
    void FflonkProver<Engine>::round5() {
        // STEP 5.1 - Compute random evaluation point y ∈ F
        // STEP 4.1 - Compute challenge alpha ∈ F
        LOG_TRACE("> Computing challenge y");
        Keccak256Transcript<Engine> *transcript = new Keccak256Transcript<Engine>(E);
        transcript->addPolCommitment(proof->getPolynomialCommitment("W1"));

        challenges["y"] = transcript->getChallenge();

        std::ostringstream ss;
        ss << "··· challenges.y: " << E.fr.toString(challenges["y"]);
        LOG_TRACE(ss);

        // STEP 5.2 - Compute L(X)
        LOG_TRACE("> Computing L polynomial");
        computeL();

        LOG_TRACE("> Computing ZTS2 polynomial");
        computeZTS2();

        FrElement ZTS2Y = polynomials["ZTS2"]->evaluate(challenges["y"]);
        E.fr.inv(ZTS2Y, ZTS2Y);
        polynomials["L"]->mulScalar(ZTS2Y);

        FrElement dividendArr[2];
        E.fr.neg(dividendArr[0], challenges["y"]);
        dividendArr[1] = E.fr.one();
        Polynomial<Engine> *polDividend = Polynomial<Engine>::fromCoefficients(E, dividendArr, 2);
        LOG_TRACE("> Computing W' = L / ZTS2 polynomial");
        Polynomial<Engine> *polRemainder = polynomials["L"]->divBy(*polDividend);

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
        LOG_TRACE("> Computing W' multi exponentiation");
        G1Point W2 = multiExponentiation(polynomials["L"]);
        proof->addPolynomialCommitment("W2", W2);
        dump->dump("[W2]_1", W2);
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeL() {
        buffers["L"] = new FrElement[zkey->domainSize * 16];

        FrElement evalR0Y = polynomials["R0"]->evaluate(challenges["y"]);
        FrElement evalR1Y = polynomials["R1"]->evaluate(challenges["y"]);
        FrElement evalR2Y = polynomials["R2"]->evaluate(challenges["y"]);
        FrElement evalZTY = polynomials["ZT"]->evaluate(challenges["y"]);

        FrElement mulL0 = E.fr.sub(challenges["y"], roots["S0h0"][0]);
        for (uint i = 1; i < 8; i++) {
            mulL0 = E.fr.mul(mulL0, E.fr.sub(challenges["y"], roots["S0h0"][i]));
        }

        FrElement mulL1 = E.fr.sub(challenges["y"], roots["S1h1"][0]);
        for (uint i = 1; i < 4; i++) {
            mulL1 = E.fr.mul(mulL1, E.fr.sub(challenges["y"], roots["S1h1"][i]));
        }

        FrElement mulL2 = E.fr.sub(challenges["y"], roots["S2h2"][0]);
        for (uint i = 1; i < 3; i++) {
            mulL2 = E.fr.mul(mulL2, E.fr.sub(challenges["y"], roots["S2h2"][i]));
        }
        for (uint i = 0; i < 3; i++) {
            mulL2 = E.fr.mul(mulL2, E.fr.sub(challenges["y"], roots["S2h3"][i]));
        }

        FrElement preL0 = E.fr.mul(mulL1, mulL2);
        FrElement preL1 = E.fr.mul(challenges["alpha"], E.fr.mul(mulL0, mulL2));
        FrElement preL2 = E.fr.mul(E.fr.square(challenges["alpha"]), E.fr.mul(mulL0, mulL1));

        toInverse["yBatch"] = preL1; //TODO review

        LOG_TRACE("··· Computing F fft");

        evaluations["F"] = new Evaluations<Engine>(E, fft, *polynomials["F"], zkey->domainSize * 16);

        LOG_TRACE("··· Computing L evaluations");

        // Set initial omega
        std::ostringstream ss;
#pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize * 16; i++) {
            if ((0 != i) && (i % 100000 == 0)) {
                ss.str("");
                ss << "    L evaluation " << i << "/" << zkey->domainSize * 16;
                //LOG_TRACE(ss);
            }

            FrElement omega = fft->root(zkeyPower + 4, i);

            FrElement c0 = evaluations["C0"]->eval[i];
            FrElement c1 = evaluations["C1"]->eval[i];
            FrElement c2 = evaluations["C2"]->eval[i];
            FrElement f = evaluations["F"]->eval[i];

            // l0 = (y-h1) (y-h1w4) (y-h1w4_2) (y-h1w4_3) (y-h2) (y-h2w3) (y-h2w3_2) (y-h3) (y-h3w3) (y-h3w3_2) (C0(X) - R0(X))
            FrElement l0 = E.fr.mul(preL0, E.fr.sub(c0, evalR0Y));

            // f1 = alpha (y-h0) (y-h0w8) (y-h0w8_2) (y-h0w8_3) (y-h0w8_4) (y-h0w8_5) (y-h0w8_6) (y-h0w8_7) (y-h0w8_8)
            //            (y-h2) (y-h2w3) (y-h2w3_2) (y-h3) (y-h3w3) (y-h3w3_2) (C1(X) - R1(X))
            FrElement l1 = E.fr.mul(preL1, E.fr.sub(c1, evalR1Y));

            // f2 = alpha^2 (y-h0) (y-h0w8) (y-h0w8_2) (y-h0w8_3) (y-h0w8_4) (y-h0w8_5) (y-h0w8_6) (y-h0w8_7) (y-h0w8_8)
            //            (y-h1) (y-h1w4) (y-h1w4_2) (y-h1w4_3) (C2(X) - R2(X))
            FrElement l2 = E.fr.mul(preL2, E.fr.sub(c2, evalR2Y));

            // l3 = ZT(y) (f(X)/ZT(X))
            // Recall f is already a f(X)/ZT(X)
            FrElement l3 = E.fr.mul(evalZTY, f);

            FrElement l = E.fr.sub(E.fr.add(E.fr.add(l0, l1), l2), l3);

            buffers["L"][i] = l;
        }

        LOG_TRACE("··· Computing L ifft");
        polynomials["L"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["L"], zkey->domainSize * 16);

        // Check degree
        if (polynomials["L"]->getDegree() >= 9 * zkey->domainSize + 18) {
            throw std::runtime_error("L Polynomial is not well calculated");
        }

        delete buffers["L"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::computeZTS2() {
        FrElement arr[10] = {roots["S1h1"][0], roots["S1h1"][1], roots["S1h1"][2], roots["S1h1"][3],
                             roots["S2h2"][0], roots["S2h2"][1], roots["S2h2"][2],
                             roots["S2h3"][0], roots["S2h3"][1], roots["S2h3"][2]};
        polynomials["ZTS2"] = Polynomial<Engine>::zerofierPolynomial(arr, 10);
    }

    template<typename Engine>
    void FflonkProver<Engine>::batchInverse(FrElement *elements, u_int64_t length) {
        // Calculate products: a, ab, abc, abcd, ...
        FrElement *products = new FrElement[length];

        products[0] = elements[0];
        for (int64_t index = 1; index < length; index++) {
            E.fr.mul(products[index], products[index - 1], elements[index]);
        }

        // Calculate inverses: 1/a, 1/ab, 1/abc, 1/abcd, ...
        FrElement *inverses = new FrElement[length];
        E.fr.inv(inverses[length - 1], products[length - 1]);
        for (int64_t index = length - 1; index > 0; index--) {
            E.fr.mul(inverses[index - 1], inverses[index], elements[index]);
        }

        elements[0] = inverses[0];
        for (int64_t index = 1; index < length; index++) {
            E.fr.mul(elements[index], inverses[index], products[index - 1]);
        }
    }

    template<typename Engine>
    typename Engine::FrElement *FflonkProver<Engine>::polynomialFromMontgomery(Polynomial<Engine> *polynomial) {
        const u_int64_t length = polynomial->getLength();

        FrElement *result = new FrElement[length];

#pragma omp parallel for
        for (u_int32_t index = 0; index < length; ++index) {
            E.fr.fromMontgomery(result[index], polynomial->coef[index]);
        }

        return result;
    }

    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getMontgomeryBatchedInverse() {
        std::ostringstream ss;
        //   · denominator needed in step 8 and 9 of the verifier to multiply by 1/Z_H(xi)
        FrElement xiN = challenges["xi"];
        for (u_int32_t i = 0; i < zkeyPower; i++) {
            xiN = E.fr.square(xiN);
        }
        toInverse["zh"] = E.fr.sub(xiN, E.fr.one());

        //   · denominator needed in step 10 and 11 of the verifier
        //     toInverse.yBatch -> Computed in round5, computeL()

        //   · denominator needed in the verifier when computing L_i^{S1}(X) and L_i^{S2}(X)
        for (uint i = 0; i < 4; i++) {
            ss.str("");
            ss << "LiS1_" << (i + 1);
            toInverse[ss.str()] = computeLiS1(i);
        }

        for (uint i = 0; i < 6; i++) {
            ss.str("");
            ss << "LiS2_" << (i + 1);
            toInverse[ss.str()] = computeLiS2(i);
        }

        //   · L_i i=1 to num public inputs, needed in step 6 and 7 of the verifier to compute L_1(xi) and PI(xi)
        u_int32_t size = std::max(1, (int) zkey->nPublic);

        FrElement w = E.fr.one();
        for (u_int32_t i = 0; i < size; i++) {
            ss.str("");
            ss << "Li_" << (i + 1);
            toInverse[ss.str()] = E.fr.mul(E.fr.set(zkey->domainSize), E.fr.sub(challenges["xi"], w));

            //w = E.fr.mul(w, zkey->w);
        }

        FrElement mulAccumulator = E.fr.one();
        for (auto &[key, value]: toInverse) {
            mulAccumulator = E.fr.mul(mulAccumulator, value);
        }

        E.fr.inv(mulAccumulator, mulAccumulator);
        return mulAccumulator;
    }

    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::computeLiS1(u_int32_t i) {
        // Compute L_i^{(S1)}(y)
        u_int32_t idx = i;
        FrElement den = E.fr.one();
        for (uint j = 0; j < 3; j++) {
            idx = (idx + 1) % 4;

            den = E.fr.mul(den, E.fr.sub(roots["S1h1"][i], roots["S1h1"][idx]));
        }
        return den;
    }

    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::computeLiS2(u_int32_t i) {
        // Compute L_i^{(S1)}(y)
        u_int32_t idx = i;
        FrElement den = E.fr.one();
        for (uint j = 0; j < 5; j++) {
            idx = (idx + 1) % 6;

            FrElement root1 = i < 3 ? roots["S2h2"][i] : roots["S2h3"][i - 3];
            FrElement root2 = idx < 3 ? roots["S2h2"][idx] : roots["S2h3"][idx - 3];
            den = E.fr.mul(den, E.fr.sub(root1, root2));
        }
        return den;
    }


    template<typename Engine>
    typename Engine::G1Point FflonkProver<Engine>::multiExponentiation(Polynomial<Engine> *polynomial) {
        G1Point value;
        FrElement *pol = this->polynomialFromMontgomery(polynomial);

        E.g1.multiMulByScalar(value, PTau, (uint8_t *) pol, sizeof(pol[0]), polynomial->getLength());

        delete[] pol;

        return value;
    }

    template<typename Engine>
    void FflonkProver<Engine>::printPol(std::string name, const Polynomial<Engine> *polynomial) {
        dump->dump(name, polynomial->coef[0]);
        dump->dump(name, polynomial->coef[polynomial->getDegree()]);
        dump->dump(name, polynomial->coef[polynomial->getLength() - 1]);
    }
}
