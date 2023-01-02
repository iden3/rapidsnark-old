#include "fflonk_prover.hpp"

#include "logger.hpp"
#include "curve_utils.hpp"
#include "zkey.hpp"
#include "zkey_fflonk.hpp"
#include "wtns_utils.hpp"

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
            auto wtns = WtnsUtils::loadHeader(fdWtns);

            LOG_TRACE("> Reading zkey file");
            zkey = Zkey::FflonkZkeyHeader::loadFflonkZkeyHeader(fdZkey);

            if (zkey->protocolId != Zkey::FFLONK_PROTOCOL_ID) {
                throw std::invalid_argument("zkey file is not fflonk");
            }

            fft = new FFT<typename Engine::Fr>(zkey->domainSize * 9 + 18);
            zkeyPower = fft->log2(zkey->domainSize);

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

            size_t sDomain = zkey->domainSize * sizeof(FrElement);

            std::ostringstream ss;
            ss << "----------------------------\n";
            ss << "  FFLONK PROVE SETTINGS\n";
            ss << "  Curve:         " << curveName << "\n";
            ss << "  Circuit power: " << zkey->power << "\n";
            ss << "  Domain size:   " << zkey->domainSize << "\n";
            ss << "  Vars:          " << zkey->nVars << "\n";
            ss << "  Public vars:   " << zkey->nPublic << "\n";
            ss << "  Constraints:   " << zkey->nConstraints << "\n";
            ss << "  Additions:     " << zkey->nAdditions << "\n";
            ss << "----------------------------\n";
            LOG_TRACE(ss);

            //Read witness data
            LOG_TRACE("> Reading witness file data");
            FrElement *buffWitness = (FrElement *) fdWtns->getSectionData(2);

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

            memcpy(polynomials["Sigma1"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION), sDomain);
            memcpy(polynomials["Sigma2"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION), sDomain);
            memcpy(polynomials["Sigma3"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION), sDomain);

            LOG_TRACE("··· Reading Sigma evaluations ");
            evaluations["Sigma1"] = new Evaluation<Engine>[zkey->domainSize];
            evaluations["Sigma2"] = new Evaluation<Engine>[zkey->domainSize];
            evaluations["Sigma3"] = new Evaluation<Engine>[zkey->domainSize];

            memcpy(polynomials["Sigma1"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION) + sDomain, sDomain * 4);
            memcpy(polynomials["Sigma2"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION) + sDomain, sDomain * 4);
            memcpy(polynomials["Sigma3"], (FrElement *) fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION) + sDomain, sDomain * 4);

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
    void FflonkProver<Engine>::calculateAdditions(BinFileUtils::BinFile *fdZkey) {
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

    //TODO !!!!!!!!!!
    template<typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getMontgomeryBatchedInverse() {
        return E.fr.zero;
    }

}