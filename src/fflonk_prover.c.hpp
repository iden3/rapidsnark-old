#include "fflonk_prover.hpp"

#include "curve_utils.hpp"
#include "zkey_fflonk.hpp"
#include "wtns_utils.hpp"
#include <sodium.h>
#include "thread_utils.hpp"
#include "polynomial/cpolynomial.hpp"

#define ELPP_NO_DEFAULT_LOG_FILE
#include "logger.hpp"
using namespace CPlusPlusLogging;

namespace Fflonk
{
    template <typename Engine>
    void FflonkProver<Engine>::initialize(void* reservedMemoryPtr, uint64_t reservedMemorySize)
    {
        zkey = NULL;
        this->reservedMemoryPtr = (FrElement *)reservedMemoryPtr;
        this->reservedMemorySize = reservedMemorySize;

        curveName = CurveUtils::getCurveNameByEngine();

        Logger::getInstance()->enableConsoleLogging();
        Logger::getInstance()->updateLogLevel(LOG_LEVEL_DEBUG);
    }

    template <typename Engine>
    FflonkProver<Engine>::FflonkProver(Engine &_E) : E(_E)
    {
        initialize(NULL);
    }

    template <typename Engine>
    FflonkProver<Engine>::FflonkProver(Engine &_E, void* reservedMemoryPtr, uint64_t reservedMemorySize) : E(_E)
    {
        initialize(reservedMemoryPtr, reservedMemorySize);
    }

    template <typename Engine>
    FflonkProver<Engine>::~FflonkProver() {
        this->removePrecomputedData();

        delete transcript;
        delete proof;
    }

    template<typename Engine>
    void FflonkProver<Engine>::removePrecomputedData() {
        // DELETE RESERVED MEMORY (if necessary)
        delete[] precomputedBigBuffer;
        delete[] mapBuffersBigBuffer;
        delete[] buffInternalWitness;

        if(NULL == reservedMemoryPtr) {
            delete[] inverses;
            delete[] products;
            delete[] nonPrecomputedBigBuffer;
        }

        delete fft;

        mapBuffers.clear();

        for (auto const &x : roots) delete[] x.second;
        roots.clear();

        for (auto const &x : omegas) delete[] x.second;
        omegas.clear();

        delete polynomials["QL"];
        delete polynomials["QR"];
        delete polynomials["QM"];
        delete polynomials["QO"];
        delete polynomials["QC"];
        delete polynomials["Sigma1"];
        delete polynomials["Sigma2"];
        delete polynomials["Sigma3"];
        delete polynomials["C0"];

        delete evaluations["QL"];
        delete evaluations["QR"];
        delete evaluations["QM"];
        delete evaluations["QO"];
        delete evaluations["QC"];
        delete evaluations["Sigma1"];
        delete evaluations["Sigma2"];
        delete evaluations["Sigma3"];
        delete evaluations["lagrange"];
    }

    template<typename Engine>
    void FflonkProver<Engine>::setZkey(BinFileUtils::BinFile *fdZkey) {
        try
        {
            if(NULL != zkey) {
                removePrecomputedData();
            }

            LOG_TRACE("> Reading zkey file");

            zkey = Zkey::FflonkZkeyHeader::loadFflonkZkeyHeader(fdZkey);

            if (zkey->protocolId != Zkey::FFLONK_PROTOCOL_ID)
            {
                throw std::invalid_argument("zkey file is not fflonk");
            }

            LOG_TRACE("> Starting fft");

            fft = new FFT<typename Engine::Fr>(zkey->domainSize * 4);
            zkeyPower = fft->log2(zkey->domainSize);

            mpz_t altBbn128r;
            mpz_init(altBbn128r);
            mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

            if (mpz_cmp(zkey->rPrime, altBbn128r) != 0)
            {
                throw std::invalid_argument("zkey curve not supported");
            }

            sDomain = zkey->domainSize * sizeof(FrElement);

            ////////////////////////////////////////////////////
            // PRECOMPUTED BIG BUFFER
            ////////////////////////////////////////////////////
            // Precomputed 1 > polynomials buffer
            uint64_t lengthPrecomputedBigBuffer = 0;
            lengthPrecomputedBigBuffer += zkey->domainSize * 1 * 8; // Polynomials QL, QR, QM, QO, QC, Sigma1, Sigma2 & Sigma3
            lengthPrecomputedBigBuffer += zkey->domainSize * 8 * 1; // Polynomial  C0
            // Precomputed 2 > evaluations buffer
            lengthPrecomputedBigBuffer += zkey->domainSize * 4 * 8; // Evaluations QL, QR, QM, QO, QC, Sigma1, Sigma2, Sigma3
            lengthPrecomputedBigBuffer += zkey->domainSize * 4 * zkey->nPublic; // Evaluations Lagrange1
            // Precomputed 3 > ptau buffer
            lengthPrecomputedBigBuffer += zkey->domainSize * 9 * sizeof(G1PointAffine) / sizeof(FrElement); // PTau buffer

            precomputedBigBuffer = new FrElement[lengthPrecomputedBigBuffer];

            polPtr["Sigma1"] = &precomputedBigBuffer[0];
            polPtr["Sigma2"] = polPtr["Sigma1"] + zkey->domainSize;
            polPtr["Sigma3"] = polPtr["Sigma2"] + zkey->domainSize;
            polPtr["QL"]     = polPtr["Sigma3"] + zkey->domainSize;
            polPtr["QR"]     = polPtr["QL"] + zkey->domainSize;
            polPtr["QM"]     = polPtr["QR"] + zkey->domainSize;
            polPtr["QO"]     = polPtr["QM"] + zkey->domainSize;
            polPtr["QC"]     = polPtr["QO"] + zkey->domainSize;
            polPtr["C0"]     = polPtr["QC"] + zkey->domainSize;

            evalPtr["Sigma1"] = polPtr["C0"] + zkey->domainSize * 8;
            evalPtr["Sigma2"] = evalPtr["Sigma1"] + zkey->domainSize * 4;
            evalPtr["Sigma3"] = evalPtr["Sigma2"] + zkey->domainSize * 4;
            evalPtr["QL"]     = evalPtr["Sigma3"] + zkey->domainSize * 4;
            evalPtr["QR"]     = evalPtr["QL"] + zkey->domainSize * 4;
            evalPtr["QM"]     = evalPtr["QR"] + zkey->domainSize * 4;
            evalPtr["QO"]     = evalPtr["QM"] + zkey->domainSize * 4;
            evalPtr["QC"]     = evalPtr["QO"] + zkey->domainSize * 4;
            evalPtr["lagrange"] = evalPtr["QC"] + zkey->domainSize * 4;

            PTau = (G1PointAffine *)(evalPtr["lagrange"] + zkey->domainSize * 4 * zkey->nPublic);

            // Read Q selectors polynomials and evaluations
            LOG_TRACE("... Loading QL, QR, QM, QO, & QC polynomial coefficients and evaluations");

            // Reserve memory for Q's polynomials
            polynomials["QL"] = new Polynomial<Engine>(E, polPtr["QL"], zkey->domainSize);
            polynomials["QR"] = new Polynomial<Engine>(E, polPtr["QR"], zkey->domainSize);
            polynomials["QM"] = new Polynomial<Engine>(E, polPtr["QM"], zkey->domainSize);
            polynomials["QO"] = new Polynomial<Engine>(E, polPtr["QO"], zkey->domainSize);
            polynomials["QC"] = new Polynomial<Engine>(E, polPtr["QC"], zkey->domainSize);

            int nThreads = omp_get_max_threads() / 2;

            // Read Q's polynomial coefficients from zkey file
            ThreadUtils::parcpy(polynomials["QL"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["QR"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QR_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["QM"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QM_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["QO"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QO_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["QC"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QC_SECTION),
                                sDomain, nThreads);

            polynomials["QL"]->fixDegree();
            polynomials["QR"]->fixDegree();
            polynomials["QM"]->fixDegree();
            polynomials["QO"]->fixDegree();
            polynomials["QC"]->fixDegree();

            std::ostringstream ss;
            ss << "... Reading Q selector evaluations ";

            // Reserve memory for Q's evaluations
            evaluations["QL"] = new Evaluations<Engine>(E, evalPtr["QL"], zkey->domainSize * 4);
            evaluations["QR"] = new Evaluations<Engine>(E, evalPtr["QR"], zkey->domainSize * 4);
            evaluations["QM"] = new Evaluations<Engine>(E, evalPtr["QM"], zkey->domainSize * 4);
            evaluations["QO"] = new Evaluations<Engine>(E, evalPtr["QO"], zkey->domainSize * 4);
            evaluations["QC"] = new Evaluations<Engine>(E, evalPtr["QC"], zkey->domainSize * 4);

            // Read Q's evaluations from zkey file
            ThreadUtils::parcpy(evaluations["QL"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QL_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["QR"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QR_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["QM"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QM_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["QO"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QO_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["QC"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_QC_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);

            // Read Sigma polynomial coefficients and evaluations from zkey file
            LOG_TRACE("... Loading Sigma1, Sigma2 & Sigma3 polynomial coefficients and evaluations");

            polynomials["Sigma1"] = new Polynomial<Engine>(E, polPtr["Sigma1"], zkey->domainSize);
            polynomials["Sigma2"] = new Polynomial<Engine>(E, polPtr["Sigma2"], zkey->domainSize);
            polynomials["Sigma3"] = new Polynomial<Engine>(E, polPtr["Sigma3"], zkey->domainSize);

            ThreadUtils::parcpy(polynomials["Sigma1"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["Sigma2"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION),
                                sDomain, nThreads);
            ThreadUtils::parcpy(polynomials["Sigma3"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION),
                                sDomain, nThreads);

            polynomials["Sigma1"]->fixDegree();
            polynomials["Sigma2"]->fixDegree();
            polynomials["Sigma3"]->fixDegree();

            evaluations["Sigma1"] = new Evaluations<Engine>(E, evalPtr["Sigma1"], zkey->domainSize * 4);
            evaluations["Sigma2"] = new Evaluations<Engine>(E, evalPtr["Sigma2"], zkey->domainSize * 4);
            evaluations["Sigma3"] = new Evaluations<Engine>(E, evalPtr["Sigma3"], zkey->domainSize * 4);

            ThreadUtils::parcpy(evaluations["Sigma1"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA1_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["Sigma2"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA2_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);
            ThreadUtils::parcpy(evaluations["Sigma3"]->eval,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_SIGMA3_SECTION) + zkey->domainSize,
                                sDomain * 4, nThreads);

            LOG_TRACE("... Loading C0 polynomial coefficients");
            polynomials["C0"] = new Polynomial<Engine>(E, polPtr["C0"], zkey->domainSize * 8);
            ThreadUtils::parcpy(polynomials["C0"]->coef,
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_C0_SECTION),
                                sDomain * 8, nThreads);
            polynomials["C0"]->fixDegree();

            // Read Lagrange polynomials & evaluations from zkey file
            LOG_TRACE("... Loading Lagrange evaluations");
            evaluations["lagrange"] = new Evaluations<Engine>(E, evalPtr["lagrange"], zkey->domainSize * 4 * zkey->nPublic);
            for(uint64_t i = 0 ; i < zkey->nPublic ; i++) {
                ThreadUtils::parcpy(evaluations["lagrange"]->eval + zkey->domainSize * 4 * i,
                                    (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_LAGRANGE_SECTION) + zkey->domainSize + zkey->domainSize * 5 * i,
                                    sDomain * 4, nThreads);
            }
            LOG_TRACE("... Loading Powers of Tau evaluations");

            ThreadUtils::parset(PTau, 0, sizeof(G1PointAffine) * zkey->domainSize * 9, nThreads);

            // domainSize * 9 = SRS length in the zkey saved in setup process.
            // it corresponds to the maximum SRS length needed, specifically to commit C2
            ThreadUtils::parcpy(this->PTau,
                                (G1PointAffine *)fdZkey->getSectionData(Zkey::ZKEY_FF_PTAU_SECTION),
                                (zkey->domainSize * 9) * sizeof(G1PointAffine), nThreads);

            // Load A, B & C map buffers
            LOG_TRACE("... Loading A, B & C map buffers");

            u_int64_t byteLength = sizeof(u_int32_t) * zkey->nConstraints;

            mapBuffersBigBuffer = new u_int32_t[zkey->nConstraints * 3];

            mapBuffers["A"] = mapBuffersBigBuffer;
            mapBuffers["B"] = mapBuffers["A"] + zkey->nConstraints;
            mapBuffers["C"] = mapBuffers["B"] + zkey->nConstraints;

            buffInternalWitness = new FrElement[zkey->nAdditions];

            LOG_TRACE("··· Loading additions");
            additionsBuff = (Zkey::Addition<Engine> *)fdZkey->getSectionData(Zkey::ZKEY_FF_ADDITIONS_SECTION);

            LOG_TRACE("··· Loading map buffers");
            ThreadUtils::parset(mapBuffers["A"], 0, byteLength * 3, nThreads);

            // Read zkey sections and fill the buffers
            ThreadUtils::parcpy(mapBuffers["A"],
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_A_MAP_SECTION),
                                byteLength, nThreads);
            ThreadUtils::parcpy(mapBuffers["B"],
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_B_MAP_SECTION),
                                byteLength, nThreads);
            ThreadUtils::parcpy(mapBuffers["C"],
                                (FrElement *)fdZkey->getSectionData(Zkey::ZKEY_FF_C_MAP_SECTION),
                                byteLength, nThreads);

            transcript = new Keccak256Transcript<Engine>(E);
            proof = new SnarkProof<Engine>(E, "fflonk");

            omegas["w8"] = new FrElement[8];
            omegas["w4"] = new FrElement[4];
            omegas["w3"] = new FrElement[3];
            roots["S0h0"] = new FrElement[8];
            roots["S1h1"] = new FrElement[4];
            roots["S2h2"] = new FrElement[6];

            lengthBatchInversesBuffer = zkey->domainSize * 2;

            if(NULL == this->reservedMemoryPtr) {
                inverses = new FrElement[zkey->domainSize];
                products = new FrElement[zkey->domainSize];
            } else {
                inverses = this->reservedMemoryPtr;
                products = inverses + zkey->domainSize;
            }

            ////////////////////////////////////////////////////
            // NON-PRECOMPUTED BIG BUFFER
            ////////////////////////////////////////////////////
            // Non-precomputed 1 > polynomials buffer
            lengthNonPrecomputedBigBuffer = 0;
            lengthNonPrecomputedBigBuffer += zkey->domainSize * 16 * 1; // Polynomial L (A, B & C will (re)use this buffer)
            lengthNonPrecomputedBigBuffer += zkey->domainSize * 8  * 1; // Polynomial C1
            lengthNonPrecomputedBigBuffer += zkey->domainSize * 16 * 1; // Polynomial C2
            lengthNonPrecomputedBigBuffer += zkey->domainSize * 16 * 1; // Polynomial F
            lengthNonPrecomputedBigBuffer += zkey->domainSize * 16 * 1; // Polynomial tmp (Z, T0, T1, T1z, T2 & T2z will (re)use this buffer)
            // Non-precomputed 2 > evaluations buffer
            lengthNonPrecomputedBigBuffer += zkey->domainSize * 4  * 3; // Evaluations A, B & C
            lengthNonPrecomputedBigBuffer += zkey->domainSize * 4  * 1; // Evaluations Z
            // Non-precomputed 3 > buffers buffer
            buffersLength = 0;
            buffersLength   += zkey->domainSize * 4  * 3; // Evaluations A, B & C
            buffersLength   += zkey->domainSize * 16 * 1; // Evaluations tmp (Z, numArr, denArr, T0, T1, T1z, T2 & T2z will (re)use this buffer)
            lengthNonPrecomputedBigBuffer += buffersLength;

            if(NULL == this->reservedMemoryPtr) {
                nonPrecomputedBigBuffer = new FrElement[lengthNonPrecomputedBigBuffer];
            } else {
                if((lengthBatchInversesBuffer + lengthNonPrecomputedBigBuffer) * sizeof(FrElement) > reservedMemorySize) {
                    ss.str("");
                    ss << "Not enough reserved memory to generate a prove. Increase reserved memory size at least to "
                        << (lengthBatchInversesBuffer + lengthNonPrecomputedBigBuffer) * sizeof(FrElement) << " bytes";
                    throw std::runtime_error(ss.str());
                }

                nonPrecomputedBigBuffer = this->reservedMemoryPtr + lengthBatchInversesBuffer;
            }
            
            polPtr["L"] = &nonPrecomputedBigBuffer[0];
            polPtr["C1"]  = polPtr["L"]  + zkey->domainSize * 16;
            polPtr["C2"]  = polPtr["C1"] + zkey->domainSize * 8;
            polPtr["F"]   = polPtr["C2"] + zkey->domainSize * 16;
            polPtr["tmp"] = polPtr["F"]  + zkey->domainSize * 16;
            // Reuses
            polPtr["A"]   = polPtr["L"];
            polPtr["B"]   = polPtr["A"]  + zkey->domainSize;
            polPtr["C"]   = polPtr["B"]  + zkey->domainSize;
            polPtr["Z"]   = polPtr["tmp"];
            polPtr["T0"]  = polPtr["Z"]  + zkey->domainSize * 2;
            polPtr["T1"]  = polPtr["T0"] + zkey->domainSize * 4;
            polPtr["T1z"] = polPtr["T1"] + zkey->domainSize * 2;
            polPtr["T2"]  = polPtr["T1"] + zkey->domainSize * 2;
            polPtr["T2z"] = polPtr["T2"] + zkey->domainSize * 4;

            evalPtr["A"] = polPtr["F"];
            evalPtr["B"] = evalPtr["A"]  + zkey->domainSize * 4;
            evalPtr["C"] = evalPtr["B"]  + zkey->domainSize * 4;
            evalPtr["Z"] = evalPtr["C"]  + zkey->domainSize * 4;

            buffers["A"]   = polPtr["tmp"]  + zkey->domainSize * 16;
            buffers["B"]   = buffers["A"]  + zkey->domainSize;
            buffers["C"]   = buffers["B"]  + zkey->domainSize;
            buffers["tmp"] = buffers["C"]  + zkey->domainSize;

            // Reuses
            buffers["Z"]      = buffers["tmp"];
            buffers["numArr"] = buffers["tmp"];
            buffers["denArr"] = buffers["numArr"] + zkey->domainSize;
            buffers["T0"]     = buffers["tmp"];
            buffers["T1"]     = buffers["tmp"];
            buffers["T1z"]    = buffers["tmp"] + zkey->domainSize * 2;
            buffers["T2"]     = buffers["tmp"];
            buffers["T2z"]    = buffers["tmp"] + zkey->domainSize * 4;
        }
        catch (const std::exception &e)
        {
            std::cerr << "EXCEPTION: " << e.what() << "\n";
            exit(EXIT_FAILURE);
        }
    }

    template<typename Engine>
    std::tuple <json, json> FflonkProver<Engine>::prove(BinFileUtils::BinFile *fdZkey, BinFileUtils::BinFile *fdWtns) {

        this->setZkey(fdZkey);
        return this->prove(fdWtns);
    }

    template <typename Engine>
    std::tuple<json, json> FflonkProver<Engine>::prove(BinFileUtils::BinFile *fdZkey, FrElement *buffWitness, WtnsUtils::Header* wtnsHeader)
    {
        this->setZkey(fdZkey);
        return this->prove(buffWitness, wtnsHeader);
    }

    template<typename Engine>
    std::tuple <json, json> FflonkProver<Engine>::prove(BinFileUtils::BinFile *fdWtns) {
        LOG_TRACE("> Reading witness file header");
        auto wtnsHeader = WtnsUtils::loadHeader(fdWtns);

        // Read witness data
        LOG_TRACE("> Reading witness file data");
        buffWitness = (FrElement *)fdWtns->getSectionData(2);

        return this->prove(buffWitness, wtnsHeader.get());
    }

    template <typename Engine>
    std::tuple<json, json> FflonkProver<Engine>::prove(FrElement *buffWitness, WtnsUtils::Header* wtnsHeader)
    {
        if(NULL == zkey) {
            throw std::runtime_error("Zkey data not set");
        }

        try
        {
            LOG_TRACE("FFLONK PROVER STARTED");

            this->buffWitness = buffWitness;

            if(NULL != wtnsHeader) {
                if (mpz_cmp(zkey->rPrime, wtnsHeader->prime) != 0)
                {
                    throw std::invalid_argument("Curve of the witness does not match the curve of the proving key");
                }

                if (wtnsHeader->nVars != zkey->nVars - zkey->nAdditions)
                {
                    std::ostringstream ss;
                    ss << "Invalid witness length. Circuit: " << zkey->nVars << ", witness: " << wtnsHeader->nVars << ", "
                    << zkey->nAdditions;
                    throw std::invalid_argument(ss.str());
                }
            }

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

            transcript->reset();
            proof->reset();

            // First element in plonk is not used and can be any value. (But always the same).
            // We set it to zero to go faster in the exponentiations.
            buffWitness[0] = E.fr.zero();

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
            // toInverse property is the variable to store the values to be inverted

            // Until this point all calculations made are circuit depending and independent from the data, so we could
            // this big function into two parts: until here circuit dependent and from here is the proof calculation

            double startTime = omp_get_wtime();

            LOG_TRACE("> Computing Additions");

            int nThreads = omp_get_max_threads() / 2;
            // Set 0's to buffers["A"], buffers["B"], buffers["C"] & buffers["Z"]
            ThreadUtils::parset(buffers["A"], 0, buffersLength * sizeof(FrElement), nThreads);

            calculateAdditions();

            // START FFLONK PROVER PROTOCOL

            // ROUND 1. Compute C1(X) polynomial
            LOG_TRACE("> ROUND 1");
            round1();

            // ROUND 2. Compute C2(X) polynomial
            LOG_TRACE("> ROUND 2");
            round2();

            // ROUND 3. Compute opening evaluations
            LOG_TRACE("> ROUND 3");
            round3();

            // ROUND 4. Compute W(X) polynomial
            LOG_TRACE("> ROUND 4");
            round4();

            // ROUND 5. Compute W'(X) polynomial
            LOG_TRACE("> ROUND 5");
            round5();

            proof->addEvaluationCommitment("inv", getMontgomeryBatchedInverse(challenges["xi"], challenges["y"], roots));

            // Prepare public inputs
            json publicSignals;
            FrElement montgomery;
            for (u_int32_t i = 1; i <= zkey->nPublic; i++)
            {
                E.fr.toMontgomery(montgomery, buffWitness[i]);
                publicSignals.push_back(E.fr.toString(montgomery).c_str());
            }

            LOG_TRACE("FFLONK PROVER FINISHED");

            ss.str("");
            ss << "Execution time: " << omp_get_wtime() - startTime << "\n";
            LOG_TRACE(ss);

            delete polynomials["A"];
            delete polynomials["B"];
            delete polynomials["C"];
            delete polynomials["C1"];
            delete polynomials["C2"];
            delete polynomials["F"];
            delete polynomials["L"];
            delete polynomials["R0"];
            delete polynomials["R1"];
            delete polynomials["R2"];
            delete polynomials["T0"];
            delete polynomials["T1"];
            delete polynomials["T1z"];
            delete polynomials["T2"];
            delete polynomials["T2z"];
            delete polynomials["Z"];
            delete polynomials["ZT"];
            delete polynomials["ZTS2"];

            delete evaluations["A"];
            delete evaluations["B"];
            delete evaluations["C"];
            delete evaluations["Z"];

            return {proof->toJson(), publicSignals};
        }
        catch (const std::exception &e)
        {
            std::cerr << "EXCEPTION: " << e.what() << "\n";
            exit(EXIT_FAILURE);
        }
    }

    template <typename Engine>
    void FflonkProver<Engine>::calculateAdditions()
    {
        for (u_int32_t i = 0; i < zkey->nAdditions; i++)
        {
            // Get witness value
            FrElement witness1 = getWitness(additionsBuff[i].signalId1);
            FrElement witness2 = getWitness(additionsBuff[i].signalId2);

            // Calculate final result
            witness1 = E.fr.mul(additionsBuff[i].factor1, witness1);
            witness2 = E.fr.mul(additionsBuff[i].factor2, witness2);
            buffInternalWitness[i] = E.fr.add(witness1, witness2);
        }
    }

    template <typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getWitness(u_int64_t idx)
    {
        u_int32_t diff = zkey->nVars - zkey->nAdditions;
        if (idx < diff)
        {
            return buffWitness[idx];
        }
        else if (idx < zkey->nVars)
        {
            return buffInternalWitness[idx - diff];
        }

        return E.fr.zero();
    }

    // ROUND 1
    template <typename Engine>
    void FflonkProver<Engine>::round1()
    {
        // STEP 1.1 - Generate random blinding scalars (b_1, ..., b9) ∈ F

        // 0 index not used, set to zero
        for (u_int32_t i = 1; i < BLINDINGFACTORSLENGTH; i++)
        {
            // memset((void *)&(blindingFactors[i].v[0]), 0, sizeof(FrElement));
            // randombytes_buf((void *)&(blindingFactors[i].v[0]), sizeof(FrElement)-1);
            memset((void *)&(blindingFactors[i].v[0]), 0, sizeof(FrElement));
            blindingFactors[i].v[0] = 1;
        }

        // STEP 1.2 - Compute wire polynomials a(X), b(X) and c(X)
        LOG_TRACE("> Computing A, B, C wire polynomials");
        computeWirePolynomials();

        // STEP 1.3 - Compute the quotient polynomial T0(X)
        LOG_TRACE("> Computing T0 polynomial");
        computeT0();

        // STEP 1.4 - Compute the FFT-style combination polynomial C1(X)
        LOG_TRACE("> Computing C1 polynomial");
        computeC1();

        // The first output of the prover is ([C1]_1)
        LOG_TRACE("> Computing C1 multi exponentiation");
        u_int64_t lengths[4] = {polynomials["A"]->getDegree() + 1,
                                polynomials["B"]->getDegree() + 1,
                                polynomials["C"]->getDegree() + 1,
                                polynomials["T0"]->getDegree() + 1};
        G1Point C1 = multiExponentiation(polynomials["C1"], 4, lengths);
        proof->addPolynomialCommitment("C1", C1);
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeWirePolynomials()
    {
        // Build A, B and C evaluations buffer from zkey and witness files
        FrElement bFactorsA[2] = {blindingFactors[2], blindingFactors[1]};
        FrElement bFactorsB[2] = {blindingFactors[4], blindingFactors[3]};
        FrElement bFactorsC[2] = {blindingFactors[6], blindingFactors[5]};

        computeWirePolynomial("A", bFactorsA);
        computeWirePolynomial("B", bFactorsB);
        computeWirePolynomial("C", bFactorsC);

        // Check degrees
        if (polynomials["A"]->getDegree() >= zkey->domainSize)
        {
            throw std::runtime_error("A Polynomial is not well calculated");
        }
        if (polynomials["B"]->getDegree() >= zkey->domainSize)
        {
            throw std::runtime_error("B Polynomial is not well calculated");
        }
        if (polynomials["C"]->getDegree() >= zkey->domainSize)
        {
            throw std::runtime_error("C Polynomial is not well calculated");
        }
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeWirePolynomial(std::string polName, FrElement blindingFactors[])
    {

        // Compute all witness from signal ids and set them to the polynomial buffers
        #pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->nConstraints; ++i)
        {
            FrElement witness = getWitness(mapBuffers[polName][i]);
            E.fr.toMontgomery(buffers[polName][i], witness);
        }

        buffers[polName][zkey->domainSize-2] = blindingFactors[1];
        buffers[polName][zkey->domainSize-1] = blindingFactors[0];

        // Create the polynomial
        // and compute the coefficients of the wire polynomials from evaluations
        std::ostringstream ss;
        ss << "··· Computing " << polName << " ifft";
        LOG_TRACE(ss);
        polynomials[polName] = Polynomial<Engine>::fromEvaluations(E, fft, buffers[polName], polPtr[polName], zkey->domainSize);

        // Compute the extended evaluations of the wire polynomials
        ss.str("");
        ss << "··· Computing " << polName << " fft";
        LOG_TRACE(ss);
        evaluations[polName] = new Evaluations<Engine>(E, fft, evalPtr[polName], *polynomials[polName], zkey->domainSize * 4);
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeT0()
    {
        #pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->domainSize * 4; i++)
        {
            //            if ((0 != i) && (i % 100000 == 0)) {
            //                ss.str("");
            //                ss << "      t0 evaluation " << i << "/" << zkey->domainSize * 4;
            ////                LOG_TRACE(ss);
            //            }

            // Get related evaluations to compute current T0 evaluation
            FrElement a = evaluations["A"]->eval[i];
            FrElement b = evaluations["B"]->eval[i];
            FrElement c = evaluations["C"]->eval[i];

            FrElement ql = evaluations["QL"]->eval[i];
            FrElement qr = evaluations["QR"]->eval[i];
            FrElement qm = evaluations["QM"]->eval[i];
            FrElement qo = evaluations["QO"]->eval[i];
            FrElement qc = evaluations["QC"]->eval[i];

            // Compute current public input
            FrElement pi = E.fr.zero();
            for (u_int32_t j = 0; j < zkey->nPublic; j++)
            {
                u_int32_t offset = (j * 4 * zkey->domainSize) + i;
                FrElement lPol = evaluations["lagrange"]->eval[offset];
                FrElement aVal = buffers["A"][j];

                pi = E.fr.sub(pi, E.fr.mul(lPol, aVal));
            }

            // T0(X) = [q_L(X)·a(X) + q_R(X)·b(X) + q_M(X)·a(X)·b(X) + q_O(X)·c(X) + q_C(X) + PI(X)] · 1/Z_H(X)
            // Compute first T0(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            // expression 1 -> q_L(X)·a(X)
            FrElement e1 = E.fr.mul(a, ql);

            // expression 2 -> q_R(X)·b(X)
            FrElement e2 = E.fr.mul(b, qr);

            // expression 3 -> q_M(X)·a(X)·b(X)
            FrElement e3 = E.fr.mul(E.fr.mul(a, b), qm);

            // expression 4 -> q_O(X)·c(X)
            FrElement e4 = E.fr.mul(c, qo);

            // t0 = expressions 1 + expression 2 + expression 3 + expression 4 + qc + pi
            FrElement t0 = E.fr.add(e1, E.fr.add(e2, E.fr.add(e3, E.fr.add(e4, E.fr.add(qc, pi)))));

            buffers["T0"][i] = t0;
        }

        // Compute the coefficients of the polynomial T0(X) from buffers.T0
        LOG_TRACE("··· Computing T0 ifft");
        polynomials["T0"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T0"], polPtr["T0"], zkey->domainSize * 4);

        // Divide the polynomial T0 by Z_H(X)
        LOG_TRACE("··· Computing T0 / ZH");
        polynomials["T0"]->divZh(zkey->domainSize);

        // Check degree
        if (polynomials["T0"]->getDegree() >= 2 * zkey->domainSize - 2)
        {
            throw std::runtime_error("T0 Polynomial is not well calculated");
        }
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeC1()
    {
        // C1(X) := a(X^4) + X · b(X^4) + X^2 · c(X^4) + X^3 · T0(X^4)
        CPolynomial<Engine> *C1 = new CPolynomial(E, 4);

        C1->addPolynomial(0, polynomials["A"]);
        C1->addPolynomial(1, polynomials["B"]);
        C1->addPolynomial(2, polynomials["C"]);
        C1->addPolynomial(3, polynomials["T0"]);

        polynomials["C1"] = C1->getPolynomial(polPtr["C1"]);

        // Check degree
        if (polynomials["C1"]->getDegree() >= 8 * zkey->domainSize - 8)
        {
            throw std::runtime_error("C1 Polynomial is not well calculated");
        }

        delete C1;
    }

    // ROUND 2
    template <typename Engine>
    void FflonkProver<Engine>::round2()
    {
        // STEP 2.1 - Compute permutation challenge beta and gamma ∈ F
        // Compute permutation challenge beta
        LOG_TRACE("> Computing challenges beta and gamma");
        transcript->reset();

        G1Point C0;
        E.g1.copy(C0, *((G1PointAffine *)zkey->C0));
        transcript->addPolCommitment(C0);

        for (u_int32_t i = 0; i < zkey->nPublic; i++)
        {
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
        u_int64_t lengths[3] = {polynomials["Z"]->getDegree() + 1,
                                polynomials["T1"]->getDegree() + 1,
                                polynomials["T2"]->getDegree() + 1};
        G1Point C2 = multiExponentiation(polynomials["C2"], 3, lengths);
        proof->addPolynomialCommitment("C2", C2);
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeZ()
    {
        FrElement *numArr = buffers["numArr"];
        FrElement *denArr = buffers["denArr"];

        LOG_TRACE("··· Computing Z evaluations");

        std::ostringstream ss;
        #pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize; i++)
        {
            //            if ((0 != i) && (i % 100000 == 0)) {
            //                ss.str("");
            //                ss << "> Computing Z evaluation " << i << "/" << zkey->domainSize;
            //                //LOG_TRACE(ss);
            //            }

            FrElement omega = fft->root(zkeyPower, i);

            // Z(X) := numArr / denArr
            // numArr := (a + beta·ω + gamma)(b + beta·ω·k1 + gamma)(c + beta·ω·k2 + gamma)
            FrElement betaw = E.fr.mul(challenges["beta"], omega);

            FrElement num1 = buffers["A"][i];
            num1 = E.fr.add(num1, betaw);
            num1 = E.fr.add(num1, challenges["gamma"]);

            FrElement num2 = buffers["B"][i];
            num2 = E.fr.add(num2, E.fr.mul(*((FrElement *)zkey->k1), betaw));
            num2 = E.fr.add(num2, challenges["gamma"]);

            FrElement num3 = buffers["C"][i];
            num3 = E.fr.add(num3, E.fr.mul(*((FrElement *)zkey->k2), betaw));
            num3 = E.fr.add(num3, challenges["gamma"]);

            numArr[i] = E.fr.mul(num1, E.fr.mul(num2, num3));

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

            denArr[i] = E.fr.mul(den1, E.fr.mul(den2, den3));
        }

        FrElement numPrev = numArr[0];
        FrElement denPrev = denArr[0];
        FrElement numCur, denCur;

        for (u_int64_t i = 0; i < zkey->domainSize - 1; i++)
        {
            numCur = numArr[i + 1];
            denCur = denArr[i + 1];

            numArr[i + 1] = numPrev;
            denArr[i + 1] = denPrev;

            numPrev = E.fr.mul(numPrev, numCur);
            denPrev = E.fr.mul(denPrev, denCur);
        }

        numArr[0] = numPrev;
        denArr[0] = denPrev;

        // Compute the inverse of denArr to compute in the next command the
        // division numArr/denArr by multiplying num · 1/denArr
        batchInverse(denArr, zkey->domainSize);

        // Multiply numArr · denArr where denArr was inverted in the previous command
        #pragma omp parallel for
        for (u_int32_t i = 0; i < zkey->domainSize; i++)
        {
            buffers["Z"][i] = E.fr.mul(numArr[i], denArr[i]);
        }

        if (!E.fr.eq(buffers["Z"][0], E.fr.one()))
        {
            throw std::runtime_error("Copy constraints does not match");
        }

        // Compute polynomial coefficients z(X) from buffers.Z
        LOG_TRACE("··· Computing Z ifft");
        polynomials["Z"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["Z"], polPtr["Z"], zkey->domainSize, 3);

        // Compute extended evaluations of z(X) polynomial
        LOG_TRACE("··· Computing Z fft");
        evaluations["Z"] = new Evaluations<Engine>(E, fft, evalPtr["Z"], *polynomials["Z"], zkey->domainSize * 4);

        // Blind z(X) polynomial coefficients with blinding scalars b
        FrElement bFactors[3] = {blindingFactors[9], blindingFactors[8], blindingFactors[7]};
        polynomials["Z"]->blindCoefficients(bFactors, 3);

        // Check degree
        if (polynomials["Z"]->getDegree() >= zkey->domainSize + 3)
        {
            throw std::runtime_error("Z Polynomial is not well calculated");
        }
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeT1()
    {
        LOG_TRACE("··· Computing T1 evaluations");

        std::ostringstream ss;

        #pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize * 2; i++)
        {
            //            if ((0 != i) && (i % 100000 == 0)) {
            //                ss.str("");
            //                ss << "    T1 evaluation " << i << "/" << zkey->domainSize * 4;
            //                //LOG_TRACE(ss);
            //            }

            FrElement omega = fft->root(zkeyPower + 1, i);
            FrElement omega2 = E.fr.square(omega);

            FrElement z = evaluations["Z"]->eval[i * 2];
            FrElement zp = E.fr.add(E.fr.add(
                                        E.fr.mul(blindingFactors[7], omega2), E.fr.mul(blindingFactors[8], omega)),
                                    blindingFactors[9]);

            // T1(X) := (z(X) - 1) · L_1(X)
            // Compute first T1(X)·Z_H(X), so divide later the resulting polynomial by Z_H(X)
            FrElement lagrange1 = evaluations["lagrange"]->eval[i * 2];
            FrElement t1 = E.fr.mul(E.fr.sub(z, E.fr.one()), lagrange1);
            FrElement t1z = E.fr.mul(zp, lagrange1);

            buffers["T1"][i] = t1;
            buffers["T1z"][i] = t1z;
        }

        // Compute the coefficients of the polynomial T1(X) from buffers.T1
        LOG_TRACE("··· Computing T1 ifft");
        polynomials["T1"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T1"], polPtr["T1"], zkey->domainSize * 2);

        // Divide the polynomial T1 by Z_H(X)
        polynomials["T1"]->divZh(zkey->domainSize, 2);

        // Compute the coefficients of the polynomial T1z(X) from buffers.T1z
        LOG_TRACE("··· Computing T1z ifft");
        polynomials["T1z"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T1z"], polPtr["T1z"], zkey->domainSize * 2);

        // Add the polynomial T1z to T1 to get the final polynomial T1
        polynomials["T1"]->add(*polynomials["T1z"]);

        // Check degree
        if (polynomials["T1"]->getDegree() >= zkey->domainSize + 2)
        {
            throw std::runtime_error("T1 Polynomial is not well calculated");
        }
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeT2()
    {
        LOG_TRACE("··· Computing T2 evaluations");

        std::ostringstream ss;
        #pragma omp parallel for
        for (u_int64_t i = 0; i < zkey->domainSize * 4; i++)
        {
            //            if ((0 != i) && (i % 100000 == 0)) {
            //                ss.str("");
            //                ss << "    T2 evaluation " << i << "/" << zkey->domainSize * 4;
            //                //LOG_TRACE(ss);
            //            }

            FrElement omega = fft->root(zkeyPower + 2, i);
            FrElement omega2 = E.fr.square(omega);
            FrElement omegaW = E.fr.mul(omega, fft->root(zkeyPower, 1));
            FrElement omegaW2 = E.fr.square(omegaW);

            FrElement a = evaluations["A"]->eval[i];
            FrElement b = evaluations["B"]->eval[i];
            FrElement c = evaluations["C"]->eval[i];
            FrElement z = evaluations["Z"]->eval[i];
            FrElement zW = evaluations["Z"]->eval[(zkey->domainSize * 4 + 4 + i) % (zkey->domainSize * 4)];

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

            FrElement e12 = E.fr.add(b, E.fr.mul(betaX, *((FrElement *)zkey->k1)));
            e12 = E.fr.add(e12, challenges["gamma"]);

            FrElement e13 = E.fr.add(c, E.fr.mul(betaX, *((FrElement *)zkey->k2)));
            e13 = E.fr.add(e13, challenges["gamma"]);

            FrElement e1 = E.fr.mul(E.fr.mul(E.fr.mul(e11, e12), e13), z);
            FrElement e1z = E.fr.mul(E.fr.mul(E.fr.mul(e11, e12), e13), zp);

            // expression 2 -> (a(X) + beta·sigma1(X) + gamma)(b(X) + beta·sigma2(X) + gamma)(c(X) + beta·sigma3(X) + gamma)z(Xω)
            FrElement e21 = E.fr.add(a, E.fr.mul(challenges["beta"], sigma1));
            e21 = E.fr.add(e21, challenges["gamma"]);

            FrElement e22 = E.fr.add(b, E.fr.mul(challenges["beta"], sigma2));
            e22 = E.fr.add(e22, challenges["gamma"]);

            FrElement e23 = E.fr.add(c, E.fr.mul(challenges["beta"], sigma3));
            e23 = E.fr.add(e23, challenges["gamma"]);

            FrElement e2 = E.fr.mul(E.fr.mul(E.fr.mul(e21, e22), e23), zW);
            FrElement e2z = E.fr.mul(E.fr.mul(E.fr.mul(e21, e22), e23), zWp);

            FrElement t2 = E.fr.sub(e1, e2);
            FrElement t2z = E.fr.sub(e1z, e2z);

            buffers["T2"][i] = t2;
            buffers["T2z"][i] = t2z;
        }

        // Compute the coefficients of the polynomial T2(X) from buffers.T2
        LOG_TRACE("··· Computing T2 ifft");
        polynomials["T2"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T2"], polPtr["T2"], zkey->domainSize * 4);

        // Divide the polynomial T2 by Z_H(X)
        polynomials["T2"]->divZh(zkey->domainSize);

        // Compute the coefficients of the polynomial T2z(X) from buffers.T2z
        LOG_TRACE("··· Computing T2z ifft");
        polynomials["T2z"] = Polynomial<Engine>::fromEvaluations(E, fft, buffers["T2z"], polPtr["T2z"], zkey->domainSize * 4);

        // Add the polynomial T2z to T2 to get the final polynomial T2
        polynomials["T2"]->add(*polynomials["T2z"]);

        // Check degree
        if (polynomials["T2"]->getDegree() >= 3 * zkey->domainSize)
        {
            throw std::runtime_error("T2 Polynomial is not well calculated");
        }
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeC2()
    {
        // C2(X) := z(X^3) + X · T1(X^3) + X^2 · T2(X^3)
        CPolynomial<Engine> *C2 = new CPolynomial(E, 3);

        C2->addPolynomial(0, polynomials["Z"]);
        C2->addPolynomial(1, polynomials["T1"]);
        C2->addPolynomial(2, polynomials["T2"]);

        polynomials["C2"] = C2->getPolynomial(polPtr["C2"]);

        // Check degree
        if (polynomials["C2"]->getDegree() >= 9 * zkey->domainSize)
        {
            throw std::runtime_error("C2 Polynomial is not well calculated");
        }

        delete C2;
    }

    // ROUND 3
    template <typename Engine>
    void FflonkProver<Engine>::round3()
    {
        LOG_TRACE("> Computing challenge xi");
        // STEP 3.1 - Compute evaluation challenge xi ∈ S
        transcript->reset();
        transcript->addScalar(challenges["gamma"]);
        transcript->addPolCommitment(proof->getPolynomialCommitment("C2"));

        // Obtain a xi_seeder from the transcript
        // To force h1^4 = xi, h2^3 = xi and h_3^2 = xiω
        // we compute xi = xi_seeder^12, h1 = xi_seeder^3, h2 = xi_seeder^4 and h3 = xi_seeder^6
        challenges["xiSeed"] = transcript->getChallenge();
        FrElement xiSeed2;
        E.fr.square(xiSeed2, challenges["xiSeed"]);

        // Compute omega8, omega4 and omega3
        omegas["w8"][0] = E.fr.one();
        for (uint i = 1; i < 8; i++)
        {
            omegas["w8"][i] = E.fr.mul(omegas["w8"][i - 1], *((FrElement *)zkey->w8));
        }

        // Compute omega3 and omega4
        omegas["w4"][0] = E.fr.one();
        for (uint i = 1; i < 4; i++)
        {
            omegas["w4"][i] = E.fr.mul(omegas["w4"][i - 1], *((FrElement *)zkey->w4));
        }

        omegas["w3"][0] = E.fr.one();
        omegas["w3"][1] = *((FrElement *)zkey->w3);
        E.fr.square(omegas["w3"][2], omegas["w3"][1]);

        // Compute h0 = xiSeeder^3
        roots["S0h0"][0] = E.fr.mul(xiSeed2, challenges["xiSeed"]);
        for (uint i = 1; i < 8; i++)
        {
            roots["S0h0"][i] = E.fr.mul(roots["S0h0"][0], omegas["w8"][i]);
        }

        // Compute h1 = xi_seeder^6
        roots["S1h1"][0] = E.fr.square(roots["S0h0"][0]);
        for (uint i = 1; i < 4; i++)
        {
            roots["S1h1"][i] = E.fr.mul(roots["S1h1"][0], omegas["w4"][i]);
        }

        // Compute h2 = xi_seeder^8
        roots["S2h2"][0] = E.fr.mul(roots["S1h1"][0], xiSeed2);
        roots["S2h2"][1] = E.fr.mul(roots["S2h2"][0], omegas["w3"][1]);
        roots["S2h2"][2] = E.fr.mul(roots["S2h2"][0], omegas["w3"][2]);

        // Multiply h3 by third-root-omega to obtain h_3^3 = xiω
        // So, h3 = xi_seeder^8 ω^{1/3}
        roots["S2h2"][3] = E.fr.mul(roots["S2h2"][0], *((FrElement *)zkey->wr));
        roots["S2h2"][4] = E.fr.mul(roots["S2h2"][3], omegas["w3"][1]);
        roots["S2h2"][5] = E.fr.mul(roots["S2h2"][3], omegas["w3"][2]);

        // Compute xi = xi_seeder^24
        challenges["xi"] = E.fr.mul(E.fr.square(roots["S2h2"][0]), roots["S2h2"][0]);

        std::ostringstream ss;
        ss << "··· challenges.xi: " << E.fr.toString(challenges["xi"]);
        LOG_TRACE(ss);

        // STEP 3.2 - Compute opening evaluations and add them to the proof (third output of the prover)
        LOG_TRACE("··· Computing evaluations");
        proof->addEvaluationCommitment("ql", polynomials["QL"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qr", polynomials["QR"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qm", polynomials["QM"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qo", polynomials["QO"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("qc", polynomials["QC"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("s1", polynomials["Sigma1"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("s2", polynomials["Sigma2"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("s3", polynomials["Sigma3"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("a", polynomials["A"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("b", polynomials["B"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("c", polynomials["C"]->fastEvaluate(challenges["xi"]));
        proof->addEvaluationCommitment("z", polynomials["Z"]->fastEvaluate(challenges["xi"]));

        challenges["xiw"] = E.fr.mul(challenges["xi"], fft->root(zkeyPower, 1));
        proof->addEvaluationCommitment("zw", polynomials["Z"]->fastEvaluate(challenges["xiw"]));
        proof->addEvaluationCommitment("t1w", polynomials["T1"]->fastEvaluate(challenges["xiw"]));
        proof->addEvaluationCommitment("t2w", polynomials["T2"]->fastEvaluate(challenges["xiw"]));
    }

    // ROUND 4
    template <typename Engine>
    void FflonkProver<Engine>::round4()
    {
        LOG_TRACE("> Computing challenge alpha");

        FrElement challengeAlpha = computeChallengeAlpha(challenges["xiSeed"]);
        challenges["alpha"] = challengeAlpha;
    
        std::ostringstream ss;
        ss << "··· challenges.alpha: " << E.fr.toString(challengeAlpha);
        LOG_TRACE(ss);

        std::map<std::string, Polynomial<Engine> *> polynomialsF;
        polynomialsF["C0"] = polynomials["C0"];
        polynomialsF["C1"] = polynomials["C1"];
        polynomialsF["C2"] = polynomials["C2"];

        // STEP 4.2 - Compute F(X)
        LOG_TRACE("> Computing R polynomials");
        std::map<std::string, Polynomial<Engine> *> polynomialsR = computeR(roots, polynomialsF);
        polynomials["R0"] = polynomialsR["R0"];
        polynomials["R1"] = polynomialsR["R1"];
        polynomials["R2"] = polynomialsR["R2"];

        LOG_TRACE("> Computing F polynomial");
        polynomials["F"] = computeW(challengeAlpha, polynomialsF, polynomialsR);

        // The fourth output of the prover is ([W1]_1), where W1:=(f/Z_t)(x)
        LOG_TRACE("> Computing W1 multi exponentiation");
        u_int64_t lengths[1] = {polynomials["F"]->getDegree() + 1};
        G1Point W1 = multiExponentiation(polynomials["F"], 1, lengths);

        proof->addPolynomialCommitment("W1", W1);
    }

    template <typename Engine>
    Polynomial<Engine> * FflonkProver<Engine>::computeRi(uint32_t degree, FrElement* rootsRi, Polynomial<Engine> * polynomialFi) {

        FrElement* evals = new FrElement[degree];

        for(uint32_t i = 0; i < degree; ++i) {
            evals[i] = polynomialFi->fastEvaluate(rootsRi[i]);
        }

        Polynomial<Engine> * polynomialRi = Polynomial<Engine>::lagrangePolynomialInterpolation(rootsRi, evals, degree);
        
        // Check the degree of R0(X) < (degree - 1)
        if (polynomialRi->getDegree() > (degree - 1))
        {
            throw std::runtime_error("Polynomial is not well calculated");
        }

        return polynomialRi;
    }

    template <typename Engine>
    std::map<std::string, Polynomial<Engine> *> FflonkProver<Engine>::computeR(std::map<std::string, FrElement *> rootsMap, std::map<std::string, Polynomial<Engine> *> polynomialsF) {
        uint32_t degree[3] = {8, 4, 6};
        
        std::map<std::string, Polynomial<Engine> *> polynomialsR;

        for(uint32_t i = 0; i < rootsMap.size(); ++i) {
            std::string log = "> Computing R" + std::to_string(i) + " polynomial"; 
            LOG_TRACE(log);
            polynomialsR["R" + std::to_string(i)] = computeRi(degree[i], roots["S" + std::to_string(i) + "h" + std::to_string(i)], polynomialsF["C" + std::to_string(i)]);
        }

        return polynomialsR;
    }

    template <typename Engine>
    Polynomial<Engine> * FflonkProver<Engine>::computeW(FrElement challengeAlpha, std::map<std::string, Polynomial<Engine> *> polynomialsF, std::map<std::string, Polynomial<Engine> *> polynomialsR)
    {

        LOG_TRACE("··· Computing W polynomial");

        FrElement alpha = E.fr.one();

        Polynomial<Engine> * polynomialF = new Polynomial<Engine>(E, polPtr["F"], zkey->domainSize * 16);
        
        for(uint32_t i = 0; i < roots.size(); ++i) {
            auto fTmp = Polynomial<Engine>::fromPolynomial(E, *polynomialsF["C" + std::to_string(i)], polPtr["tmp"]);
            fTmp->sub(*polynomialsR["R" + std::to_string(i)]);
            fTmp->mulScalar(alpha);
            if(i == 0) {
                fTmp->divByZerofier(8, challenges["xi"]);
            } else if(i == 1) {
                fTmp->divByZerofier(4, challenges["xi"]);
            } else if (i == 2) {
                fTmp->divByZerofier(3, challenges["xi"]);
                fTmp->divByZerofier(3, challenges["xiw"]);
            }

            polynomialF->add(*fTmp);

            alpha = E.fr.mul(alpha, challengeAlpha);
           
        }

        // Check degree
        if (polynomialF->getDegree() >= 9 * zkey->domainSize - 6)
        {
            throw std::runtime_error("W Polynomial is not well calculated");
        }

        return polynomialF;
        
    }

    template <typename Engine>
    Polynomial<Engine>* FflonkProver<Engine>::computeZT(std::map<std::string, FrElement *> rootsMap)
    {
        uint32_t degree[3] = {8, 4, 6};
        
        uint32_t nRoots = 0;
        for(uint32_t i = 0; i < rootsMap.size(); ++i) {
            nRoots += degree[i];
        }

        FrElement arr[nRoots];

        uint32_t index = 0;
        for(uint32_t i = 0; i < rootsMap.size(); ++i) {
            for(uint32_t j = 0; j < degree[i]; ++j) {
                arr[index++] = rootsMap["S" + std::to_string(i) + "h" + std::to_string(i)][j];
            }

        }

        return Polynomial<Engine>::zerofierPolynomial(arr, nRoots);
    }

    // ROUND 5
    template <typename Engine>
    void FflonkProver<Engine>::round5()
    {
        FrElement challengeY = computeChallengeY(challenges["alpha"], proof->getPolynomialCommitment("W1"));
        challenges["y"] = challengeY;

        std::ostringstream ss;
        ss << "··· challenges.y: " << E.fr.toString(challengeY);
        LOG_TRACE(ss);

        LOG_TRACE("> Compute Wp polynomial");
        std::map<std::string, Polynomial<Engine> *> polynomialsR;
        polynomialsR["R0"] = polynomials["R0"];
        polynomialsR["R1"] = polynomials["R1"];
        polynomialsR["R2"] = polynomials["R2"];

        std::map<std::string, Polynomial<Engine> *> polynomialsF;
        polynomialsF["C0"] = polynomials["C0"];
        polynomialsF["C1"] = polynomials["C1"];
        polynomialsF["C2"] = polynomials["C2"];


        polynomials["L"] = computeWp(challenges["y"], challenges["alpha"], roots, polynomialsF, polynomialsR, polynomials["F"]);

        // The fifth output of the prover is ([W2]_1), where W2:=(f/Z_t)(x)
        LOG_TRACE("> Computing W' multi exponentiation");
        u_int64_t lengths[1] = {polynomials["L"]->getDegree() + 1};
        G1Point W2 = multiExponentiation(polynomials["L"], 1, lengths);

        proof->addPolynomialCommitment("W2", W2);
    }

    template <typename Engine>
    Polynomial<Engine> * FflonkProver<Engine>::computeWp(FrElement challengeY, FrElement challengeAlpha, std::map<std::string, FrElement *> rootsMap, std::map<std::string, Polynomial<Engine> *> polynomialsF, std::map<std::string, Polynomial<Engine> *> polynomialsR, Polynomial<Engine> * polynomialW) {
        LOG_TRACE("> Computing L polynomial");
        Polynomial<Engine>* polynomialL = computeL(challengeY, challengeAlpha, rootsMap, polynomialsF, polynomialsR, polynomialW);

        LOG_TRACE("> Computing ZTS2 polynomial");
        Polynomial<Engine>* polynomialZTS2 = computeZTS2(rootsMap);

        FrElement ZTS2Y = polynomialZTS2->fastEvaluate(challengeY);
        E.fr.inv(ZTS2Y, ZTS2Y);
        polynomialL->mulScalar(ZTS2Y);

        LOG_TRACE("> Computing W' = L / ZTS2 polynomial");
        polynomialL->divByZerofier(1, challengeY);

        if (polynomialL->getDegree() >= 9 * zkey->domainSize - 1)
        {
            throw std::runtime_error("Degree of L(X)/(ZTS2(y)(X-y)) is not correct");
        }

        return polynomialL;
    }

    template <typename Engine>
    Polynomial<Engine> * FflonkProver<Engine>::computeL(FrElement challengeY, FrElement challengeAlpha, std::map<std::string, FrElement *> rootsMap, std::map<std::string, Polynomial<Engine> *> polynomialsF, std::map<std::string, Polynomial<Engine> *> polynomialsR, Polynomial<Engine> * polynomialW)
    {
        LOG_TRACE("··· Computing L polynomial");

        FrElement mulL[rootsMap.size()];
        FrElement preL[rootsMap.size()];
        FrElement evalRiY[rootsMap.size()];

        uint32_t degree[3] = {8, 4, 6};

        for(uint32_t i = 0; i < rootsMap.size(); ++i) {
            mulL[i] = E.fr.sub(challengeY, rootsMap["S" + std::to_string(i) + "h" + std::to_string(i)][0]);
            for (uint32_t j = 1; j < degree[i]; j++)
            {
                mulL[i] = E.fr.mul(mulL[i], E.fr.sub(challengeY, rootsMap["S" + std::to_string(i) + "h" + std::to_string(i)][j]));
            }
        }

       
        FrElement alpha = E.fr.one();
        for(uint32_t i = 0; i < rootsMap.size(); ++i) {
            preL[i] = alpha;
            for(uint32_t j = 0; j < rootsMap.size(); ++j) {
                if(i != j) {
                    preL[i] = E.fr.mul(preL[i], mulL[j]);
                }
            }
            alpha = E.fr.mul(alpha,challengeAlpha);
        }

        for(uint32_t i = 0; i < rootsMap.size(); ++i) {
           evalRiY[i] = polynomialsR["R" + std::to_string(i)]->fastEvaluate(challengeY);
        }

        // COMPUTE L(X)
        Polynomial<Engine>* polynomialL = new Polynomial<Engine>(E, polPtr["L"], zkey->domainSize * 16);
        for(uint32_t i = 0; i < rootsMap.size(); ++i) {
            auto fTmp = Polynomial<Engine>::fromPolynomial(E, *polynomialsF["C" + std::to_string(i)], polPtr["tmp"]);
            fTmp->subScalar(evalRiY[i]);
            fTmp->mulScalar(preL[i]);
            polynomialL->add(*fTmp);

        }
       
        LOG_TRACE("> Computing ZT polynomial");
        Polynomial<Engine>* polynomialZT = computeZT(rootsMap);

        FrElement evalZTY = polynomialZT->fastEvaluate(challengeY);
        polynomialW->mulScalar(evalZTY);
        polynomialL->sub(*polynomialW);

        // Check degree
        if (polynomialL->getDegree() >= 9 * zkey->domainSize)
        {
            throw std::runtime_error("L Polynomial is not well calculated");
        }

        toInverse["denH1"] = mulL[1];
        toInverse["denH2"] = mulL[2];

        return polynomialL;

        
    }

    template <typename Engine>
    Polynomial<Engine>* FflonkProver<Engine>::computeZTS2(std::map<std::string, FrElement *> rootsMap)
    {
         uint32_t degree[3] = {8, 4, 6};
        
        uint32_t nRoots = 0;
        for(uint32_t i = 1; i < roots.size(); ++i) {
            nRoots += degree[i];
        }

        FrElement arr[nRoots];

        uint32_t index = 0;
        for(uint32_t i = 1; i < roots.size(); ++i) {
            for(uint32_t j = 0; j < degree[i]; ++j) {
                arr[index++] = roots["S" + std::to_string(i) + "h" + std::to_string(i)][j];
            }

        }

        return Polynomial<Engine>::zerofierPolynomial(arr, nRoots);
    }

    template <typename Engine>
    void FflonkProver<Engine>::batchInverse(FrElement *elements, u_int64_t length)
    {
        // Calculate products: a, ab, abc, abcd, ...
        products[0] = elements[0];
        for (u_int64_t index = 1; index < length; index++)
        {
            E.fr.mul(products[index], products[index - 1], elements[index]);
        }

        // Calculate inverses: 1/a, 1/ab, 1/abc, 1/abcd, ...
        E.fr.inv(inverses[length - 1], products[length - 1]);
        for (uint64_t index = length - 1; index > 0; index--)
        {
            E.fr.mul(inverses[index - 1], inverses[index], elements[index]);
        }

        elements[0] = inverses[0];
        for (u_int64_t index = 1; index < length; index++)
        {
            E.fr.mul(elements[index], inverses[index], products[index - 1]);
        }
    }

    template <typename Engine>
    void FflonkProver<Engine>::computeLiSinglePoint(uint32_t index, uint32_t degree, FrElement challengesY, std::map<std::string, FrElement *> rootsMap)
    {    
        std::ostringstream ss;

        FrElement den1 = E.fr.set(degree);
        for (u_int64_t i = 0; i < (degree - 2); i++)
        {
            den1 = E.fr.mul(den1, rootsMap["S" + std::to_string(index) + "h" + std::to_string(index)][0]);
        }

        // Compute L_i^{(S0)}(y)
        for (uint j = 0; j < degree; j++) {

            FrElement den2 = rootsMap["S" + std::to_string(index) + "h" + std::to_string(index)][((degree - 1) * j) % degree];
            FrElement den3 = E.fr.sub(challengesY, rootsMap["S" + std::to_string(index) + "h" + std::to_string(index)][j]);

            ss.str("");
            ss << "LiS" << index << "_" << (j + 1) << " ";
            toInverse[ss.str()] = E.fr.mul(E.fr.mul(den1, den2), den3);
        }
        return;
    }

    
    template <typename Engine>
    void FflonkProver<Engine>::computeLiS2()
    {    
        std::ostringstream ss;

        // Compute L_i^{(S2)}(y)
        FrElement den1 = E.fr.mul(E.fr.mul(E.fr.set(3), roots["S2h2"][0]), E.fr.sub(challenges["xi"], challenges["xiw"]));
        for (uint j = 0; j < 3; j++) {
            FrElement den2 = roots["S2h2"][2*j % 3];
            FrElement den3 = E.fr.sub(challenges["y"], roots["S2h2"][j]);
             
            ss.str("");
            ss << "LiS2_" << (j + 1) << " ";
            toInverse[ss.str()] = E.fr.mul(E.fr.mul(den1, den2), den3);
        }

        den1 = E.fr.mul(E.fr.mul(E.fr.set(3), roots["S2h2"][3]), E.fr.sub(challenges["xiw"], challenges["xi"]));
        for (uint j = 0; j < 3; j++) {
            FrElement den2 = roots["S2h2"][(2*j % 3) + 3];
            FrElement den3 = E.fr.sub(challenges["y"], roots["S2h2"][j + 3]);
             
            ss.str("");
            ss << "LiS2_" << (j + 1 + 3) << " ";
            toInverse[ss.str()] = E.fr.mul(E.fr.mul(den1, den2), den3);
        }
        return;
    }

    template <typename Engine>
    typename Engine::FrElement *FflonkProver<Engine>::polynomialFromMontgomery(Polynomial<Engine> *polynomial)
    {
        const u_int64_t length = polynomial->getLength();

        FrElement *result = buffers["tmp"];
        int nThreads = omp_get_max_threads() / 2;
        ThreadUtils::parset(result, 0, length * sizeof(FrElement), nThreads);

#pragma omp parallel for
        for (u_int32_t index = 0; index < length; ++index)
        {
            E.fr.fromMontgomery(result[index], polynomial->coef[index]);
        }

        return result;
    }

    template <typename Engine>
    typename Engine::G1Point FflonkProver<Engine>::multiExponentiation(Polynomial<Engine> *polynomial)
    {
        G1Point value;
        FrElement *pol = this->polynomialFromMontgomery(polynomial);

        E.g1.multiMulByScalar(value, PTau, (uint8_t *)pol, sizeof(pol[0]), polynomial->getDegree() + 1);

        return value;
    }

    template <typename Engine>
    typename Engine::G1Point
    FflonkProver<Engine>::multiExponentiation(Polynomial<Engine> *polynomial, u_int32_t nx, u_int64_t x[])
    {
        G1Point value;
        FrElement *pol = this->polynomialFromMontgomery(polynomial);

        E.g1.multiMulByScalar(value, PTau, (uint8_t *)pol, sizeof(pol[0]), polynomial->getDegree() + 1, nx, x);

        return value;
    }

    //HELPERS
    template <typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::computeChallengeXiSeed(FrElement previousChallenge)
    {    
        uint32_t indexes[1] = {2};

        transcript->reset();
        if(previousChallenge != NULL) {
            transcript->addScalar(previousChallenge);
        }

        for(uint32_t i = 0; i < sizeof(indexes); ++i) {
            transcript->addPolCommitment(proof->getPolynomialCommitment("C" + std::to_string(indexes[i])));
        }

        return transcript->getChallenge();
    }


    template <typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::computeChallengeAlpha(FrElement challengeXiSeed)
    {    
        uint32_t nEvaluations = proof->getNumberEvaluations();

        // TODO: SORT !!!
        //string evaluations[nEvaluations] = proof->getEvaluationsNames();
        string evaluations[nEvaluations] = {"ql", "qr", "qm", "qo", "qc", "s1", "s2", "s3", "a", "b", "c", "z", "zw", "t1w", "t2w"};

        transcript->reset();
        transcript->addScalar(challengeXiSeed);

        for (uint32_t i = 0; i < nEvaluations; ++i) {
            transcript->addScalar(proof->getEvaluationCommitment(evaluations[i]));
        }
        return transcript->getChallenge();
    }

    template <typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::computeChallengeY(FrElement challengeAlpha, G1Point W)
    {    
        transcript->reset();
        transcript->addScalar(challengeAlpha);
        transcript->addPolCommitment(W);
        return transcript->getChallenge();
    }

    template <typename Engine>
    std::map<std::string, typename Engine::FrElement*> FflonkProver<Engine>::calculateRoots(FrElement challengeXiSeed) {
        uint32_t degrees[3] = {8, 4, 3};

        std::map<std::string, FrElement *> rootsMap;
        std::map<std::string, FrElement *> omegasMap;

        for(uint32_t i = 0; i < sizeof(degrees); i++) {

            omegasMap["w" + std::to_string(degrees[i])][0] = E.fr.one();
            for (uint i = 1; i < degrees[i]; i++)
            {   
                //TODO: GENERALISE THIS
                FrElement val;
                if(degrees[i] == 8) {
                    val = *((FrElement *)zkey->w8);
                } else if(degrees[i] == 4) {
                    val = *((FrElement *)zkey->w4);
                } else if(degrees[i] == 3) {
                    val = *((FrElement *)zkey->w3);
                }
                omegasMap["w" + std::to_string(degrees[i])][i] = E.fr.mul(omegasMap["w" + std::to_string(degrees[i])][i - 1], val);
            }

            rootsMap["S" + std::to_string(i) + "h" + std::to_string(i)] = new FrElement[degrees[i]];
        }

        FrElement xiSeed2;
        E.fr.square(xiSeed2, challengeXiSeed);

        // Compute h0 = xiSeeder^3
        rootsMap["S0h0"][0] = E.fr.mul(xiSeed2, challengeXiSeed);
        for (uint i = 1; i < 8; i++)
        {
            rootsMap["S0h0"][i] = E.fr.mul(rootsMap["S0h0"][0], omegasMap["w8"][i]);
        }

        // Compute h1 = xi_seeder^6
        rootsMap["S1h1"][0] = E.fr.square(rootsMap["S0h0"][0]);
        for (uint i = 1; i < 4; i++)
        {
            rootsMap["S1h1"][i] = E.fr.mul(rootsMap["S1h1"][0], omegasMap["w4"][i]);
        }

        // Compute h2 = xi_seeder^8
        rootsMap["S2h2"][0] = E.fr.mul(rootsMap["S1h1"][0], xiSeed2);
        rootsMap["S2h2"][1] = E.fr.mul(rootsMap["S2h2"][0], omegasMap["w3"][1]);
        rootsMap["S2h2"][2] = E.fr.mul(rootsMap["S2h2"][0], omegasMap["w3"][2]);

        // Multiply h3 by third-root-omega to obtain h_3^3 = xiω
        // So, h3 = xi_seeder^8 ω^{1/3}
        rootsMap["S2h2"][3] = E.fr.mul(rootsMap["S2h2"][0], *((FrElement *)zkey->wr));
        rootsMap["S2h2"][4] = E.fr.mul(rootsMap["S2h2"][3], omegasMap["w3"][1]);
        rootsMap["S2h2"][5] = E.fr.mul(rootsMap["S2h2"][3], omegasMap["w3"][2]);


        return rootsMap;
    }

    template <typename Engine>
    void FflonkProver<Engine>::open(FrElement previousChallenge) {
        FrElement challengeXiSeed = computeChallengeXiSeed(previousChallenge);

        // Calculate roots
        std::map<std::string, FrElement *> rootsMap = calculateRoots(challengeXiSeed);

        // TODO: 24 is equal to powerW
        FrElement challengeXi = E.fr.one();
        for(uint32_t i = 0; i < 24; ++i) {
            challengeXi = E.fr.mul(challengeXi, challengeXiSeed);
        }

        //TODO: Calculate Evaluations
        

        FrElement challengeAlpha = computeChallengeAlpha(challengeXiSeed);

        std::map<std::string, Polynomial<Engine> *> polynomialsF;
        polynomialsF["C0"] = polynomials["C0"];
        polynomialsF["C1"] = polynomials["C1"];
        polynomialsF["C2"] = polynomials["C2"];

        std::map<std::string, Polynomial<Engine> *> polynomialsR = computeR(rootsMap, polynomialsF);

        

        Polynomial<Engine>* polynomialW = computeW(challengeAlpha, polynomialsF, polynomialsR);

        LOG_TRACE("> Computing W multi exponentiation");
        G1Point W = multiExponentiation(polynomialW, 1, {polynomialW->getDegree() + 1});
        proof->addPolynomialCommitment("W1", W);

        FrElement challengeY = computeChallengeY(challengeAlpha, W);

        Polynomial<Engine>* polynomialWp = computeWp(challengeY, challengeAlpha, rootsMap, polynomialsR, polynomialW);
        G1Point Wp = multiExponentiation(polynomialWp, 1, {polynomialWp->getDegree() + 1});
        proof->addPolynomialCommitment("W2", Wp);

        proof->addEvaluationCommitment("inv", getMontgomeryBatchedInverse(challengeXi, challengeY, roots));
    }

     template <typename Engine>
    typename Engine::FrElement FflonkProver<Engine>::getMontgomeryBatchedInverse(FrElement challengeXi, FrElement challengeY, std::map<std::string, FrElement *> rootsMap)
    {
        std::map <std::string, FrElement> inverseElements;

        std::ostringstream ss;
        //   · denominator needed in step 8 and 9 of the verifier to multiply by 1/Z_H(xi)
        FrElement xiN = challengeXi;
        for (u_int64_t i = 0; i < zkeyPower; i++)
        {
            xiN = E.fr.square(xiN);
        }

        inverseElements["zh"] = E.fr.sub(xiN, E.fr.one());

        FrElement mulL[rootsMap.size() - 1];

        uint32_t degree[3] = {8, 4, 6};

        for(uint32_t i = 0; i < rootsMap.size() - 1; ++i) {
            mulL[i] = E.fr.sub(challengeY, rootsMap["S" + std::to_string(i+1) + "h" + std::to_string(i+1)][0]);
            for (uint32_t j = 1; j < degree[i]; j++)
            {
                mulL[i] = E.fr.mul(mulL[i], E.fr.sub(challengeY, rootsMap["S" + std::to_string(i+1) + "h" + std::to_string(i+1)][j]));
            }
            inverseElements["mulL" + std::to_string(i+1)] = mulL[i];
        }

        computeLiSinglePoint(0, 8, challenges["y"], roots);

        computeLiSinglePoint(1, 4, challenges["y"], roots);

        // TODO
        computeLiS2();

        FrElement mulAccumulator = E.fr.one();
        for (auto &[key, value] : toInverse)
        {
            mulAccumulator = E.fr.mul(mulAccumulator, value);
        }

        E.fr.inv(mulAccumulator, mulAccumulator);
        return mulAccumulator;
    }

    
}
