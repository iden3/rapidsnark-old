#ifndef FFLONK_PROVER_HPP
#define FFLONK_PROVER_HPP

#include <string>
#include <map>
#include "snark_proof.hpp"
#include "binfile_utils.hpp"
#include <gmp.h>
#include "fft.hpp"
#include "zkey_fflonk.hpp"
#include "polynomial/polynomial.hpp"
#include "polynomial/evaluations.hpp"
#include <nlohmann/json.hpp>
#include "mul_z.hpp"

using json = nlohmann::json;

namespace Fflonk {

    template<typename Engine>
    class FflonkProver {
        using FrElement = typename Engine::FrElement;
        using G1Point = typename Engine::G1Point;

        Engine &E;
        FFT<typename Engine::Fr> *fft;
        MulZ<Engine> *mulZ;

        // To compute in parallel we will precompute all the omegas
        // We're using Fr.w[zkey,power+2] and Fr.w[zkey,power+4], we'll only compute the last one and
        // when Fr.w[zkey,power+2] must be used we'll take one out of four.
        FrElement *omegaBuffer;

        BinFileUtils::BinFile *fdZkey;
        BinFileUtils::BinFile *fdWtns;

        Zkey::FflonkZkeyHeader *zkey;
        u_int32_t zkeyPower;
        std::string curveName;
        size_t sDomain;

        G1Point *PTau;

        FrElement *buffWitness;
        FrElement *buffInternalWitness;

        std::map<std::string, FrElement[]> buffers;
        std::map <std::string, Polynomial<Engine>*> polynomials;
        std::map <std::string, Evaluations<Engine>*> evaluations;

        std::map <std::string, FrElement> toInverse;
        std::map <std::string, FrElement> challenges;
        std::map<std::string, FrElement[]> roots;
        FrElement blindingFactors[10];


        SnarkProof<Engine> *proof;
    public:
        FflonkProver(Engine &E);

        ~FflonkProver();

        std::tuple <json, json> prove(BinFileUtils::BinFile *fdZkey, BinFileUtils::BinFile *fdWtns);

        void calculateAdditions(BinFileUtils::BinFile *fdZkey);

        FrElement getWitness(u_int64_t idx);

        void round1();

        void round2();

        void round3();

        void round4();

        void round5();

        //ROUND 1 functions
        void computeWirePolynomials();

        void computeWirePolynomial(std::string polName, FrElement blindingFactors[]);

        void computeT0();

        void computeC1();

        //ROUND 2 functions
        void computeZ();

        void computeT1();

        void computeT2();

        void computeC2();

        //ROUND 4 functions
        void computeR1();

        void computeR2();

        void computeF();

        void computeZT();

        //ROUND 5 functions
        void computeL();

        void computeZTS2();

        FrElement getMontgomeryBatchedInverse();

        G1Point expTau(const FrElement *polynomial, int64_t from, int64_t count);
    };
}

#include "fflonk_prover.cpp"

#endif
