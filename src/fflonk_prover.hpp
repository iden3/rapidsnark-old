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
#include "dump.hpp"

using json = nlohmann::json;

namespace Fflonk {

    template<typename Engine>
    class FflonkProver {
        using FrElement = typename Engine::FrElement;
        using G1Point = typename Engine::G1Point;
        using G1PointAffine = typename Engine::G1PointAffine;

        Dump::Dump<Engine> *dump;

        Engine &E;
        FFT<typename Engine::Fr> *fft;
        MulZ<Engine> *mulZ;

        BinFileUtils::BinFile *fdZkey;
        BinFileUtils::BinFile *fdWtns;

        Zkey::FflonkZkeyHeader *zkey;
        u_int32_t zkeyPower;
        std::string curveName;
        size_t sDomain;

        G1PointAffine *PTau;

        FrElement *buffWitness;
        FrElement *buffInternalWitness;

        std::map<std::string, u_int32_t*> mapBuffers;
        std::map<std::string, FrElement*> buffers;
        std::map <std::string, Polynomial<Engine>*> polynomials;
        std::map <std::string, Evaluations<Engine>*> evaluations;

        std::map <std::string, FrElement> toInverse;
        std::map <std::string, FrElement> challenges;
        std::map<std::string, FrElement*> roots;
        FrElement blindingFactors[10];
        G1Point C1, C2, W1, W2;

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

        void batchInverse(FrElement *elements, u_int64_t length);

        FrElement* polynomialFromMontgomery(Polynomial<Engine> *polynomial);

        G1Point multiExponentiation(Polynomial<Engine> *polynomial);

        void printPol(std::string name, const Polynomial<Engine> *polynomial);
    };
}

#include "fflonk_prover.cpp"

#endif
