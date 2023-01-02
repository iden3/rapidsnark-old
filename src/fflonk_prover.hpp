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
#include "polynomial/evaluation.hpp"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Fflonk {

    template<typename Engine>
    class FflonkProver {
        using FrElement = typename Engine::FrElement;
        using FrElements = typename Engine::FrElement[];
        using G1Point = typename Engine::G1Point;

        Engine &E;
        FFT<typename Engine::Fr> *fft;

        Zkey::FflonkZkeyHeader *zkey;
        u_int32_t zkeyPower;
        std::string curveName;

        G1Point *PTau;

        FrElement *buffWitness;
        FrElements *buffInternalWitness;

        std::map<std::string, FrElement[]> buffers;
        std::map<std::string, Polynomial<Engine>> polynomials;
        std::map<std::string, Evaluation<Engine>> evaluations;

        std::map <std::string, FrElement> toInverse;
        std::map <std::string, FrElement> challenges;
        FrElements *blindingFactors;
        std::map<std::string, FrElement[]> roots;


        SnarkProof<Engine> *proof;
    public:
        FflonkProver();

        ~FflonkProver();

        std::tuple<json, json> prove(BinFileUtils::BinFile *fdZkey, BinFileUtils::BinFile *fdWtns);

        void calculateAdditions(BinFileUtils::BinFile *fdZkey);

        FrElement getWitness(u_int64_t idx);

        void round1();
        void round2();
        void round3();
        void round4();
        void round5();

        FrElement getMontgomeryBatchedInverse();
    };
}

#endif
