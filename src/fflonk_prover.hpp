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
#include "keccak_256_transcript.hpp"
#include "wtns_utils.hpp"
#include "zkey.hpp"

using json = nlohmann::json;
using namespace std::chrono;

#define BLINDINGFACTORSLENGTH 10

namespace Fflonk {

    template<typename Engine>
    class FflonkProver {
        using FrElement = typename Engine::FrElement;
        using G1Point = typename Engine::G1Point;
        using G1PointAffine = typename Engine::G1PointAffine;

        Engine &E;
        FFT<typename Engine::Fr> *fft = NULL;

        Zkey::FflonkZkeyHeader *zkey;
        u_int32_t zkeyPower;
        std::string curveName;
        size_t sDomain;

        FrElement *reservedMemoryPtr;
        uint64_t reservedMemorySize;

        FrElement *precomputedBigBuffer;
        G1PointAffine *PTau;

        u_int64_t lengthNonPrecomputedBigBuffer;
        FrElement *nonPrecomputedBigBuffer;

        u_int32_t *mapBuffersBigBuffer;

        FrElement *buffInternalWitness;
        FrElement *buffWitness;

        Zkey::Addition<Engine> *additionsBuff;

        u_int64_t lengthBatchInversesBuffer;

        FrElement *inverses;
        FrElement *products;

        // This is the length of the buffer that must be zeroed after each proof (starting from buffers["A"] pointer)
        u_int64_t buffersLength;

        std::map<std::string, FrElement *> polPtr;
        std::map<std::string, FrElement *> evalPtr;

        std::map<std::string, u_int32_t *> mapBuffers;
        std::map<std::string, FrElement *> buffers;
        std::map<std::string, Polynomial<Engine> *> polynomials;
        std::map<std::string, Evaluations<Engine> *> evaluations;

        std::map <std::string, FrElement> toInverse;
        std::map <std::string, FrElement> challenges;
        std::map<std::string, FrElement *> roots;
        std::map<std::string, FrElement *> omegas;

        FrElement blindingFactors[BLINDINGFACTORSLENGTH];

        Keccak256Transcript<Engine> *transcript;
        SnarkProof<Engine> *proof;
    public:
        FflonkProver(Engine &E);
        FflonkProver(Engine &E, void* reservedMemoryPtr, uint64_t reservedMemorySize);

        ~FflonkProver();

        void setZkey(BinFileUtils::BinFile *fdZkey);

        std::tuple <json, json> prove(BinFileUtils::BinFile *fdZkey, BinFileUtils::BinFile *fdWtns);
        std::tuple <json, json> prove(BinFileUtils::BinFile *fdZkey, FrElement *wtns, WtnsUtils::Header* wtnsHeader = NULL);

        std::tuple <json, json> prove(BinFileUtils::BinFile *fdWtns);
        std::tuple <json, json> prove(FrElement *wtns, WtnsUtils::Header* wtnsHeader = NULL);

    protected:
        void initialize(void* reservedMemoryPtr, uint64_t reservedMemorySize = 0);

        void removePrecomputedData();

        void calculateAdditions();

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
        Polynomial<Engine>* computeRi(uint32_t degree, FrElement* rootsRi, Polynomial<Engine> * polynomialsF);

        std::map<std::string, Polynomial<Engine> *> computeR(std::map<std::string, FrElement *> rootsMap, std::map<std::string, Polynomial<Engine> *> polynomialsF);

        Polynomial<Engine>* computeW(FrElement challengeAlpha, std::map<std::string, Polynomial<Engine> *> polynomialsF, std::map<std::string, Polynomial<Engine> *> polynomialsR);

        Polynomial<Engine>* computeZT(std::map<std::string, FrElement *> rootsMap);

        //ROUND 5 functions
        Polynomial<Engine>* computeWp(FrElement challengeY, FrElement challengeAlpha, std::map<std::string, FrElement *> rootsMap, std::map<std::string, Polynomial<Engine> *> polynomialsF, std::map<std::string, Polynomial<Engine> *> polynomialsR, Polynomial<Engine> * polynomialW);

        Polynomial<Engine>* computeL(FrElement challengeY, FrElement challengeAlpha, std::map<std::string, FrElement *> rootsMap, std::map<std::string, Polynomial<Engine> *> polynomialsF, std::map<std::string, Polynomial<Engine> *> polynomialsR, Polynomial<Engine> * polynomialW);

        Polynomial<Engine>* computeZTS2(std::map<std::string, FrElement *> rootsMap);

        void batchInverse(FrElement *elements, u_int64_t length);

        FrElement *polynomialFromMontgomery(Polynomial<Engine> *polynomial);

        FrElement getMontgomeryBatchedInverse(FrElement challengeXi, FrElement challengeY, std::map<std::string, FrElement *> rootsMap);

        void computeLiSinglePoint(uint32_t index, uint32_t degree, FrElement challengesY, std::map<std::string, FrElement *> rootsMap);

        void computeLiS2();

        G1Point multiExponentiation(Polynomial<Engine> *polynomial);

        G1Point multiExponentiation(Polynomial<Engine> *polynomial, u_int32_t nx, u_int64_t x[]);

        // Helpers

        FrElement computeChallengeXiSeed(FrElement previousChallenge);

        std::map<std::string, FrElement *> calculateRoots(FrElement challengeXiSeed);

        FrElement computeChallengeAlpha(FrElement challengeXiSeed);
        
        FrElement computeChallengeY(FrElement challengeAlpha, G1Point W);
        
        void open(FrElement previousChallenge);

    };
}

#include "fflonk_prover.c.hpp"

#endif