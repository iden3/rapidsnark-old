#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <gmp.h>
#include <memory>
#include <stdexcept>
#include <nlohmann/json.hpp>

#include <alt_bn128.hpp>
#include "binfile_utils.hpp"
#include "zkey_utils.hpp"
#include "wtns_utils.hpp"
#include "groth16.hpp"
#include "fflonk_prover.hpp"
#include "zkey.hpp"
#include "logger.hpp"
#include "benchmark.hpp"

using namespace CPlusPlusLogging;

using json = nlohmann::json;

#define handle_error(msg) \
           do { perror(msg); exit(EXIT_FAILURE); } while (0)

mpz_t altBbn128r;

void groth16Prove(BinFileUtils::BinFile *zkey, BinFileUtils::BinFile *wtns, std::string proofFilename, std::string publicFilename);
void fflonkProve(BinFileUtils::BinFile *zkey, BinFileUtils::BinFile *wtns, std::string proofFilename, std::string publicFilename);

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "Invalid number of parameters:\n";
        std::cerr << "Usage: prove <circuit.zkey> <witness.wtns> <proof.json> <public.json>\n";
        return -1;
    }

    Logger::getInstance()->enableConsoleLogging();
    Logger::getInstance()->updateLogLevel(LOG_LEVEL_DEBUG);

    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    try {
        std::string zkeyFilename = argv[1];
        std::string wtnsFilename = argv[2];
        std::string proofFilename = argv[3];
        std::string publicFilename = argv[4];

        LOG_TRACE("> Opening zkey file");
        auto zkey = BinFileUtils::openExisting(zkeyFilename, "zkey", 1);

        LOG_TRACE("> Opening wtns file");
        auto wtns = BinFileUtils::openExisting(wtnsFilename, "wtns", 2);

        const int protocolId = Zkey::getProtocolIdFromZkey(zkey.get());

        switch(protocolId) {
            case Zkey::FFLONK_PROTOCOL_ID:
                LOG_TRACE("> FFLONK protocol detected");
                fflonkProve(zkey.get(), wtns.get(), proofFilename, publicFilename);
                break;
            case Zkey::GROTH16_PROTOCOL_ID:
                groth16Prove(zkey.get(), wtns.get(), proofFilename, publicFilename);
                LOG_TRACE("> GROTH16 protocol detected");
                break;
            default:
                throw std::invalid_argument("protocol not supported");
        }
    } catch (std::exception& e) {
        mpz_clear(altBbn128r);
        std::cerr << e.what() << '\n';
        return -1;
    }
    mpz_clear(altBbn128r);
    exit(EXIT_SUCCESS);
}

void fflonkProve(BinFileUtils::BinFile *zkey, BinFileUtils::BinFile *wtns, std::string proofFilename, std::string publicFilename) {
    auto prover = new Fflonk::FflonkProver<AltBn128::Engine>(AltBn128::Engine::engine);

    auto [proofJson, publicSignalsJson] = prover->prove(zkey, wtns);

    std::ofstream file;
    file.open(proofFilename);
    file << proofJson;
    file.close();

    file.open(publicFilename);
    file << publicSignalsJson;
    file.close();
}

void groth16Prove(BinFileUtils::BinFile *zkey, BinFileUtils::BinFile *wtns, std::string proofFilename, std::string publicFilename) {
    auto zkeyHeader = ZKeyUtils::loadHeader(zkey);

    std::string proofStr;
    if (mpz_cmp(zkeyHeader->rPrime, altBbn128r) != 0) {
        throw std::invalid_argument("zkey curve not supported");
    }

    auto wtnsHeader = WtnsUtils::loadHeader(wtns);

    if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
        throw std::invalid_argument("different wtns curve");
    }

    auto prover = Groth16::makeProver<AltBn128::Engine>(
            zkeyHeader->nVars,
            zkeyHeader->nPublic,
            zkeyHeader->domainSize,
            zkeyHeader->nCoefs,
            zkeyHeader->vk_alpha1,
            zkeyHeader->vk_beta1,
            zkeyHeader->vk_beta2,
            zkeyHeader->vk_delta1,
            zkeyHeader->vk_delta2,
            zkey->getSectionData(4),    // Coefs
            zkey->getSectionData(5),    // pointsA
            zkey->getSectionData(6),    // pointsB1
            zkey->getSectionData(7),    // pointsB2
            zkey->getSectionData(8),    // pointsC
            zkey->getSectionData(9)     // pointsH1
    );
    AltBn128::FrElement *wtnsData = (AltBn128::FrElement *) wtns->getSectionData(2);
    auto proof = prover->prove(wtnsData);

    std::ofstream proofFile;
    proofFile.open(proofFilename);
    proofFile << proof->toJson();
    proofFile.close();

    std::ofstream publicFile;
    publicFile.open(publicFilename);

    json jsonPublic;
    AltBn128::FrElement aux;
    for (u_int32_t i = 1; i <= zkeyHeader->nPublic; i++) {
        AltBn128::Fr.toMontgomery(aux, wtnsData[i]);
        jsonPublic.push_back(AltBn128::Fr.toString(aux));
    }

    publicFile << jsonPublic;
    publicFile.close();
}
