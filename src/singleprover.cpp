#include <singleprover.hpp>
#include "fr.hpp"

#include "logger.hpp"
#include "wtns_utils.hpp"

SingleProver::SingleProver(std::string zkeyFileName) {

    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    zKey = BinFileUtils::openExisting(zkeyFileName, "zkey", 1);
    zkHeader = ZKeyUtils::loadHeader(zKey.get());

    std::string proofStr;
    if (mpz_cmp(zkHeader->rPrime, altBbn128r) != 0) {
        throw std::invalid_argument( "zkey curve not supported" );
    }
    
    prover = Groth16::makeProver<AltBn128::Engine>(
        zkHeader->nVars,
        zkHeader->nPublic,
        zkHeader->domainSize,
        zkHeader->nCoefs,
        zkHeader->vk_alpha1,
        zkHeader->vk_beta1,
        zkHeader->vk_beta2,
        zkHeader->vk_delta1,
        zkHeader->vk_delta2,
        zKey->getSectionData(4),    // Coefs
        zKey->getSectionData(5),    // pointsA
        zKey->getSectionData(6),    // pointsB1
        zKey->getSectionData(7),    // pointsB2
        zKey->getSectionData(8),    // pointsC
        zKey->getSectionData(9)     // pointsH1
    );

    LOG_INFO("SingleProver::SingleProver initialized from zkey");
}

SingleProver::~SingleProver()
{
    mpz_clear(altBbn128r);
}

json SingleProver::startProve(std::string input)
{
    LOG_INFO("SingleProver::startProve begin");

    // TODO: Do not log PII. prover->prove also logs the ZK proof, which we might not want to log.
    LOG_DEBUG(input);

    json j = json::parse(input);
    std::ofstream file("./build/input.json");
    file << j;
    file.close();

    std::string witnessFile("./build/witness.wtns");
    std::string command("./build/zkLogin ./build/input.json " + witnessFile);
    LOG_INFO(command);
    std::array<char, 128> buffer;
    std::string result;

    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe)
    {
        throw std::runtime_error("Couldn't start command");
    }
    while (fgets(buffer.data(), 128, pipe) != NULL) {
        result += buffer.data();
    }
    auto returnCode = pclose(pipe);

    if (result != "") {
        LOG_INFO("Unexpected result");
        LOG_INFO(result);
    }

    if (returnCode != 0) {
        LOG_INFO("Unexpected return code");
        auto str = std::to_string(returnCode);
        LOG_INFO(str);
    }

    auto wtns = BinFileUtils::openExisting(witnessFile, "wtns", 2);
    auto wtnsHeader = WtnsUtils::loadHeader(wtns.get());
    if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
        throw std::invalid_argument( "different wtns curve" );
    }

    AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns->getSectionData(2);
    auto proof = prover->prove(wtnsData);

    return proof->toJson();
}
