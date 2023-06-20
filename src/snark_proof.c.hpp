#include <sstream>

#include "snark_proof.hpp"
#include "logger.hpp"
#include "curve_utils.hpp"

using namespace CPlusPlusLogging;

template<typename Engine>
SnarkProof<Engine>::SnarkProof(Engine &_E, const std::string &protocol) : E(_E) {
    this->protocol = protocol;
    this->curve = CurveUtils::getCurveNameByEngine();
    this->reset();
}

template<typename Engine>
void SnarkProof<Engine>::reset() {
    this->polynomialCommitments.clear();
    this->evaluationCommitments.clear();
}

template<typename Engine>
void SnarkProof<Engine>::addPolynomialCommitment(const std::string &key, G1Point &polynomialCommmitment) {
    if (0 != polynomialCommitments.count(key)) {
        std::ostringstream ss;
        ss << "!!! SnarkProof::addPolynomialCommitment. '" << key << "' already exist in proof";
        LOG_ALARM(ss);
    }
    this->polynomialCommitments[key] = polynomialCommmitment;
}

template<typename Engine>
typename Engine::G1Point SnarkProof<Engine>::getPolynomialCommitment(const std::string &key) {
    if (0 == polynomialCommitments.count(key)) {
        std::ostringstream ss;
        ss << "!!! SnarkProof::addPolynomialCommitment. '" << key << "' does not exist in proof";
        LOG_ALARM(ss);
    }
    return this->polynomialCommitments[key];
}

template<typename Engine>
void SnarkProof<Engine>::addEvaluationCommitment(const std::string &key, FrElement evaluationCommitment) {
    if (0 != evaluationCommitments.count(key)) {
        std::ostringstream ss;
        ss << "!!! SnarkProof::addPolynomialCommitment. '" << key << "' already exist in proof";
        LOG_ALARM(ss);
    }
    this->evaluationCommitments[key] = evaluationCommitment;
}

template<typename Engine>
typename Engine::FrElement SnarkProof<Engine>::getEvaluationCommitment(const std::string &key) {
    if (0 == evaluationCommitments.count(key)) {
        std::ostringstream ss;
        ss << "!!! SnarkProof::addPolynomialCommitment. '" << key << "' does not exist in proof";
        LOG_ALARM(ss);
    }
    return this->evaluationCommitments[key];
}

template<typename Engine>
json SnarkProof<Engine>::toJson() {
    json jsonProof;

    jsonProof["polynomials"] = {};
    jsonProof["evaluations"] = {};

    for (auto &[key, point]: this->polynomialCommitments) {
        G1PointAffine tmp;
        E.g1.copy(tmp, point);

        jsonProof["polynomials"][key] = {};

        std::string x = E.f1.toString(tmp.x);
        std::string y = E.f1.toString(tmp.y);

        jsonProof["polynomials"][key].push_back(x);
        jsonProof["polynomials"][key].push_back(y);
        jsonProof["polynomials"][key].push_back("1");
    }

    for (auto &[key, element]: this->evaluationCommitments) {
        jsonProof["evaluations"][key] = E.fr.toString(element);
    }

    jsonProof["protocol"] = this->protocol;
    jsonProof["curve"] = this->curve;

    return jsonProof;
}

template<typename Engine>
uint32_t SnarkProof<Engine>::getNumberEvaluations() {
    return this->evaluationCommitments.size();
}


template<typename Engine>
std::string* SnarkProof<Engine>::getEvaluationsNames() {
    std::cout << this->evaluationCommitments.size() << std::endl;
    std::string* evaluationsNames = new std::string[this->evaluationCommitments.size()];

    uint32_t index = 0;
    for (auto &[key, element]: this->evaluationCommitments) {
        evaluationsNames[index++] = key;
    }

    std::cout << evaluationsNames->size() << std::endl;
    for (uint32_t i = 0; i < this->evaluationCommitments.size(); ++i) {
        std::cout << evaluationsNames[i] << std::endl;
    }

    return evaluationsNames;
}