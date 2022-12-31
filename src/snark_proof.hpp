#ifndef SNARK_PROOF_HPP
#define SNARK_PROOF_HPP

#include <map>
#include <string>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

template<typename Engine>
class SnarkProof {
    using FrElement = typename Engine::FrElement;
    using G1PointAffine = typename Engine::G1PointAffine;

    Engine &E;
    std::string protocol;
    std::string curve;

    std::map <std::string, G1PointAffine> polynomialCommitments;
    std::map <std::string, FrElement> evaluationCommitments;
public :
    SnarkProof(const std::string protocol);

    void resetProof();

    void addPolynomialCommitment(std::string &key, G1PointAffine &polynomialCommitment);

    G1PointAffine getPolynomialCommitment(std::string &key);

    void addEvaluationCommitment(std::string &key, FrElement &evaluationCommitment);

    FrElement getEvaluationCommitment(std::string &key);

    std::string toJsonStr();

    json toJson();
};

#endif //SNARK_PROOF_HPP
