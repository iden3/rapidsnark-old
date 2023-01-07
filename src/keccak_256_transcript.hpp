#ifndef KECCAK_256_TRANSCRIPT_HPP
#define KECCAK_256_TRANSCRIPT_HPP

#include <any>
#include <vector>
#include <sstream>

template<typename Engine>
class Keccak256Transcript {
    using FrElement = typename Engine::FrElement;
    using G1Point = typename Engine::G1Point;

    int fieldElements;
    int groupElements;

    std::vector<std::any> elements;
public:
    Keccak256Transcript();

    void addScalar(FrElement value);

    void addPolCommitment(G1Point value);

    void reset();

    typename Engine::FrElement getChallenge();
};

#endif
