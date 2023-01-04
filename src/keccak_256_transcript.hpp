#ifndef KECCAK_256_TRANSCRIPT_HPP
#define KECCAK_256_TRANSCRIPT_HPP

template<typename Engine>
class Keccak256Transcript {
public:
    using FrElement = typename Engine::FrElement;
    using G1Point = typename Engine::G1Point;

    Keccak256Transcript();

    void addScalar(FrElement value);

    void addPolCommitment(G1Point value);

    void reset();

    typename Engine::FrElement getChallenge();
};


#endif
