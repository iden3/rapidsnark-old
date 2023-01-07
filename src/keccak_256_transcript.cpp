#include "keccak_256_transcript.hpp"


template<typename Engine>
Keccak256Transcript<Engine>::Keccak256Transcript() {
    reset();
}

template<typename Engine>
void Keccak256Transcript<Engine>::addScalar(FrElement value) {
    this->elements.push_back(value);
    fieldElements++;
}

template<typename Engine>
void Keccak256Transcript<Engine>::addPolCommitment(G1Point value) {
    this->elements.push_back(value);
    groupElements++;
}

template<typename Engine>
void Keccak256Transcript<Engine>::reset() {
    fieldElements = 0;
    groupElements = 0;

    this->elements.clear();
}

template<typename Engine>
typename Engine::FrElement Keccak256Transcript<Engine>::getChallenge() {
    u_int8_t data[Engine::engine.fr.bytes(2) * fieldElements + Engine::engine.g1.F.bytes() * groupElements * 2 * 3];
    u_int64_t bytes = 0;

    memset(data, 0, sizeof(data));

    for (int i = 0; i < this->elements.size(); i++) {
        bytes += toRprBE(elements[i], data, bytes, sizeof(elements[i]));
    }

    FrElement res;
    hashToFr(res, data, bytes);

    return res;
}
