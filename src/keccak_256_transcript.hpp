#ifndef KECCAK_256_TRANSCRIPT_HPP
#define KECCAK_256_TRANSCRIPT_HPP

#include <any>
#include <vector>
#include <sstream>

template<typename Engine>
class Keccak256Transcript {
    using FrElement = typename Engine::FrElement;
    using G1Point = typename Engine::G1Point;
    using G1PointAffine = typename Engine::G1PointAffine;

    enum ElementTypeEnum {
        FrType, G1Type, G1AffineType
    };
    struct ElementTypeStruct {
        ElementTypeEnum type;
        std::any element;
    };

    Engine &E;

    int fieldElements;
    int groupElements;

    std::vector<ElementTypeStruct> elements;

    u_int64_t toRprBE(G1Point &point, uint8_t *data, int64_t seek, int64_t size);
    u_int64_t toRprBE(G1PointAffine &point, uint8_t *data, int64_t seek, int64_t size);

    void hashToFr(FrElement &element, u_int8_t *data, int64_t size);

public:
    Keccak256Transcript(Engine &E);

    void addScalar(FrElement value);

    void addPolCommitment(G1Point value);

    void addPolCommitment(G1PointAffine value);

    void reset();

    typename Engine::FrElement getChallenge();
};

#include "keccak_256_transcript.c.hpp"

#endif
