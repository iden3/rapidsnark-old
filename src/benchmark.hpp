#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include <string>
#include <map>
#include "binfile_utils.hpp"
#include <gmp.h>
#include "fft.hpp"
#include "zkey_fflonk.hpp"
#include "polynomial/polynomial.hpp"
#include "polynomial/evaluations.hpp"
#include "dump.hpp"

using json = nlohmann::json;

namespace Benchmark {

    template<typename Engine>
    class Benchmark {
        using FrElement = typename Engine::FrElement;
        using G1Point = typename Engine::G1Point;
        using G1PointAffine = typename Engine::G1PointAffine;

        Dump::Dump<Engine> *dump;

        Engine &E;
        FFT<typename Engine::Fr> *fft = NULL;

        typename Engine::FrElement* polynomialFromMontgomery(Polynomial<Engine> *polynomial);

        typename Engine::G1Point multiExponentiation(G1PointAffine* ptau, Polynomial<Engine> *polynomial);

        void createBuffer(FrElement *buffer, u_int32_t domainSize);

        double runIfft(u_int32_t domainSize, int iterations);

        double runFft(u_int32_t domainSize, int iterations);

        double runMultiexp(G1PointAffine* ptau, u_int32_t domainSize, int iterations);

    public:
        Benchmark(Engine &E);

        ~Benchmark();

        void benchmarkMultiply(uint64_t n);

        void benchmarkFft(int initialPower, int finalPower, int iterations);

        void benchmarkMultiexp(std::string ptauFile, int initialPower, int finalPower, int iterations);
    };
}

#include "benchmark.cpp"

#endif
