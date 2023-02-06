#include "benchmark.hpp"

#include "curve_utils.hpp"
#include "zkey.hpp"
#include "zkey_fflonk.hpp"
#include "wtns_utils.hpp"
#include <sodium.h>
#include "thread_utils.hpp"
#include "polynomial/cpolynomial.hpp"
#include "polynomial/cpolynomial.hpp"

namespace Benchmark {

    template<typename Engine>
    Benchmark<Engine>::Benchmark(Engine &_E) :E(_E) {
    }

    template<typename Engine>
    Benchmark<Engine>::~Benchmark() {
        if (NULL != fft) {
            delete fft;
        }
    }

    template<typename Engine>
    void Benchmark<Engine>::run(int initialPower, int finalPower) {
        std::cout << "BENCHMARK STARTED";

        dump = new Dump::Dump<Engine>(E);
        u_int64_t finalDomainSize = std::pow(2, finalPower);
        fft = new FFT<typename Engine::Fr>(finalDomainSize);

        std::ostringstream ss;

        //Benchmark IFFT's
        for (int power = initialPower; power <= finalPower; power++) {
            u_int64_t domainSize = std::pow(2, power);

            ss << "BENCHMARK ifft's domain size: " << domainSize << "\n";
            std::cout << ss.str() << "\n";

            double time = benchmarkIfft(domainSize, 1);
            ss << "BENCHMARK result: " << time << "\n";
            std::cout << ss.str() << "\n";

            for (int extension = 1; extension <= 16; extension = extension * 2) {
                ss << "BENCHMARK fft's domain size: " << domainSize << "   extension: " << extension << "\n";
                std::cout << ss.str() << "\n";

                double time = benchmarkFft(domainSize, extension, 1);
                ss << "BENCHMARK result: " << time << "\n";
                std::cout << ss.str() << "\n";
            }
        }

        std::cout << "BENCHMARK FINISHED";
    }

    template<typename Engine>
    void Benchmark<Engine>::createBuffer(FrElement *buffer, u_int32_t domainSize) {
        buffer[0] = E.fr.one();
        buffer[1] = E.fr.one();
        for (int i = 2; i < domainSize; i++) {
            buffer[i] = E.fr.add(buffer[i - 1], buffer[i - 2]);
        }
    }

    template<typename Engine>
    double Benchmark<Engine>::benchmarkIfft(u_int32_t domainSize, int iterations) {
        FrElement *buffer = new FrElement[domainSize];
        FrElement *bufferCoef = new FrElement[domainSize];

        createBuffer(buffer, domainSize);

        double bestTime = numeric_limits<double>::max();
        double time;
        for (int i = 0; i < iterations; i++) {
            time = omp_get_wtime();
            Polynomial<Engine> *polynomial = Polynomial<Engine>::fromEvaluations(E, fft, buffer, bufferCoef,
                                                                                 domainSize);
            time = omp_get_wtime() - time;
            bestTime = std::min(bestTime, time);
        }

        delete[] buffer;
        delete[] bufferCoef;

        return bestTime;
    }

    template<typename Engine>
    double Benchmark<Engine>::benchmarkFft(u_int32_t domainSize, int extension, int iterations) {
        FrElement *buffer = new FrElement[domainSize];
        FrElement *bufferCoef = new FrElement[domainSize];
        FrElement *bufferEval = new FrElement[domainSize * extension];

        createBuffer(buffer, domainSize);

        double bestTime = numeric_limits<double>::max();
        double time;
        for (int i = 0; i < iterations; i++) {
            Polynomial<Engine> *polynomial = Polynomial<Engine>::fromEvaluations(E, fft, buffer, bufferCoef,
                                                                                 domainSize);

            time = omp_get_wtime();
            Evaluations<Engine> *evaluations = new Evaluations<Engine>(E, fft, bufferEval, *polynomial,
                                                                       domainSize * extension);
            time = omp_get_wtime() - time;
            bestTime = std::min(bestTime, time);
        }

        delete[] buffer;
        delete[] bufferCoef;
        delete[] bufferEval;

        return bestTime;
    }
}
