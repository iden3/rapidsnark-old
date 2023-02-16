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
    void Benchmark<Engine>::benchmarkMultiply(uint64_t n) {
        std::cout << "BENCHMARK MULTIPLY STARTED\n";
        FrElement op1 = E.fr.set(123);
        FrElement op2 = E.fr.set(247);

        double time = omp_get_wtime();
        for (uint64_t i = 0; i < n; i++) {
            op1 = E.fr.mul(op1, op2);
        }
        time = omp_get_wtime() - time;
        cout << time << "\n";
        std::cout << "BENCHMARK MULTIPLY FINISHED\n";
    }

    template<typename Engine>
    void Benchmark<Engine>::benchmarkFft(int initialPower, int finalPower, int iterations) {
        std::cout << "BENCHMARK FFT STARTED\n";

        dump = new Dump::Dump<Engine>(E);
        u_int64_t finalDomainSize = std::pow(2, finalPower);
        fft = new FFT<typename Engine::Fr>(finalDomainSize);

        double resultsIfft[finalPower - initialPower];
        double resultsFft[finalPower - initialPower];

        //Benchmark IFFT's
        for (int power = initialPower; power <= finalPower; power++) {
            cout << "power: " << power << "\n";
            u_int64_t domainSize = std::pow(2, power);
            cout << "ifft " << power << "\n";

            resultsIfft[power - initialPower] = runIfft(domainSize, iterations);

            cout << "fft " << power << "\n";
            resultsFft[power - initialPower] = runFft(domainSize, iterations);
        }
        std::ostringstream ss;

        // PRINT IFFT RESULTS
        ss << "\nBENCHMARK IFFT\n";
        for(int power = initialPower; power <= finalPower; power++) {
            ss << "2^" << power << "\t" << resultsIfft[power - initialPower] << "\n";
        }
        std::cout << ss.str();

        ss.str("");
        ss << "\nBENCHMARK FFT\n";
        for (int power = initialPower; power <= finalPower; power++) {
            ss << "2^" << power << "\t" << resultsFft[power - initialPower] << "\n";
        }
        std::cout << ss.str();

        std::cout << "\nBENCHMARK FFT FINISHED\n";
    }

    template<typename Engine>
    void Benchmark<Engine>::benchmarkMultiexp(std::string ptauFile, int initialPower, int finalPower, int iterations) {
        std::cout << "BENCHMARK MULTIEXP STARTED\n";

        dump = new Dump::Dump<Engine>(E);
        u_int64_t finalDomainSize = std::pow(2, finalPower);
        fft = new FFT<typename Engine::Fr>(finalDomainSize);

        double resultMexp[finalPower - initialPower];
        auto fdZkey = BinFileUtils::openExisting(ptauFile, "zkey", 2);

        G1PointAffine* PTau = new G1PointAffine[finalDomainSize];
        memset(PTau, 0, sizeof(G1PointAffine) * finalDomainSize);

        int nThreads = omp_get_max_threads() / 2;
        ThreadUtils::parcpy(PTau,
                            (G1PointAffine *) fdZkey->getSectionData(Zkey::ZKEY_FF_PTAU_SECTION),
                            (finalDomainSize) * sizeof(G1PointAffine), nThreads);

        //Benchmark Multiexponentiations
        for (int power = initialPower; power <= finalPower; power++) {
            cout << "power: " << power << "\n";
            u_int64_t domainSize = std::pow(2, power);
            cout << "ifft " << power << "\n";

            resultMexp[power - initialPower] = runMultiexp(&PTau[0], domainSize, iterations);
        }
        std::ostringstream ss;

        ss << "\nBENCHMARK MULTIEXP\n";
        for(int power = initialPower; power <= finalPower; power++) {
            ss << "2^" << power << "\t" << resultMexp[power - initialPower] << "\n";
        }

        std::cout << ss.str();

        std::cout << "\nBENCHMARK MULTIEXP FINISHED\n";
    }

    template<typename Engine>
    void Benchmark<Engine>::createBuffer(FrElement *buffer, u_int32_t domainSize) {
        buffer[0] = E.fr.one();
        buffer[1] = E.fr.one();
        for (uint64_t i = 2; i < domainSize; i++) {
            buffer[i] = E.fr.add(buffer[i - 1], buffer[i - 2]);
        }
    }

    template<typename Engine>
    double Benchmark<Engine>::runIfft(u_int32_t domainSize, int iterations) {
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
            delete polynomial;
            bestTime = std::min(bestTime, time);
        }

        delete[] buffer;
        delete[] bufferCoef;

        return bestTime;
    }

    template<typename Engine>
    double Benchmark<Engine>::runFft(u_int32_t domainSize, int iterations) {
        FrElement *buffer = new FrElement[domainSize];
        FrElement *bufferCoef = new FrElement[domainSize];
        FrElement *bufferEval = new FrElement[domainSize];

        createBuffer(buffer, domainSize);

        double bestTime = numeric_limits<double>::max();
        double time;

        for (int i = 0; i < iterations; i++) {
            Polynomial<Engine> *polynomial = Polynomial<Engine>::fromEvaluations(E, fft, buffer, bufferCoef,
                                                                                 domainSize);
            time = omp_get_wtime();
            Evaluations<Engine> *evaluations = new Evaluations<Engine>(E, fft, bufferEval, *polynomial, domainSize);
            time = omp_get_wtime() - time;
            delete evaluations;
            bestTime = std::min(bestTime, time);
        }

        delete[] buffer;
        delete[] bufferCoef;
        delete[] bufferEval;

        return bestTime;
    }

    template<typename Engine>
    double Benchmark<Engine>::runMultiexp(G1PointAffine* ptau, u_int32_t domainSize, int iterations) {
        FrElement *buffer = new FrElement[domainSize];
        FrElement *bufferCoef = new FrElement[domainSize];

        createBuffer(buffer, domainSize);

        double bestTime = numeric_limits<double>::max();
        double time;
        for (int i = 0; i < iterations; i++) {
            Polynomial<Engine> *polynomial = Polynomial<Engine>::fromEvaluations(E, fft, buffer, bufferCoef,
                                                                                 domainSize);
            time = omp_get_wtime();
            typename Engine::G1Point commit = multiExponentiation(ptau, polynomial);
            dump->dump("...", commit);
            time = omp_get_wtime() - time;
            bestTime = std::min(bestTime, time);
        }

        delete[] buffer;
        delete[] bufferCoef;

        return bestTime;
    }

    template<typename Engine>
    typename Engine::FrElement *Benchmark<Engine>::polynomialFromMontgomery(Polynomial<Engine> *polynomial) {
        const u_int64_t length = polynomial->getLength();

#pragma omp parallel for
        for (u_int32_t index = 0; index < length; ++index) {
            E.fr.fromMontgomery(polynomial->coef[index], polynomial->coef[index]);
        }

        return polynomial->coef;
    }

    template<typename Engine>
    typename Engine::G1Point Benchmark<Engine>::multiExponentiation(G1PointAffine* ptau, Polynomial<Engine> *polynomial) {
        G1Point value;

        FrElement *pol = this->polynomialFromMontgomery(polynomial);

        E.g1.multiMulByScalar(value, ptau, (uint8_t *) pol, sizeof(pol[0]), polynomial->getDegree() + 1);

        return value;
    }

}
