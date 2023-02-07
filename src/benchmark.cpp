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
        FrElement op1 = E.fr.set(123);
        FrElement op2 = E.fr.set(247);

        double time = omp_get_wtime();
        for (uint64_t i = 0; i < n; i++) {
            op1 = E.fr.mul(op1, op2);
        }
        time = omp_get_wtime() - time;
        cout << time << "\n";
    }

    template<typename Engine>
    void Benchmark<Engine>::run(std::string ptauFile, int initialPower, int finalPower, int iterations) {
        std::cout << "BENCHMARK STARTED\n";

        dump = new Dump::Dump<Engine>(E);
        u_int64_t finalDomainSize = std::pow(2, finalPower);
        fft = new FFT<typename Engine::Fr>(finalDomainSize);

        double resultsIfft[finalPower - initialPower];
        double resultsFft[finalPower - initialPower];
        double resultsMultiexp[finalPower - initialPower];

//        auto zkey = BinFileUtils::openExisting(zkeyFilename, "zkey", 1);
//
//        PTau = new G1PointAffine[finalDomainSize];
//        memset(PTau, 0, sizeof(PTau));
//
//        ThreadUtils::parcpy(this->PTau,
//                            (G1PointAffine *) fdZkey->getSectionData(Zkey::ZKEY_FF_PTAU_SECTION),
//                            (finalDomainSize) * sizeof(G1PointAffine), nThreads);

        //Benchmark IFFT's
        for (int power = initialPower; power <= finalPower; power++) {
            cout << "power: " << power << "\n";
            u_int64_t domainSize = std::pow(2, power);
            cout << "ifft " << power << "\n";

            resultsIfft[power - initialPower] = benchmarkIfft(domainSize, iterations);

            cout << "fft " << power << "\n";
            resultsFft[power - initialPower] = benchmarkFft(domainSize, iterations);

//            cout << "mexp " << power << "\n";
//            resultsMultiexp[power - initialPower] = benchmarkMultiexp(domainSize, iterations);
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

//        ss.str("");
//        ss << "\nBENCHMARK MULTIEXP\n";
//        for(int power = initialPower; power <= finalPower; power++) {
//            ss << "2^" << power << "\t" << resultsMultiexp[power - initialPower] << "\n";
//        }
//
//        std::cout << ss.str();

        std::cout << "\nBENCHMARK FINISHED\n";
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
    double Benchmark<Engine>::benchmarkFft(u_int32_t domainSize, int iterations) {
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
            bestTime = std::min(bestTime, time);
        }

        delete[] buffer;
        delete[] bufferCoef;
        delete[] bufferEval;

        return bestTime;
    }

    template<typename Engine>
    double Benchmark<Engine>::benchmarkMultiexp(u_int32_t domainSize, int iterations) {
        FrElement *buffer = new FrElement[domainSize];
        FrElement *bufferCoef = new FrElement[domainSize];

        createBuffer(buffer, domainSize);

        double bestTime = numeric_limits<double>::max();
        double time;
        for (int i = 0; i < iterations; i++) {
            Polynomial<Engine> *polynomial = Polynomial<Engine>::fromEvaluations(E, fft, buffer, bufferCoef,
                                                                                 domainSize);
            time = omp_get_wtime();
            typename Engine::G1Point commit = multiExponentiation(polynomial);
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
    typename Engine::G1Point Benchmark<Engine>::multiExponentiation(Polynomial<Engine> *polynomial) {
        G1Point value;
        cout << "1\n";
        this->polynomialFromMontgomery(polynomial);
        cout << "2\n";

        E.g1.multiMulByScalar(value, PTau, (uint8_t *) polynomial, sizeof(polynomial[0]), polynomial->getDegree() + 1);

        cout << "3\n";
        return value;
    }

}
