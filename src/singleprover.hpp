#ifndef SINGLEPROVER_H
#define SINGLEPROVER_H

#include <string>

#include "alt_bn128.hpp"
#include "groth16.hpp"
#include "zkey_utils.hpp"

/** 
 * One-shot prover server
 * 
 * 1. At server initialization, it reads a zkey file and instantiates a prover
 * 2. During its run, it reads a witness file and returns a proof
 * 
 */
class SingleProver {

    mpz_t altBbn128r;
    std::unique_ptr<Groth16::Prover<AltBn128::Engine> > prover;
    std::unique_ptr<ZKeyUtils::Header> zkHeader;
    std::unique_ptr<BinFileUtils::BinFile> zKey;

public: 
    SingleProver(std::string zkeyFileName);
    ~SingleProver();
    json startProve(std::string wtns);
};

#endif // SINGLEPROVER_H
