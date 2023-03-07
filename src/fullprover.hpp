#ifndef FULLPROVER_H
#define FULLPROVER_H

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <mutex>
#include "alt_bn128.hpp"
#include "groth16.hpp"
#include "binfile_utils.hpp"
#include "zkey_utils.hpp"

class FullProver {
    enum Status {aborted = -2, busy = -1,  failed = 0, success = 1, unverified =2, uninitialized=3, initializing=5, ready=6 };
    Status status;
    std::mutex mtx;

    std::string pendingInput;
    std::string executingInput;
    std::string pendingCircuit;
    std::string executingCircuit;

    std::map<std::string, std::unique_ptr<Groth16::Prover<AltBn128::Engine>>> provers;
    std::map<std::string, std::unique_ptr<ZKeyUtils::Header>> zkHeaders;
    std::map<std::string, std::unique_ptr<BinFileUtils::BinFile>> zKeys;

    mpz_t altBbn128r;

    json proof;
    json pubData;
    std::string errString;

    bool canceled;

    bool isCanceled();
    void calcFinished();
    void thread_calculateProve();
    void checkPending();



public: 
    FullProver(std::string zkeyFileNames[], int size);
    ~FullProver();
    void startProve(std::string input, std::string circuit);
    void abort();
    json getStatus();
    std::string &getErrString() { return errString; };


};

#endif // FULLPROVER_H