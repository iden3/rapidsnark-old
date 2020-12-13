#ifndef FULLPROVER_H
#define FULLPROVER_H

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "alt_bn128.hpp"
#include "groth16.hpp"
#include "calcwit.hpp"
#include "binfile_utils.hpp"

class FullProver {
    enum Status {aborted = -2, busy = -1,  failed = 0, success = 1, unverified =2, uninitialized=3, initializing=5, ready=6 };
    Status status;
    std::mutex mtx;

    std::string pendingInput;
    std::string executingInput;

    json proof;
    json pubData;
    std::string errString;

    Circom_Circuit *circuit;
    Circom_CalcWit *calcWit;

    AltBn128::FrElement *wtns;
    bool canceled;

    std::unique_ptr<BinFileUtils::BinFile> zkey;

    std::unique_ptr<Groth16::Prover<AltBn128::Engine>> prover;

    Circom_Circuit *loadCircuit(std::string const &datFileName);
    bool isCanceled();
    void calcFinished();
    void thread_calculateProve();
    void checkPending();



public: 
    FullProver(std::string datFileName, std::string zkeyFileName);
    ~FullProver();
    void startProve(std::string input);
    void abort();
    json getStatus();
    std::string &getErrString() { return errString; };


};

#endif // FULLPROVER_H