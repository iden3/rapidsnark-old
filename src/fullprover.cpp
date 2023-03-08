#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "fullprover.hpp"
#include "fr.hpp"

#include "logger.hpp"
#include "wtns_utils.hpp"

using namespace CPlusPlusLogging;

std::string getfilename(std::string path)
{
    path = path.substr(path.find_last_of("/\\") + 1);
    size_t dot_i = path.find_last_of('.');
    return path.substr(0, dot_i);
}

FullProver::FullProver(std::string zkeyFileNames[], int size) {
    pendingInput="";
    pendingCircuit="";
    canceled = false;

    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    for(int i = 0; i < size; i++) {
        std::string circuit = getfilename(zkeyFileNames[i]);
        zKeys[circuit] = BinFileUtils::openExisting(zkeyFileNames[i], "zkey", 1);
        zkHeaders[circuit] = ZKeyUtils::loadHeader(zKeys[circuit].get());

        std::string proofStr;
        if (mpz_cmp(zkHeaders[circuit]->rPrime, altBbn128r) != 0) {
            throw std::invalid_argument( "zkey curve not supported" );
        }
        
        std::ostringstream ss1;
        ss1 << "circuit: " << circuit;
        LOG_DEBUG(ss1);

        provers[circuit] = Groth16::makeProver<AltBn128::Engine>(
            zkHeaders[circuit]->nVars,
            zkHeaders[circuit]->nPublic,
            zkHeaders[circuit]->domainSize,
            zkHeaders[circuit]->nCoefs,
            zkHeaders[circuit]->vk_alpha1,
            zkHeaders[circuit]->vk_beta1,
            zkHeaders[circuit]->vk_beta2,
            zkHeaders[circuit]->vk_delta1,
            zkHeaders[circuit]->vk_delta2,
            zKeys[circuit]->getSectionData(4),    // Coefs
            zKeys[circuit]->getSectionData(5),    // pointsA
            zKeys[circuit]->getSectionData(6),    // pointsB1
            zKeys[circuit]->getSectionData(7),    // pointsB2
            zKeys[circuit]->getSectionData(8),    // pointsC
            zKeys[circuit]->getSectionData(9)     // pointsH1
        );
    }

    status = ready;
}

FullProver::~FullProver() {
    mpz_clear(altBbn128r);
}

void FullProver::startProve(std::string input, std::string circuit) {
    LOG_TRACE("FullProver::startProve begin");
    LOG_DEBUG(input);
    std::lock_guard<std::mutex> guard(mtx);
    pendingInput = input;
    pendingCircuit = circuit;
    if (status == busy) {
        abort();
    }
    checkPending();
    LOG_TRACE("FullProver::startProve end");
}

void FullProver::checkPending() {
    LOG_TRACE("FullProver::checkPending begin");
    if (status != busy) {
        std::string input = pendingInput;
        std::string circuit = pendingCircuit;
        if (input != "" && circuit != "") {
            status = busy;
            executingInput = pendingInput;
            executingCircuit = pendingCircuit;
            pendingInput = "";
            pendingCircuit = "";
            errString = "";
            canceled = false;
            proof = nlohmann::detail::value_t::null;
            std::thread th(&FullProver::thread_calculateProve, this);
            th.detach();
        }
    }
    LOG_TRACE("FullProver::checkPending end");
}

void FullProver::thread_calculateProve() {
    LOG_TRACE("FullProver::thread_calculateProve start");
    
    try {
        LOG_TRACE(executingInput);
        // Generate witness
        json j = json::parse(executingInput);
        std::string circuit = executingCircuit;
        
        std::ofstream file("./build/input_"+ circuit +".json");
        file << j;
        file.close();

        std::string witnessFile("./build/" + circuit + ".wtns");
        std::string command("./build/" + circuit + " ./build/input_"+ circuit +".json " + witnessFile);
        LOG_TRACE(command);
        std::array<char, 128> buffer;
        std::string result;

        // std::cout << "Opening reading pipe" << std::endl;
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe)
        {
            std::cerr << "Couldn't start command." << std::endl;
        }
        while (fgets(buffer.data(), 128, pipe) != NULL) {
            // std::cout << "Reading..." << std::endl;
            result += buffer.data();
        }
        auto returnCode = pclose(pipe);

        std::cout << result << std::endl;
        std::cout << returnCode << std::endl;
        
        // Load witness
        auto wtns = BinFileUtils::openExisting(witnessFile, "wtns", 2);
        auto wtnsHeader = WtnsUtils::loadHeader(wtns.get());
                
        if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
            throw std::invalid_argument( "different wtns curve" );
        }

        AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns->getSectionData(2);

        pubData.clear();
        AltBn128::FrElement aux;
        for (int i=1; i<=zkHeaders[circuit]->nPublic; i++) {
            AltBn128::Fr.toMontgomery(aux, wtnsData[i]);
            pubData.push_back(AltBn128::Fr.toString(aux));
        }
        
        if (!isCanceled()) {
            proof = provers[circuit]->prove(wtnsData)->toJson();
        } else {
            LOG_TRACE("AVOIDING prove");
            proof = {};
        }
       

        calcFinished();
    } catch (std::runtime_error e) {
        if (!isCanceled()) {
            errString = e.what();
        }
        calcFinished();
    } 

    LOG_TRACE("FullProver::thread_calculateProve end");
}


void FullProver::calcFinished() {
    std::lock_guard<std::mutex> guard(mtx);
    LOG_TRACE("FullProver::calcFinished start");
    if (canceled) {
        LOG_TRACE("FullProver::calcFinished aborted");
        status = aborted;
    } else if (errString != "") {
        LOG_TRACE("FullProver::calcFinished failed");
        status = failed;
    } else {
        LOG_TRACE("FullProver::calcFinished success");
        status = success;
    }
    canceled = false;
    executingInput = "";
    checkPending();
    LOG_TRACE("FullProver::calcFinished end");
}


bool FullProver::isCanceled() {
    std::lock_guard<std::mutex> guard(mtx);
    LOG_TRACE("FullProver::isCanceled start");
    if (canceled) {
        LOG_TRACE("FullProver::isCanceled canceled==true");
    }
    LOG_TRACE("FullProver::isCanceled end");
    return canceled;
}

void FullProver::abort() {
    std::lock_guard<std::mutex> guard(mtx);
    LOG_TRACE("FullProver::abort start");
    if (status!= busy) {
        LOG_TRACE("FullProver::abort end -> not usy");
        return;
    }
    canceled = true;
    LOG_TRACE("FullProver::abort end -> canceled=true");
}


json FullProver::getStatus() {
    LOG_TRACE("FullProver::getStatus start");
    json st;
    if (status == ready) {
        LOG_TRACE("ready");
        st["status"] = "ready";
    } else if (status == aborted) {
        LOG_TRACE("aborted");
        st["status"] = "aborted";
    } else if (status == failed) {
        LOG_TRACE("failed");
        st["status"] = "failed";
        st["error"] = errString;
    } else if (status == success) {
        LOG_TRACE("success");
        st["status"] = "success";
        st["proof"] = proof.dump();
        st["pubData"] = pubData.dump();
    } else if (status == busy) {
        LOG_TRACE("busy");
        st["status"] = "busy";
    }
    LOG_TRACE("FullProver::getStatus end");
    return st;
}
