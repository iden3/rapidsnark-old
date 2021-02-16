#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "zkey_utils.hpp"

#include "fullprover.hpp"
#include "fr.hpp"

#include "logger.hpp"
using namespace CPlusPlusLogging;

FullProver::FullProver(std::string datFileName, std::string zkeyFileName) {
    pendingInput="";
    wtns=NULL;
    canceled = false;

    circuit = loadCircuit(datFileName);

    // open output
    calcWit = new Circom_CalcWit(circuit);

    zkey = BinFileUtils::openExisting(zkeyFileName, "zkey", 1);
    auto zkeyHeader = ZKeyUtils::loadHeader(zkey.get());

    prover = Groth16::makeProver<AltBn128::Engine>(
        zkeyHeader->nVars,
        zkeyHeader->nPublic,
        zkeyHeader->domainSize,
        zkeyHeader->nCoefs,
        zkeyHeader->vk_alpha1,
        zkeyHeader->vk_beta1,
        zkeyHeader->vk_beta2,
        zkeyHeader->vk_delta1,
        zkeyHeader->vk_delta2,
        zkey->getSectionData(4),    // Coefs
        zkey->getSectionData(5),    // pointsA
        zkey->getSectionData(6),    // pointsB1
        zkey->getSectionData(7),    // pointsB2
        zkey->getSectionData(8),    // pointsC
        zkey->getSectionData(9)     // pointsH1
    );

    wtns = new AltBn128::FrElement[circuit->NVars];

    status = ready;
}

FullProver::~FullProver() {
    delete calcWit;
    delete wtns;
}

void FullProver::startProve(std::string input) {
    LOG_TRACE("FullProver::startProve begin");
    LOG_DEBUG(input);
    std::lock_guard<std::mutex> guard(mtx);
    pendingInput = input;
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
        if (input != "") {
            status = busy;
            executingInput = pendingInput;
            pendingInput = "";
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
        calcWit->calculateProve(wtns, executingInput, [this](){ return /* TODO isCanceled() */ false; });
        pubData.clear();
        LOG_TRACE("FullProver::thread_calculateProve calculating prove");
        for (int i=1; i<=circuit->NPublic; i++) {
            AltBn128::FrElement aux;
            AltBn128::Fr.toMontgomery(aux, wtns[i]);
            pubData.push_back(AltBn128::Fr.toString(aux));
        }
        if (!isCanceled()) {
            proof = prover->prove(wtns)->toJson();
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


#define ADJ_P(a) *((void **)&a) = (void *)(((char *)circuit)+ (uint64_t)(a))

Circom_Circuit *FullProver::loadCircuit(std::string const &datFileName) {
    LOG_TRACE("FullProver::loadCircuit start");
    Circom_Circuit *circuitF;
    Circom_Circuit *circuit;

    int fd;
    struct stat sb;

    fd = open(datFileName.c_str(), O_RDONLY);
    if (fd == -1) {
        std::ostringstream ss;
        ss << ".dat file not found: " << datFileName << "\n";
        LOG_ERROR(ss);
        throw std::system_error(errno, std::generic_category(), "open");
    }

    if (fstat(fd, &sb) == -1) {         /* To obtain file size */
        throw std::system_error(errno, std::generic_category(), "fstat");
    }

    circuitF = (Circom_Circuit *)mmap(NULL, sb.st_size, PROT_READ , MAP_PRIVATE, fd, 0);
    close(fd);

    circuit = (Circom_Circuit *)malloc(sb.st_size);
    memcpy((void *)circuit, (void *)circuitF, sb.st_size);

    munmap(circuitF, sb.st_size);

    ADJ_P(circuit->wit2sig);
    ADJ_P(circuit->components);
    ADJ_P(circuit->mapIsInput);
    ADJ_P(circuit->constants);
    ADJ_P(circuit->P);
    ADJ_P(circuit->componentEntries);

    for (int i=0; i<circuit->NComponents; i++) {
        ADJ_P(circuit->components[i].hashTable);
        ADJ_P(circuit->components[i].entries);
        circuit->components[i].fn = _functionTable[  (uint64_t)circuit->components[i].fn];
    }

    for (int i=0; i<circuit->NComponentEntries; i++) {
        ADJ_P(circuit->componentEntries[i].sizes);
    }

    LOG_TRACE("FullProver::loadCircuit end");
    return circuit;
}



