#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "zkey_utils.hpp"

#include "fullprover.hpp"
#include "fr.hpp"

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
    std::lock_guard<std::mutex> guard(mtx);
    pendingInput = input;
    if (status == busy) {
        abort();
    }
    checkPending();
}


void FullProver::checkPending() {
    if (status != busy) {
        std::string input = pendingInput;
        if (input != "") {
            executingInput = pendingInput;
            pendingInput = "";
            status = busy;
            errString = "";
            canceled = false;
            proof = nlohmann::detail::value_t::null;
            std::thread th(&FullProver::thread_calculateProve, this);
            th.detach();
        }
    }
}

void FullProver::thread_calculateProve() {
    try {
        calcWit->calculateProve(wtns, executingInput, [this](){ return isCanceled(); });
        pubData.clear();
        for (int i=1; i<=circuit->NPublic; i++) {
            AltBn128::FrElement aux;
            AltBn128::Fr.toMontgomery(aux, wtns[i]);
            pubData.push_back(AltBn128::Fr.toString(aux));
        }
        proof = prover->prove(wtns)->toJson();

        calcFinished();
    } catch (std::runtime_error e) {
        if (!isCanceled()) {
            errString = e.what();
        }
        calcFinished();
    }
}


void FullProver::calcFinished() {
    std::lock_guard<std::mutex> guard(mtx);
    if (canceled) {
        status = aborted;
    } else if (errString != "") {
        status = failed;
    } else {
        status = success;
    }
    canceled = false;
    executingInput = "";
    checkPending();
}


bool FullProver::isCanceled() {
    std::lock_guard<std::mutex> guard(mtx);
    return canceled;
}

void FullProver::abort() {
    std::lock_guard<std::mutex> guard(mtx);
    if (status!= busy) return;
    canceled = true;
}


json FullProver::getStatus() {
    json st;
    if (status == ready) {
        st["status"] = "ready";
    } else if (status == aborted) {
        st["status"] = "aborted";
    } else if (status == failed) {
        st["status"] = "failed";
        st["error"] = errString;
    } else if (status == success) {
        st["status"] = "success";
        st["proof"] = proof;
        st["pubData"] = pubData;
    } else if (status == busy) {
        st["status"] = "busy";
    }
    return st;
}


#define ADJ_P(a) *((void **)&a) = (void *)(((char *)circuit)+ (uint64_t)(a))

Circom_Circuit *FullProver::loadCircuit(std::string const &datFileName) {
    Circom_Circuit *circuitF;
    Circom_Circuit *circuit;

    int fd;
    struct stat sb;

    fd = open(datFileName.c_str(), O_RDONLY);
    if (fd == -1) {
        std::cout << ".dat file not found: " << datFileName << "\n";
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

    return circuit;
}



