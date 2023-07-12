const { sh, cli } = require("tasksfile");

function cleanAll() {
    sh("rm -rf build");
}

function createFieldSources() {
    sh("mkdir -p build");
    sh("npm install", {cwd: "depends/ffiasm"});
    sh("node ../depends/ffiasm/src/buildzqfield.js -q 21888242871839275222246405745257275088696311157297823662689037894645226208583 -n Fq", {cwd: "build"});
    sh("node ../depends/ffiasm/src/buildzqfield.js -q 21888242871839275222246405745257275088548364400416034343698204186575808495617 -n Fr", {cwd: "build"});

    if (process.platform === "darwin") {
        sh("nasm -fmacho64 --prefix _ fq.asm", {cwd: "build"});
    }  else if (process.platform === "linux") {
        sh("nasm -felf64 fq.asm", {cwd: "build"});
    } else throw("Unsupported platform");

    if (process.platform === "darwin") {
        sh("nasm -fmacho64 --prefix _ fr.asm", {cwd: "build"});
    }  else if (process.platform === "linux") {
        sh("nasm -felf64 fr.asm", {cwd: "build"});
    } else throw("Unsupported platform");
}

function buildPistche() {
    sh("git submodule init && git submodule update");
    sh("mkdir -p build", {cwd: "depends/pistache"});
    sh("cmake -G \"Unix Makefiles\" -DCMAKE_BUILD_TYPE=Release ..", {cwd: "depends/pistache/build"});
    sh("make", {cwd: "depends/pistache/build"});
}


function buildProverServer() {
    sh("cp " + process.argv[3] + " build/circuit.cpp", {cwd: ".", nopipe: true});
    sh("g++" +
        " -I."+
        " -I../src"+
        " -I../depends/pistache/include"+
        " -I../depends/json/single_include"+
        " -I../depends/ffiasm/c"+
        " -I../depends/circom_runtime/c"+
        " ../src/main_proofserver.cpp"+
        " ../src/curve_utils.cpp"+
        " ../src/proverapi.cpp"+
        " ../src/fullprover.cpp"+
        " ../src/binfile_utils.cpp"+
        " ../src/zkey_utils.cpp"+
        " ../src/logger.cpp"+
        " ../depends/circom_runtime/c/calcwit.cpp"+
        " ../depends/circom_runtime/c/utils.cpp"+
        " ../depends/ffiasm/c/misc.cpp"+
        " ../depends/ffiasm/c/naf.cpp"+
        " ../depends/ffiasm/c/splitparstr.cpp"+
        " ../depends/ffiasm/c/alt_bn128.cpp"+
        " fq.cpp"+
        " fq.o"+
        " fr.cpp"+
        " fr.o"+
        " circuit.cpp"+
        " -L../depends/pistache/build/src -lpistache"+
        " -o proverServer"+
        " -fmax-errors=5 -pthread -std=c++17 -fopenmp -lgmp -lsodium -g -DSANITY_CHECK", {cwd: "build", nopipe: true}
    );
}


function buildProver() {
    sh("g++" +
        " -I."+
        " -I../src"+
        " -I../depends/ffiasm/c"+
        " -I../depends/json/single_include"+
        " ../src/main_prover.cpp"+
        " ../src/binfile_utils.cpp"+
        " ../src/zkey_utils.cpp"+
        " ../src/zkey.cpp"+
        " ../src/zkey_fflonk.cpp"+
        " ../src/curve_utils.cpp"+
        " ../src/wtns_utils.cpp"+
        " ../src/keccak_wrapper.cpp"+
        " ../src/logger.cpp"+
        " ../depends/ffiasm/c/misc.cpp"+
        " ../depends/ffiasm/c/naf.cpp"+
        " ../depends/ffiasm/c/splitparstr.cpp"+
        " ../depends/ffiasm/c/alt_bn128.cpp"+
        " fq.cpp"+
        " fq.o"+
        " fr.cpp"+
        " fr.o"+
        " -o prover" +
        " -Wall"+
        " -fmax-errors=5 -std=c++17 -pthread -lgmp -lsodium -O3 -fopenmp", {cwd: "build", nopipe: true}
    );
}

function buildTest() {
    sh("g++" +
        " -I."+
        " -I../src"+
        " -I../depends/ffiasm/c"+
        " -I../depends/json/single_include"+
        " ../test/main_test.cpp"+
        " ../test/polynomial.test.cpp"+
        " ../src/binfile_utils.cpp"+
        " ../src/zkey_utils.cpp"+
        " ../src/zkey.cpp"+
        " ../src/zkey_fflonk.cpp"+
        " ../src/curve_utils.cpp"+
        " ../src/wtns_utils.cpp"+
        " ../src/keccak_wrapper.cpp"+
        " ../src/logger.cpp"+
        " ../depends/ffiasm/c/misc.cpp"+
        " ../depends/ffiasm/c/naf.cpp"+
        " ../depends/ffiasm/c/splitparstr.cpp"+
        " ../depends/ffiasm/c/alt_bn128.cpp"+
        " fq.cpp"+
        " fq.o"+
        " fr.cpp"+
        " fr.o"+
        " -o proverTest" +
        " -Wall"+
        " -fmax-errors=5 -std=c++17 -pthread -lgmp -lsodium -O3 -fopenmp -lgtest", {cwd: "build", nopipe: true}
    );
}

function compile() {
    sh("g++ -c" +
        " -I."+
        " -I../src"+
        " -I/opt/homebrew/Cellar/gmp/6.2.1_1/include/"+
        " -I/opt/homebrew/opt/libomp/include"+
        " -I/opt/homebrew/Cellar/libsodium/1.0.18_1/include"+
        " -I../depends/ffiasm/c"+
        " -I../depends/json/single_include"+
        " ../src/main_prover.cpp"+
        " ../src/binfile_utils.cpp"+
        " ../src/zkey_utils.cpp"+
        " ../src/zkey.cpp"+
        " ../src/zkey_fflonk.cpp"+
        " ../src/curve_utils.cpp"+
        " ../src/wtns_utils.cpp"+
        " ../src/keccak_wrapper.cpp"+
        " ../src/logger.cpp"+
        " ../depends/ffiasm/c/misc.cpp"+
        " ../depends/ffiasm/c/naf.cpp"+
        " ../depends/ffiasm/c/splitparstr.cpp"+
        " ../depends/ffiasm/c/alt_bn128.cpp"+
        " fq.cpp"+
        " fr.cpp"+
        " -Wall"+
        " -fmax-errors=20 -std=c++17 -pthread -O3 -lssl -lcrypto -fopenmp", {cwd: "build", nopipe: true}
    );
}

cli({
    cleanAll,
    createFieldSources,
    buildPistche,
    buildProverServer,
    buildProver,
    buildTest,
    compile
})
