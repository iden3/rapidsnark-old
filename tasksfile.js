const { sh, cli } = require("tasksfile");

function cleanAll() {
    sh("rm -rf build");
}

function downloadGoogleTest() {
    sh("mkdir -p build");
    sh("wget https://github.com/google/googletest/archive/release-1.10.0.tar.gz", {cwd: "build"});
    sh("tar xzf release-1.10.0.tar.gz", {cwd: "build"});
    sh("rm  release-1.10.0.tar.gz", {cwd: "build"});
}

function compileGoogleTest() {
    sh("g++ -Igoogletest -Igoogletest/include -c googletest/src/gtest-all.cc", {cwd: "build/googletest-release-1.10.0"});
    sh("ar -rv libgtest.a gtest-all.o",{cwd: "build/googletest-release-1.10.0"});
}

function createFieldSources() {
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
    sh("cp " + process.argv[3] + " circuit.cpp", {cwd: "build", nopipe: true});
    sh("g++" +
        " -Igoogletest-release-1.10.0/googletest/include"+
        " -I."+
        " -I../src"+
        " -I../depends/pistache/include"+
        " -I../depends/json/single_include"+
        " ../src/main_proofserver.cpp"+
        " ../src/proverapi.cpp"+
        " ../src/fullprover.cpp"+
        " ../src/calcwit.cpp"+
        " ../src/utils.cpp"+
        " ../src/misc.cpp"+
        " ../src/naf.cpp"+
        " ../src/splitparstr.cpp"+
        " fq.cpp"+
        " fq.o"+
        " fr.cpp"+
        " fr.o"+
        " ../src/alt_bn128.cpp"+
        " ../src/binfile_utils.cpp"+
        " ../src/zkey_utils.cpp"+
        " circuit.cpp"+
        " -L../depends/pistache/build/src -lpistache"+
        " -o proverServer"+
        " -fmax-errors=5 -pthread -std=c++17 -fopenmp -lgmp -g -DSANITY_CHECK", {cwd: "build", nopipe: true}
    );
}


function buildProver() {
    sh("g++" +
        " -Igoogletest-release-1.10.0/googletest/include"+
        " -I."+
        " -I../c"+
        " ../src/misc.cpp"+
        " ../src/naf.cpp"+
        " ../src/splitparstr.cpp"+
        " ../src/alt_bn128.cpp"+
        " ../src/main_prover.cpp"+
        " ../src/binfile_utils.cpp"+
        " ../src/zkey_utils.cpp"+
        " ../src/wtns_utils.cpp"+
        " fq.cpp"+
        " fq.o"+
        " fr.cpp"+
        " fr.o"+
        " -o prover" +
        " -fmax-errors=5 -std=c++17 -pthread -lgmp -O3 -fopenmp", {cwd: "build", nopipe: true}
    );
}


cli({
    cleanAll,
    downloadGoogleTest,
    compileGoogleTest,
    createFieldSources,
    buildPistche,
    buildProverServer,
    buildProver
});
