# rapidsnark

rapid snark is a zkSnark proof generation written in C++ and intel assembly. That generates proofs created in [circom](https://github.com/iden3/snarkjs) and [snarkjs](https://github.com/iden3/circom) very fast.

## dependencies

You should have installed gcc, libsodium, and gmp (development)

In ubuntu:

````
sudo apt install build-essential
sudo apt-get install libgmp-dev
sudo apt-get install libsodium-dev
sudo apt-get install nasm
````

## compile prover

````sh
npm install
git submodule init
git submodule update
npx task createFieldSources
npx task buildPistche
````

## compile prover in stand alone mode
````
npx task buildProver <circuit_cpp location>
````
## compile prover in server mode
````
npx task buildProverServer <circuit_cpp location>
````

## Lunch prover in server mode
````
./rapidsnark/build/proverServer  <circuit_dat> <circuit_zkey>
````

## Building proof

You have a full prover compiled in the build directory.

So you can replace snarkjs command:

````sh
snarkjs groth16 prove <circuit.zkey> <witness.wtns> <proof.json> <public.json>
````

by this one
````sh
./build/prove <circuit.zkey> <witness.wtns> <proof.json> <public.json>
````

## Benchmark

This prover uses intel assembly with ADX extensions and parallelizes as much as it can the proof generation. 

The prover is much faster that snarkjs and faster than bellman.

[TODO] Some comparation tests should be done.


## License

rapidsnark is part of the iden3 project copyright 2021 0KIMS association and published with GPL-3 license. Please check the COPYING file for more details.
