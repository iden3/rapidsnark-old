# rapidsnark


## compile prover

````sh
npm install
git submodule init
git submodule update
npx task createFieldSources
npx task buildProver
````

## Building proof

You have a full prover compiled in the build directory.

So you can replace snarkjs command:

````sh
snarkjs groth16 prove <circuit.zkey> <witness.wtns> <proof.json>
````

by this one
````sh
./build/prove <circuit.zkey> <witness.wtns> <proof.json>
````

## Benchmark

This prover uses intel assembly with ADX extensions and parallelizes as much as it can the proof generation. 

The prover is much faster that snarkjs and faster than bellman.

[TODO] Some comparation tests should be done.


## License

rapidsnark is part of the iden3 project copyright 2021 0KIMS association and published with GPL-3 license. Please check the COPYING file for more details.