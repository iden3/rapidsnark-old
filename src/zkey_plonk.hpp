#ifndef ZKEY_PLONK_H
#define ZKEY_PLONK_H

#include <gmp.h>

namespace Zkey {
    const int ZKEY_PL_NSECTIONS = 14;

    const int ZKEY_PL_HEADER_SECTION = 2;
    const int ZKEY_PL_ADDITIONS_SECTION = 3;
    const int ZKEY_PL_A_MAP_SECTION = 4;
    const int ZKEY_PL_B_MAP_SECTION = 5;
    const int ZKEY_PL_C_MAP_SECTION = 6;
    const int ZKEY_PL_QM_SECTION = 7;
    const int ZKEY_PL_QL_SECTION = 8;
    const int ZKEY_PL_QR_SECTION = 9;
    const int ZKEY_PL_QO_SECTION = 10;
    const int ZKEY_PL_QC_SECTION = 11;
    const int ZKEY_PL_SIGMA_SECTION = 12;
    const int ZKEY_PL_LAGRANGE_SECTION = 13;
    const int ZKEY_PL_PTAU_SECTION = 14;

    class PlonkZkeyHeader  {
    public:
        int protocolId;

        u_int32_t n8q;
        mpz_t qPrime;
        u_int32_t n8r;
        mpz_t rPrime;

        u_int32_t nVars;
        u_int32_t nPublic;
        u_int32_t domainSize;
        u_int32_t nAdditions;
        u_int32_t nConstraints;

        void *k1;
        void *k2;
        void *QM;
        void *QL;
        void *QR;
        void *QO;
        void *QC;
        void *S1;
        void *S2;
        void *S3;
        void *X2;

        PlonkZkeyHeader();

        ~PlonkZkeyHeader();

        static PlonkZkeyHeader* loadPlonkZkeyHeader(BinFileUtils::BinFile *f);
    };
}

#endif
