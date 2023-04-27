#include <stdexcept>

#include "zkey.hpp"
#include "zkey_plonk.hpp"

namespace Zkey {
    PlonkZkeyHeader::PlonkZkeyHeader() {
        this->protocolId = Zkey::PLONK_PROTOCOL_ID;
    }

    PlonkZkeyHeader::~PlonkZkeyHeader() {
        mpz_clear(qPrime);
        mpz_clear(rPrime);
    }

    PlonkZkeyHeader* PlonkZkeyHeader::loadPlonkZkeyHeader(BinFileUtils::BinFile *f) {
        auto plonkZkeyHeader = new PlonkZkeyHeader();

        f->startReadSection(Zkey::ZKEY_PL_HEADER_SECTION);

        plonkZkeyHeader->n8q = f->readU32LE();
        mpz_init(plonkZkeyHeader->qPrime);
        mpz_import(plonkZkeyHeader->qPrime, plonkZkeyHeader->n8q, -1, 1, -1, 0, f->read(plonkZkeyHeader->n8q));

        plonkZkeyHeader->n8r = f->readU32LE();
        mpz_init(plonkZkeyHeader->rPrime);
        mpz_import(plonkZkeyHeader->rPrime, plonkZkeyHeader->n8r, -1, 1, -1, 0, f->read(plonkZkeyHeader->n8r));

        plonkZkeyHeader->nVars = f->readU32LE();
        plonkZkeyHeader->nPublic = f->readU32LE();
        plonkZkeyHeader->domainSize = f->readU32LE();
        plonkZkeyHeader->nAdditions = f->readU32LE();
        plonkZkeyHeader->nConstraints = f->readU32LE();

        plonkZkeyHeader->k1 = f->read(plonkZkeyHeader->n8r);
        plonkZkeyHeader->k2 = f->read(plonkZkeyHeader->n8r);

        plonkZkeyHeader->QM = f->read(plonkZkeyHeader->n8q * 2);
        plonkZkeyHeader->QL = f->read(plonkZkeyHeader->n8q * 2);
        plonkZkeyHeader->QR = f->read(plonkZkeyHeader->n8q * 2);
        plonkZkeyHeader->QO = f->read(plonkZkeyHeader->n8q * 2);
        plonkZkeyHeader->QC = f->read(plonkZkeyHeader->n8q * 2);

        plonkZkeyHeader->S1 = f->read(plonkZkeyHeader->n8q * 2);
        plonkZkeyHeader->S2 = f->read(plonkZkeyHeader->n8q * 2);
        plonkZkeyHeader->S3 = f->read(plonkZkeyHeader->n8q * 2);

        plonkZkeyHeader->X2 = f->read(plonkZkeyHeader->n8q * 4);

        f->endReadSection();

        return plonkZkeyHeader;
    }
}

