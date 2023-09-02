#include <pistache/router.h>
#include <pistache/endpoint.h>
#include "singleprover.hpp"

using namespace Pistache;

class ProverAPI {
    SingleProver &prover;
public:
    ProverAPI(SingleProver &_prover) : prover(_prover) {};
    void postInput(const Rest::Request& request, Http::ResponseWriter response);
};