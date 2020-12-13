#include <pistache/router.h>
#include <pistache/endpoint.h>
#include "fullprover.hpp"

using namespace Pistache;

class ProverAPI {
    FullProver &fullProver;
public:
    ProverAPI(FullProver &_fullProver) : fullProver(_fullProver) {};
    void postStart(const Rest::Request& request, Http::ResponseWriter response);
    void postStop(const Rest::Request& request, Http::ResponseWriter response);
    void postInput(const Rest::Request& request, Http::ResponseWriter response);
    void postCancel(const Rest::Request& request, Http::ResponseWriter response);
    void getStatus(const Rest::Request& request, Http::ResponseWriter response);

    void getConfig(const Rest::Request& request, Http::ResponseWriter response);
    void postConfig(const Rest::Request& request, Http::ResponseWriter response);

};