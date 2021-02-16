#include <pistache/router.h>
#include <pistache/endpoint.h>

#include "proverapi.hpp"
#include "fullprover.hpp"
#include "logger.hpp"

using namespace CPlusPlusLogging;
using namespace Pistache;
using namespace Pistache::Rest;

int main(int argc, char **argv) {

    Logger::getInstance()->enableConsoleLogging();
    Logger::getInstance()->updateLogLevel(LOG_LEVEL_DEBUG);
    LOG_INFO("Initializing server...");

    FullProver fullProver(argv[1], argv[2]);
    ProverAPI proverAPI(fullProver);
    Address addr(Ipv4::any(), Port(9080));

    auto opts = Http::Endpoint::options().threads(1).maxRequestSize(128000000);
    Http::Endpoint server(addr);
    server.init(opts);
    Router router;
    Routes::Get(router, "/status", Routes::bind(&ProverAPI::getStatus, &proverAPI));
    Routes::Post(router, "/start", Routes::bind(&ProverAPI::postStart, &proverAPI));
    Routes::Post(router, "/stop", Routes::bind(&ProverAPI::postStop, &proverAPI));
    Routes::Post(router, "/input", Routes::bind(&ProverAPI::postInput, &proverAPI));
    Routes::Post(router, "/cancel", Routes::bind(&ProverAPI::postCancel, &proverAPI));
    server.setHandler(router.handler());
    LOG_INFO("Server ready on port 9080...");
    server.serve();
}
