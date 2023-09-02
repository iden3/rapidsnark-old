#include <pistache/router.h>
#include <pistache/endpoint.h>
#include "proverapi.hpp"
#include "singleprover.hpp"
#include "logger.hpp"

using namespace CPlusPlusLogging;
using namespace Pistache;
using namespace Pistache::Rest;

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Invalid number of parameters:\n";
        std::cerr << "Usage: proverServer <port> <circuit.zkey> \n";
        return -1;
    }

    Logger::getInstance()->enableConsoleLogging();
    Logger::getInstance()->updateLogLevel(LOG_LEVEL_DEBUG);
    LOG_INFO("Initializing server...");
    int port = std::stoi(argv[1]); // parse port
    // parse the zkeys
    std::string zkeyFileName = argv[2];

    SingleProver prover(zkeyFileName);
    ProverAPI proverAPI(prover);
    Address addr(Ipv4::any(), Port(port));

    auto opts = Http::Endpoint::options().threads(1).maxRequestSize(128000000);
    Http::Endpoint server(addr);
    server.init(opts);
    Router router;
    Routes::Post(router, "/input", Routes::bind(&ProverAPI::postInput, &proverAPI));
    server.setHandler(router.handler());
    std::string serverReady("Server ready on port " + std::to_string(port) + "...");
    LOG_INFO(serverReady);
    server.serve();
}
