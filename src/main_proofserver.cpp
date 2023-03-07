#include <pistache/router.h>
#include <pistache/endpoint.h>
#include "proverapi.hpp"
#include "fullprover.hpp"
#include "logger.hpp"

using namespace CPlusPlusLogging;
using namespace Pistache;
using namespace Pistache::Rest;

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Invalid number of parameters:\n";
        std::cerr << "Usage: proverServer <port> <circuit1.zkey> <circuit2.zkey> ... <circuitN.zkey> \n";
        return -1;
    }

    Logger::getInstance()->enableConsoleLogging();
    Logger::getInstance()->updateLogLevel(LOG_LEVEL_DEBUG);
    LOG_INFO("Initializing server...");
    int port = std::stoi(argv[1]); // parse port
    // parse the zkeys
    std::string zkeyFileNames[argc - 2];
    for (int i = 0; i < argc - 2; i++) {
        zkeyFileNames[i] = argv[i + 2];
    }

    FullProver fullProver(zkeyFileNames, argc - 2);
    ProverAPI proverAPI(fullProver);
    Address addr(Ipv4::any(), Port(port));

    auto opts = Http::Endpoint::options().threads(1).maxRequestSize(128000000);
    Http::Endpoint server(addr);
    server.init(opts);
    Router router;
    Routes::Get(router, "/status", Routes::bind(&ProverAPI::getStatus, &proverAPI));
    Routes::Post(router, "/start", Routes::bind(&ProverAPI::postStart, &proverAPI));
    Routes::Post(router, "/stop", Routes::bind(&ProverAPI::postStop, &proverAPI));
    Routes::Post(router, "/input/:circuit", Routes::bind(&ProverAPI::postInput, &proverAPI));
    Routes::Post(router, "/cancel", Routes::bind(&ProverAPI::postCancel, &proverAPI));
    server.setHandler(router.handler());
    std::string serverReady("Server ready on port " + std::to_string(port) + "...");
    LOG_INFO(serverReady);
    server.serve();
}
