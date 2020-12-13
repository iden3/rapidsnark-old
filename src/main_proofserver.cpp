#include <pistache/router.h>
#include <pistache/endpoint.h>

#include "proverapi.hpp"
#include "fullprover.hpp"

using namespace Pistache;
using namespace Pistache::Rest;

int main(int argc, char **argv) {
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
    cout << "Server ready...\n";
    server.serve();
}
