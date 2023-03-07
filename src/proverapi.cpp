#include "proverapi.hpp"
#include "nlohmann/json.hpp"
#include "logger.hpp"

using namespace Pistache;
using json = nlohmann::json;


void ProverAPI::postInput(const Rest::Request& request, Http::ResponseWriter response) {
    std::string circuit(request.param(":circuit").as<std::string>());
    LOG_TRACE(circuit);
    fullProver.startProve(request.body(), circuit);
    response.send(Http::Code::Ok);
}

void ProverAPI::postCancel(const Rest::Request& request, Http::ResponseWriter response) {
    fullProver.abort();
    response.send(Http::Code::Ok);
}

void ProverAPI::getStatus(const Rest::Request& request, Http::ResponseWriter response) {
    json j = fullProver.getStatus();
    LOG_DEBUG(j.dump().c_str());
    response.send(Http::Code::Ok, j.dump(), MIME(Application, Json));
}

void ProverAPI::postStart(const Rest::Request& request, Http::ResponseWriter response) {
    response.send(Http::Code::Ok);
}

void ProverAPI::postStop(const Rest::Request& request, Http::ResponseWriter response) {
    response.send(Http::Code::Ok);
}

void ProverAPI::getConfig(const Rest::Request& request, Http::ResponseWriter response) {
    response.send(Http::Code::Ok);
}

void ProverAPI::postConfig(const Rest::Request& request, Http::ResponseWriter response) {
    response.send(Http::Code::Ok);
}

