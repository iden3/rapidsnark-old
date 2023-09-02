#include "proverapi.hpp"
#include "nlohmann/json.hpp"
#include "logger.hpp"

using namespace Pistache;
using json = nlohmann::json;

void ProverAPI::postInput(const Rest::Request& request, Http::ResponseWriter response) {
    try {
        json j = prover.startProve(request.body());
        LOG_DEBUG(j.dump().c_str());
        response.send(Http::Code::Ok, j.dump(), MIME(Application, Json));
    } catch (const std::exception &e) {
        auto errString = e.what();
        LOG_ERROR(errString);
        response.send(Http::Code::Bad_Request, errString);
    }
}
