#include "tsd_seasonal.hpp"
#include "parameters.hpp"

int main(int argc, char *argv[])
{
    Parameters params;

    params.file_name = argv[1];
    params.max_simulation_time = std::stoul(argv[2]);
    params.simulation_time_change = std::stoul(argv[3]);
    params.survival_prob[male] = std::stod(argv[4]);
    params.survival_prob[female] = std::stod(argv[5]);
    params.omega_t[male] = std::stod(argv[6]);
    params.omega_t[female] = std::stod(argv[7]);
    params.t_opt[male] = std::stod(argv[8]);
    params.t_opt[female] = std::stod(argv[9]);
    params.skip_output = static_cast<unsigned>(std::stoi(argv[10]));
    params.temperature_intercept_change = std::stod(argv[11]);
    params.mu_t = std::stod(argv[12]);
    params.mu_tb = std::stod(argv[13]);
    params.mu_time_threshold = std::stod(argv[14]);
    params.temp_error_sd = std::stod(argv[15]);
    params.cue_error = std::stod(argv[16]);

    TSDSeasonal sim(params);
}
