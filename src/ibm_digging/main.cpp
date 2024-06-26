#include "tsd_seasonal.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;

    params.file_name = argv[1];
    params.max_simulation_time = std::stod(argv[2]);
    params.survival_prob[male] = std::stod(argv[3]);
    params.survival_prob[female] = std::stod(argv[4]);
    params.omega_t[male] = std::stod(argv[5]);
    params.omega_t[female] = std::stod(argv[6]);
    params.t_opt[male] = std::stod(argv[7]);
    params.t_opt[female] = std::stod(argv[8]);
    params.skip_output = std::stoi(argv[9]);
    params.temperature_intercept_change = std::stod(argv[10]);
    params.slope = std::stod(argv[11]);
    params.mu_t = std::stod(argv[12]);
    params.mu_depth = std::stod(argv[13]);

    TSDSeasonal sim(params);
}
