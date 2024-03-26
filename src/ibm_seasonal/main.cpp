#include "tsd_seasonal.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;

    params.file_name = argv[1];
    params.max_simulation_time = std::stod(argv[2]);
    params.survival_prob[male] = std::stod(argv[3]);
    params.survival_prob[female] = std::stod(argv[4]);
    params.t_opt[female] = std::stod(argv[5]);
    params.frequency = std::stod(argv[6]);

    TSDSeasonal sim(params);
}
