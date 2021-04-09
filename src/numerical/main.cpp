#include "tsd_multivariate.hpp"

int main(int argc, char **argv)
{
    TSD_Multivariate simulation{};

    simulation.init_arguments(argc, argv);

    simulation.run();
}
