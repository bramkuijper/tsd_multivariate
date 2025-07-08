#include "tsd_multivariate.hpp"

int main(int argc, char **argv)
{
    
    parameters.surv[Male][0] = std::atof(argv[1]);
    parameters.surv[Male][1] = std::atof(argv[2]);
    parameters.surv[Female][0] = std::atof(argv[3]);
    parameters.surv[Female][1] = std::atof(argv[4]);
    parameters.init_d[Female][0] = std::atof(argv[5]);
    parameters.init_d[Male][0] = std::atof(argv[6]);
    parameters.init_d[Female][1] = std::atof(argv[7]);
    parameters.init_d[Male][1] = std::atof(argv[8]);
    parameters.init_b = std::atof(argv[9]);
    parameters.init_s[0] = std::atof(argv[10]);
    parameters.init_s[1] = std::atof(argv[11]);

    parameters.sigma[0][1] = std::atof(argv[12]);
    parameters.sigma[0][0] = 1.0 - sigma[0][1];

    parameters.sigma[1][0] = std::atof(argv[13]);
    parameters.sigma[1][1] = 1.0 - sigma[1][0];
    
    parameters.n[Female] = std::atof(argv[14]);
    parameters.n[Male] = std::atof(argv[15]);

    parameters.eul_d = atof(argv[16]);
    parameters.eul_sr = atof(argv[17]);
    parameters.eul_b = atof(argv[18]);
    parameters.file_name = argv[19];

    parameters.p[0] = sigma[1][0] / (sigma[1][0] + sigma[0][1]);
    parameters.p[1] = 1.0 - p[0];

    parameters.delta_surv = false;

    TSD_Multivariate simulation{parameters};



} // end init_arguments()

    simulation.init_arguments(argc, argv);

    simulation.run();
}
