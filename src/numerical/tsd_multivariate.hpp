#ifndef TSD_MULTIVARIATE_HPP
#include <fstream>
#include "parameters.hpp"

class TSD_Multivariate
{
    private:
        Parameters pars;

        // the actual traits of interest
        double d[2][2]; // dispersal first index is sex, second index

        // eigenvalue
        double lambda;

        // reproductive values
        // first index sex
        // second index environment
        double v[2][2];

        // stable class frequencies
        // first index sex
        // second index environment
        double u[2][2];

        // relatedness coefficients between focal adult
        // and random offspring in local patch
        // first index is sex of focal adult
        // second index is environment of focal adult
        double rj[2][2];

        // coefficients of consanguinity between two randomly
        // sampled focal adults
        // index is environment
        double Qff[2];
        double Qfm[2];
        double Qmm[2];

        // private functions

        // calculate eigenvectors
        void eigenvectors(bool const output);

        // total number of competing offspring in the local site
        double C(Sex const sex_t1, bool const envt_t1);

        // derivative of the local competition function wrt to 
        // the mutant dispersal trait d[sex_d]
        double dCddlocal(Sex const sex_t1, 
                bool const envt_t1, 
                Sex const sex_d);
        
        // derivative of the local competition function wrt to 
        // the mutant sex ratio trait s[envt_d]
        double dCdsrlocal(Sex const sex_t1, 
                bool const envt_t1, 
                bool const envt_d);
        
        // derivative of the local competition function wrt to 
        // the burrowing trait b
        double dCdblocal(Sex const sex_t1, 
                bool const envt_t1);

        // total fecundity * survival of offsprign born in envt_t1 and
        // of sex sex_t1
        double fecundity_survival(bool const envt_t1, Sex const sex_t1);

        // derivative of fecundity functin wrt to trait s[envt_dsr]
        double dfecundity_survival_ds(
                bool const envt_t1
                ,Sex const sex_t1
                ,bool const envt_dsr
                );

        // derivative of fecundity functin wrt to trait b
        double dfecundity_survival_db(bool const envt_t1,Sex const sex_t1);

        void write_data(std::ofstream &output_file, int const generation);

        void write_parameters(std::ofstream &output_file);

        void write_data_headers(std::ofstream &output_file);

    public: 
        
        TSD_Multivariate();
        
        TSD_Multivariate(parstruct const &parstruct);

        // run the simulation
        void run();

        // initialize using command line arguments
        void init_arguments(int argc, char **argv);
        
        // element a((envt_t,sex_t) -> (envt_t1,sex_t1)) of the resident
        // transition matrix
        double A_resident(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1);

        double dB_dsr(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1
                ,bool const envt_dsr);

        double dB_db(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1);
        
        double dB_dd(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1
                ,Sex const sex_d);

        bool update_traits();

        void iterate_relatedness();
};

#endif
