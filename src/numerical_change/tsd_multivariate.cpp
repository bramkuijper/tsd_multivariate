#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


#include "tsd_multivariate.hpp"

TSD_Multivariate::TSD_Multivariate() :
    surv{{0.0,0.0},{0.0,0.0}}
    ,d{0.0,0.0}
    ,b{0.0}
    ,s{0.0,0.0}
    ,sigma{{0.0,0.0},{0.0,0.0}}
    ,p{0.0,0.0}
    ,lambda{0.0}
    ,v{{0.0,0.0},{0.0,0.0}}
    ,u{{0.0,0.0},{0.0,0.0}}
    ,n{0.0,0.0}
    ,delta_surv{false}
    ,eul_b{0.0}
    ,eul_sr{0.0}
    ,eul_d{0.0}
    ,base{}
{
} // end void TSD_Multivariate::TSD_Multivariate

// constructor with arguments
TSD_Multivariate::TSD_Multivariate(
        parstruct const &pstruct) :
            surv{{pstruct.surv[0][0],
                    pstruct.surv[0][1]},
                {pstruct.surv[1][0],
                    pstruct.surv[1][1]}}
            ,d{pstruct.d[0],d[1]}
            ,b{pstruct.b}
            ,s{pstruct.s[0],pstruct.s[1]}
            ,sigma{
                {pstruct.sigma[0][0]
                    ,pstruct.sigma[0][1]},
                {pstruct.sigma[1][0]
                    ,pstruct.sigma[1][1]}}
        ,p{pstruct.p[0],pstruct.p[1]}
        ,lambda{1.0}
        ,v{
            {pstruct.v[0][0],pstruct.v[0][1]}
            ,{pstruct.v[1][0],pstruct.v[1][1]}}
        ,u{
            {pstruct.u[0][0],pstruct.u[0][1]}
            ,{pstruct.u[1][0],pstruct.u[1][1]}}
        ,n{pstruct.n[0],pstruct.n[1]}
        ,delta_surv{pstruct.delta_surv}
        ,eul_d{pstruct.eul_d}
        ,eul_sr{pstruct.eul_sr}
        ,eul_b{pstruct.eul_b}
        ,base{pstruct.base}
{
}

// initialize arguments
void TSD_Multivariate::init_arguments(int argc, char **argv)
{
    surv_t0[Male][0] = std::atof(argv[1]);
    surv_t0[Male][1] = std::atof(argv[2]);
    surv_t0[Female][0] = std::atof(argv[3]);
    surv_t0[Female][1] = std::atof(argv[4]);
    surv_tend[Male][0] = std::atof(argv[5]);
    surv_tend[Male][1] = std::atof(argv[6]);
    surv_tend[Female][0] = std::atof(argv[7]);
    surv_tend[Female][1] = std::atof(argv[8]);

    d[Female] = std::atof(argv[9]);
    d[Male] = std::atof(argv[10]);

    b = std::atof(argv[11]);

    s[0] = std::atof(argv[12]);
    s[1] = std::atof(argv[13]);

    sigma_t0[0][1] = std::atof(argv[14]);
    sigma_t0[0][0] = 1.0 - sigma[0][1];

    sigma_t0[1][0] = std::atof(argv[15]);
    sigma_t0[1][1] = 1.0 - sigma[1][0];
    
    sigma_tend[0][1] = std::atof(argv[16]);
    sigma_tend[0][0] = 1.0 - sigma[0][1];

    sigma_tend[1][0] = std::atof(argv[17]);
    sigma_tend[1][1] = 1.0 - sigma[1][0];
   
    // set sigma to initial value
    sigma[0][1] = sigma_t0[0][1];
    sigma[0][0] = sigma_t0[0][0];
    sigma[1][0] = sigma_t0[1][0];
    sigma[1][1] = sigma_t0[1][1];

    lambda = std::atof(argv[18]);
    delta_surv = std::atoi(argv[19]);

    n[Female] = std::atof(argv[20]);
    n[Male] = std::atof(argv[21]);

    eul_d = atof(argv[22]);
    eul_sr = atof(argv[23]);
    eul_b = atof(argv[24]);
    base = argv[25];

    p[0] = sigma[1][0] / (sigma[1][0] + sigma[0][1]);
    p[1] = 1.0 - p[0];

} // end init_arguments()

// run the iteration
void TSD_Multivariate::run()
{
    // keep track of whether 
    // the envt has changed or not
    bool envt_changed = false;


    base = base + ".csv";

//    int skip_rows = 10;

    // initialize the file only when you run the thing
    std::ofstream output_file{base};

    write_parameters(output_file);

    write_data_headers(output_file);

    bool converged;

    for (int time_step = 0; time_step < max_time; ++time_step)
    {
        // update the eigenvectors during every time step
        eigenvectors(true);

        // update relatedness coefficients
        iterate_relatedness();

        converged = update_traits();

        if (converged)
        {
            if (envt_changed)
            {
                // if there is convergence and the envt has already changed
                // then end the simulation, we are done as we have converged
                // traits to their novel environment
                write_data(output_file, time_step);
                break;
            }
            else
            {
                write_data(output_file, time_step);
                change_envt();
            }
        }

        if (time_step % skip_rows == 0)
        {
            write_data(output_file, time_step);
        }
    }
} // end run()


// changes characteristics of the environment
void TSD_Multivariate::change_envt()
{
    for (unsigned idx1 = 0; idx1 < 2; ++idx1)
    {
        for (unsigned idx2 = 0; idx2 < 2; ++idx2)
        {
            sigma[idx1][idx2] = sigma_tend[idx1][idx2];
            surv[idx1][idx2] = surv_tend[idx1][idx2];
        }
    }
} // end change_envt()

void TSD_Multivariate::write_data_headers(std::ofstream &output_file)
{
    output_file << "generation;b;lambda;";

    for (int iter_i = 0; iter_i < 2; ++iter_i)
    {
        output_file << "d" << (static_cast<Sex>(iter_i) == Male ? "m" : "f") << ";";
        output_file << "sr" << iter_i + 1 << ";";
        output_file << "Qff" << iter_i + 1 << ";";
        output_file << "Qmm" << iter_i + 1 << ";";
        output_file << "Qfm" << iter_i + 1 << ";";

        for (int iter_j = 0; iter_j < 2; ++iter_j)
        {
            output_file << "v" << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << iter_j + 1 << ";"; 

            output_file << "u" << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << iter_j + 1 << ";";
            
            output_file << "r" << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << iter_j + 1 << ";";
        }
    }

    output_file << std::endl;
}

void TSD_Multivariate::write_data(std::ofstream &output_file, int const generation)
{
    output_file << generation << ";" << b << ";" << lambda << ";";

    for (int iter_i = 0; iter_i < 2; ++iter_i)
    {
        output_file << d[iter_i] << ";";
        output_file << s[iter_i] << ";";
        output_file << Qff[iter_i] << ";";
        output_file << Qmm[iter_i] << ";";
        output_file << Qfm[iter_i] << ";";

        for (int iter_j = 0; iter_j < 2; ++iter_j)
        {
            output_file << v[iter_i][iter_j] 
                << ";" << u[iter_i][iter_j] 
                << ";" << rj[iter_i][iter_j] << ";";
        }
    }

    output_file << std::endl;
}

void TSD_Multivariate::write_parameters(std::ofstream &output_file)
{
    output_file << std::endl
        << std::endl
        << "burrow_surv;" << delta_surv << std::endl
        << "eul_d;" << eul_d << std::endl
        << "eul_b;" << eul_b << std::endl
        << "eul_sr;" << eul_sr << std::endl;

    for (int iter_i = 0; iter_i < 2; ++iter_i)
    {
        // output p_i
        output_file << "p" << (iter_i + 1) << ";" << p[iter_i] << std::endl;

        // output nf, nm
        output_file << "n" << 
                (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << ";" << n[iter_i] << std::endl;
        output_file << "d_init" << 
                (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << ";" << d[iter_i] << std::endl;

        for (int iter_j = 0; iter_j < 2; ++iter_j)
        {
            output_file << "surv" 
                << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << (iter_j + 1) 
                << ";"
                << surv[iter_i][iter_j] 
                << std::endl;
            
            output_file << "surv_tend" 
                << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << (iter_j + 1) 
                << ";"
                << surv_tend[iter_i][iter_j] 
                << std::endl;

            output_file << "surv_t0" 
                << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << (iter_j + 1) 
                << ";"
                << surv_t0[iter_i][iter_j] 
                << std::endl;

            output_file << "sigma" 
                << (iter_i + 1) 
                << (iter_j + 1) 
                << ";"
                << sigma[iter_i][iter_j]
                << std::endl;

            << output_file << "sigma_tend"
                << (iter_i + 1) 
                << (iter_j + 1) 
                << ";"
                << sigma_tend[iter_i][iter_j]
                << std::endl;
            
            << output_file << "sigma_t0"
                << (iter_i + 1) 
                << (iter_j + 1) 
                << ";"
                << sigma_t0[iter_i][iter_j]
                << std::endl;
        } // end for iter_i
    } // end for iter_i
} // end write_parameters

/*
 * derivative of fecundity/survival function relative to b
 *
 * \param envt_t1  a boolean specifying environment in which offspring is born
 *
 * \param sex_t1 Sex enum specifying sex of offspring
 */

double TSD_Multivariate::dfecundity_survival_db(
        bool const envt_t1
        ,Sex const sex_t1
        ) 
{
    if (sex_t1 == Male) 
    {
        return(envt_t1 == 0 ? 
                0.0
                :
                -1.0 * s[envt_t1] * surv[sex_t1][envt_t1] + 
                    1.0 * s[!envt_t1] * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
    }

    // fecundity / survival in terms of 
    // expected numbers of daughters produced
    return(envt_t1 == 0 ?
            0.0
            :
            -1.0 * (1.0 - s[envt_t1]) * surv[sex_t1][envt_t1] + 
                    1.0 * (1.0 - s[!envt_t1]) * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
} // end d fecundity_survival dsr

/*
 * derivative of fecundity/survival function relative to sex ratio trait i
 *
 * \param envt_t1  a boolean specifying environment in which offspring is born
 *
 * \param sex_t1 Sex enum specifying sex of offspring
 *
 * \param envt_dsr a boolean specifying over which two of the sex ratio traits 
 * the derivative is taken
 */

double TSD_Multivariate::dfecundity_survival_ds(
        bool const envt_t1
        ,Sex const sex_t1
        ,bool const envt_dsr // 
        ) 
{
    // set up the derivatives
    double ds[2] = {0.0,0.0};
  
    // set the 'correct' derivative to 1
    ds[envt_dsr] = 1.0;

    if (sex_t1 == Male) 
    {
        return(envt_t1 == 0 ? 
                ds[envt_t1] * surv[sex_t1][envt_t1]
                :
                (1.0 - b) * ds[envt_t1] * surv[sex_t1][envt_t1] + 
                    b * ds[!envt_t1] * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
    }

    // fecundity / survival in terms of 
    // expected numbers of daughters produced
    return(envt_t1 == 0 ?
            - ds[envt_t1] * surv[sex_t1][envt_t1]
            :
            (1.0 - b) * (- ds[envt_t1]) * surv[sex_t1][envt_t1] + 
                    b * (- ds[!envt_t1]) * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
} // end d fecundity_survival dsr

// total fecundity * survival of offsprign born in envt_t1 and
// of sex sex_t1
double TSD_Multivariate::fecundity_survival(bool const envt_t1, Sex const sex_t1) 
{
    if (sex_t1 == Male) 
    {
        return(envt_t1 == 0 ? 
                s[envt_t1] * surv[sex_t1][envt_t1]
                :
                (1.0 - b) * s[envt_t1] * surv[sex_t1][envt_t1] + 
                    b * s[!envt_t1] * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
    }

    return(envt_t1 == 0 ?
            (1.0 - s[envt_t1]) * surv[sex_t1][envt_t1]
            :
            (1.0 - b) * (1.0 - s[envt_t1]) * surv[sex_t1][envt_t1] + 
                    b * (1.0 - s[!envt_t1]) * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
            
} // end fecundity_survival


// derivative of the local competition function wrt to 
// the mutant dispersal trait d[sex_d]
double TSD_Multivariate::dCdsrlocal(Sex const sex_t1
        ,bool const envt_t1
        ,bool const envt_d
        )
{
    return((1.0 - d[sex_t1]) * dfecundity_survival_ds(envt_t1, sex_t1, envt_d));
} // end TSD_Multivariate::dCdsrlocal()

// derivative of the local competition function wrt to 
// the mutant burrowing trait 
double TSD_Multivariate::dCdblocal(Sex const sex_t1
        ,bool const envt_t1
        )
{
    return((1.0 - d[sex_t1]) * dfecundity_survival_db(envt_t1, sex_t1));
} // end TSD_Multivariate::dCddlocal()
// derivative of the local competition function wrt to 
// the mutant dispersal trait d[sex_d]
double TSD_Multivariate::dCddlocal(Sex const sex_t1
        ,bool const envt_t1
        ,Sex const sex_d
        )
{
    if (sex_t1 != sex_d)
    {
        return(0.0);
    }

    return(-1.0 * fecundity_survival(envt_t1, sex_t1)); 
} // end TSD_Multivariate::dCddlocal()

// total number of competing offspring
double TSD_Multivariate::C(Sex const sex_t1, bool const envt_t1)
{
    return((1.0 - d[sex_t1]) * fecundity_survival(envt_t1, sex_t1) 
            + d[sex_t1] * 
                    (p[0] * fecundity_survival(0, sex_t1)
                     +  p[1] * fecundity_survival(1, sex_t1)));
} // endl TSD_Multivariate::C()

/**
 * iterate relatedness coefficients
 * and coefficients of consanguinity
 **/
void TSD_Multivariate::iterate_relatedness()
{
    double Qmmtplus1[2] = {0.0,0.0};
    double Qfftplus1[2] = {0.0,0.0};
    double Qfmtplus1[2] = {0.0,0.0};

    // aux variable to store inner sum of coeffs
    // of consanguinity
    double Qxksum = 0.0;

    // probability that individual establishes
    // breeding position in natal patch of envt_k
    double g_male[2] = {
        (1.0 - d[Male]) * fecundity_survival(0, Male) / C(Male, 0),
        (1.0 - d[Male]) * fecundity_survival(1, Male) / C(Male, 1)
    };

    double g_female[2] = {
        (1.0 - d[Female]) * fecundity_survival(0, Female) / C(Female, 0),
        (1.0 - d[Female]) * fecundity_survival(1, Female) / C(Female, 1)
    };
            
    // total probability environment is in state j
    double sum_lj[2] = {
        p[0] * sigma[0][0] + p[1] * sigma[1][0],
        p[0] * sigma[0][1] + p[1] * sigma[1][1]
    };

    bool converged = false;

    // variable to store value of iterator
    int iter_t;

    for (iter_t = 0; iter_t <= max_time; ++iter_t)
    {
        converged = true;

        for (int envt_j = 0; envt_j < 2; ++envt_j)
        {
            Qmmtplus1[envt_j] = 0.0;
            Qfftplus1[envt_j] = 0.0;
            Qfmtplus1[envt_j] = 0.0;

            for (int envt_k = 0; envt_k < 2; ++envt_k)
            {
                Qxksum = 0.25 *(1.0/n[Female] + (n[Female] - 1.0) / n[Female] * Qff[envt_k])
                         + 0.5 * Qfm[envt_k]
                         + 0.25 * (1.0 / n[Male] + (n[Male] - 1.0) / n[Male]  * Qmm[envt_k]);

                Qmmtplus1[envt_j] += p[envt_k] * sigma[envt_k][envt_j] * 
                    pow(g_male[envt_k],2.0) / sum_lj[envt_j] * Qxksum;

                Qfftplus1[envt_j] += p[envt_k] * sigma[envt_k][envt_j] * 
                    pow(g_female[envt_k],2.0) / sum_lj[envt_j] * Qxksum;

                Qfmtplus1[envt_j] += p[envt_k] * sigma[envt_k][envt_j] * 
                    g_female[envt_k] * g_male[envt_k] / sum_lj[envt_j] * Qxksum;

                assert(std::isnan(Qmmtplus1[envt_j]) == 0);
                assert(std::isinf(Qmmtplus1[envt_j]) == 0);

                assert(std::isnan(Qfftplus1[envt_j]) == 0);
                assert(std::isinf(Qfftplus1[envt_j]) == 0);

                assert(std::isnan(Qfmtplus1[envt_j]) == 0);
                assert(std::isinf(Qfmtplus1[envt_j]) == 0);
            } // end for int envt_k
        } // end for (int envt_j = 0; envt_j < 2; ++envt_j)

        for (int envt_j = 0; envt_j < 2; ++envt_j)
        {
            if (fabs(Qmmtplus1[envt_j] - Qmm[envt_j]) > vanish_bound)
            {
                converged = false;
            }
            
            if (fabs(Qfmtplus1[envt_j] - Qfm[envt_j]) > vanish_bound)
            {
                converged = false;
            }
            
            if (fabs(Qfftplus1[envt_j] - Qff[envt_j]) > vanish_bound)
            {
                converged = false;
            }

            Qmm[envt_j] = Qmmtplus1[envt_j];
            Qff[envt_j] = Qfftplus1[envt_j];
            Qfm[envt_j] = Qfmtplus1[envt_j];
        }

        if (converged)
        {
            break;
        }
    } // end for (int iter_t = 0; iter_t <= max_time; ++iter_t)

    if (iter_t == max_time - 1)
    {
        std::cout << "relatedness coefficients did not converge" << std::endl;
    }

    // Wrights equilibrium inbreeding coefficient
    double Qj;

    for (int envt_j = 0; envt_j < 2; ++envt_j)
    {
        Qj = 0.5 * (1 + Qfm[envt_j]);

        // relatedness between a focal adult female and any local offspring of 
        // a particular sex
        rj[Male][envt_j] = 0.5 * 1.0 / n[Male] * (Qj + Qfm[envt_j]) +
                0.5 * (n[Male] - 1.0) / n[Male] * (Qmm[envt_j] + Qfm[envt_j]);
        
        rj[Female][envt_j] = 0.5 * 1.0 / n[Female] * (Qj + Qfm[envt_j]) +
            0.5 * (n[Female] - 1.0)/n[Female] * (Qff[envt_j] + Qfm[envt_j]);
    }
} // end iterate_relatedness()

// partial derivative of mutant transition matrix
// wrt burrowing depth
double TSD_Multivariate::dB_db(
        bool const envt_t
        ,Sex const sex_t
        ,bool const envt_t1
        ,Sex const sex_t1)
{
    double dbij_db_focal = 0.0;
    double dbij_db_local = 0.0;

    double rx_own = sex_t == Male ? Qfm[envt_t] : 0.5 * (1.0 + Qfm[envt_t]);
    double rx_local = sex_t == Male ? 
        Qfm[envt_t] 
        : 
        1.0 / n[Female] * 0.5 * (1.0 + Qfm[envt_t]) + 
            (n[Female] - 1.0) / n[Female] * Qff[envt_t];

    if (sex_t == Female)
    {
        dbij_db_focal = sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) *
            dfecundity_survival_db(envt_t, sex_t1) / C(sex_t1, envt_t)
            + (
                p[envt_t1] * (1.0 - sigma[envt_t1][!envt_t1]) / C(sex_t1, envt_t1) 
                + 
                p[!envt_t1] * sigma[!envt_t1][envt_t1] / C(sex_t1, !envt_t1)
            ) * d[sex_t1] *
                        dfecundity_survival_db(envt_t, sex_t1); 
    }
    else
    {
        // if focal is a male, only local derivatives matter
        dbij_db_local = sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) *
            dfecundity_survival_db(envt_t, sex_t1) / C(sex_t1, envt_t)
            + (
                p[envt_t1] * (1.0 - sigma[envt_t1][!envt_t1]) / C(sex_t1, envt_t1) 
                + 
                p[!envt_t1] * sigma[!envt_t1][envt_t1] / C(sex_t1, !envt_t1)
            ) * d[sex_t1] *
                        dfecundity_survival_db(envt_t, sex_t1); 
    }

    // derivative of the local traits in the denominator
    dbij_db_local += sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) * 
                                -fecundity_survival(envt_t, sex_t1) * 
                                    dCdblocal(sex_t1, envt_t) / 
                                        (C(sex_t1, envt_t) * C(sex_t1, envt_t));
    // remote fitness component 0 for bij_local
    return(rx_own * dbij_db_focal + rx_local * dbij_db_local);
} // end TSD_Multivariate::dB_db


// partial derivative of mutant transition matrix
// wrt to the mutant dispersal trait of sex_d
double TSD_Multivariate::dB_dd(
        bool const envt_t
        ,Sex const sex_t
        ,bool const envt_t1
        ,Sex const sex_t1
        ,Sex const sex_d)
{
    double dbij_dd_focal = 0.0;
    double dbij_dd_local = 0.0;

    double dsexd[2] = {0.0,0.0};

    dsexd[sex_d] = 1.0;

    dbij_dd_focal = sigma[envt_t][envt_t1] * (-dsexd[sex_t1]) *
        fecundity_survival(envt_t, sex_t1) / C(sex_t1, envt_t)
        + (
            p[envt_t1] * (1.0 - sigma[envt_t1][!envt_t1]) / C(sex_t1, envt_t1) 
            + 
            p[!envt_t1] * sigma[!envt_t1][envt_t1] / C(sex_t1, !envt_t1)
        ) * dsexd[sex_t1] * fecundity_survival(envt_t, sex_t1); 
        
    dbij_dd_local = -sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) *
            fecundity_survival(envt_t, sex_t1) * dCddlocal(sex_t1, envt_t, sex_d) / 
                (C(sex_t1, envt_t) * C(sex_t1, envt_t));

    return((0.5 + 0.5 * Qfm[envt_t]) * dbij_dd_focal 
            + rj[sex_t][envt_t] * dbij_dd_local);
} // end TSD_Multivariate::dB_db

// partial derivatives of the mutant transition matrix
// wrt to mutant sex ratio trait used in environment envt_dsr in {0,1}
double TSD_Multivariate::dB_dsr(
        bool const envt_t
        ,Sex const sex_t
        ,bool const envt_t1
        ,Sex const sex_t1
        ,bool const envt_dsr)
{
    double dbij_dsr_focal = 0.0;
    double dbij_dsr_local = 0.0;

    if (sex_t == Female)
    {
         dbij_dsr_focal = sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) *
            dfecundity_survival_ds(envt_t, sex_t1, envt_dsr) / C(sex_t1, envt_t)
            + (
                p[envt_t1] * (1.0 - sigma[envt_t1][!envt_t1]) / C(sex_t1, envt_t1) 
                + 
                p[!envt_t1] * sigma[!envt_t1][envt_t1] / C(sex_t1, !envt_t1)
            ) * d[sex_t1] *
                        dfecundity_survival_ds(envt_t, sex_t1, envt_dsr); 
    }
    else // sex_t == Male
    {
        // if focal is a male, only local derivatives matter
        dbij_dsr_local = sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) *
            dfecundity_survival_ds(envt_t, sex_t1, envt_dsr) / C(sex_t1, envt_t)
            + (
                p[envt_t1] * (1.0 - sigma[envt_t1][!envt_t1]) / C(sex_t1, envt_t1) 
                + 
                p[!envt_t1] * sigma[!envt_t1][envt_t1] / C(sex_t1, !envt_t1)
            ) * d[sex_t1] *
                        dfecundity_survival_ds(envt_t, sex_t1, envt_dsr); 
    }

    // derivative of the local traits in the denominator
    dbij_dsr_local += sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) * 
                                -fecundity_survival(envt_t, sex_t1) * 
                                    dCdsrlocal(sex_t1, envt_t, envt_dsr) / 
                                        (C(sex_t1, envt_t) * C(sex_t1, envt_t));
               // remote fitness component 0 for bij_local

    return((0.5 + 0.5 * Qfm[envt_t]) * dbij_dsr_focal 
            + rj[sex_t][envt_t] * dbij_dsr_local);
} // end dB_dsr

// update the trait values
// due to successive mutations of small effect
bool TSD_Multivariate::update_traits() 
{
    // total selection gradient on sex ratio loci
    double delta_sr[2] = {0.0,0.0};
    // total selection gradient on philopatry loci
    double delta_d[2] = {0.0,0.0};
    // total selection gradient on burrowing depth locus
    double delta_b = 0.0;

    bool converged = true;

    Sex Sex_t1;
    Sex Sex_t;

    // go through all entries of mutant transition matrix
    for (int envt_t = 0; envt_t < 2; ++envt_t)
    {
        for (int sex_t = 0; sex_t < 2; ++sex_t)
        {
            Sex_t = static_cast<Sex>(sex_t);
            for (int envt_t1 = 0; envt_t1 < 2; ++envt_t1)
            {
                for (int sex_t1 = 0; sex_t1 < 2; ++sex_t1)
                {
                    Sex_t1 = static_cast<Sex>(sex_t1);

                    for (int envt_dsr = 0; envt_dsr < 2; ++envt_dsr)
                    {
                        delta_sr[envt_dsr] += v[Sex_t1][envt_t1] * u[Sex_t][envt_t] *
                            dB_dsr(envt_t, Sex_t, envt_t1, Sex_t1, envt_dsr);

//                        std::cout  
//                           << " v_" << (Sex_t1 == Male ? "m" : "f")
//                           << "_" << (envt_t1 + 1) 
//                           << " " << v[Sex_t1][envt_t1] << "   "
//                           << " u_" << (Sex_t == Male ? "m" : "f")
//                           << "_" << (envt_t + 1) 
//                            << " " <<  u[Sex_t][envt_t]
//                            << " deriv over "  << (envt_dsr + 1) << " "
//                            << dB_dsr(envt_t, Sex_t, envt_t1, Sex_t1, envt_dsr) << std::endl;


                        assert(std::isnan(delta_sr[envt_dsr]) == 0); 
                        assert(std::isinf(delta_sr[envt_dsr]) == 0); 
                    }


                    for (int sex_d = 0; sex_d < 2; ++sex_d)
                    {
                        delta_d[sex_d] += v[Sex_t1][envt_t1] * u[Sex_t][envt_t] *
                            dB_dd(envt_t, Sex_t, envt_t1, Sex_t1, static_cast<Sex>(sex_d));
                        
                        assert(std::isnan(delta_d[sex_d]) == 0); 
                        assert(std::isinf(delta_d[sex_d]) == 0); 
                    }

//                    std::cout  
//                       << " u_" << (Sex_t == Male ? "m" : "f")
//                       << "_" << (envt_t + 1) 
//                        << ": " <<  u[Sex_t][envt_t]
//                       << " v_" << (Sex_t1 == Male ? "m" : "f")
//                       << "_" << (envt_t1 + 1) 
//                       << ": " << v[Sex_t1][envt_t1] << "   "
//                        << " dB_db: "
//                        << dB_db(envt_t, Sex_t, envt_t1, Sex_t1) << std::endl;

                    delta_b += v[Sex_t1][envt_t1] * u[Sex_t][envt_t] *
                            dB_db(envt_t, Sex_t, envt_t1, Sex_t1); 
                        
                    assert(std::isinf(delta_b) == 0); 
                    assert(std::isnan(delta_b) == 0); 
                }
            }
        }
    }

    double ztplus1;

    for (int envt_dsr = 0; envt_dsr < 2; ++envt_dsr)
    {
        ztplus1 = std::clamp(s[envt_dsr] + eul_sr * delta_sr[envt_dsr]/lambda,0.0,1.0);


        if (fabs(ztplus1 - s[envt_dsr]) > vanish_bound)
        {
            converged = false;
        }

        s[envt_dsr] = ztplus1;
    }
    
    for (int sex_d = 0; sex_d < 2; ++sex_d)
    {
        ztplus1 = std::clamp(d[sex_d] + eul_d * delta_d[sex_d]/lambda,0.0,1.0);

        if (fabs(ztplus1 - d[sex_d]) > vanish_bound)
        {
            converged = false;
        }

        d[sex_d] = ztplus1;
    }

    ztplus1 = std::clamp(b + eul_b * delta_b/lambda,0.0,1.0);

    if (fabs(ztplus1 - b) > vanish_bound)
    {
        converged = false;
    }

    b = ztplus1;

    return(converged);
} // end TSD_Multivariate::update_traits()

// get entries of the resident transition matrix
double TSD_Multivariate::A_resident(
        bool const envt_t
        ,Sex const sex_t
        ,bool const envt_t1
        ,Sex const sex_t1)
{
    // philopatric fitness component
    double aij = sigma[envt_t][envt_t1] * (1.0 - d[sex_t1]) * 
                                fecundity_survival(envt_t, sex_t1) / C(sex_t1, envt_t)
               // remote fitness component 
                    + (
                        p[envt_t1] * (1.0 - sigma[envt_t1][!envt_t1]) / C(sex_t1, envt_t1)
                        + 
                        p[!envt_t1] * sigma[!envt_t1][envt_t1] / C(sex_t1, !envt_t1)
                    ) * d[sex_t1] *
                            fecundity_survival(envt_t, sex_t1); 

    return(aij);                                        
} // end A_resident

// calculate left and right eigenvectors
void TSD_Multivariate::eigenvectors(bool const output)
{
    size_t ndim = 4;

    // allocate GSL things to work on the eigenvectors
    //
    // allocate a vector to calculate eigenvalues
    gsl_vector_complex *eval = gsl_vector_complex_alloc(ndim);

    // allocate a matrix to calculate eigenvectors
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(ndim, ndim);
     
    // allocate a vector to calculate eigenvalues for the tranpose
    gsl_vector_complex *evalT = gsl_vector_complex_alloc(ndim);
    
    // allocate a matrix to calculate eigenvectors for the tranpose
    gsl_matrix_complex *evecT = gsl_matrix_complex_alloc(ndim, ndim);

    gsl_matrix *m = gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *mT = gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *m_stats = gsl_matrix_alloc(ndim,ndim);
     
    gsl_eigen_nonsymmv_workspace * workspace_w = gsl_eigen_nonsymmv_alloc(ndim);
    gsl_eigen_nonsymmv_workspace * workspace_wT = gsl_eigen_nonsymmv_alloc(ndim);

    // auxiliary variable to see whether stuff has been converged
    bool converged;

    // auxiliary variables to store the lhs
    double v_i_tplus1, u_i_tplus1, ev_tplus1;

    double Atmp;

    size_t row_i, col_j;
        
    // now calculate solutions of the eigenvectors
    // until convergence
    for (long int iter_time_step = 0; 
            iter_time_step < max_time;
            ++iter_time_step)
    {
        row_i = 0;
        col_j = 0;

        for (size_t envt_i = 0; envt_i < 2; ++envt_i)
        {
            for (size_t sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (size_t envt_j = 0; envt_j < 2; ++envt_j)
                {
                    for (size_t sex_j = 0; sex_j < 2; ++sex_j)
                    {
                        Atmp = A_resident(
                                envt_i
                                ,static_cast<Sex>(sex_i)
                                ,envt_j
                                ,static_cast<Sex>(sex_j));

//                        std::cout << "Atmp_" <<  
//                                (static_cast<Sex>(sex_i) == Male ? "m" : "f")
//                                << "_" << (envt_i + 1) << "_" <<
//                                (static_cast<Sex>(sex_j) == Male ? "m" : "f")
//                                << "_" << (envt_j + 1) << " " << Atmp << " " << std::endl;
                        
                        gsl_matrix_set(m, row_i, col_j, Atmp);
                        gsl_matrix_set(mT, col_j, row_i, Atmp);
                        gsl_matrix_set(m_stats, row_i, col_j, Atmp);

                        // increment column after 
                        // each (sex_i, envt_i) combination
                        ++col_j;

                        // if maximum column number has been reached
                        // reset to 0
                        if (col_j > 3)
                        {
                            col_j = 0;
                        }
                        
                    }// end for sex_j
                } // end for envt_j

                // increment row after each (sex_i, envt_i) combination
                ++row_i;

                // if max row number has been reached
                // reset to 0
                if (row_i > 3)
                {
                    row_i = 0;
                }
            } // end for envt_i
        } //end for sex_i
    
        // now calculate evs
        // set up the workspaces
        gsl_eigen_nonsymmv(m, eval, evec, workspace_w);
        gsl_eigen_nonsymmv(mT, evalT, evecT, workspace_wT);
                
        

        gsl_eigen_nonsymmv_sort(eval, evec, 
                GSL_EIGEN_SORT_ABS_DESC);
            
        // get the first column from the eigenvector matrix
        // giving us the right eigenvector
        gsl_vector_complex_view evec_i 
                   = gsl_matrix_complex_column(evec, 0);

        // get the first column from the eigenvector matrix
        // of the transpose of the left eigenvector
        // giving us the left eigenvector
        gsl_eigen_nonsymmv_sort (evalT, evecT, 
                                GSL_EIGEN_SORT_ABS_DESC);

        // get the last column from the eigenvector
        gsl_vector_complex_view evecT_i 
                   = gsl_matrix_complex_column(evecT, 0);

        double product_vu = 0;
        double total_u = 0;

        // make variables that normalize both eigenvectors
        for (size_t iter1 = 0; iter1 < ndim; ++iter1)
        {
            // norm for right ev
            total_u += fabs(GSL_REAL(
                   gsl_vector_complex_get(&evec_i.vector, iter1)));

            // norm for left ev
            product_vu += fabs(GSL_REAL(
                   gsl_vector_complex_get(&evec_i.vector, iter1)))
                   * fabs(GSL_REAL(
                   gsl_vector_complex_get(&evecT_i.vector, iter1)));
        }

        product_vu /= total_u;
       
        // by default we assume convergence until proved otherwise
        converged = true;

        bool envt_idx;
        Sex sex_idx;

        // get the results from the eigenvector
        for (size_t dim_i = 0; dim_i < ndim; ++dim_i)
        {
            // get the right eigenvalue out of the eigenvector of A
            u_i_tplus1 = fabs(GSL_REAL(
                        gsl_vector_complex_get(&evec_i.vector,dim_i)))/
                total_u;
            
            // get the left eigenvalue out of the eigenvector of A^T
            v_i_tplus1 = fabs(GSL_REAL(
                        gsl_vector_complex_get(&evecT_i.vector,dim_i)))/
                product_vu;

            // translate dimensions into
            // envt sex, with
            // 0: envt = 0, sex = 0 
            sex_idx = static_cast<Sex>(dim_i % 2);
            envt_idx = static_cast<Sex>(std::floor(dim_i/2.0));

            if (fabs(u[sex_idx][envt_idx] - u_i_tplus1) > vanish_bound)
            {
                converged = false;
            }

            if (fabs(v[sex_idx][envt_idx] - v_i_tplus1) > vanish_bound)
            {
                converged = false;
            }

            v[sex_idx][envt_idx] = v_i_tplus1;
            u[sex_idx][envt_idx] = u_i_tplus1;

            assert(std::isnan(v[sex_idx][envt_idx]) == 0);
            assert(std::isnan(u[sex_idx][envt_idx]) == 0);
        } // for dim_i
        
        // now get the eigenvalue
        ev_tplus1 = GSL_REAL(gsl_vector_complex_get(eval, 0));

        if (fabs(lambda - ev_tplus1) > vanish_bound)
        {
            converged = false;
        }

        lambda = ev_tplus1;

        if (converged)
        {
            break;
        }

    } // end iteration

//    for (size_t envt_idx = 0; envt_idx < 2; ++envt_idx)
//    {
//        for (size_t sex_idx = 0; sex_idx < 2; ++sex_idx)
//        {
//            std::cout << "v_" <<  
//                   (static_cast<Sex>(sex_idx) == Male ? "m" : "f")
//                   << "_" << (envt_idx + 1) << " " << v[sex_idx][envt_idx] << std::endl;
//            
//            std::cout << "u_" <<  
//                   (static_cast<Sex>(sex_idx) == Male ? "m" : "f")
//                   << "_" << (envt_idx + 1) << " " << u[sex_idx][envt_idx] << std::endl;
//        }
//    }

    // free all the stuff we've used
    gsl_eigen_nonsymmv_free (workspace_w);
    gsl_eigen_nonsymmv_free (workspace_wT);
    gsl_vector_complex_free(evalT);
    gsl_matrix_complex_free(evecT);
    gsl_matrix_free(mT);
    gsl_matrix_free(m);
    gsl_matrix_free(m_stats);

}

