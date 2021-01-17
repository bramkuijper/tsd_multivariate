#include "individual.hpp"

// default constructor
Individual::Individual() :
    t{0.0,0.0}, // member initializers, for more info see http://www.icce.rug.nl/documents/cplusplus/cplusplus07.html#l143 
    tprime{0.0,0.0},
    p{0.0,0.0},
    p_phen{0.0},
    s{0.0},
    envt_quality_high{false}
{
}

// copy constructor
Individual::Individual(Individual const &other) :
    t{other.t[0],other.t[1]},
    tprime{other.tprime[0],other.tprime[1]},
    p{other.p[0],other.p[1]},
    p_phen{other.p_phen},
    s{other.s},
    envt_quality_high{other.envt_quality_high}
{
}

void Individual::operator=(Individual const &other)
{
    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        t[allele_idx] = other.t[allele_idx];
        tprime[allele_idx] = other.tprime[allele_idx];
        p[allele_idx] = other.p[allele_idx];
    }

    p_phen = other.p_phen;
    s = other.s;
    envt_quality_high = other.envt_quality_high;
} // end Individual::operator=()
// 
