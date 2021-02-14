#include "individual.hpp"

// default constructor
Individual::Individual() :
    sr{{0.0,0.0},{0.0,0.0}}, // member initializers, for more info see http://www.icce.rug.nl/documents/cplusplus/cplusplus07.html#l143 
    d{{0.0,0.0},{0.0,0.0}},
    b{0.0,0.0},
    viability{0.0},
    sex{0}
{
}

// copy constructor
Individual::Individual(Individual const &other) :
    sr{{other.sr[0][0],other.sr[0][1]},{other.sr[1][0],other.sr[1][1]}}, 
    d{{other.d[0][0],other.d[0][1]},{other.d[1][0],other.d[1][1]}},
    b{other.b[0],other.b[1]},
    viability{other.viability},
    sex{other.sex}
{
}

void Individual::operator=(Individual const &other)
{
    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        sr[0][allele_idx] = other.sr[0][allele_idx];
        sr[1][allele_idx] = other.sr[1][allele_idx];
        d[0][allele_idx] = other.d[0][allele_idx];
        d[1][allele_idx] = other.d[1][allele_idx];

        b[allele_idx] = other.b[allele_idx];
    }


    viability = other.viability;
    sex = other.sex;
} // end Individual::operator=()
// 
