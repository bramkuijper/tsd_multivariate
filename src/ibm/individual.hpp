#ifndef INDIVIDUAL_HPP_
#define INDIVIDUAL_HPP_

enum Sex {
    Female = 0,
    Male = 1
};

class Individual
{
    public:
        // sex ratio in envt 1 and 2, each coded by diploid loci
        double sr[2][2];

        // male and female dispersal, each coded by diploid loci
        double d[2][2];

        double b[2];

        double viability;

        bool sex;

        // constructor (i.e., build this individual)
        Individual();

        // copy constructor (i.e., copy this individual from another one)
        Individual(Individual const &other);

        // change the assignment operator so that also individuals
        // can be assigned, i..e., Individual1 = Individual2;
        void operator=(Individual const &other);

}; // end class individual





#endif
