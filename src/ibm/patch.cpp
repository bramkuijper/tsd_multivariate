#include "individual.hpp"
#include "patch.hpp"
// an individual patch of the metapopulation

Patch::Patch()
{
}//end Patch::Patch
    
Patch::Patch(Patch const &other)
{
    breedersF = other.breedersF;
    breedersM = other.breedersM;
    phil_juvsF = other.phil_juvsF;
    phil_juvsM = other.phil_juvsM;
}

void Patch::operator=(Patch const &other)
{
    breedersF = other.breedersF;
    breedersM = other.breedersM;
    phil_juvsF = other.phil_juvsF;
    phil_juvsM = other.phil_juvsM;
}
