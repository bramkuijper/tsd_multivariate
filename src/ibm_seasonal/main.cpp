#include "tsd_seasonal.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;

    params.file_name = argv[1];

    TSDSeasonal sim(params);
}
