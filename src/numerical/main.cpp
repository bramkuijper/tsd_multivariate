#include "tsd_turtle.hpp"

int main(int argc, char **argv)
{
    TSDTurtle simulation{};

    simulation.init_arguments(argc, argv);

    simulation.run();

}
