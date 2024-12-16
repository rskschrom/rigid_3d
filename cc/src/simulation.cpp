#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include "io.h"
#include "simulation.h"

void Simulation::incrSimStep()
{
    simStep += 1;
}

int Simulation::getSimStep()
{
    return simStep;
}
