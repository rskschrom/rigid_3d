#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <math.h>
#include <vector>
#include "particle.h"

int main(void){

    float pmass = 920.;
    Particle par = Particle("crystal_0000_1.5.txt", pmass);
    par.write();
}
