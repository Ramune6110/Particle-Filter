#include <iostream>
#include <Eigen/Dense>
#include <cstdio>
#include "particle_filter.h"

using namespace std;
using namespace Eigen;

int main() 
{
    Particle_Filter PF;

    PF.simulation();

    return 0;
}
