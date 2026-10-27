#include "fluidsim.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* Build with gcc -std=c11 -O3 -Wall -Wextra -o fluid3d fluidsim.c util.c advection.c diffusion.c projection.c diagnostics.c setup.c -lm */
double u[kUN], v[kVN], w[kWN], T[kCN], p[kCN];
double tu[kUN], tv[kVN], tw[kWN], tc[kCN];
double rhs[kCN], cgR[kCN], cgZ[kCN], cgD[kCN], cgQ[kCN];
int heaterOn = 1, buoyancyOn = 1;

void Step(double dt)
{
    AdvectAll(dt);
    DiffuseVelocity(dt);
    DiffuseScalar(T, kAlpha, dt);
    ApplyBuoyancy(dt);
    ApplyHeater();
    Project(dt);
}

int main(int argc, char **argv)
{
    if (argc > 1 && strcmp(argv[1], "--test") == 0)
        return RunTests();

    long steps = (argc > 1) ? atol(argv[1]) : 1000;
    for (int c = 0; c < kCN; c++) T[c] = kTAmb;

    for (long s = 0; s < steps; s++) {
        double dt = PickDt();
        Step(dt);
    }
    return 0;
}
