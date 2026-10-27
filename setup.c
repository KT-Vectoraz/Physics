#include "fluidsim.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ApplyBuoyancy(double dt)
{
    if (!buoyancyOn) return;
    for (int i = 0; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 0; k < N; k++) {
                double Tf = 0.5 * (T[kIC(i, j - 1, k)] + T[kIC(i, j, k)]);
                v[kIV(i, j, k)] += dt * kGrav * kBeta * (Tf - kTAmb);
            }
}

/* heats a square patch on the floor, centered in x and z */
void ApplyHeater(void)
{
    if (!heaterOn) return;
    for (int i = 0; i < N; i++) {
        double x = (i + 0.5) * H;
        if (x < 0.08 || x > 0.12) continue;
        for (int k = 0; k < N; k++) {
            double z = (k + 0.5) * H;
            if (z < 0.08 || z > 0.12) continue;
            T[kIC(i, 0, k)] = kTHeat;
        }
    }
}

void InitNoise(void)
{
    srand(1);
    memset(u, 0, sizeof u);
    memset(v, 0, sizeof v);
    memset(w, 0, sizeof w);
    for (int i = 1; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                u[kIU(i, j, k)] = 0.5 * (rand() / (double)RAND_MAX - 0.5);
    for (int i = 0; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 0; k < N; k++)
                v[kIV(i, j, k)] = 0.5 * (rand() / (double)RAND_MAX - 0.5);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 1; k < N; k++)
                w[kIW(i, j, k)] = 0.5 * (rand() / (double)RAND_MAX - 0.5);
    Project(1e-3);
}

int RunTests(void)
{
    int fail = 0;
    heaterOn = 0;
    buoyancyOn = 0;
    for (int c = 0; c < kCN; c++) T[c] = kTAmb;
    InitNoise();

    double ke0 = Kinetic(), kePrev = ke0, tTotal = 0.0;
    int bumps = 0;
    int steps = 200;
    for (int s = 0; s < steps; s++) {
        double dt = PickDt();
        Step(dt);
        tTotal += dt;
        double ke = Kinetic();
        if (ke > kePrev * (1.0 + 1e-6)) bumps++;
        kePrev = ke;
    }

    double um = Umax();
    double rel = Maxdiv() * H / (um > 1e-12 ? um : 1e-12);

    printf("Closed box (%d^3), no slip walls, no heating, %d steps %.4f s simulated\n",
           N, steps, tTotal);
    printf("  Kinetic energy: %.4e -> %.4e J  dissipated %.2f%%\n",
           ke0, kePrev, 100.0 * (1.0 - kePrev / ke0));
    printf("  steps where energy rose: %d  expect 0\n", bumps);
    printf("  relative divergence after projection: %.2e  expect < 1e-2\n", rel);

    if (kePrev >= ke0 || bumps > 0 || rel > 1e-2) fail = 1;
    printf("%s\n", fail ? "FAIL" : "PASS");
    return fail;
}
