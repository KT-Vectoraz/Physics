#include "fluidsim.h"
#include <string.h>

void AdvectScalar(double *a, double dt)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                double x = (i + 0.5) * H, y = (j + 0.5) * H, z = (k + 0.5) * H;
                double bx = x - dt * SampU(x, y, z);
                double by = y - dt * SampV(x, y, z);
                double bz = z - dt * SampW(x, y, z);
                tc[kIC(i, j, k)] = SampC(a, bx, by, bz);
            }
    memcpy(a, tc, sizeof tc);
}

void AdvectAll(double dt)
{
    memcpy(tu, u, sizeof u);
    memcpy(tv, v, sizeof v);
    memcpy(tw, w, sizeof w);

    for (int i = 1; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                double x = i * H, y = (j + 0.5) * H, z = (k + 0.5) * H;
                double bx = x - dt * SampU(x, y, z);
                double by = y - dt * SampV(x, y, z);
                double bz = z - dt * SampW(x, y, z);
                tu[kIU(i, j, k)] = SampU(bx, by, bz);
            }

    for (int i = 0; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 0; k < N; k++) {
                double x = (i + 0.5) * H, y = j * H, z = (k + 0.5) * H;
                double bx = x - dt * SampU(x, y, z);
                double by = y - dt * SampV(x, y, z);
                double bz = z - dt * SampW(x, y, z);
                tv[kIV(i, j, k)] = SampV(bx, by, bz);
            }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 1; k < N; k++) {
                double x = (i + 0.5) * H, y = (j + 0.5) * H, z = k * H;
                double bx = x - dt * SampU(x, y, z);
                double by = y - dt * SampV(x, y, z);
                double bz = z - dt * SampW(x, y, z);
                tw[kIW(i, j, k)] = SampW(bx, by, bz);
            }

    AdvectScalar(T, dt);

    memcpy(u, tu, sizeof u);
    memcpy(v, tv, sizeof v);
    memcpy(w, tw, sizeof w);
}
