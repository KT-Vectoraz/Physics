#include "fluidsim.h"
#include <math.h>

double Kinetic(void)
{
    double s = 0.0;
    for (int q = 0; q < kUN; q++) s += u[q] * u[q];
    for (int q = 0; q < kVN; q++) s += v[q] * v[q];
    for (int q = 0; q < kWN; q++) s += w[q] * w[q];
    return 0.5 * kRho * s * H * H * H;
}

double Maxdiv(void)
{
    double m = 0.0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                double d = (u[kIU(i + 1, j, k)] - u[kIU(i, j, k)]
                          + v[kIV(i, j + 1, k)] - v[kIV(i, j, k)]
                          + w[kIW(i, j, k + 1)] - w[kIW(i, j, k)]) / H;
                m = fmax(m, fabs(d));
            }
    return m;
}
