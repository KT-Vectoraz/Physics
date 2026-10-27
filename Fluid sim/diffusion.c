#include "fluidsim.h"
#include <string.h>

void DiffuseVelocity(double dt)
{
    double a = kNu * dt / (H * H);
    memcpy(tu, u, sizeof u);
    memcpy(tv, v, sizeof v);
    memcpy(tw, w, sizeof w);

    /* x is the normal direction for u, while y and z get no slip ghost cells */
    for (int i = 1; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                double c = u[kIU(i, j, k)];
                double l = u[kIU(i - 1, j, k)], r = u[kIU(i + 1, j, k)];
                double d = (j == 0)     ? -c : u[kIU(i, j - 1, k)];
                double t = (j == N - 1) ? -c : u[kIU(i, j + 1, k)];
                double f = (k == 0)     ? -c : u[kIU(i, j, k - 1)];
                double b = (k == N - 1) ? -c : u[kIU(i, j, k + 1)];
                tu[kIU(i, j, k)] = c + a * (l + r + d + t + f + b - 6.0 * c);
            }

    /* y is the normal direction for v, while x and z get no slip ghost cells */
    for (int i = 0; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 0; k < N; k++) {
                double c = v[kIV(i, j, k)];
                double d = v[kIV(i, j - 1, k)], t = v[kIV(i, j + 1, k)];
                double l = (i == 0)     ? -c : v[kIV(i - 1, j, k)];
                double r = (i == N - 1) ? -c : v[kIV(i + 1, j, k)];
                double f = (k == 0)     ? -c : v[kIV(i, j, k - 1)];
                double b = (k == N - 1) ? -c : v[kIV(i, j, k + 1)];
                tv[kIV(i, j, k)] = c + a * (l + r + d + t + f + b - 6.0 * c);
            }

    /* z is the normal direction for w, while x and y get no slip ghost cells */
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 1; k < N; k++) {
                double c = w[kIW(i, j, k)];
                double f = w[kIW(i, j, k - 1)], b = w[kIW(i, j, k + 1)];
                double l = (i == 0)     ? -c : w[kIW(i - 1, j, k)];
                double r = (i == N - 1) ? -c : w[kIW(i + 1, j, k)];
                double d = (j == 0)     ? -c : w[kIW(i, j - 1, k)];
                double t = (j == N - 1) ? -c : w[kIW(i, j + 1, k)];
                tw[kIW(i, j, k)] = c + a * (l + r + d + t + f + b - 6.0 * c);
            }

    memcpy(u, tu, sizeof u);
    memcpy(v, tv, sizeof v);
    memcpy(w, tw, sizeof w);
}

void DiffuseScalar(double *a, double D, double dt)
{
    double k0 = D * dt / (H * H);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                double c = a[kIC(i, j, k)];
                double l = (i > 0)     ? a[kIC(i - 1, j, k)] : c;
                double r = (i < N - 1) ? a[kIC(i + 1, j, k)] : c;
                double d = (j > 0)     ? a[kIC(i, j - 1, k)] : c;
                double t = (j < N - 1) ? a[kIC(i, j + 1, k)] : c;
                double f = (k > 0)     ? a[kIC(i, j, k - 1)] : c;
                double b = (k < N - 1) ? a[kIC(i, j, k + 1)] : c;
                tc[kIC(i, j, k)] = c + k0 * (l + r + d + t + f + b - 6.0 * c);
            }
    memcpy(a, tc, sizeof tc);
}
