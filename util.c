#include "fluidsim.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double Clampd(double x, double lo, double hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}

/* samples an array at continuous coordinates using trilinear interpolation */
static double Trilerp(const double *a, int nx, int ny, int nz,
                      double fx, double fy, double fz)
{
    fx = Clampd(fx, 0.0, nx - 1);
    fy = Clampd(fy, 0.0, ny - 1);
    fz = Clampd(fz, 0.0, nz - 1);
    int i0 = (int)fx, j0 = (int)fy, k0 = (int)fz;
    if (i0 > nx - 2) i0 = nx - 2;
    if (j0 > ny - 2) j0 = ny - 2;
    if (k0 > nz - 2) k0 = nz - 2;
    double sx = fx - i0, sy = fy - j0, sz = fz - k0;

#define kAT(di, dj, dk) a[(((i0 + (di)) * ny) + (j0 + (dj))) * nz + (k0 + (dk))]
    double x00 = (1 - sz) * kAT(0, 0, 0) + sz * kAT(0, 0, 1);
    double x01 = (1 - sz) * kAT(0, 1, 0) + sz * kAT(0, 1, 1);
    double x10 = (1 - sz) * kAT(1, 0, 0) + sz * kAT(1, 0, 1);
    double x11 = (1 - sz) * kAT(1, 1, 0) + sz * kAT(1, 1, 1);
#undef kAT
    double y0 = (1 - sy) * x00 + sy * x01;
    double y1 = (1 - sy) * x10 + sy * x11;
    return (1 - sx) * y0 + sx * y1;
}

/* these turn world coordinates into staggered grid lookups */
double SampU(double x, double y, double z)
{
    return Trilerp(u, N + 1, N, N, x / H, y / H - 0.5, z / H - 0.5);
}
double SampV(double x, double y, double z)
{
    return Trilerp(v, N, N + 1, N, x / H - 0.5, y / H, z / H - 0.5);
}
double SampW(double x, double y, double z)
{
    return Trilerp(w, N, N, N + 1, x / H - 0.5, y / H - 0.5, z / H);
}
double SampC(const double *a, double x, double y, double z)
{
    return Trilerp(a, N, N, N, x / H - 0.5, y / H - 0.5, z / H - 0.5);
}

double Umax(void)
{
    double m = 0.0;
    for (int q = 0; q < kUN; q++) m = fmax(m, fabs(u[q]));
    for (int q = 0; q < kVN; q++) m = fmax(m, fabs(v[q]));
    for (int q = 0; q < kWN; q++) m = fmax(m, fabs(w[q]));
    return m;
}

double PickDt(void)
{
    double um = Umax();
    if (um < 1e-9) um = 1e-9;
    double dt = kCfl * H / um;
    double dv = 0.2 * H * H / fmax(kNu, kAlpha);   /* keep the time step small enough for diffusion to stay stable */
    if (dt > dv) dt = dv;
    if (dt > kDtMax) dt = kDtMax;
    return dt;
}
