#include "fluidsim.h"
#include <math.h>

static double Dot(const double *a, const double *b, int n)
{
    double s = 0.0;
    for (int q = 0; q < n; q++) s += a[q] * b[q];
    return s;
}

int Neighbours(int i, int j, int k)
{
    return (i > 0) + (i < N - 1) + (j > 0) + (j < N - 1) + (k > 0) + (k < N - 1);
}

/* computes y = A x, where A is the negative Laplacian with Neumann walls */
void ApplyA(const double *x, double *y)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                int c = kIC(i, j, k);
                double s = 0.0;
                if (i > 0)     s += x[kIC(i - 1, j, k)];
                if (i < N - 1) s += x[kIC(i + 1, j, k)];
                if (j > 0)     s += x[kIC(i, j - 1, k)];
                if (j < N - 1) s += x[kIC(i, j + 1, k)];
                if (k > 0)     s += x[kIC(i, j, k - 1)];
                if (k < N - 1) s += x[kIC(i, j, k + 1)];
                y[c] = (Neighbours(i, j, k) * x[c] - s) / (H * H);
            }
}

void Project(double dt)
{
    /* solve lap(p) = (rho/dt) div(u*); write it as A p = b where b is neg(rho/dt) div */
    double bMean = 0.0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                double div = (u[kIU(i + 1, j, k)] - u[kIU(i, j, k)]
                            + v[kIV(i, j + 1, k)] - v[kIV(i, j, k)]
                            + w[kIW(i, j, k + 1)] - w[kIW(i, j, k)]) / H;
                rhs[kIC(i, j, k)] = -(kRho / dt) * div;
                bMean += rhs[kIC(i, j, k)];
            }
    bMean /= kCN;
    for (int c = 0; c < kCN; c++) rhs[c] -= bMean;

    /* preconditioned conjugate gradient, warm started from the previous pressure */
    double bNorm = sqrt(Dot(rhs, rhs, kCN));
    if (bNorm > 0.0) {
        ApplyA(p, cgQ);
        for (int c = 0; c < kCN; c++) cgR[c] = rhs[c] - cgQ[c];
        for (int c = 0; c < kCN; c++) {
            int i = c / (N * N), j = (c / N) % N, k = c % N;
            cgZ[c] = cgR[c] * H * H / Neighbours(i, j, k);
            cgD[c] = cgZ[c];
        }
        double rz = Dot(cgR, cgZ, kCN);

        for (int it = 0; it < kCgMax; it++) {
            ApplyA(cgD, cgQ);
            double alpha = rz / Dot(cgD, cgQ, kCN);
            for (int c = 0; c < kCN; c++) {
                p[c] += alpha * cgD[c];
                cgR[c] -= alpha * cgQ[c];
            }
            if (sqrt(Dot(cgR, cgR, kCN)) < kCgTol * bNorm) break;
            for (int c = 0; c < kCN; c++) {
                int i = c / (N * N), j = (c / N) % N, k = c % N;
                cgZ[c] = cgR[c] * H * H / Neighbours(i, j, k);
            }
            double rzNew = Dot(cgR, cgZ, kCN);
            double beta = rzNew / rz;
            rz = rzNew;
            for (int c = 0; c < kCN; c++) cgD[c] = cgZ[c] + beta * cgD[c];
        }
    }

    /* correct the velocities by subtracting the pressure gradient on interior faces */
    for (int i = 1; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                u[kIU(i, j, k)] -= dt / kRho * (p[kIC(i, j, k)] - p[kIC(i - 1, j, k)]) / H;

    for (int i = 0; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 0; k < N; k++)
                v[kIV(i, j, k)] -= dt / kRho * (p[kIC(i, j, k)] - p[kIC(i, j - 1, k)]) / H;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 1; k < N; k++)
                w[kIW(i, j, k)] -= dt / kRho * (p[kIC(i, j, k)] - p[kIC(i, j, k - 1)]) / H;
}
