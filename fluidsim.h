#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#ifndef N
#define N 48                      /* how many grid cells on each side */
#endif
#define kDom 0.20                  /* how long one side of the box is, in metres */
#define H   (kDom / N)             /* size of one cell, in metres */

/* air properties, roughly room temperature at one atmosphere */
#define kRho     1.204
#define kNu      1.516e-5
#define kAlpha   2.2e-5
#define kGrav    9.80665
#define kTAmb   293.15
#define kBeta    (1.0 / kTAmb)
#define kTHeat  313.15

/* time stepping constants */
#define kCfl     0.5
#define kDtMax  1.0e-3
#define kCgTol  1.0e-7
#define kCgMax  500

#define kUN ((N + 1) * N * N)     /* u lives on x faces, so i runs from 0 to N while j and k run from 0 to N-1 */
#define kVN (N * (N + 1) * N)     /* v lives on y faces, so j runs from 0 to N */
#define kWN (N * N * (N + 1))     /* w lives on z faces, so k runs from 0 to N */
#define kCN (N * N * N)           /* p and T live at cell centers */

#define kIC(i, j, k) ((((i) * N) + (j)) * N + (k))        /* index into a cell centered array */
#define kIU(i, j, k) ((((i) * N) + (j)) * N + (k))        /* index into the u array */
#define kIV(i, j, k) ((((i) * (N + 1)) + (j)) * N + (k))  /* index into the v array */
#define kIW(i, j, k) ((((i) * N) + (j)) * (N + 1) + (k))  /* index into the w array */

/* global state */
extern double u[kUN], v[kVN], w[kWN], T[kCN], p[kCN];
extern double tu[kUN], tv[kVN], tw[kWN], tc[kCN];
extern double rhs[kCN], cgR[kCN], cgZ[kCN], cgD[kCN], cgQ[kCN];
extern int heaterOn, buoyancyOn;

/* Util */
double SampU(double x, double y, double z);
double SampV(double x, double y, double z);
double SampW(double x, double y, double z);
double SampC(const double *a, double x, double y, double z);
double Umax(void);
double PickDt(void);

/* Advection */
void AdvectScalar(double *a, double dt);
void AdvectAll(double dt);

/* Diffusion */
void DiffuseVelocity(double dt);
void DiffuseScalar(double *a, double D, double dt);

/* Forcing */
void ApplyBuoyancy(double dt);
void ApplyHeater(void);

/* Projection */
int Neighbours(int i, int j, int k);
void ApplyA(const double *x, double *y);
void Project(double dt);

/* Diagnostics */
double Kinetic(void);
double Maxdiv(void);

/* Setup / Tests */
void InitNoise(void);
int RunTests(void);

/* Time step */
void Step(double dt);

#endif
