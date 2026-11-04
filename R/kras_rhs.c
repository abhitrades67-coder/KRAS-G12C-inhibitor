/* kras_rhs.c — compiled deSolve right-hand side for the KRAS G12C QSP model.
 * Mirrors qsp_rhs() in R/model.R exactly; ~20x faster. States and parameters are passed
 * in the fixed orders below (see R/model.R make_p/make_cpar). Build: R CMD SHLIB kras_rhs.c
 */
#include <R.h>
#include <math.h>

/* parameter vector (order fixed; filled by initmod via odeparms) */
static double p[20];
#define kact_base  p[0]
#define kgap       p[1]
#define Kd_sig     p[2]
#define P0         p[3]
#define fb         p[4]
#define kfb        p[5]
#define kprog      p[6]
#define kprog_doce p[7]
#define kdeg       p[8]
#define ksyn       p[9]
#define kcov_eff   p[10]
#define cn         p[11]
#define kprol      p[12]
#define kdeath0    p[13]
#define kkill      p[14]
#define dep        p[15]
#define bypass     p[16]
#define kres       p[17]
#define doce_res   p[18]
#define edir_amp   p[19]

static double clipd(double x, double lo, double hi){ return x<lo?lo:(x>hi?hi:x); }

void initmod(void (*odeparms)(int *, double *)) {
  int n = 20; odeparms(&n, p);
}

/* y = {Rnc, Rb, U, res, Ts, Tr}; ydot same order; yout = {S, Ss, a, edir, occ, res} */
void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
  double Rnc = y[0] < 0 ? 0 : y[0];
  double Rb  = y[1];
  double U   = y[2] < 0 ? 0 : y[2];
  double res = clipd(y[3], 0.0, 1.0);
  double Ts  = y[4] < 0 ? 0 : y[4];
  double Tr  = y[5] < 0 ? 0 : y[5];

  double a  = (kact_base * U) / (kact_base * U + kgap);
  double Rg = (1.0 - a) * Rnc;
  double Ra = a * Rnc;

  double P = Ra / (Ra + Kd_sig);
  double S = clipd(P / P0, 0.0, 1.6);

  double Utarget = 1.0 + fb * (1.0 - S > 0 ? 1.0 - S : 0.0);
  double dU = kfb * (Utarget - U);

  double capture = kcov_eff * cn * Rg;
  double occ = Rb / (Rnc + Rb + 1e-12);
  double dRnc = ksyn - kdeg * Rnc - capture;
  double dRb  = capture - kdeg * Rb;
  double dres = kprog * occ * (1.0 - res);

  double edir = edir_amp * exp(-kprog_doce * (*t));

  double Ss = clipd(S + (1.0 - S) * res,    0.0, 1.6);
  double Sr = clipd(S + (1.0 - S) * bypass, 0.0, 1.6);
  double gS = kprol * (dep * Ss + (1.0 - dep)) - (kdeath0 + kkill * dep * (1.0 - Ss) + edir);
  double gR = kprol * (dep * Sr + (1.0 - dep)) - (kdeath0 + kkill * dep * (1.0 - Sr) + edir * doce_res);
  double dTs = gS * Ts - kres * Ts;
  double dTr = gR * Tr + kres * Ts;

  ydot[0] = dRnc; ydot[1] = dRb; ydot[2] = dU; ydot[3] = dres; ydot[4] = dTs; ydot[5] = dTr;
  yout[0] = S; yout[1] = Ss; yout[2] = a; yout[3] = edir; yout[4] = occ; yout[5] = res;
}
