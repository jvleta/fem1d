#pragma once

int main();
void assemble(double adiag[], double aleft[], double arite[], double f[],
              double h[], int indx[], int nl, int node[], int nu, int nquad,
              int nsub, double ul, double ur, double xn[], double xquad[]);
double ff(double x);
void geometry(double h[], int ibc, int indx[], int nl, int node[], int nsub,
              int *nu, double xl, double xn[], double xquad[], double xr);
void init(int *ibc, int *nquad, double *ul, double *ur, double *xl, double *xr);
void output(double f[], int ibc, int indx[], int nsub, int nu, double ul,
            double ur, double xn[]);
void phi(int il, double x, double *phii, double *phiix, double xleft,
         double xrite);
double pp(double x);
void prsys(double adiag[], double aleft[], double arite[], double f[], int nu);
double qq(double x);
void solve(double adiag[], double aleft[], double arite[], double f[], int nu);
void timestamp();