#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "fem1d.h"

using namespace std;

//****************************************************************************80

int main() {
#define NSUB 5
#define NL 2
  double adiag[NSUB + 1];
  double aleft[NSUB + 1];
  double arite[NSUB + 1];
  double f[NSUB + 1];
  double h[NSUB];
  int ibc;
  int indx[NSUB + 1];
  int node[NL * NSUB];
  int nquad;
  int nu;
  double ul;
  double ur;
  double xl;
  double xn[NSUB + 1];
  double xquad[NSUB];
  double xr;

  timestamp();

  cout << "\n";
  cout << "FEM1D\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Solve the two-point boundary value problem\n";
  cout << "\n";
  cout << "  - d/dX (P dU/dX) + Q U  =  F\n";
  cout << "\n";
  cout << "  on the interval [XL,XR], specifying\n";
  cout << "  the value of U or U' at each end.\n";
  cout << "\n";
  cout << "  The interval [XL,XR] is broken into NSUB = " << NSUB
       << " subintervals\n";
  cout << "  Number of basis functions per element is NL = " << NL << "\n";
  //
  //  Initialize the data that defines the problem.
  //
  init(&ibc, &nquad, &ul, &ur, &xl, &xr);
  //
  //  Compute the quantities which define the geometry of the
  //  problem.
  //
  geometry(h, ibc, indx, NL, node, NSUB, &nu, xl, xn, xquad, xr);
  //
  //  Assemble the linear system.
  //
  assemble(adiag, aleft, arite, f, h, indx, NL, node, nu, nquad, NSUB, ul, ur,
           xn, xquad);
  //
  //  Print out the linear system.
  //
  prsys(adiag, aleft, arite, f, nu);
  //
  //  Solve the linear system.
  //
  solve(adiag, aleft, arite, f, nu);
  //
  //  Print out the solution.
  //
  output(f, ibc, indx, NSUB, nu, ul, ur, xn);
  //
  //  Terminate.
  //
  cout << "\n";
  cout << "FEM1D:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp();

  return 0;
#undef NL
#undef NSUB
}
//****************************************************************************80

void assemble(double adiag[], double aleft[], double arite[], double f[],
              double h[], int indx[], int nl, int node[], int nu, int nquad,
              int nsub, double ul, double ur, double xn[], double xquad[]) {
  double aij;
  double he;
  int i;
  int ie;
  int ig;
  int il;
  int iq;
  int iu;
  int jg;
  int jl;
  int ju;
  double phii;
  double phiix;
  double phij;
  double phijx;
  double x;
  double xleft;
  double xquade;
  double xrite;
  //
  //  Zero out the arrays that hold the coefficients of the matrix
  //  and the right hand side.
  //
  for (i = 0; i < nu; i++) {
    f[i] = 0.0;
  }
  for (i = 0; i < nu; i++) {
    adiag[i] = 0.0;
  }
  for (i = 0; i < nu; i++) {
    aleft[i] = 0.0;
  }
  for (i = 0; i < nu; i++) {
    arite[i] = 0.0;
  }
  //
  //  For interval number IE,
  //
  for (ie = 0; ie < nsub; ie++) {
    he = h[ie];
    xleft = xn[node[0 + ie * 2]];
    xrite = xn[node[1 + ie * 2]];
    //
    //  consider each quadrature point IQ,
    //
    for (iq = 0; iq < nquad; iq++) {
      xquade = xquad[ie];
      //
      //  and evaluate the integrals associated with the basis functions
      //  for the left, and for the right nodes.
      //
      for (il = 1; il <= nl; il++) {
        ig = node[il - 1 + ie * 2];
        iu = indx[ig] - 1;

        if (0 <= iu) {
          phi(il, xquade, &phii, &phiix, xleft, xrite);
          f[iu] = f[iu] + he * ff(xquade) * phii;
          //
          //  Take care of boundary nodes at which U' was specified.
          //
          if (ig == 0) {
            x = 0.0;
            f[iu] = f[iu] - pp(x) * ul;
          } else if (ig == nsub) {
            x = 1.0;
            f[iu] = f[iu] + pp(x) * ur;
          }
          //
          //  Evaluate the integrals that take a product of the basis
          //  function times itself, or times the other basis function
          //  that is nonzero in this interval.
          //
          for (jl = 1; jl <= nl; jl++) {
            jg = node[jl - 1 + ie * 2];
            ju = indx[jg] - 1;

            phi(jl, xquade, &phij, &phijx, xleft, xrite);

            aij = he * (pp(xquade) * phiix * phijx + qq(xquade) * phii * phij);
            //
            //  If there is no variable associated with the node, then it's
            //  a specified boundary value, so we multiply the coefficient
            //  times the specified boundary value and subtract it from the
            //  right hand side.
            //
            if (ju < 0) {
              if (jg == 0) {
                f[iu] = f[iu] - aij * ul;
              } else if (jg == nsub) {
                f[iu] = f[iu] - aij * ur;
              }
            }
            //
            //  Otherwise, we add the coefficient we've just computed to the
            //  diagonal, or left or right entries of row IU of the matrix.
            //
            else {
              if (iu == ju) {
                adiag[iu] = adiag[iu] + aij;
              } else if (ju < iu) {
                aleft[iu] = aleft[iu] + aij;
              } else {
                arite[iu] = arite[iu] + aij;
              }
            }
          }
        }
      }
    }
  }
  return;
}

double ff(double x) {
  double value;

  value = 0.0;

  return value;
}

void geometry(double h[], int ibc, int indx[], int nl, int node[], int nsub,
              int *nu, double xl, double xn[], double xquad[], double xr) {
  int i;
  //
  //  Set the value of XN, the locations of the nodes.
  //
  cout << "\n";
  cout << "  Node      Location\n";
  cout << "\n";
  for (i = 0; i <= nsub; i++) {
    xn[i] = ((double)(nsub - i) * xl + (double)i * xr) / (double)(nsub);
    cout << "  " << setw(8) << i << "  " << setw(14) << xn[i] << "\n";
  }
  //
  //  Set the lengths of each subinterval.
  //
  cout << "\n";
  cout << "Subint    Length\n";
  cout << "\n";
  for (i = 0; i < nsub; i++) {
    h[i] = xn[i + 1] - xn[i];
    cout << "  " << setw(8) << i + 1 << "  " << setw(14) << h[i] << "\n";
  }
  //
  //  Set the quadrature points, each of which is the midpoint
  //  of its subinterval.
  //
  cout << "\n";
  cout << "Subint    Quadrature point\n";
  cout << "\n";
  for (i = 0; i < nsub; i++) {
    xquad[i] = 0.5 * (xn[i] + xn[i + 1]);
    cout << "  " << setw(8) << i + 1 << "  " << setw(14) << xquad[i] << "\n";
  }
  //
  //  Set the value of NODE, which records, for each interval,
  //  the node numbers at the left and right.
  //
  cout << "\n";
  cout << "Subint  Left Node  Right Node\n";
  cout << "\n";
  for (i = 0; i < nsub; i++) {
    node[0 + i * 2] = i;
    node[1 + i * 2] = i + 1;
    cout << "  " << setw(8) << i + 1 << "  " << setw(8) << node[0 + i * 2]
         << "  " << setw(8) << node[1 + i * 2] << "\n";
  }
  //
  //  Starting with node 0, see if an unknown is associated with
  //  the node.  If so, give it an index.
  //
  *nu = 0;
  //
  //  Handle first node.
  //
  i = 0;
  if (ibc == 1 || ibc == 3) {
    indx[i] = -1;
  } else {
    *nu = *nu + 1;
    indx[i] = *nu;
  }
  //
  //  Handle nodes 1 through nsub-1
  //
  for (i = 1; i < nsub; i++) {
    *nu = *nu + 1;
    indx[i] = *nu;
  }
  //
  //  Handle the last node.
  //
  i = nsub;

  if (ibc == 2 || ibc == 3) {
    indx[i] = -1;
  } else {
    *nu = *nu + 1;
    indx[i] = *nu;
  }

  cout << "\n";
  cout << "  Number of unknowns NU = " << *nu << "\n";
  cout << "\n";
  cout << "  Node  Unknown\n";
  cout << "\n";
  for (i = 0; i <= nsub; i++) {
    cout << "  " << setw(8) << i << "  " << setw(8) << indx[i] << "\n";
  }

  return;
}

void init(int *ibc, int *nquad, double *ul, double *ur, double *xl,
          double *xr) {
  //
  //  IBC declares what the boundary conditions are.
  //
  *ibc = 1;
  //
  //  NQUAD is the number of quadrature points per subinterval.
  //  The program as currently written cannot handle any value for
  //  NQUAD except 1//
  //
  *nquad = 1;
  //
  //  Set the values of U or U' at the endpoints.
  //
  *ul = 0.0;
  *ur = 1.0;
  //
  //  Define the location of the endpoints of the interval.
  //
  *xl = 0.0;
  *xr = 1.0;
  //
  //  Print out the values that have been set.
  //
  cout << "\n";
  cout << "  The equation is to be solved for\n";
  cout << "  X greater than XL = " << *xl << "\n";
  cout << "  and less than XR = " << *xr << "\n";
  cout << "\n";
  cout << "  The boundary conditions are:\n";
  cout << "\n";

  if (*ibc == 1 || *ibc == 3) {
    cout << "  At X = XL, U = " << *ul << "\n";
  } else {
    cout << "  At X = XL, U' = " << *ul << "\n";
  }

  if (*ibc == 2 || *ibc == 3) {
    cout << "  At X = XR, U = " << *ur << "\n";
  } else {
    cout << "  At X = XR, U' = " << *ur << "\n";
  }

  cout << "\n";
  cout << "  Number of quadrature points per element is " << *nquad << "\n";

  return;
}

void output(double f[], int ibc, int indx[], int nsub, int nu, double ul,
            double ur, double xn[]) {
  int i;
  double u;

  cout << "\n";
  cout << "  Computed solution coefficients:\n";
  cout << "\n";
  cout << "  Node    X(I)        U(X(I))\n";
  cout << "\n";

  for (i = 0; i <= nsub; i++) {
    //
    //  If we're at the first node, check the boundary condition.
    //
    if (i == 0) {
      if (ibc == 1 || ibc == 3) {
        u = ul;
      } else {
        u = f[indx[i] - 1];
      }
    }
    //
    //  If we're at the last node, check the boundary condition.
    //
    else if (i == nsub) {
      if (ibc == 2 || ibc == 3) {
        u = ur;
      } else {
        u = f[indx[i] - 1];
      }
    }
    //
    //  Any other node, we're sure the value is stored in F.
    //
    else {
      u = f[indx[i] - 1];
    }

    cout << "  " << setw(8) << i << "  " << setw(8) << xn[i] << "  " << setw(14)
         << u << "\n";
  }

  return;
}

void phi(int il, double x, double *phii, double *phiix, double xleft,
         double xrite) {
  if (xleft <= x && x <= xrite) {
    if (il == 1) {
      *phii = (xrite - x) / (xrite - xleft);
      *phiix = -1.0 / (xrite - xleft);
    } else {
      *phii = (x - xleft) / (xrite - xleft);
      *phiix = 1.0 / (xrite - xleft);
    }
  }
  //
  //  If X is outside of the interval, just set everything to 0.
  //
  else {
    *phii = 0.0;
    *phiix = 0.0;
  }

  return;
}

double pp(double x) {
  double value;

  value = 1.0;

  return value;
}

void prsys(double adiag[], double aleft[], double arite[], double f[], int nu) {
  int i;

  cout << "\n";
  cout << "Printout of tridiagonal linear system:\n";
  cout << "\n";
  cout << "Equation  ALEFT  ADIAG  ARITE  RHS\n";
  cout << "\n";

  for (i = 0; i < nu; i++) {
    cout << "  " << setw(8) << i + 1 << "  " << setw(14) << aleft[i] << "  "
         << setw(14) << adiag[i] << "  " << setw(14) << arite[i] << "  "
         << setw(14) << f[i] << "\n";
  }

  return;
}

double qq(double x) {
  double value;

  value = 0.0;

  return value;
}

void solve(double adiag[], double aleft[], double arite[], double f[], int nu) {
  int i;
  //
  //  Carry out Gauss elimination on the matrix, saving information
  //  needed for the backsolve.
  //
  arite[0] = arite[0] / adiag[0];

  for (i = 1; i < nu - 1; i++) {
    adiag[i] = adiag[i] - aleft[i] * arite[i - 1];
    arite[i] = arite[i] / adiag[i];
  }
  adiag[nu - 1] = adiag[nu - 1] - aleft[nu - 1] * arite[nu - 2];
  //
  //  Carry out the same elimination steps on F that were done to the
  //  matrix.
  //
  f[0] = f[0] / adiag[0];
  for (i = 1; i < nu; i++) {
    f[i] = (f[i] - aleft[i] * f[i - 1]) / adiag[i];
  }
  //
  //  And now carry out the steps of "back substitution".
  //
  for (i = nu - 2; 0 <= i; i--) {
    f[i] = f[i] - arite[i] * f[i + 1];
  }
  return;
}

void timestamp() {
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  cout << time_buffer << "\n";

  return;
#undef TIME_SIZE
}
