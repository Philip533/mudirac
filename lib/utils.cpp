/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * utils.cpp
 *
 * Utility functions of various sorts
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "utils.hpp"
#include "constants.hpp"

/**
 * @brief  Compute the effective mass of a two-body system
 * @note   Compute the effective mass of a two-body system.
 * Treat negative masses as infinite.
 *
 * @param  m1:  First mass
 * @param  m2:  Second mass
 * @retval      Effective mass
 */
double effectiveMass(double m1, double m2) {
  double im1 = m1 > 0 ? 1.0 / m1 : 0.0;
  double im2 = m2 > 0 ? 1.0 / m2 : 0.0;
  return 1.0 / (im1 + im2);
}

/**
 * @brief  Factorial of n
 *
 * @param  n:   Argument
 * @param force_zero: If n< 0, then we return zero. This is used for CG coefficients
 * @retval      n!
 */
int factorial(int n, bool force_zero) {
  int fn = 1;

  if (n < 0) {

    if (force_zero){
      return 0;
    }else{
      throw "Factorial of a negative number";
    }
  }

  for (int i = 1; i <= n; ++i) {
    fn *= i;
  }

  return fn;
}

/**
 * @brief Sinc function, or sin(x)/x
 *
 * @param x         Argument
 * @retval          sin(x)/x
 */
double sinc(double x) {
  return x <= 0 ? 1 : sin(x) / x;
}

/**
  * @brief  Linear interpolation
  * @note   Lineral interpolation between two values with a
  * parameter t
  *
  * @param  a:      First value
  * @param  b:      Second value
  * @param  t:      Interpolation parameter (t=0 -> a, t=1 -> b)
  * @retval         Interpolated value
 */
double lerp(double a, double b, double t) {
  return a + (b - a) * t;
}

/**
 * @brief  Generate a linear grid
 * @note   Generate a linearly spaced grid between two given points.
 *
 * @param  x0: Starting point.
 * @param  x1: End point.
 * @param  n: Number of points (default = 100).
 * @retval x: Vector containing the grid.
 */
vector<double> linGrid(double x0, double x1, int n) {
  double dx = (x1 - x0) / (n - 1.0);
  vector<double> x(n);

  for (int i = 0; i < n; ++i) {
    x[i] = i * dx + x0;
  }

  return x;
}

/**
 * @brief  Generate a logarithmic grid
 * @note   Generate a logarithmically spaced grid between two given points.
 *
 * @param  x0: Starting point.
 * @param  x1: End point.
 * @param  n: Number of points (default = 100).
 * @retval grid: Vectors containing the grids, [x, exp(x)].
 */
vector<vector<double>> logGrid(double x0, double x1, int n) {
  double dlog = log(x1 / x0) / (n - 1.0);
  vector<double> x(n), ex(n);

  for (int i = 0; i < n; ++i) {
    x[i] = i * dlog;
    ex[i] = x0 * exp(x[i]);
  }

  return vector<vector<double>>({x, ex});
}

/**
 * @brief  Generate a logarithmic grid
 * @note   Generate a logarithmically spaced grid from a central point
 * and a step, ranging between the given indices.
 *
 * @param  xc:  Central point
 * @param  dx:  Logarithmic step
 * @param  i0:  Starting index
 * @param  i1:  End index
 * @retval grid: Vectors containing the grids, [x, exp(x)].
 */
vector<vector<double>> logGrid(double xc, double dx, int i0, int i1) {
  int n = i1 - i0 + 1;
  if (n <= 0) {
    throw invalid_argument("i1 must be greater or equal than i0 in logGrid");
  }
  vector<double> x(n), ex(n);

  for (int i = i0; i <= i1; ++i) {
    x[i - i0] = i * dx;
    ex[i - i0] = xc * exp(x[i - i0]);
  }

  return vector<vector<double>>({x, ex});
}

/**
 * @brief  Compute the value of a generalised Laguerre polynomial at x
 *
 * @param  x:       Argument of the polynomial
 * @param  n:       Parameter n of the polynomial
 * @param  alpha:   Parameter alpha of the polynomial
 * @retval          Value of the polynomial
 */
double genLaguerrePoly(double x, int n, double alpha) {
  // We do this recursively. Not the most efficient way, but since we only need
  // small numbers, it should be sufficient.

  double L = 1.0;

  if (n < 0) {
    throw "Parameter k of generalised Laguerre polynomial must be >= 0.";
  }

  switch (n) {
    case 0:
      return L;
      break;
    case 1:
      return 1.0 + alpha - x;
      break;
    case 2:
      return 0.5 * x * x - (alpha + 2.0) * x + 0.5 * alpha * alpha + 1.5 * alpha + 1.0;
      break;
    case 3:
      return -(x * x * x) / 6.0 + (alpha + 3.0) * x * x / 2.0 - (0.5 * alpha * alpha + 2.5 * alpha + 3) * x + (alpha + 1) * (alpha + 2) * (alpha + 3) / 6.0;
      break;
    default:
      return ((2 * n - 1 + alpha - x) * genLaguerrePoly(x, n - 1, alpha) -
              (n - 1 + alpha) * genLaguerrePoly(x, n - 2, alpha)) /
             n;
      break;
  }
}

/**
 * @brief  Compute the value of a generalised Laguerre polynomial at multiple x
 *
 * @param  x:       Arguments of the polynomial
 * @param  n:       Parameter n of the polynomial
 * @param  alpha:   Parameter alpha of the polynomial
 * @retval          Value of the polynomial
 */
vector<double> genLaguerrePoly(vector<double> x, int n, int alpha) {
  // Just iteratively go over all x
  vector<double> L(x.size());

  for (int i = 0; i < x.size(); ++i) {
    L[i] = genLaguerrePoly(x[i], n, alpha);
  }

  return L;
}

/**
 * @brief  Count the nodes in a given function
 * @note   Count the nodes (times the value crosses zero) in a
 * function. Includes only an interval comprised between two extremes that are
 * above a given tolerance.
 *
 * @param  v:   Vector containing the function
 * @param  tol: Tolerance
 * @retval      Number of nodes
 */
int countNodes(vector<double> v, double tol) {
  int nc = 0;
  int i0, i1;
  double max = 0;

  for (int i = 0; i < v.size(); ++i) {
    double y = abs(v[i]);
    if (y > max) {
      max = y;
    }
  }

  LOG(TRACE) << "Counting nodes with tolerance " << tol << ", max = " << max << ", thr = " << max * tol << "\n";

  for (i0 = 0; i0 < v.size(); ++i0) {
    if (abs(v[i0]) / max > tol) {
      break;
    }
  }
  for (i1 = v.size() - 1; i1 > i0; --i1) {
    if (abs(v[i1]) / max > tol) {
      break;
    }
  }

  LOG(TRACE) << "Counting nodes between indices " << i0 << " and " << i1 << "\n";

  for (int i = i0 + 1; i <= i1; ++i) {
    bool node = ((v[i] * v[i - 1]) < 0);
    nc += node;
    if (node) {
      LOG(TRACE) << "Node found at index " << i << "\n";
    }
  }

  return nc;
}

/**
 * @brief  Convert quantum numbers l and s to k
 * @note   Convert quantum numbers l and s as used in
 * the radial Schroedinger equation to their equivalent k
 * in the Dirac one
 *
 * @param  l:       Orbital quantum number
 * @param  s:       Spin quantum number
 * @param  &k:      Dirac quantum number
 * @retval None
 */
void qnumSchro2Dirac(int l, bool s, int &k) {
  k = ((!s && l > 0) ? l : -l - 1);
}

/**
 * @brief  Convert quantum number k to l and s
 * @note   Convert quantum number k as used in
 * the radial Dirac equation to its equivalent l
 * and s in the Schroedinger one.
 *
 * @param  k:        Dirac quantum number
 * @param  &l:       Orbital quantum number
 * @param  &s:       Spin quantum number
 * @retval None
 */
void qnumDirac2Schro(int k, int &l, bool &s) {
  s = (k < 0);
  l = s ? -(k + 1) : k;
}

/**
 * @brief  Convert from nodes and l to n
 * @note   Convert from number of nodes in the wavefunction and
 * orbital quantum number l to principal quantum number
 *
 * @param  nodes:   Number of nodes
 * @param  l:       Orbital quantum number
 * @param  &n:      Principal quantum number
 * @retval None
 */
void qnumNodes2Principal(int nodes, int l, int &n) {
  n = nodes + l + 1;
}

/**
 * @brief  Convert from n and l to nodes
 * @note   Convert from principal quantum number
 * and orbital quantum number l to number of nodes in the wavefunction
 *
 * @param  n:       Principal quantum number
 * @param  l:       Orbital quantum number
 * @param  &nodes:  Number of nodes
 * @retval None
 */
void qnumPrincipal2Nodes(int n, int l, int &nodes) {
  nodes = n - l - 1;
}

/**
 * @brief  Clebsch-Gordan coefficients for spin 1/2 particle
 * @note   Clebsch-Gordan coefficients for a spin 1/2 particle,
 * as a function of k, m (magnetic quantum number for the total spin j)
 * and spin orientation. Returns:
 *
 * c(l, 1/2; m-s, s)
 *
 * (see Development Guide for more details).
 *
 * @param  k:       Dirac quantum number
 * @param  m:       Magnetic quantum number
 * @param  s:       Spin
 * @retval Clebsch-Gordan coefficient
 */
double cgCoeff(int k, double m, bool s) {
  if (s) {
    return -k / abs(k) * sqrt((k - m + 0.5) / (2.0 * k + 1.0));
  } else {
    return sqrt((k + m + 0.5) / (2 * k + 1.0));
  }
}

/**
  * @brief Calculate a general CG coefficient
  * @brief This is a wrapper around GSL's Wigner3j coefficient routine
  * @note The arguments are multiplied by 2 as this is specified by GSL's interface
 */
double generalCgCoeff(double j1, double j2, double J, double m1, double m2, double M){

  double wigner3j = gsl_sf_coupling_3j(int(2*j1), int(2*j2), int(2*J), int(2*m1), int(2*m2), -int(2*M));
  // return first_term * second_term * third_term;
  return std::pow(-1,j1-j2+M) * std::sqrt(2*J + 1.0) * wigner3j;

}

/**
  * These matrix elements are complex valued
  * It is < psi1 | Y_1^m \alpha_x | psi2 >
  * @param k1: Initial state k value
  * @param k2: Final state k value
  * @param mu1: Initial state magnetic value
  * @param mu2: Final state magnetic value
  * @param m:    Order of the spherical harmonic
  * @param J12:  Radial integral \int P_b Q_a r
  * @param J21:  Radial integral \int P_a Q_b r
 */
std::complex<double> Y1mAlphaX(int k1, int k2, double mu1, double mu2, int m, double J12, double J21){

  int l1, l2;
  bool s1, s2;

  // We have the selection rules  mu1 + m + 1 = mu2 
  // and                          mu1 + m - 1 = mu2
  // Convert k value to l and s
  qnumDirac2Schro(k1, l1, s1);
  qnumDirac2Schro(k2, l2, s2);

  // Sign of the k values
  int sgk1 = (k1 < 0 ? -1 : 1);
  int sgk2 = (k2 < 0 ? -1 : 1);

  double prefactor1 = std::sqrt(3.0*(2.0*l1 - 2.0*sgk1 + 1.0)/((4.0*Physical::pi)*(2.0*l2 + 1.0)));

  double prefactor2 = std::sqrt(3.0*(2.0*l1 + 1.0)/((4.0*Physical::pi)*(2.0*l2 - 2.0*sgk2 + 1.0)));

  double val1 = l1-sgk1;
  double val2 = l2;
  double val3 = l1;
  double val4 = l2-sgk2;

  // These are the initial prefactors CG coeffs for each of the two terms in our expression
  double c1 = generalCgCoeff(val1,1.0,val2,0.0,0.0,0.0);
  double c2 = generalCgCoeff(val3,1.0,val4,0.0,0.0,0.0);

  // There are 8 unique CG coeffs that take the k value as an argument
  double u1 = cgCoeff(k2,mu2,true);
  double u2 = cgCoeff(-k2,mu2,false);
  double u3 = cgCoeff(k2,mu2,false);
  double u4 = cgCoeff(-k2,mu2,true);
  double v1 = cgCoeff(k1,mu1,true);
  double v2 = cgCoeff(-k1,mu1,false);
  double v3 = cgCoeff(k1,mu1,false);
  double v4 = cgCoeff(-k1,mu1,true);

  // Now we need the 4 magnetic CG coeffs
  double c3 = generalCgCoeff(l1-sgk1,1.0,l2,mu1+0.5,m,mu2-0.5);
  double c4 = generalCgCoeff(l1-sgk1,1.0,l2,mu1-0.5,m,mu2+0.5);
  double c5 = generalCgCoeff(l1,1.0,l2-sgk2,mu1+0.5,m,mu2-0.5);
  double c6 = generalCgCoeff(l1,1.0,l2-sgk2,mu1-0.5,m,mu2+0.5);

  // Multiply everything together
  double term1 = J21 * prefactor1 * c1 * (u1 * v2 * c3 + u3 * v4 * c4);

  double term2 = J12 * prefactor2 * c2 * (u4 * v3 * c5 + u2 * v1 * c6);

  // Final complex result
  std::complex<double> matel(0.0,term1-term2);

  return matel;

}

std::complex<double> Y1mAlphaY(int k1, int k2, double mu1, double mu2, int m, double J12, double J21){

  int l1, l2;
  bool s1, s2;

  // We have the selection rules  mu1 + m + 1 = mu2 
  // and                          mu1 + m - 1 = mu2
  // Convert k value to l and s
  qnumDirac2Schro(k1, l1, s1);
  qnumDirac2Schro(k2, l2, s2);

  // Sign of the k values
  int sgk1 = (k1 < 0 ? -1 : 1);
  int sgk2 = (k2 < 0 ? -1 : 1);

  double prefactor1 = std::sqrt(3.0*(2.0*l1 - 2.0*sgk1 + 1.0)/((4.0*Physical::pi)*(2.0*l2 + 1.0)));

  double prefactor2 = std::sqrt(3.0*(2.0*l1 + 1.0)/((4.0*Physical::pi)*(2.0*l2 - 2.0*sgk2 + 1.0)));

  double val1 = l1-sgk1;
  double val2 = l2;
  double val3 = l1;
  double val4 = l2-sgk2;

  // These are the initial prefactors CG coeffs for each of the two terms in our expression
  double c1 = generalCgCoeff(val1,1.0,val2,0.0,0.0,0.0);
  double c2 = generalCgCoeff(val3,1.0,val4,0.0,0.0,0.0);

  // There are 8 unique CG coeffs that take the k value as an argument
  double u1 = cgCoeff(k2,mu2,true);
  double u2 = cgCoeff(-k2,mu2,false);
  double u3 = cgCoeff(k2,mu2,false);
  double u4 = cgCoeff(-k2,mu2,true);
  double v1 = cgCoeff(k1,mu1,true);
  double v2 = cgCoeff(-k1,mu1,false);
  double v3 = cgCoeff(k1,mu1,false);
  double v4 = cgCoeff(-k1,mu1,true);

  // Now we need the 4 magnetic CG coeffs
  double c3 = generalCgCoeff(l1-sgk1,1.0,l2,mu1+0.5,m,mu2-0.5);
  double c4 = generalCgCoeff(l1-sgk1,1.0,l2,mu1-0.5,m,mu2+0.5);
  double c5 = generalCgCoeff(l1,1.0,l2-sgk2,mu1+0.5,m,mu2-0.5);
  double c6 = generalCgCoeff(l1,1.0,l2-sgk2,mu1-0.5,m,mu2+0.5);

  // Multiply everything together
  double term1 = J21 * prefactor1 * c1 * (u1 * v2 * c3 - u3 * v4 * c4);

  double term2 = J12 * prefactor2 * c2 * (u2 * v1 * c6 - u4 * v3 * c5);

  // Final complex result
  std::complex<double> matel(term1+term2, 0.0);

  return matel;

}
std::complex<double> Y1mAlphaZ(int k1, int k2, double mu1, double mu2, int m, double J12, double J21){

  int l1, l2;
  bool s1, s2;

  // Convert k value to l and s
  qnumDirac2Schro(k1, l1, s1);
  qnumDirac2Schro(k2, l2, s2);

  // Sign of the k values
  int sgk1 = (k1 < 0 ? -1 : 1);
  int sgk2 = (k2 < 0 ? -1 : 1);

  double prefactor1 = std::sqrt(3.0*(2.0*l1 - 2.0*sgk1 + 1.0)/((4.0*Physical::pi)*(2.0*l2 + 1.0)));

  double prefactor2 = std::sqrt(3.0*(2.0*l1 + 1.0)/((4.0*Physical::pi)*(2.0*l2 - 2.0*sgk2 + 1.0)));

  double val1 = l1-sgk1;
  double val2 = l2;
  double val3 = l1;
  double val4 = l2-sgk2;

  // These are the initial prefactors CG coeffs for each of the two terms in our expression
  double c1 = generalCgCoeff(val1,1.0,val2,0.0,0.0,0.0);
  double c2 = generalCgCoeff(val3,1.0,val4,0.0,0.0,0.0);

  // There are 8 unique CG coeffs that take the k value as an argument
  double u1 = cgCoeff(k2,mu2,true);
  double u2 = cgCoeff(-k2,mu2,false);
  double u3 = cgCoeff(k2,mu2,false);
  double u4 = cgCoeff(-k2,mu2,true);
  double v1 = cgCoeff(k1,mu1,true);
  double v2 = cgCoeff(-k1,mu1,false);
  double v3 = cgCoeff(k1,mu1,false);
  double v4 = cgCoeff(-k1,mu1,true);

  // Now we need the 4 magnetic CG coeffs
  double c3 = generalCgCoeff(l1-sgk1,1.0,l2,mu1-0.5,m,mu2-0.5);
  double c4 = generalCgCoeff(l1-sgk1,1.0,l2,mu1+0.5,m,mu2+0.5);
  double c5 = generalCgCoeff(l1,1.0,l2-sgk2,mu1-0.5,m,mu2-0.5);
  double c6 = generalCgCoeff(l1,1.0,l2-sgk2,mu1+0.5,m,mu2+0.5);

  // Multiply everything together
  double term1 = J21 * prefactor1 * c1 * (u1 * v4 * c3 - u3 * v2 * c4);

  double term2 = J12 * prefactor2 * c2 * (u4 * v1 * c5 - u2 * v3 * c6);

  // std::cout << "Alpha Z terms " << term1 << " " << term2 << std::endl;
  // Final complex result
  std::complex<double> matel(0.0, term1-term2);

  return matel;

}
/**
  * @brief  Parse an atomic state's quantum numbers from IUPAC notation
  * @note   Parse an atomic state's quantum numbers from IUPAC X-ray notation.
  *
  * @param  istate:     IUPAC notation string
  * @param  &n:         Principal quantum number
  * @param  &l:         Orbital quantum number
  * @param  &s:         Spin quantum number
  * @retval None
 */
void parseIupacState(string istate, int &n, int &l, bool &s) {
  int len = istate.size();

  if (len < 2) {
    if (istate == "K") {
      // A special case...
      n = 1;
      l = 0;
      s = false;
      return;
    } else {
      throw invalid_argument("istate is not a valid IUPAC notation state designation");
    }
  }
  // The string is interpreted as a letter + a number
  char shell = istate[0];
  int orbit = stoi(istate.substr(1));

  n = shell - 'J';
  l = orbit / 2;
  s = (orbit % 2) == 1;

  if (n < 1 || l > n - 1) {
    throw invalid_argument("istate is not a valid IUPAC notation state designation");
  }
}

/**
 * @brief  Return a IUPAC string describing the given state
 * @note   Return a string description of a state of given n,
 * l and s in standard IUPAC X-ray notation.
 *
 * @param  n:   Principal quantum number
 * @param  l:   Orbital quantum number
 * @param  s:   Spin quantum number
 * @retval
 */
string printIupacState(int n, int l, bool s) {
  string ans = "";

  if (n < 1 || l < 0 || l >= n) {
    throw invalid_argument("Invalid quantum numbers passed to printIupacState");
  }

  ans += 'J' + n;
  ans += to_string(l == 0 ? 1 : (2 * l + s));

  return ans;
}

/**
  * @brief  Parse a range of atomic states' quantum numbers from IUPAC notation
  * @note   Parse a range of atomic states' quantum numbers from IUPAC notation.
  *
  * @param  irange:          IUPAC notation range string (e.g.: K1:L3 will give all states of the first two shells)
  * @param  &nrange:         Principal quantum numbers
  * @param  &lrange:         Orbital quantum numbers
  * @param  &srange:         Spin quantum numbers
  * @retval None
 */
void parseIupacRange(string irange, vector<int> &nrange, vector<int> &lrange, vector<bool> &srange) {
  nrange.clear();
  lrange.clear();
  srange.clear();

  int n1, l1, n2, l2;
  bool s1, s2;

  vector<string> limits = splitString(irange, ":");
  if (limits.size() == 1) {
    // It's not a range
    parseIupacState(limits[0], n1, l1, s1);
    nrange.push_back(n1);
    lrange.push_back(l1);
    srange.push_back(s1);
  } else if (limits.size() == 2) {
    // It's a range
    parseIupacState(limits[0], n1, l1, s1);
    parseIupacState(limits[1], n2, l2, s2);

    for (int ni = n1; ni <= n2; ++ni) {
      int omin = ni == n1 ? 2 * l1 + s1 : 1;
      int omax = ni == n2 ? 2 * l2 + s2 : 2 * ni - 1;
      for (int oi = omin; oi <= omax; ++oi) {
        nrange.push_back(ni);
        lrange.push_back(oi / 2);
        srange.push_back(oi % 2 == 1 && oi != 1);
      }
    }
  } else {
    throw invalid_argument("Invalid range passed to parseIupacRange");
  }
}

/**
 * @brief  Split a string with a given separator
 * @note   Split a string s using each occurrence of a separator sep
 *
 * @param  s:       String to split
 * @param  sep:     Separator
 * @param  merge:   If true, merge multiple consecutive occurrences of sep
 * @param  maxn:    Maximum number of splits (counting from the beginning)
 * @retval          Vector of string fragments
 */
vector<string> splitString(string s, string sep, bool merge, int maxn) {
  int pos = -1;
  vector<string> fragments;

  maxn = maxn < 0 ? s.size() : maxn;

  while (s.size() > 0) {
    if (maxn == 0) {
      fragments.push_back(s);
      break;
    }
    pos = s.find_first_of(sep);
    while (pos == 0 && merge) {
      s = s.substr(1);
      pos = s.find_first_of(sep);
    }
    if (s.size() == 0)
      break;
    fragments.push_back(s.substr(0, pos));
    s = s.substr(pos + 1);
    if (pos == string::npos)
      break;
    maxn--;
  }

  return fragments;
}

/**
 * @brief   Strip a string of the given characters on both ends
 * @note    Strip a string of any instances of the given characters
 * on both ends. If no characters are passed, spaces are used.
 *
 * @param   s:          String to strip
 * @param   strip:      Characters to strip
 * @retval              Stripped string
 */
string stripString(string s, string strip) {
  int pos;
  pos = s.find_first_not_of(strip);
  if (pos == -1) {
    // All to strip
    return "";
  }
  s = s.substr(pos);
  pos = s.find_last_not_of(strip);
  s = s.substr(0, pos + 1);

  return s;
}

/**
 * @brief   Make a string uppercase
 * @note    Make a string uppercase
 *
 * @param   s:      String to modify
 * @retval          Uppercase string
 */
string upperString(string s) {
  string s2;
  for (int i = 0; i < s.size(); ++i) {
    s2.push_back(toupper(s[i]));
  }

  return s2;
}

/**
 * @brief Write a generic tabulated file with two columns
 * @note Write a generic tabulated file with two columns of data. Used mostly
 * for debug purposes.
 *
 * @param  col1:            Vector with first column of data
 * @param  col2:            Vector with second column of data
 * @param  fname:           Name of file to save
 * @retval None
 */
void writeTabulated2ColFile(vector<double> col1, vector<double> col2, string fname) {

  ofstream out(fname);

  int N = min(col1.size(), col2.size());

  for (int i = 0; i < N; ++i) {
    out << col1[i] << "\t" << col2[i] << "\n";
  }

  out.close();
}