/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * auger.cpp
 *
 * Routines for calculating Auger electron emission rates 
 *
 * @author Philip Jones
 * @version 1.0 14/04/25
 */

#include "auger.hpp"

/**
 * Calculate the electron integral in the Auger rate in
 * Equation A.18 of Akylas' thesis for L = 1
 *
 * @param x_axis: radial grid values
 * @param bound_elec: bound electron values on grid
 * @paran unbound_elec: unbound electron values on grid
 *
 */
double electronAugerDipole(vector<double> x_axis, vector<double> bound_elec, vector<double> unbound_elec){

  vector<double> product_state = productWvfn(bound_elec, unbound_elec);

  double dipole_integral;

  for(int i = 0; i< product_state.size();i++){
    // std::cout << x_axis[i] << " " << bound_elec[i] << " "<< unbound_elec[i] << " " << product_state[i] << std::endl;
  }
  dipole_integral = trapzInt(x_axis, product_state);

  return dipole_integral;
}

/**
 * Calculate the electron integral in the Auger rate in
 * Equation A.18 of Akylas' thesis for L = 2
 *
 * @param x_axis: radial grid values
 * @param bound_elec: bound electron values on grid
 * @paran unbound_elec: unbound electron values on grid
 *
 */
double electronAugerQuadrupole(vector<double> x_axis, vector<double> bound_elec, vector<double> unbound_elec){

  int N = x_axis.size();

  vector<double> product_state = productWvfn(bound_elec, unbound_elec);
  vector<double> integrand(N,0.0);

  double quadrupole_rate;

  // We need to divide by r for the quadrupole rate
  for(int i = 0; i < N; i++){
    integrand[i] = product_state[i] / x_axis[i];
  }

  quadrupole_rate = trapzInt(x_axis, integrand);

  return quadrupole_rate;
}

/**
 * Multiply any two initial and final states 
 * in the same coordinate i.e both muon or both electron
 * They have to be on the same grid.
 */
vector<double> productWvfn(vector<double> initial_state, vector<double> final_state){

  int N = initial_state.size();
  if (final_state.size() != N) {
    throw invalid_argument("Invalid size for arrays passed to ProductWvfn");
  }

  vector<double> product(N,0.0);

  for(int i = 0; i < N; i++){
    product[i] = initial_state[i] * final_state[i];
  }

  return product;
}

/** Calculates the angular integrals needed for the Auger rate.
 * These include spinor spherical harmonics for the muon, so we retain 
 * a k-dependence in the CG coefficients.
 * @param L : Multipolarity of transition
 * @param li: Initial orbital quantum number of muon
 * @param lf: Final orbital quantum number of muon
 * @param l : Angular momentum of ejected electron
 * @param lp: Angular momentum of bound electron
 * @param m : Ejected electron magnetic number
 * @param mp: Bound electron magnetic number
 * @param mi: Initial muon magnetic number
 * @param mf: Final muon magnetic number
 * @param si: Initial muon spin
 * @param sf: Final muon spin 
 */
double augerAngularIntegrals(int L, int li, int lf, int l, int lp, int m, int mp, double mi, double mf, bool si, bool sf){


  using namespace Physical;

  int ki, kf;

  // Move to the Dirac numbers
  qnumSchro2Dirac(li, si, ki);
  qnumSchro2Dirac(lf, sf, kf);

  double total_sum = 0.0;
  double coeff1 = generalCgCoeff(L, l, lp, 0,0,0);

  // This is the prefactor that's independent of M,
  // and accounts for both the electrons and the muon
  double prefactor = coeff1*std::sqrt((2*li+1)/(2*lf+1))*std::sqrt((2*l+1)/((2*lp+1)));

  // M independent CG coefficient for the muon
  double coeff5 = generalCgCoeff(L,li,lf,0.0,0.0,0.0);
  double electron_sum = 0.0;

  // First we need to loop over M
  for(int M = -L; M <= L; M++){
    
    // Parity 
    double parity = std::pow(-1, M);

    // M dependent electron coupling
    double coeff2 = generalCgCoeff(L,l,lp,M,m,mp);

    // Muon couplings
    double coeff3 = cgCoeff(kf,mf,true);
    double coeff4 = cgCoeff(ki,mi,true);
    double coeff6 = generalCgCoeff(L,li,lf,-M,mi-0.5,mf-0.5);
    double coeff7 = cgCoeff(kf,mf,false);
    double coeff8 = cgCoeff(ki,mi,false);
    double coeff10 = generalCgCoeff(L,li,lf,-M,mi+0.5,mf+0.5);

    // Combine it all together
    double term1 = coeff3 * coeff4 * coeff6;
    double term2 = coeff7 * coeff8 * coeff10;

    electron_sum += parity * coeff2;
    // Perform the sum
    total_sum += parity * coeff2 * coeff5*(term1 + term2);
  
  }

  return total_sum*prefactor; 

}
