/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * integrate.hpp
 *
 * Routines for integrating differential equations with the shooting method - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <math.h>
#include <vector>
#include <stdexcept>
#include "utils.hpp"
#include "hydrogenic.hpp"
#include "state.hpp"
#include "integrate.hpp"
#include "constants.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

double electronAugerDipole(vector<double> x_axis, vector<double> bound_elec, vector<double> unbound_elec);
double electronAugerQuadrupole(vector<double> x_axis, vector<double> bound_elec, vector<double> unbound_elec);
vector<double> productWvfn(vector<double> initial_state, vector<double> final_state);
double augerAngularIntegrals(int L, int li, int lf, int l, int lp, int m, int mp, double mi, double mf, bool s1, bool s2);
