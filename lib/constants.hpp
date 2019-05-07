/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * constants.hpp
 * 
 * Useful physical constants, atomic units
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

// All values from CODATA 2014 unless otherwise specified

#ifndef PHYSICAL_DEF
#define PHYSICAL_DEF
namespace Physical
{
const double alpha = 7.2973525664e-3; // Fine structure constant
const double c = 137.035999139;       // Speed of light

// Particle constants
const double m_e = 1.0;                    // Electron mass
const double m_mu = 1.0 / 4.83633170e-3;   // Muon mass
const double m_p = 1.0 / 5.44617021352e-4; // Proton mass

// Units
const double m = 1.0 / 5.2917721067e-11;  // Metre
const double ang = 1.0 / 5.2917721067e-1; // Angstrom
const double fm = 1.0 / 5.2917721067e4;   // Femtometre
const double eV = 1.0 / 27.211385;        // Electron volt
const double amu = 1822.888486192;        // Atomic mass unit
} // namespace Physical
#endif
