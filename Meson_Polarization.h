#ifndef __MESON_POLARIZATION_H
#define __MESON_POLARIZATION_H

#include <limits> 
#include <string>
#include <cmath>  // For sqrt, exp, log
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <omp.h>
// Constants
const double pi = 3.141592653589793;



// Modified Heaviside step function
double heaviside(double x);

// Fermi-Dirac distribution function

double n(double E, double T);


double ImPi_Pion(double w, double q, double T, double p, double m, double mu_star, double Delta, double thermal_cutoff, double Nf = 2.0, double Lambda_CutOff = 651.0);


double ImPi_Diquark(double w, double q, double T, double p, double m, double mu_star, double Delta, double thermal_cutoff, double Nf = 2.0, double Lambda_CutOff = 651.0);


double RePi_Pion(double w, double q, double T, double mu, double m,
                                        double mu_star, double delta, double m0,
                                        double mu_star0, double delta0,
                                        double thermal_cutoff,
                                        double Nf = 2.0, double Lambda_CutOff = 651.0);
                                        
                                        
                                        
double RePi_Diquark(double w, double q, double T, double mu, double m,
                                        double mu_star, double delta, double m0,
                                        double mu_star0, double delta0,
                                        double thermal_cutoff,
                                        double Nf = 2.0, double Lambda_CutOff = 651.0);
                                        
                  

#endif // REPI_DIQUARK_H
