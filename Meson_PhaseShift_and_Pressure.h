#ifndef __MESON_PHASESHIFT_AND_PRESSURE_H
#define __MESON_PHASESHIFT_AND_PRESSURE_H


#include "Meson_Polarization.h"
#include <string>
#include <cmath>  // For sqrt, exp, log
#include <iostream>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include <vector>

#include <complex>



double phaseshift_Diquark(double Gd, double w, double q, double T, double mu, double m, 
                          double mu_star, double delta, double m0, double mu_star0, 
                          double delta0, double thermal_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);
                          
      
double phaseshift_Pion(double Gs, double w, double q, double T, double mu, double m,
                       double mu_star, double delta, double m0, double mu_star0,
                       double delta0,  double thermal_cutoff, double Nf = 2.0, double Lambda_CutOff = 651.0);




double phaseshift_approx_Diquark(double Gd, double w, double q, double T, double mu, double m, 
                          double mu_star, double delta, double m0, double mu_star0, 
                          double delta0, double thermal_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);
                          
                          

                          
double Pressure_Integrand_Diquark(double Gd, double q, double T, double mu,
                               double m, double mu_star, double delta,
                              double m0, double mu_star0, double delta0,
                               double thermal_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);
                          
                          
double pressure_diquark(double Gd, double T, double mu, double m, double mu_star, 
                       double delta, double m0, double mu_star0, double delta0, double thermal_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);





double Pressure_Integrand_Diquark_QP(double Gd, double q, double T, double mu,
                               double m, double mu_star, double delta,
                              double m0, double mu_star0, double delta0,
                               double thermal_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);
            

double pressure_diquark_QP(double Gd, double T, double mu, double m, double mu_star, 
                       double delta, double m0, double mu_star0, double delta0, double thermal_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);
                          
                       
                       
                       
                       
                       
double Pressure_Integrand_Diquark_LD(double Gd, double q, double T, double mu,
                               double m, double mu_star, double delta,
                              double m0, double mu_star0, double delta0,
                               double thermal_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);
            
            
double pressure_diquark_LD(double Gd, double T, double mu, double m, double mu_star, 
                       double delta, double m0, double mu_star0, double delta0, double thermal_cutoff, double landau_cutoff, 
                          double Nf = 2.0, double Lambda_CutOff = 651.0);


#endif // REPI_DIQUARK_H
