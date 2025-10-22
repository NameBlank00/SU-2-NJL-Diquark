#include "Meson_PhaseShift_and_Pressure.h"



#include <limits>

#include <string>
#include <algorithm>
#include <atomic>

#include <map>

#include "Meson_Polarization.h"



template <typename T> int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}



double phaseshift_Diquark(double Gd, double w, double q, double T, double mu, double m, double mu_star, double delta, double m0, double mu_star0, double delta0,  double thermal_cutoff, double Nf, double Lambda_CutOff) {



    double rePi = RePi_Diquark(w, q, T, mu, m, mu_star, delta,
                               m0, mu_star0, delta0,
                                thermal_cutoff, Nf, Lambda_CutOff);

    double imPi = 0.0;

    
    if (thermal_cutoff == 1) {
        imPi = ImPi_Diquark(w, q, T, 1, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff);
    } else if (thermal_cutoff == 0){
        imPi = ImPi_Diquark(w, q, 0, 0, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff) * heaviside(4*Lambda_CutOff*Lambda_CutOff + 4*m*m - (w + 2*mu_star)*(w + 2*mu_star) + q * q )
                                + ImPi_Diquark(w, q, T, 1, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff) - ImPi_Diquark(w, q, 0, 0, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff);  // ImPi(0)*cutoff + ImPi(T) - ImPi(0)
    }


    const double reDen = 1/(2*Gd) - rePi;
    const double imDen = -imPi ;


    double phi = std::atan( -imDen / reDen );

    if (reDen <= 0.0 && w > 0) phi += pi;
    if (reDen <= 0.0 && w < 0) phi -= pi;

    return phi;
}





double phaseshift_Pion(double Gs, double w, double q, double T, double mu, double m, double mu_star, double delta, double m0, double mu_star0, double delta0,  double thermal_cutoff, double Nf, double Lambda_CutOff) {



    double rePi = RePi_Pion(w, q, T, mu, m, mu_star, delta,
                               m0, mu_star0, delta0,
                                thermal_cutoff, Nf, Lambda_CutOff);

    double imPi = 0.0;

    
    if (thermal_cutoff == 1) {
        imPi = ImPi_Pion(w, q, T, 1, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff);
    } else if (thermal_cutoff == 0){
        imPi = ImPi_Pion(w, q, 0, 0, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff) * heaviside(4*Lambda_CutOff*Lambda_CutOff + 4*m*m - w*w + q * q )
                                + ImPi_Pion(w, q, T, 1, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff) - ImPi_Pion(w, q, 0, 0, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff);
    }


    const double reDen = 1/(2*Gs) - rePi;
    const double imDen = -imPi ;


    double phi = std::atan( -imDen / reDen );

    if (reDen <= 0.0 && w > 0) phi += pi;
    if (reDen <= 0.0 && w < 0) phi -= pi;

    return phi;
}







// Lorentz Boosted Approximation
double phaseshift_approx_Diquark(double Gd, double w, double q, double T, double mu, double m, double mu_star, double delta, double m0, double mu_star0, double delta0,  double thermal_cutoff, double Nf, double Lambda_CutOff) {

    if (w*w - q*q <= 0) return 0.0;

    double E = std::sqrt(w * w - q * q);

    double rePi = RePi_Diquark(E, 0.0, T, mu, m, mu_star, delta,
                               m0, mu_star0, delta0,
                                thermal_cutoff, Nf, Lambda_CutOff);
    double imPi = ImPi_Diquark(E, 0.0, T, 1, m, mu_star, delta,
                                thermal_cutoff, Nf, Lambda_CutOff);

    // Denominator pieces (include small regulator)
    const double reDen = 1/(2*Gd) - rePi;
    const double imDen = -imPi ;

    // EXACTLY: phi = atan(imD/reD) with your objects, but without denomSq
    double phi = std::atan( -imDen / reDen );

    // EXACTLY: if (reD >= 0) <=> if (reDen <= 0)
    if (reDen <= 0.0) phi += pi;

    return phi;
}





struct Params {
    double Gd, T, mu, m, mu_star, delta, m0, mu_star0, delta0, Nf, Lambda_CutOff;
    double thermal_cutoff; // pointer to avoid copying the string per-iteration
    double q;
};



// double phaseshift_diquark_integrand_w(double w, void *params_void) {
//     Params* p = static_cast<Params*>(params_void);

//     // small-T guard
//     const double T_zero_thresh = 1e-12;
//     if (p->T <= T_zero_thresh) return 0.0;

//     // stable Bose factor handling (does not involve thermal_cutoff)
//     double x = w / p->T;
//     if (x > 700.0) return 0.0;                // exp overflow guard -> integrand ~ 0
//     double Bose = std::exp(x) - 1.0;         

//     // Now call phaseshift with thermal_cutoff (only relevant for phase shift, not Bose)
//     double shift = phaseshift_Diquark(p->Gd, w, p->q, p->T, p->mu, p->m,
//                                       p->mu_star, p->delta, p->m0, p->mu_star0, p->delta0,
//                                       p->thermal_cutoff, p->Nf, p->Lambda_CutOff);

//     const double prefactor = 3.0 / (2.0 * pi * pi * pi);
//     return prefactor * p->q * p->q * shift / Bose; // Bose factor unaffected by thermal_cutoff

// }
















double phaseshift_diquark_integrand_w(double w, void *params_void) {
    Params* p = static_cast<Params*>(params_void);

    // small-T guard
    const double T_zero_thresh = 1e-12;
    if (p->T <= T_zero_thresh) return 0.0;

    // stable Bose factor handling (does not involve thermal_cutoff)
    double x = abs(w) / p->T;
    if (x > 700.0) return 0.0;                // exp overflow guard -> integrand ~ 0
    double Bose = sgn(w)*(std::exp(x) - 1.0);   
      

    // Now call phaseshift with thermal_cutoff (only relevant for phase shift, not Bose)
    double shift = phaseshift_Diquark(p->Gd, w, p->q, p->T, p->mu, p->m,
                                      p->mu_star, p->delta, p->m0, p->mu_star0, p->delta0,
                                      p->thermal_cutoff, p->Nf, p->Lambda_CutOff);

    const double prefactor = 1/(2*pi);
    return prefactor * shift / Bose; // Bose factor unaffected by thermal_cutoff

}








// --- helper for coarse scanning (no Bose factor) ---
double coarse_probe(double w, Params* p) {
    return std::abs(
        phaseshift_Diquark(p->Gd, w, p->q, p->T, p->mu, p->m,
                           p->mu_star, p->delta,
                           p->m0, p->mu_star0, p->delta0,
                           p->thermal_cutoff, p->Nf, p->Lambda_CutOff)
    );

}



double coarse_probe_with_Bose(double w, Params* p) {
    // Small-T guard
    const double T_zero_thresh = 1e-12;
    if (p->T <= T_zero_thresh) return 0.0;

    // Bose factor calculation
    double x = w / p->T;
    if (x > 700.0) return 0.0; // exp overflow guard -> integrand ~ 0
    double Bose = std::exp(x) - 1.0;

    // Now call phaseshift with thermal_cutoff
    double shift = phaseshift_Diquark(p->Gd, w, p->q, p->T, p->mu, p->m,
                                      p->mu_star, p->delta, p->m0, p->mu_star0, p->delta0,
                                      p->thermal_cutoff, p->Nf, p->Lambda_CutOff);

    return std::abs(shift / Bose); // Include the Bose factor in the result
}


double coarse_probe_with_switch(double w, Params* p) {
    // Check if thermal_cutoff is non-zero
    if (p->thermal_cutoff > 0) {
        // Call the original coarse_probe (without Bose factor)
        return coarse_probe(w, p);  // Using the coarse_probe without Bose
    } else {
        // Call the modified coarse_probe that includes the Bose factor
        return coarse_probe_with_Bose(w, p);  // Using the coarse_probe with Bose factor
    }
}


double Pressure_Integrand_Diquark(double Gd, double q, double T, double mu,
                                         double m, double mu_star, double delta,
                                         double m0, double mu_star0, double delta0,
                                         double thermal_cutoff, double Nf, double Lambda_CutOff)
                                         
{
    const double epsabs = 1e-6;
    const double epsrel = 1e-6;
    const size_t limit  = 100;

    Params p { Gd, T, mu, m, mu_star, delta,
               m0, mu_star0, delta0,
               Nf, Lambda_CutOff,
               thermal_cutoff, q };

    gsl_function F;
    F.function = &phaseshift_diquark_integrand_w;
    F.params   = &p;

    gsl_integration_workspace *work = gsl_integration_workspace_alloc(limit);

    double total_result = 0.0;

    // --- seeds near known nonzero regions ---
    std::vector<double> seeds;
    if (std::abs(mu_star) < 1e-12) {
        seeds = { std::sqrt(4*Lambda_CutOff*Lambda_CutOff + 4*m*m + q*q) };
    } else {
        seeds = {
            -std::sqrt(4*Lambda_CutOff*Lambda_CutOff + 4*m*m + q*q) - 2*mu_star,
             std::sqrt(4*Lambda_CutOff*Lambda_CutOff + 4*m*m + q*q) - 2*mu_star
        };
    }

    // --- coarse scan parameters ---
    const double dw        = 10 - 5*(heaviside(100 - q) + heaviside(mu_star - 5) - heaviside(100 - q)*heaviside(mu_star - 5));
    const double threshold = 1e-6;   
    const int N_consec     = 5;

    std::vector<std::pair<double,double>> active_regions;
    std::vector<std::pair<double,double>> explored;

    auto already_scanned = [&](double x) {
        for (auto [a,b] : explored) {
            if (x >= a && x <= b) return true;
        }
        return false;
    };

    auto expand_from_seed = [&](double seed) {
        if (already_scanned(seed)) return;

        int zeros = 0;
        double a = seed;
        double region_start = a;

        // expand in + direction
        while (zeros < N_consec) {
            double mid = a + 0.5*dw;
            if (already_scanned(mid)) break;

            double fmid = coarse_probe(mid, &p);

            if (fmid < threshold) zeros++;
            else { zeros = 0; }

            if (zeros == N_consec && region_start < a) {
                a += 2*dw; // include 2 bins past zeros
                active_regions.push_back({region_start, a});
                explored.push_back({region_start, a});
                break;
            }
            a += dw;
        }

        // expand in – direction
        zeros = 0;
        a = seed;
        region_start = a;

        while (zeros < N_consec) {
            double mid = a - 0.5*dw;
            if (already_scanned(mid)) break;

            double fmid = coarse_probe(mid, &p);

            if (fmid < threshold) zeros++;
            else { zeros = 0; }

            if (zeros == N_consec && region_start > a) {
                a -= 2*dw; // include 2 bins past zeros
                active_regions.push_back({a, region_start});
                explored.push_back({a, region_start});
                break;
            }
            a -= dw;
        }
    };

    for (double seed : seeds)
        expand_from_seed(seed);

    // --- integrate active regions ---
    const double w_cutoff = 1e-8;

    for (auto [a,b] : active_regions) {
        if (a < 0 && b > 0) {
            double r1=0, e1=0, r2=0, e2=0;
            gsl_integration_qag(&F, a, -w_cutoff, epsabs, epsrel, limit,
                                GSL_INTEG_GAUSS21, work, &r1, &e1);
            gsl_integration_qag(&F, w_cutoff, b, epsabs, epsrel, limit,
                                GSL_INTEG_GAUSS21, work, &r2, &e2);
            total_result += r1 + r2;
        } else if (b <= -w_cutoff || a >= w_cutoff) {
            double r=0, e=0;
            gsl_integration_qag(&F, a, b, epsabs, epsrel, limit,
                                GSL_INTEG_GAUSS21, work, &r, &e);
            total_result += r;
        }
    }

    gsl_integration_workspace_free(work);

    // symmetry doubling
    if (std::abs(mu_star) < 1e-12) {
        total_result *= 2.0;
    }

    return total_result;
}




double pressure_diquark_Test(double Gd, double T, double mu, double m, double mu_star, double delta,
                        double m0, double mu_star0, double delta0,
                        double thermal_cutoff, double Nf, double Lambda_CutOff) {

    auto integrand_q = [&](double q) -> double {
        return 3.0 / (2.0 * pi * pi) * q * q * Pressure_Integrand_Diquark(
            Gd, q, T, mu, m, mu_star, delta,
            m0, mu_star0, delta0,
            thermal_cutoff, Nf, Lambda_CutOff
        );
    };

    gsl_integration_workspace* work = gsl_integration_workspace_alloc(100);

    gsl_function F;
    F.function = [](double q, void* params) -> double {
        auto* f = static_cast<decltype(integrand_q)*>(params);
        return (*f)(q);
    };
    F.params = &integrand_q;

    double result = 0.0, error = 0.0;

    int status = gsl_integration_qag(&F, 0.0, 5e3, 1e-3, 1e-3, 10,
                                      GSL_INTEG_GAUSS15, work, &result, &error);

    if (status != GSL_SUCCESS || std::isnan(result)) {
        std::cerr << "[WARNING] GSL integration failed in pressure_diquark_Test:\n";
        std::cerr << "    status = " << gsl_strerror(status) << "\n";
        std::cerr << "    result = " << result << ", error = " << error << "\n";
        std::cerr << "    parameters: "
                  << " Gd=" << Gd << " T=" << T << " mu=" << mu
                  << " m=" << m << " mu*=" << mu_star << " delta=" << delta
                  << " thermal_cutoff=" << thermal_cutoff
                  << " Nf=" << Nf << " Lambda=" << Lambda_CutOff
                  << std::endl;
    }

    gsl_integration_workspace_free(work);
    return result;
}



//--------------------------    Quasi-Particle Pressure

double Pressure_Integrand_Diquark_QP(double Gd, double q, double T, double mu,
                                         double m, double mu_star, double delta,
                                         double m0, double mu_star0, double delta0,
                                         double thermal_cutoff, double Nf, double Lambda_CutOff)
{
    const double epsabs = 1e-9;
    const double epsrel = 1e-9;
    const size_t limit  = 5000;

    double A = -q - 2*mu_star;  //A < B
    double B = q - 2*mu_star;

    Params p { Gd, T, mu, m, mu_star, delta,
               m0, mu_star0, delta0,
               Nf, Lambda_CutOff,
               thermal_cutoff, q };

    gsl_function F;
    F.function = &phaseshift_diquark_integrand_w;
    F.params   = &p;

    gsl_integration_workspace *work = gsl_integration_workspace_alloc(limit);

    double total_result = 0.0;

    // --- seeds near known nonzero regions ---
    std::vector<double> seeds;

    seeds = {
        -std::sqrt(0*Lambda_CutOff*Lambda_CutOff + 4*m*m + q*q) - 2*mu_star,
         std::sqrt(0*Lambda_CutOff*Lambda_CutOff + 4*m*m + q*q) - 2*mu_star
        };

    // --- coarse scan parameters ---
    const double dw        = 5;
    const double threshold = 1e-6;   
    const int N_consec     = 5;

    std::vector<std::pair<double,double>> active_regions;
    std::vector<std::pair<double,double>> explored;

    auto already_scanned = [&](double x) {
        for (auto [a,b] : explored) {
            if (x >= a && x <= b) return true;
        }
        return false;
    };

    auto expand_from_seed = [&](double seed) {
        if (already_scanned(seed)) return;

        int zeros = 0;
        double a = seed;
        double region_start = a;

        // expand in + direction
        while (zeros < N_consec) {
            double mid = a + 0.5*dw;
            if (already_scanned(mid)) break;

            double fmid = coarse_probe_with_switch(mid, &p);

            if (fmid < threshold) zeros++;
            else { zeros = 0; }

            if (zeros == N_consec && region_start < a) {
                a += 2*dw; // include 2 bins past zeros
                active_regions.push_back({region_start, a});
                explored.push_back({region_start, a});
                break;
            }
            a += dw;
        }

        // expand in – direction
        zeros = 0;
        a = seed;
        region_start = a;

        while (zeros < N_consec) {
            double mid = a - 0.5*dw;
            if (already_scanned(mid)) break;

            double fmid = coarse_probe_with_switch(mid, &p);

            if (fmid < threshold) zeros++;
            else { zeros = 0; }

            if (zeros == N_consec && region_start > a) {
                a -= 2*dw; // include 2 bins past zeros
                active_regions.push_back({a, region_start});
                explored.push_back({a, region_start});
                break;
            }
            a -= dw;
        }
    };

    for (double seed : seeds)
        expand_from_seed(seed);

    // --- integrate active regions, skipping (A, B) ---
    const double w_cutoff = 1e-8;

    for (auto [a,b] : active_regions) {
        // Skip regions that overlap with (A, B)
        if (b <= A || a >= B) { 
            // If the region is entirely outside of (A, B), integrate normally
            if (a < 0 && b > 0) {
                double r1=0, e1=0, r2=0, e2=0;
                gsl_integration_qag(&F, a, -w_cutoff, epsabs, epsrel, limit,
                                    GSL_INTEG_GAUSS21, work, &r1, &e1);
                gsl_integration_qag(&F, w_cutoff, b, epsabs, epsrel, limit,
                                    GSL_INTEG_GAUSS21, work, &r2, &e2);
                total_result += r1 + r2;
            } else if (b <= -w_cutoff || a >= w_cutoff) {
                double r=0, e=0;
                gsl_integration_qag(&F, a, b, epsabs, epsrel, limit,
                                    GSL_INTEG_GAUSS21, work, &r, &e);
                total_result += r;
            }
        } else {
            // If the region overlaps (A, B), split it into two separate integrations
            if (a < A) {
                double r1=0, e1=0;
                gsl_integration_qag(&F, a, A, epsabs, epsrel, limit,
                                    GSL_INTEG_GAUSS21, work, &r1, &e1);
                total_result += r1;
            }
            if (b > B) {
                double r2=0, e2=0;
                gsl_integration_qag(&F, B, b, epsabs, epsrel, limit,
                                    GSL_INTEG_GAUSS21, work, &r2, &e2);
                total_result += r2;
            }
        }
    }

    gsl_integration_workspace_free(work);


    return total_result;
}




// double Pressure_Integrand_Diquark_QP(double Gd, double q, double T, double mu,
//                                       double m, double mu_star, double delta,
//                                       double m0, double mu_star0, double delta0,
//                                       double thermal_cutoff, double Nf, double Lambda_CutOff)
// {
//     const double epsabs = 1e-9;
//     const double epsrel = 1e-9;
//     const size_t limit  = 500;

//     double A = -q - 2*mu_star;  // A < B
//     double B = q - 2*mu_star;

//     Params p { Gd, T, mu, m, mu_star, delta,
//                m0, mu_star0, delta0,
//                Nf, Lambda_CutOff,
//                thermal_cutoff, q };

//     gsl_function F;
//     F.function = &phaseshift_diquark_integrand_w;
//     F.params   = &p;

//     gsl_integration_workspace *work = gsl_integration_workspace_alloc(limit);

//     double total_result = 0.0;

//     // --- seeds near known nonzero regions ---
//     std::vector<double> seeds;
//     if (std::abs(mu_star) < 1e-12) {
//         seeds = { std::sqrt(4*Lambda_CutOff*Lambda_CutOff + 4*m*m + q*q) };
//     } else {
//         seeds = {
//             -std::sqrt(4*Lambda_CutOff*Lambda_CutOff + 0*m*m + q*q) - 2*mu_star,
//             std::sqrt(4*Lambda_CutOff*Lambda_CutOff + 0*m*m + q*q) - 2*mu_star
//         };
//     }

//     // --- coarse scan parameters ---
//     const double dw_initial  = 5;
//     const double dwmin = 1e-9;  // minimum dw allowed
//     const double threshold = 1e-9;   // looser threshold for detecting small values
//     const int N_consec     = 4;

//     std::vector<std::pair<double,double>> active_regions;
//     std::vector<std::pair<double,double>> explored;

//     auto already_scanned = [&](double x) {
//         for (auto [a,b] : explored) {
//             if (x >= a && x <= b) return true;
//         }
//         return false;
//     };

//     // Function to detect zeros or small values and refine dw
//     auto refine_dw = [&](double a, double dw) {
//         // Check at midpoint between a and a + dw
//         double midpoint = a + 0.5 * dw;
//         double fmid = coarse_probe(midpoint, &p);

//         // If function goes to zero, halve dw and check again
//         if (fmid < threshold) {
//             dw = std::max(dw * 0.5, dwmin);
//            // std::cout << "Refined dw to: " << dw << " at position: " << midpoint << std::endl;
//         }

//         return dw;
//     };

//     // Expand from seeds using the refined dw strategy
//     auto expand_from_seed = [&](double seed) {
//         if (already_scanned(seed)) return;

//         double dw = dw_initial;

//         // Expand from the seed using an adaptive refinement strategy
//         int zeros = 0;
//         double a = seed;
//         double region_start = a;

//         // Expand in + direction with refinement
//         while (zeros < N_consec) {
//             double mid = a + 0.5*dw;
//             if (already_scanned(mid)) break;

//             double fmid = coarse_probe(mid, &p);

//             if (fmid < threshold) zeros++;
//             else { zeros = 0; }

//             // If we found a zero, halve the step size and continue searching
//             if (zeros == N_consec) {
//                 dw = refine_dw(mid, dw);  // dynamically refine dw
//             }

//             a += dw;
//         }

//         // Expand in - direction with refinement
//         zeros = 0;
//         a = seed;
//         region_start = a;

//         while (zeros < N_consec) {
//             double mid = a - 0.5*dw;
//             if (already_scanned(mid)) break;

//             double fmid = coarse_probe(mid, &p);

//             if (fmid < threshold) zeros++;
//             else { zeros = 0; }

//             // If we found a zero, halve the step size and continue searching
//             if (zeros == N_consec) {
//                 dw = refine_dw(mid, dw);  // dynamically refine dw
//             }

//             a -= dw;
//         }

//         // Mark the region as explored
//         active_regions.push_back({region_start, a});
//         explored.push_back({region_start, a});
//     };

//     // Iterate through each seed and expand to find regions
//     for (double seed : seeds)
//         expand_from_seed(seed);

//     // --- integrate active regions, skipping (A, B) ---
//     const double w_cutoff = 1e-8;

//     for (auto [a,b] : active_regions) {
//         // Skip regions that overlap with (A, B)
//         if (b <= A || a >= B) {
//             // If the region is entirely outside of (A, B), integrate normally
//             if (a < 0 && b > 0) {
//                 double r1=0, e1=0, r2=0, e2=0;
//                 gsl_integration_qag(&F, a, -w_cutoff, epsabs, epsrel, limit,
//                                     GSL_INTEG_GAUSS21, work, &r1, &e1);
//                 gsl_integration_qag(&F, w_cutoff, b, epsabs, epsrel, limit,
//                                     GSL_INTEG_GAUSS21, work, &r2, &e2);
//                 total_result += r1 + r2;
//             } else if (b <= -w_cutoff || a >= w_cutoff) {
//                 double r=0, e=0;
//                 gsl_integration_qag(&F, a, b, epsabs, epsrel, limit,
//                                     GSL_INTEG_GAUSS21, work, &r, &e);
//                 total_result += r;
//             }
//         } else {
//             // If the region overlaps (A, B), split it into two separate integrations
//             if (a < A) {
//                 double r1=0, e1=0;
//                 gsl_integration_qag(&F, a, A, epsabs, epsrel, limit,
//                                     GSL_INTEG_GAUSS21, work, &r1, &e1);
//                 total_result += r1;
//             }
//             if (b > B) {
//                 double r2=0, e2=0;
//                 gsl_integration_qag(&F, B, b, epsabs, epsrel, limit,
//                                     GSL_INTEG_GAUSS21, work, &r2, &e2);
//                 total_result += r2;
//             }
//         }
//     }

//     gsl_integration_workspace_free(work);

//     // symmetry doubling
//     if (std::abs(mu_star) < 1e-12) {
//         total_result *= 2.0;
//     }

//     return total_result;
// }


// Pressure That Does Not Consider the Landau Term
double pressure_diquark_QP(double Gd, double T, double mu, double m, double mu_star, double delta,
                        double m0, double mu_star0, double delta0,
                        double thermal_cutoff, double Nf, double Lambda_CutOff) {

    auto integrand_q = [&](double q) -> double {
        return 3.0 / (2.0 * pi * pi) * q * q * Pressure_Integrand_Diquark_QP(
            Gd, q, T, mu, m, mu_star, delta,
            m0, mu_star0, delta0,
            thermal_cutoff, Nf, Lambda_CutOff
        );
    };

    gsl_integration_workspace* work = gsl_integration_workspace_alloc(100);

    gsl_function F;
    F.function = [](double q, void* params) -> double {
        auto* f = static_cast<decltype(integrand_q)*>(params);
        return (*f)(q);
    };
    F.params = &integrand_q;

    double result = 0.0, error = 0.0;

    int status = gsl_integration_qag(&F, 0.0, 5e3, 1e-4, 1e-4, 5,
                                      GSL_INTEG_GAUSS15, work, &result, &error);

    if (status != GSL_SUCCESS || std::isnan(result)) {
        std::cerr << "[WARNING] GSL integration failed in pressure_diquark_Test:\n";
        std::cerr << "    status = " << gsl_strerror(status) << "\n";
        std::cerr << "    result = " << result << ", error = " << error << "\n";
        std::cerr << "    parameters: "
                  << " Gd=" << Gd << " T=" << T << " mu=" << mu
                  << " m=" << m << " mu*=" << mu_star << " delta=" << delta
                  << " thermal_cutoff=" << thermal_cutoff
                  << " Nf=" << Nf << " Lambda=" << Lambda_CutOff
                  << std::endl;
    }

    gsl_integration_workspace_free(work);
    return result;
}





//--------------------------    Landau Pressure



// double Pressure_Integrand_Diquark_LD(double Gd, double q, double T, double mu,
//                                                      double m, double mu_star, double delta,
//                                                      double m0, double mu_star0, double delta0,
//                                                      double thermal_cutoff,
//                                                      double Nf, double Lambda_CutOff)
// {
//     const double epsabs = 1e-6;
//     const double epsrel = 1e-6;
//     const size_t limit  = 200;

//     double A = -q - 2*mu_star;
//     double B = q - 2*mu_star;

//     Params p { Gd, T, mu, m, mu_star, delta,
//                m0, mu_star0, delta0,
//                Nf, Lambda_CutOff,
//                thermal_cutoff, q };

//     gsl_function F;
//     F.function = &phaseshift_diquark_integrand_w;
//     F.params   = &p;

//     gsl_integration_workspace *work = gsl_integration_workspace_alloc(limit);

//     double total_result = 0.0;

//     // --- seeds near known nonzero regions ---
//     std::vector<double> seeds;
//     if (std::abs(mu_star) < 1e-12) {
//         seeds = { std::sqrt(4*Lambda_CutOff*Lambda_CutOff + 4*m*m + q*q) };
//     } else {
//         seeds = {
//             - 2*mu_star
//         };
//     }

//     // --- coarse scan parameters ---
//     const double dw_default = 5.0;  // default step size
//     const double dw_near_zero = 0.2; // smaller step size near zero
//     const double threshold = 1e-6;   // looser, since we only detect activity
//     const int N_consec     = 5;

//     std::vector<std::pair<double,double>> active_regions;
//     std::vector<std::pair<double,double>> explored;

//     auto already_scanned = [&](double x) {
//         for (auto [a,b] : explored) {
//             if (x >= a && x <= b) return true;
//         }
//         return false;
//     };

//     auto expand_from_seed = [&](double seed) {
//         if (already_scanned(seed)) return;

//         int zeros = 0;
//         double a = seed;
//         double region_start = a;

//         // Function to dynamically adjust dw based on proximity to zero
//         auto get_dynamic_dw = [&](double x) {
//             if (std::abs(x) < 1e-3) {  // Adjust threshold as necessary
//                 return dw_near_zero;
//             } else {
//                 return dw_default;
//             }
//         };

//         // expand in + direction
//         while (zeros < N_consec) {
//             double dw = get_dynamic_dw(a);
//             double mid = a + 0.5 * dw;
//             if (already_scanned(mid)) break;

//             double fmid = coarse_probe(mid, &p);

//             if (fmid < threshold) zeros++;
//             else { zeros = 0; }

//             if (zeros == N_consec && region_start < a) {
//                 a += 2 * dw;  // include 2 bins past zeros
//                 active_regions.push_back({region_start, a});
//                 explored.push_back({region_start, a});
//                 break;
//             }
//             a += dw;
//         }

//         // expand in - direction
//         zeros = 0;
//         a = seed;
//         region_start = a;

//         while (zeros < N_consec) {
//             double dw = get_dynamic_dw(a);
//             double mid = a - 0.5 * dw;
//             if (already_scanned(mid)) break;

//             double fmid = coarse_probe(mid, &p);

//             if (fmid < threshold) zeros++;
//             else { zeros = 0; }

//             if (zeros == N_consec && region_start > a) {
//                 a -= 2 * dw;  // include 2 bins past zeros
//                 active_regions.push_back({a, region_start});
//                 explored.push_back({a, region_start});
//                 break;
//             }
//             a -= dw;
//         }
//     };

//     for (double seed : seeds)
//         expand_from_seed(seed);

//     // --- integrate active regions within [A, B] ---
//     const double w_cutoff = 1e-8;

//     for (auto [a, b] : active_regions) {
//         // Skip regions that are completely outside [A, B]
//         if (b <= A || a >= B) {
//             continue;
//         }

//         // Ensure that the region is within the bounds [A, B]
//         double region_start = std::max(a, A);
//         double region_end = std::min(b, B);

//         if (region_start < region_end) {
//             double r = 0, e = 0;
//             gsl_integration_qag(&F, region_start, region_end, epsabs, epsrel, limit,
//                                 GSL_INTEG_GAUSS21, work, &r, &e);
//             total_result += r;
//         }
//     }

//     gsl_integration_workspace_free(work);

//     // symmetry doubling
//     if (std::abs(mu_star) < 1e-12) {
//         total_result *= 2.0;
//     }

//     return total_result;
// }




double Pressure_Integrand_Diquark_LD(double Gd, double q, double T, double mu,
                                     double m, double mu_star, double delta,
                                     double m0, double mu_star0, double delta0,
                                     double thermal_cutoff,
                                     double Nf, double Lambda_CutOff)
{
    const double epsabs = 1e-4;  // Absolute tolerance for integration
    const double epsrel = 1e-4;  // Relative tolerance for integration
    const size_t limit  = 50;   // Maximum number of iterations for integration

    double A = -q - 2 * mu_star;  // Lower limit of integration
    double B = q - 2 * mu_star;   // Upper limit of integration

    Params p { Gd, T, mu, m, mu_star, delta,
               m0, mu_star0, delta0,
               Nf, Lambda_CutOff,
               thermal_cutoff, q };

    gsl_function F;
    F.function = &phaseshift_diquark_integrand_w;
    F.params   = &p;

    gsl_integration_workspace *work = gsl_integration_workspace_alloc(limit);

    double total_result = 0.0;

    // --- Define the points array (integration limits + singularity point) ---
    // If the singularity is at w = 0, we'll insert it in the array appropriately
    double pts[3] = {A, 0.0, B};  // pts[0] = A, pts[1] = x_1, pts[2] = 0.0 (singularity), pts[3] = B
    
    // --- Perform the integration ---
    // Since you are defining the singularity at pts[2] = 0, GSL will take care of it as a special point.
    double r = 0, e = 0;

    // GSL integration: `0.0` is the singularity point
    int status = gsl_integration_qagp(&F, pts, 3, epsabs, epsrel, limit, work, &r, &e);

    if (status != 0) {
        std::cerr << "Integration failed with error code: " << status << std::endl;
        gsl_integration_workspace_free(work);
        return 0.0;
    }

    total_result = r;

    gsl_integration_workspace_free(work);

    // Handle symmetry doubling if mu_star is close to zero (optional)

    return total_result;
}



double pressure_diquark_LD(double Gd, double T, double mu, double m, double mu_star, double delta,
                        double m0, double mu_star0, double delta0,
                        double thermal_cutoff, double landau_cutoff, double Nf, double Lambda_CutOff) {

    auto integrand_q = [&](double q) -> double {
        return 3.0 / (2.0 * pi * pi) * q * q * Pressure_Integrand_Diquark_LD(
            Gd, q, T, mu, m, mu_star, delta,
            m0, mu_star0, delta0,
            thermal_cutoff, Nf, Lambda_CutOff
        );
    };

    gsl_integration_workspace* work = gsl_integration_workspace_alloc(100);

    gsl_function F;
    F.function = [](double q, void* params) -> double {
        auto* f = static_cast<decltype(integrand_q)*>(params);
        return (*f)(q);
    };
    F.params = &integrand_q;

    double result = 0.0, error = 0.0;

    int status = gsl_integration_qag(&F, 0.0, landau_cutoff, 1e-6, 1e-6, 50,
                                      GSL_INTEG_GAUSS15, work, &result, &error);

    if (status != GSL_SUCCESS || std::isnan(result)) {
        std::cerr << "[WARNING] GSL integration failed in pressure_diquark_Test:\n";
        std::cerr << "    status = " << gsl_strerror(status) << "\n";
        std::cerr << "    result = " << result << ", error = " << error << "\n";
        std::cerr << "    parameters: "
                  << " Gd=" << Gd << " T=" << T << " mu=" << mu
                  << " m=" << m << " mu*=" << mu_star << " delta=" << delta
                  << " thermal_cutoff=" << thermal_cutoff
                  << " Nf=" << Nf << " Lambda=" << Lambda_CutOff
                  << std::endl;
    }

    gsl_integration_workspace_free(work);
    return result;
}



















double integrate_diquark_w_for_q(double q, const Params& baseParams, gsl_integration_workspace* work_w) {
    Params local = baseParams;
    local.q = q;

    gsl_function F;
    F.function = &phaseshift_diquark_integrand_w;
    F.params = &local;

    double result = 0.0, error = 0.0;
    const size_t limit = 500;

    // lower limit for w integration
    const double w_min_raw = (std::sqrt(4.0 * local.m * local.m + q * q) - 2.0 * local.mu_star)/(1 + (local.mu + 100)/(500+q));
    const double w_min = (w_min_raw >= 1e-9) ? w_min_raw : 1e-9;

    int status = GSL_SUCCESS;

    if (local.thermal_cutoff == 0) {
        // integrate from w_min to ∞
        // status = gsl_integration_qagiu(&F, w_min, 1e-8, 1e-8, limit, work_w, &result, &error);
        status  =  gsl_integration_qagi(&F, 0, 1e-8, limit, work_w, &result, &error);
    } else {
        // finite cutoff
        const double w_max = std::sqrt(4.0 * local.Lambda_CutOff * local.Lambda_CutOff +
                                       4.0 * local.m * local.m + q * q) - 
                             2.0 * local.mu_star;

        if (w_max <= w_min) {
            return 0.0;  // no valid interval
        }

        status = gsl_integration_qag(&F, w_min, w_max, 1e-8, 1e-8, 3, limit, work_w, &result, &error);
        //status  =  gsl_integration_qagi(&F, 0, 1e-8, limit, work_w, &result, &error);
    }


    return result;
}







// --- top-level pressure computation with parallel q-integration ---
// Uses 8 OpenMP threads (explicit). Each thread gets its own GSL workspace.
double pressure_diquark(double Gd, double T, double mu, double m, double mu_star, double delta,
                        double m0, double mu_star0, double delta0,
                        double thermal_cutoff, double Nf, double Lambda_CutOff) {

    // === integration grid parameters ===
    const int N_q = 1000;                // tune for accuracy; increase for better resolution
    const double q_min = 1e-8;
    const double q_max = 5e3;
    const double dq = (q_max - q_min) / static_cast<double>(N_q);

    // Prepare base params (string passed by pointer)
    Params base;
    base.Gd = Gd; base.T = T; base.mu = mu; base.m = m; base.mu_star = mu_star; base.delta = delta;
    base.m0 = m0; base.mu_star0 = mu_star0; base.delta0 = delta0;
    base.Nf = Nf; base.Lambda_CutOff = Lambda_CutOff;
    base.thermal_cutoff = thermal_cutoff;
    base.q = 0.0;


    double total = 0.0;

    // Force 8 threads for this parallel region (user-specified)
    const int nt = 5;

    #pragma omp parallel num_threads(nt)
    {
        // Each thread: its own workspace for inner integration
        gsl_integration_workspace* work_w = gsl_integration_workspace_alloc(2000);

        double local_sum = 0.0;

        #pragma omp for schedule(static)
            for (int i = 0; i < N_q; ++i) {
                double q_i = q_min + (i + 0.5) * dq;

                double inner = integrate_diquark_w_for_q(q_i, base, work_w);
                if (std::isfinite(inner)) {
                    local_sum += inner * dq;
            }
        }

        // Accumulate into total (atomic or critical). Using critical to avoid potential
        // subtle performance/precision issues with floating-point atomics on some compilers.
        #pragma omp critical
        { total += local_sum; }

        gsl_integration_workspace_free(work_w);
    } // omp parallel

    return total;
}






