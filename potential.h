#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_
using namespace std;

class EAM
{
public:
    double F_Pd(double rho);
    double F_Pd_deriv(double rho);
    double F_H(double rho);
    double F_H_deriv(double rho);
    double f_Pd(double r);
    double f_Pd_deriv(double r);
    double f_H(double r);
    double f_H_deriv(double r);
    double phi_Pd(double r);
    double phi_Pd_deriv(double r);
    double phi_PdH(double r);
    double phi_PdH_deriv(double r);
    double phi_H(double r);
    double phi_H_deriv(double r);
};
#endif
