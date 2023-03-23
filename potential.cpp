/*
- EAM potential for palladium hydride
- data obtained from the appendix of X.W. Zhou et al. 2007
*/

#include <iostream>
#include <math.h>
#include <input.h>
#include <potential.h>
using namespace std;

Input input;
double eps = input.file.eps;
const double pi = 3.141592653589793;

double EAM::F_Pd(double rho) // embedding energy for Pd 
{
    double a1  = 295878.9003038662;
    double a2  = 0.20581955357385892;
    double a3  = 0.081228755904399;
    double a4  = 0.05298811034615951;
    double b11 = 2.4242616904962846;
    double b12 = 1.4791899886249564;
    double b21 = 2.1376274623740064;
    double b22 = 1.2169215689822592;
    double b31 = 1.6486007989726832;
    double b32 = 0.8159825255339774;
    double b41 = 1.0749204110338482;
    double b42 = 0.42007491336688396;
    double b51 = 0.5128056047933808;
    double b52 = 0.12468685331167456;
    rho /= 50.0;
    double y = a1 * (rho-a2) * (rho-a3) * rho * (rho+a4) * (rho*(rho-b11)+b12) * (rho*(rho-b21)+b22) * (rho*(rho-b31)+b32) * (rho*(rho-b41)+b42) * (rho*(rho-b51)+b52);

    double Fprime = 0.251546;
    y -= 50.0*Fprime*rho; 
    return y;
}

double EAM::F_Pd_deriv(double rho)
{
    return (F_Pd(rho+eps) - F_Pd(rho-eps)) / (2.0*eps);
}

double EAM::F_H(double rho) // embedding energy for H
{
    double c = 0.000197047;
    double d = 1.18860;
    double a = 9.99780;
    double b = 60.0155;
    double eps1 = 0.0540638;
    double y = -c*(1.0/(2.0+d)*pow(rho+eps1,2.0+d) - (a+b)/(1.0+d)*pow(rho+eps1,1.0+d) + a*b/d*pow(rho+eps1,d));

    double Fprime = -0.0296604;
    y -= Fprime*rho;
    return y;
}

double EAM::F_H_deriv(double rho)
{
    return (F_H(rho+eps) - F_H(rho-eps)) / (2.0*eps);
}

double EAM::f_Pd(double r) // electron density for Pd 
{
    double funRc = 0.001510569378297; // at r = 5.0;
    double funDerivRc = -0.077976566126292; // at r = 5.0;

    double a11 = 2.7267629107325706;
    double a12 = 1.8716766113599643;
    double a21 = 2.50290548851635;
    double a22 = 1.668549182690922;
    double a31 = 2.0924467509943674;
    double a32 = 1.3150372774478005;
    double a41 = 1.564328475106985;
    double a42 = 0.8987511149780485;
    double a51 = 1.009780903403673;
    double a52 = 0.5124363774128722;
    double a61 = 0.5304054524800665;
    double a62 = 0.2169886022464641;
    double a71 = 0.1356566408715063;
    double a72 = 0.035852347523891395;
    r /= 5.0;
    double rho_a = r * (r*(r-a11)+a12) * (r*(r-a21)+a22) * (r*(r-a31)+a32) * (r*(r-a41)+a42) * (r*(r-a51)+a52) * (r*(r-a61)+a62) * (r*(r-a71)+a72);
    
    double b1 = -0.02972698211669922;
    double b2 = 0.6676807403564453;
    double b3 = -255.8965835571289;
    double b4 = 14673.409149169922;
    double b5 = -2.597301181336601e7;
    double y = b1 + r*(b2 + r*(b3 + r*(b4 + b5*rho_a)));
    
    if (r*5.0<=0.0)
        y = 0.0;
    else if (r*5.0>5.0)
        y = funRc * exp(funDerivRc/funRc*(r*5.0-5.0));
    return y;
}

double EAM::f_Pd_deriv(double r)
{
    return (f_Pd(r+eps) - f_Pd(r-eps)) / (2.0*eps);
}

double EAM::f_H(double r) // electron density for H 
{
    double funRc = 0.009986990526578; // at r = 5.35;
    double funDerivRc = -0.013075667086408; // at r = 5.35;
    
    double C = 11.0025;
    double delta = 1.30927;
    double y = C*exp(-delta*r);
   
    if (r<=0.0)
        y = 0.0;
    else if (r>5.35)
        y = funRc * exp(funDerivRc/funRc*(r-5.35));
    return y;
}

double EAM::f_H_deriv(double r)
{
    return (f_H(r+eps) - f_H(r-eps)) / (2.0*eps);
}

double EAM::phi_Pd(double r) // pair energy for Pd
{
    double funRc = -0.016958780164980; // at r = 5.05;
    double funDerivRc = 0.089019679316328; // at r = 5.05;

    double a1 = -79415.24035137112;
    double a2 = 1.0699996145674568;
    double a3 = 1.06015072612581;
    double a4 = 0.42433991011376526;
    double a5 = 0.06169160085238687;
    double b11 = 2.0586473420376348;
    double b12 = 1.0683922574015199;
    double b21 = 1.6696359816422877;
    double b22 = 0.7337878627470482;
    double b31 = 1.1690370066230809;
    double b32 = 0.3909805777737639;
    double b41 = 0.2635598721249787;
    double b42 = 0.033551116514910245;
    r /= 5.0;
    double y = a1 * (r-a2) * (r-a3) * (r-a4) * (r+a5) * (r*(r-b11)+b12) * (r*(r-b21)+b22) * (r*(r-b31)+b32) * (r*(r-b41)+b42);

    r *= 5.0;
    double Fprime = 0.251546;
    y += 2*Fprime*f_Pd(r);
    if (r<=0.0)
        y = 0.0;
    else if (r>5.05)
        y = funRc * exp(funDerivRc/funRc*(r-5.05));
    return y;
}

double EAM::phi_Pd_deriv(double r)
{
    return (phi_Pd(r+eps) - phi_Pd(r-eps)) / (2.0*eps);
}

double EAM::phi_PdH(double r) // pair energy for PdH
{
    double funRc = -3.353001113951400e-04; // at r = 5.35;
    double funDerivRc = 7.147062196813336e-04; // at r = 5.35;

    double D = 0.2494540;
    double alpha = 4.82613;
    double beta = 2.13158;
    double r0 = 1.50964;
    double y = D*(beta*exp(-alpha*(r-r0)) - alpha*exp(-beta*(r-r0)));
    if (r<=0.0)
        y = 0.0;
    else if (r>5.35)
        y = funRc * exp(funDerivRc/funRc*(r-5.35));
    return y;
}

double EAM::phi_PdH_deriv(double r)
{
    return (phi_PdH(r+eps) - phi_PdH(r-eps)) / (2.0*eps);
}


double EAM::phi_H(double r) // pair energy for H
{
    double funRc = -0.003702722457147; // at r = 5.35;
    double funDerivRc = 0.005465945471933; // at r = 5.35

    double D = 0.0661496;
    double alpha = 3.67263;
    double beta = 1.47797;
    double r0 = 2.51980;
    double y = D*(beta*exp(-alpha*(r-r0)) - alpha*exp(-beta*(r-r0)));
    if (r<=0.0)
        y = 0.0;
    else if (r>5.35)
        y = funRc * exp(funDerivRc/funRc*(r-5.35));
    return y;
}

double EAM::phi_H_deriv(double r)
{
    return (phi_H(r+eps) - phi_H(r-eps)) / (2.0*eps);
}

