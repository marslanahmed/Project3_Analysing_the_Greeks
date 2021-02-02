// Author: Muhammad Arslan Ahmed
//   Date: 2020-12-11
#include <iostream>
#include <iomanip>
#include <cmath>
#include "EuropeanOption.hpp"
using namespace std;

// CONSTRUCTOR
EuropeanOption::EuropeanOption(double S, double K, double r, double sigma, double D, double T, char typ)
{
    if ( ( typ != OPTION_CALL ) &&
         ( typ != OPTION_PUT ) )
    {
        throw INVALID_OPTION_TYPE;
    }
    Set_Spot_Price(S);
    Set_Strike_Price(K);
    Set_Interest_Rate(r);
    Set_Volatility(sigma);
    Set_Dividend_Yield(D);
    Set_Time_to_Maturity(T);
    Set_Option_Type(typ);
    
    // Calculate d1 and d2
    d1 = ( log(Get_Spot_Price()/Get_Strike_Price()) + (Get_Interest_Rate()+0.5*pow(Get_Volatility(),2.0))*Get_Time_to_Maturity() )/( Get_Volatility()*sqrt(Get_Time_to_Maturity()) );
    
    d2 = d1 - Get_Volatility()*sqrt( Get_Time_to_Maturity() );
}


// DESTRUCTOR
EuropeanOption::~EuropeanOption()
{
    // nothing to output
}

// STATISTICAL UTILITY
double EuropeanOption::N1(const double x)
{ // Cumulative Density Function (integral of PDF)
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
    
    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - N1(-x);
    }
}

double EuropeanOption::N2(const double x)
{ // Probability Density Function
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

// MEMBER FUCNTIONS
void EuropeanOption::Get_Information()
{
    if (optionType==OPTION_CALL) {
        cout<<"\n================== CALL " << endl;
        cout << fixed << setprecision(4) << setfill(' ');
        
        
    } else if(optionType==OPTION_PUT) {
        
    }
    
}

void EuropeanOption::Get_Greeks()
{
    int width=8;
    cout << fixed << setprecision(4) << setfill(' ');
    cout <<"---------- THE GREEKS" << endl;
    cout<< "BSM Value: " << setw(width) << Get_BSM_Value() << endl;
    cout<< "    Delta: " << setw(width) << Get_Delta() << endl;
    cout<< "    Gamma: " << setw(width) << Get_Gamma() << endl;
    cout<< "    Theta: " << setw(width) << Get_Theta() << endl;
    cout<< "    Speed: " << setw(width) << Get_Speed() << endl;
    cout<< "     Vega: " << setw(width) << Get_Vega()  << endl;
    cout<< "    Rho_r: " << setw(width) << Get_Rho_r() << endl;
    cout<< "    Rho_D: " << setw(width) << Get_Rho_D() << endl;
}


void EuropeanOption::Get_Greeks_FDM()
{
    int width=8;
    cout << fixed << setprecision(4) << setfill(' ');
    cout << "================== THE GREEKS" << endl;
    cout<< "BSM Value: " << setw(width) << Get_BSM_Value() << endl;
    cout<< "    Delta: " << setw(width) << Get_Delta() << endl;
    cout<< "    Gamma: " << setw(width) << Get_Gamma() << endl;
    cout<< "    Theta: " << setw(width) << Get_Theta() << endl;
    cout<< "    Speed: " << setw(width) << Get_Speed() << endl;
    cout<< "     Vega: " << setw(width) << Get_Vega()  << endl;
    cout<< "    Rho_r: " << setw(width) << Get_Rho_r() << endl;
    cout<< "    Rho_D: " << setw(width) << Get_Rho_D() << endl;
}

//================
// BSM Value
//================
double EuropeanOption::Get_BSM_Value()
{// Black Scholes Merton equation is used to calculate the analytical value
    double value=0.0;
    if (optionType=='C') {
        value = S*N1(d1)*exp(-D*T) - K*N1(d2)*exp(-r*T);
    } else if (optionType=='P') {
        value = -S*N1(-d1)*exp(-D*T) + K*N1(-d2)*exp(-r*T);
    }
    return value;
}

//================
// DELTA
//================
double EuropeanOption::Get_Delta()
{   /*
     This is the sensitivity to the underlying asset.
     
     dV/dS
     
     */
    double delta=0.0;
    if (optionType=='C') {
        delta = exp(-D*T)*N1(d1);
    } else if (optionType=='P') {
        delta = exp(-D*T)*(N1(d1) - 1);
    }
    return delta;
}

//================
// GAMMA
//================
double EuropeanOption::Get_Gamma()
{   /*
     This is the sensitivity of Delta to the underlying asset.
     The rate of change of delta with respect to the underlying
     asset's price.
     
     d^2V/dS^2
     
     */
    return (exp(-D*T)*N2(d1))/(sigma*S*sqrt(T));
}

//================
// THETA
//================
double EuropeanOption::Get_Theta()
{   /*
     This is the sensitivity of the Option value to time.
     
     dV/dT
     
     */
    double theta=0.0;
    if (optionType=='C') {
        theta = -(sigma*S*exp(-D*T)*N2(d1))/(2.0*sqrt(T));
        theta+= (D*S*N1(d1));
        theta-= (r*K*exp(-D*T)*N1(d2));
    } else if (optionType=='P') {
        theta = -(sigma*S*exp(-D*T)*N2(-d1))/(2.0*sqrt(T));
        theta+= (D*S*N1(-d1));
        theta-= (r*K*exp(-D*T)*N1(-d2));
    }
    return theta;
}

//================
// SPEED
//================
double EuropeanOption::Get_Speed()
{   /*
     This is the sensitivity of Gamma to the underlying asset.
     
     d^3V/dS^3
     
     */
    return -((exp(-D*T)*N2(d1))/(pow(sigma,2.0)*pow(S,2.0)*T))*(d1 + sigma*sqrt(T));
}

//================
// VEGA
//================
double EuropeanOption::Get_Vega()
{   /*
     This is the senstivity of the value of the option
     with respect to volatility.
     
     dV/do
     
     */
    return S*sqrt(T)*exp(-D*T)*N2(d1);
}

//================
// RHO_r
//================
double EuropeanOption::Get_Rho_r()
{   /*
     This is the sensitivity to the interest rate.
     
     dV/dr
     
     */
    double rho_r=0.0;
    if (optionType=='C') {
        rho_r = K*T*exp(-r*T)*N1(d2);
    } else if (optionType=='P') {
        rho_r = -K*T*exp(-r*T)*N1(-d2);
    }
    return rho_r;
}

//================
// RHO_D
//================
double EuropeanOption::Get_Rho_D()
{   /*
     This is the sensitivity to the dividend yield.
     
     dV/dD
     
     */
    double rho_D=0.0;
    if (optionType=='C') {
        rho_D = -T*S*exp(-D*T)*N1(d1);
    } else if (optionType=='P') {
        rho_D = T*S*exp(-D*T)*N1(-d1);
    }
    return rho_D;
}


//================
// MONTE CARLO
//================
double EuropeanOption::gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance"
  // is less than unity
  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);
    
    return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}


double EuropeanOption::Get_MC_Value( const int& N )
{
    /*
     Calculate the value of a European Option using
     N(input) simulations.
     */
    
    // The adjustment to the spot price
    double S_adjust = S * exp(T*(r-0.5*sigma*sigma));
    
    // Our current asset price ("spot")
    double S_cur = 0.0;
    
    // Holds the sum of all of the final option pay-offs
    double payoff_sum = 0.0;
    
    for (int i=0; i<N; i++)
    {
        // Generate a Gaussian random number via Box-Muller
        double gauss_bm = gaussian_box_muller();
        
        // Adjust the spot price via the Brownian motion final distribution
        S_cur = S_adjust * exp(sqrt( pow(sigma,2)*T )*gauss_bm);
        
        // Take the option pay-off, then add it to the rest of the pay-off values
        if (optionType == 'C') {
            payoff_sum += max(S_cur - K, 0.0); // max(S - K, 0)
        } else if (optionType == 'P') {
            payoff_sum += max(K - S_cur, 0.0); // max(K - S, 0)
        }
        
    }
    
    // Average the pay-off sum via the number of paths and then
    // discount the risk-free rate from the price
    return (payoff_sum / static_cast<double>(N)) * exp(-r*T);
}
