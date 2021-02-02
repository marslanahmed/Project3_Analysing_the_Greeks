// Author: Muhammad Arslan Ahmed
//   Date: 2020-12-11
#ifndef EuropeanOption_hpp
#define EuropeanOption_hpp

#include <stdio.h>
#include <iostream>
using namespace std;

class EuropeanOption {
public:
    // Constants
    const char OPTION_CALL = 'C';
    const char OPTION_PUT  = 'P';
    const int INVALID_OPTION_TYPE = -1;
    
    // Constructor
    EuropeanOption(double, double, double, double, double, double, char);
    // Destructor
    ~EuropeanOption();
    // Accessor Methods
    double Get_Spot_Price() { return S; }
    double Get_Strike_Price() { return K; }
    double Get_Interest_Rate() { return r; }
    double Get_Volatility() { return sigma; }
    double Get_Dividend_Yield() { return D; }
    double Get_Time_to_Maturity() { return T; }
    char  Get_Option_Type() { return optionType; }
    
    void Set_Spot_Price(double a) { S = a; }
    void Set_Strike_Price(double a) { K = a; }
    void Set_Interest_Rate(double a) { r = a; }
    void Set_Volatility(double a) { sigma = a; }
    void Set_Dividend_Yield(double a) { D = a; }
    void Set_Time_to_Maturity(double a) { T = a; }
    void Set_Option_Type(char a) { optionType = a; }
    
    // Member functions
    void Get_Information();
    void Get_Greeks();
    void Get_Greeks_FDM();
    double Get_BSM_Value();
    
    // Member functions (The Greeks)
    double Get_Delta();
    double Get_Gamma();
    double Get_Theta();
    double Get_Speed();
    double Get_Vega();
    double Get_Rho_r();
    double Get_Rho_D();
    
    // Member functions (Finite Difference Methods)
    
    // Member functions (Monte Carlo)
    double Get_MC_Value( const int& N );
    
private:
    double S;// spot price
    double K;// strike price
    double r;// interest rate
    double sigma;// volatility
    double D;// dividend yield
    double T;// time to maturity
    char optionType;// (C)all or (P)ut?
    
    double d1;
    double d2;
    
    double N1(const double);
    double N2(const double);
    
    double gaussian_box_muller();
};


#endif /* EuropeanOption_hpp */
