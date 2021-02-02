/*
 This code is designed to return option values
 and the greeks. It's very simplistic and was
 written by me for educational purposes.
 
 In here, I implement the Monte Carlo and
 Finite Difference Methods to calculate option
 values and the greeks too.
 
 Author: Muhammad Arslan Ahmed, PhD
 Date: 2020-12-11
 */

#include <iostream>
#include <time.h>
#include <fstream> // creates, reads, and writes to files
#include <sstream>
#include <vector>
#include <string>

#include "EuropeanOption.hpp"

using namespace std;


// FUNCTION PROTOTYPES
double percentage_error(double& x, double& y);


int main(int argc, const char * argv[]) {
    
    double S,K,r,sig,D,T;
    
    
    cout<<"============================="<<endl;
    cout<<" PRICING OPTIONS " << endl;
    cout<<"============================="<<endl;
    cout<<"Please input the "<<endl;
    cout<<"       Spot price: ";cin>>S;
    cout<<"     Strike price: ";cin>>K;
    cout<<"       Volatility: ";cin>>sig;
    cout<<"    Interest rate: ";cin>>r;
    cout<<"    Dividend rate: ";cin>>D;
    cout<<"Time to matrurity: ";cin>>T; // in years
    cout<<endl;
    
    
    
    // Number of Monte Carlo simulations
    vector<int> N = {1000, 10000, 100000, 1000000};
    vector<double> MC_call_vals;
    vector<double>  MC_put_vals;
    
    
    
    
    // =============================
    // EUROPEAN OPTION
    // =============================
    EuropeanOption call(S,K,r,sig,D,T,'C');
    EuropeanOption put( S,K,r,sig,D,T,'P');
    cout<<"=========================================================="<<endl;
    cout<<"\t\t\t\t\tEUROPEAN OPTION"<<endl;
    cout<<" CALL"<<endl;call.Get_Greeks();cout<<endl;
    cout<<" PUT"<<endl;put.Get_Greeks();
    
    double analytic_CALL_Val = call.Get_BSM_Value();
    double  analytic_PUT_Val =  put.Get_BSM_Value();
    
    // Analytical option value
    cout<<endl;
    cout<< "> ANALYTICAL VALUES "<<endl;
    cout<< "--------------------"<<endl;
    cout << fixed << setprecision(4) << setfill(' ');
    cout<< "Value of CALL option: " << analytic_CALL_Val << endl;
    cout<< "Value of  PUT option: " <<  analytic_PUT_Val << endl;
    
    // Monte Carlo method for option value
    cout<<endl;
    cout<<"> MONTE CARLO SIMS "<<endl;
    cout<<"--------------------"<<endl;
    cout<<"Sims \t\tCALL VALUE (%error)\t\t\tPUT VALUE (%error)"<<endl;
    for (int j=0;
         j<N.size();
         j++)
    {
        MC_call_vals.push_back(call.Get_MC_Value(N[j]) );
        MC_put_vals.push_back(  put.Get_MC_Value(N[j]) );
        
        // print
        cout << fixed << setprecision(4) << setfill(' ');
        cout<< N[j] << "\t\t\t" << MC_call_vals[j] << "\t("<< percentage_error(MC_call_vals[j], analytic_CALL_Val) << ")\t\t";
        cout<< MC_put_vals[j] << "\t("<< percentage_error(MC_put_vals[j], analytic_PUT_Val) << ")" << endl;
    }
    //cout<<"---------------------------------------------------------"<<endl;
    
    
    // Finite Difference method for option value
    cout<<endl;
    cout<<"> FINITE DIFF. METHOD"<<endl;
    cout<<"--------------------"<<endl;
    
    
    /*
     WORKING ON IT 
     */
    
    cout<<endl;
    return 0;
}


// HELPER FUNCTIONS
double percentage_error(double& x, double& y) {
    return abs(x-y)/x * 100.0;
}
