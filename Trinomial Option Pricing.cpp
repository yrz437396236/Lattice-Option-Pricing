#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>

using namespace std;


class backward_Trinomial_option_pricing {
public:

    string option_type;
    double dividend_rate, risk_free_rate, strike_price;
    double initial_stock_price, expiration_time, volatility;
    int no_of_divisions;

    double option_pricing() {

        clock_t start = clock();
        double* V = new double[2 * no_of_divisions + 1];
        double* S = new double[2 * no_of_divisions + 1];

        //vector <double> initial(25000, 0.0);
        double dt = expiration_time / no_of_divisions;
        double result;
        double u_ = exp(volatility * pow(dt, 0.5) * pow(1.5, 0.5));
        double d_ = 1 / u_;
        double v_ = risk_free_rate - dividend_rate - 0.5 * pow(volatility, 2);
        double pu_ = 1.0 / 3.0 + (v_ * pow(dt, 0.5)) / (2 * pow(1.5, 0.5) * volatility);
        double pm_ = 1.0 / 3.0;
        double pd_ = 1.0 - pm_ - pu_;
        
        //parametrization first

        //for (int i = 0; i <= 2 * no_of_divisions; i++) {
            //S[i] = initial_stock_price * pow(u_, no_of_divisions - i);
        //}

        S[0] = initial_stock_price * pow(u_, no_of_divisions);
        for (int i = 1; i <= 2 * no_of_divisions; i++) {
           S[i] = S[i-1]* d_;
        }

        for (int i = 0; i <= 2 * no_of_divisions; i++) {
            V[i] = max(S[i] - strike_price, 0.0);
        }

        //backward
        for (int j = no_of_divisions - 1; j >= 0; j--) {
            for (int i = 0; i <= 2 * j; i++) {
                V[i] = (pu_ * V[i] + pm_ * V[i + 1] + pd_ * V[i + 2]) * exp(-risk_free_rate * dt);
                V[i] = max(V[i], S[no_of_divisions - j + i] - strike_price);
            }
        }
        clock_t end = clock();
        cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
        cout << "--------------------------------------" << endl;

        double option_price = V[0];
        delete[]V;
        delete[]S;

        cout << option_price << endl;
        return option_price;
        
    }
};

