#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

inline double max(double a, double b) {
	return (b < a) ? a : b;
}

class recursive_binomial_option_pricing_brute_force {
public:

	string option_type;
	double dividend_rate, risk_free_rate, strike_price;
	double initial_stock_price, expiration_time, volatility;
	double u, p, d, R;
	int no_of_divisions;
	
	double parameter() {
		cout << "Brute Force Binomial Option Pricing(CRR)" << endl;
		cout << "Expiration Time (Years) = " << expiration_time << endl;
		cout << "Number of Divisions = " << no_of_divisions << endl;
		cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
		cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
		cout << "Initial Stock Price = " << initial_stock_price << endl;
		cout << "Strike Price = " << strike_price << endl;
		cout << "Up Factor (CRR) = " << u << endl;
		cout << "Down Factor (CRR) = " << d << endl;
		cout << "Uptick Probability (CRR) = " << p << endl;
		cout << "Dividend Rate = " << dividend_rate << endl;
		cout << "--------------------------------------" << endl;
		return 0;
	}

	double option_pricing() {

		clock_t start = clock();
		initialization();

		if (option_type == "AC") {

			//American call option
			double A_call_price_CRR = american_call_option_CRR(0, 0, initial_stock_price);

			cout << "CRR Binomial Price of an American Call Option = " << A_call_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return A_call_price_CRR;
		}
		else if (option_type == "AP") {

			//American put option
			double A_put_price_CRR = american_put_option_CRR(0, 0, initial_stock_price);

			cout << "CRR Binomial Price of an American Put Option = " << A_put_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return A_put_price_CRR;
		}
		else if (option_type == "EC") {

			// European call option
			double E_call_price_CRR = european_call_option_CRR(0, 0, initial_stock_price);

			cout << "CRR Binomial Price of an European Call Option = " << E_call_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return E_call_price_CRR;

		}
		else if (option_type == "EP") {

			//  European put option
			double E_put_price_CRR = european_put_option_CRR(0, 0, initial_stock_price);

			cout << "CRR Binomial Price of an European Put Option = " << E_put_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return E_put_price_CRR;
		}
		else {
			cout << "Unknown/Unsupport option type" << endl;
		}
	}

private:
	double initialization() {

		R = exp((risk_free_rate) * expiration_time / ((double)no_of_divisions));
		u = exp(volatility * sqrt(expiration_time / ((double)no_of_divisions)));
		d = exp((-volatility) * sqrt(expiration_time / ((double)no_of_divisions)));
		p = ((exp((risk_free_rate - dividend_rate) * expiration_time / ((double)no_of_divisions))) - d) / (u - d);
		return 0;
	}

	double american_put_option_CRR(int k, int i, double current_stock_price) {
		if (k == no_of_divisions)
			return max(0.0, (strike_price - current_stock_price));
		else{
			return max((strike_price - current_stock_price),
				(p * american_put_option_CRR(k + 1, i + 1, current_stock_price * u) +
					(1 - p) * american_put_option_CRR(k + 1, i - 1, current_stock_price * d)) / R);
		}
	}

	double american_call_option_CRR(int k, int i, double current_stock_price) {

		if (k == no_of_divisions)
			return max(0.0, (current_stock_price - strike_price));
		else{
			return max((current_stock_price - strike_price),
				(p * american_call_option_CRR(k + 1, i + 1, current_stock_price * u) +
					(1 - p) * american_call_option_CRR(k + 1, i - 1, current_stock_price * d)) / R);
		}
	}

	double european_put_option_CRR(int k, int i, double current_stock_price) {
		if (k == no_of_divisions)
			return max(0.0, (strike_price - current_stock_price));
		else{
			return (p * european_put_option_CRR(k + 1, i + 1, current_stock_price * u) +
				(1 - p) * european_put_option_CRR(k + 1, i - 1, current_stock_price * d)) / R;
		}
	}

	double european_call_option_CRR(int k, int i, double current_stock_price) {
		if (k == no_of_divisions)
			return max(0.0, (current_stock_price - strike_price));
		else{
			return (p * european_call_option_CRR(k + 1, i + 1, current_stock_price * u) +
				(1 - p) * european_call_option_CRR(k + 1, i - 1, current_stock_price * d)) / R;
		}
	}
};

class recursive_binomial_option_pricing_memorized {
public:

	string option_type;
	double dividend_rate, risk_free_rate, strike_price;
	double initial_stock_price, expiration_time, volatility;
	double u, p, d, R;
	int no_of_divisions;

	double parameter() {
		cout << "Memorized Binomial Option Pricing(CRR)" << endl;
		cout << "Expiration Time (Years) = " << expiration_time << endl;
		cout << "Number of Divisions = " << no_of_divisions << endl;
		cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
		cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
		cout << "Initial Stock Price = " << initial_stock_price << endl;
		cout << "Strike Price = " << strike_price << endl;
		cout << "Up Factor (CRR) = " << u << endl;
		cout << "Down Factor (CRR) = " << d << endl;
		cout << "Uptick Probability (CRR) = " << p << endl;
		cout << "Dividend Rate = " << dividend_rate << endl;
		cout << "--------------------------------------" << endl;
		return 0;
	}

	double option_pricing() {

		clock_t start = clock();
		initialization();

		//memorization
		double** matrix = new double* [no_of_divisions];
		for (int i = 0; i <= no_of_divisions; i++)
			matrix[i] = new double[2 * no_of_divisions];
		for (int i = 0; i <= no_of_divisions; i++)
			for (int j = 0; j <= 2 * no_of_divisions; j++)
			{
				matrix[i][j] = -1;
			}

		if (option_type == "AC") {

			//American call option
			double A_call_price_CRR = american_call_option_CRR(0, 0, initial_stock_price, matrix);

			cout << "CRR Binomial Price of an American Call Option = " << A_call_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return A_call_price_CRR;
		}
		else if (option_type == "AP") {

			//American put option
			double A_put_price_CRR = american_put_option_CRR(0, 0, initial_stock_price, matrix);

			cout << "CRR Binomial Price of an American Put Option = " << A_put_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return A_put_price_CRR;
		}
		else if (option_type == "EC") {

			// European call option
			double E_call_price_CRR = european_call_option_CRR(0, 0, initial_stock_price, matrix);

			cout << "CRR Binomial Price of an European Call Option = " << E_call_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return E_call_price_CRR;

		}
		else if (option_type == "EP") {

			//  European put option
			double E_put_price_CRR = european_put_option_CRR(0, 0, initial_stock_price, matrix);

			cout << "CRR Binomial Price of an European Put Option = " << E_put_price_CRR << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
			return E_put_price_CRR;
		}
		else {
			cout << "Unknown/Unsupport option type" << endl;
		}
	}

private:
	double initialization() {

		R = exp((risk_free_rate) * expiration_time / ((double)no_of_divisions));
		u = exp(volatility * sqrt(expiration_time / ((double)no_of_divisions)));
		d = exp((-volatility) * sqrt(expiration_time / ((double)no_of_divisions)));
		p = ((exp((risk_free_rate - dividend_rate) * expiration_time / ((double)no_of_divisions))) - d) / (u - d);
		return 0;
	}

	double american_put_option_CRR(int k, int i, double current_stock_price, double** matrix) {
		if (matrix[k][i + k] != -1)
		{
			return matrix[k][i + k];
		}
		else
		{
			if (k == no_of_divisions)
				return max(0.0, (strike_price - current_stock_price));
			else
				matrix[k][i + k] = max((strike_price - current_stock_price),
					(p * american_put_option_CRR(k + 1, i + 1, current_stock_price * u, matrix) +
						(1 - p) * american_put_option_CRR(k + 1, i - 1, current_stock_price * d, matrix)) / R);
			return matrix[k][i + k];
		}
	}

	double american_call_option_CRR(int k, int i, double current_stock_price, double** matrix) {
		if (matrix[k][i + k] != -1)
		{
			return matrix[k][i + k];
		}
		else
		{
			if (k == no_of_divisions)
				return max(0.0, (current_stock_price - strike_price));
			else
				matrix[k][i + k] = max((current_stock_price - strike_price),
					(p * american_call_option_CRR(k + 1, i + 1, current_stock_price * u, matrix) +
						(1 - p) * american_call_option_CRR(k + 1, i - 1, current_stock_price * d, matrix)) / R);
			return matrix[k][i + k];
		}
	}
	
	double european_put_option_CRR(int k, int i, double current_stock_price, double** matrix) {
		if (matrix[k][i + k] != -1)
		{
			return matrix[k][i + k];
		}
		else
		{
			if (k == no_of_divisions)
				return max(0.0, (strike_price - current_stock_price));
			else
				matrix[k][i + k] = (p * european_put_option_CRR(k + 1, i + 1, current_stock_price * u, matrix) +
						(1 - p) * european_put_option_CRR(k + 1, i - 1, current_stock_price * d, matrix)) / R;
			return matrix[k][i + k];
		}
	}

	double european_call_option_CRR(int k, int i, double current_stock_price, double** matrix) {
		if (matrix[k][i + k] != -1)
		{
			return matrix[k][i + k];
		}
		else
		{
			if (k == no_of_divisions)
				return max(0.0, (current_stock_price - strike_price));
			else
				matrix[k][i + k] = (p * european_call_option_CRR(k + 1, i + 1, current_stock_price * u, matrix) +
						(1 - p) * european_call_option_CRR(k + 1, i - 1, current_stock_price * d, matrix)) / R;
			return matrix[k][i + k];
		}
	}
};

class backward_binomial_option_pricing_optimized {
public:

	string option_type;
	double dividend_rate, risk_free_rate, strike_price;
	double initial_stock_price, expiration_time, volatility;
	double u, d, p;
	int no_of_divisions;
	bool is_simplified_u;

	double option_pricing() {

		clock_t start = clock();

		double delta_t = expiration_time / (double)no_of_divisions;
		double r_inv = exp(-risk_free_rate * delta_t);

		if (is_simplified_u) {
			//simplified U
			u = exp(volatility * sqrt(delta_t));
			d = exp(-volatility * sqrt(delta_t));
			p = (exp((risk_free_rate - dividend_rate) * delta_t) - d) / (u - d);
		}
		else {
			// complex u
			double a = exp((risk_free_rate - dividend_rate) * delta_t);
			double b2 = a * a * (exp(volatility * volatility * delta_t) - 1);
			double temp = a * a + b2 + 1;
			u = (temp + sqrt(temp * temp - 4 * a * a)) / (2 * a);
			d = 1 / u;
			p = (a - d) / (u - d);
		}

		double q = 1 - p;
		double p_ = p * r_inv;
		double q_ = q * r_inv;

		//memorization
		double* V = new double[2 * no_of_divisions + 1];
		double* S = new double[2 * no_of_divisions + 1];
		S[no_of_divisions] = initial_stock_price;

		// terminal stock value
		for (int i = 1; i <= no_of_divisions; i++) {
			S[no_of_divisions - i] = S[no_of_divisions - i + 1] * d;
			S[no_of_divisions + i] = S[no_of_divisions + i - 1] * u;
		}

		// terminal option value
		if ((option_type == "AC") || (option_type == "EC")) {
			for (int i = 0; i <= 2 * no_of_divisions; i += 2) {
				V[i] = max((S[i] - strike_price), 0);
			}
		}
		else if ((option_type == "AP") || (option_type == "EP")) {
			for (int i = 0; i <= 2 * no_of_divisions; i += 2) {
				V[i] = max((strike_price - S[i]), 0);
			}
		}
		else {
			cout << "Unknown/Unsupport option type" << endl;
		}

		//backward pricing
		if (option_type == "AC") {

			//American call option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					V[j + no_of_divisions] = max(p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1], S[j + no_of_divisions] - strike_price);
				}
			}

			cout << "CRR Binomial Price of an American Call Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else if (option_type == "AP") {

			//American put option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					V[j + no_of_divisions] = max(p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1], strike_price - S[j + no_of_divisions]);
				}
			}

			cout << "CRR Binomial Price of an American Put Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else if (option_type == "EC") {

			// European call option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					V[j + no_of_divisions] = p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1];
				}
			}

			cout << "CRR Binomial Price of an European Call Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else if (option_type == "EP") {

			//  European put option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					V[j + no_of_divisions] = p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1];
				}
			}

			cout << "CRR Binomial Price of an European Put Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else {
			cout << "Unknown/Unsupport option type" << endl;
		}

		double option_price = V[no_of_divisions];
		delete[]V;
		delete[]S;
		return option_price;
	}

	double parameter() {
		cout << "Opmitized Binomial Option Pricing(CRR)" << endl;
		cout << "Expiration Time (Years) = " << expiration_time << endl;
		cout << "Number of Divisions = " << no_of_divisions << endl;
		cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
		cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
		cout << "Initial Stock Price = " << initial_stock_price << endl;
		cout << "Strike Price = " << strike_price << endl;
		cout << "Up Factor (CRR) = " << u << endl;
		cout << "Down Factor (CRR) = " << d << endl;
		cout << "Uptick Probability (CRR) = " << p << endl;
		cout << "Dividend Rate = " << dividend_rate << endl;
		cout << "--------------------------------------" << endl;
		return 0;
	}
};

class backward_binomial_option_pricing_BBS {
public:

	string option_type;
	double dividend_rate, risk_free_rate, strike_price;
	double initial_stock_price, expiration_time, volatility;
	double u, d, p;
	int no_of_divisions;
	bool is_simplified_u;

	double option_price_call_black_scholes(const double& S,       // spot (underlying) price
		const double& K,       // strike (exercise) price,
		const double& r,       // interest rate
		const double& d,       // dividend rate
		const double& sigma,   // volatility 
		const double& time) {  // time to maturity 
		double time_sqrt = sqrt(time);
		double d1 = (log(S / K) + (r - d) * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
		double d2 = d1 - (sigma * time_sqrt);
		return S * exp(-d * time) * N(d1) - K * exp(-r * time) * N(d2);
	}

	double option_price_put_black_scholes(const double& S,      // spot price
		const double& K,      // Strike (exercise) price,
		const double& r,      // interest rate
		const double& d,       // dividend rate
		const double& sigma,  // volatility
		const double& time) {
		double time_sqrt = sqrt(time);
		double d1 = (log(S / K) + (r - d) * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
		double d2 = d1 - (sigma * time_sqrt);
		return K * exp(-r * time) * N(-d2) - S * exp(-d * time) * N(-d1);
	}

	double N(const double& z) {
		if (z > 6.0) { return 1.0; }; // this guards against overflow 
		if (z < -6.0) { return 0.0; };
		double b1 = 0.31938153;
		double b2 = -0.356563782;
		double b3 = 1.781477937;
		double b4 = -1.821255978;
		double b5 = 1.330274429;
		double p = 0.2316419;
		double c2 = 0.3989423;
		double a = fabs(z);
		double t = 1.0 / (1.0 + a * p);
		double b = c2 * exp((-z) * (z / 2.0));
		double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
		n = 1.0 - b * n;
		if (z < 0.0) n = 1.0 - n;
		return n;
	}

	double option_pricing() {

		clock_t start = clock();

		double delta_t = expiration_time / (double)no_of_divisions;
		double r_inv = exp(-risk_free_rate * delta_t);

		if (is_simplified_u) {
			//simplified U
			u = exp(volatility * sqrt(delta_t));
			d = exp(-volatility * sqrt(delta_t));
			p = (exp((risk_free_rate - dividend_rate) * delta_t) - d) / (u - d);
		}
		else {
			// complex u
			double a = exp((risk_free_rate - dividend_rate) * delta_t);
			double b2 = a * a * (exp(volatility * volatility * delta_t) - 1);
			double temp = a * a + b2 + 1;
			u = (temp + sqrt(temp * temp - 4 * a * a)) / (2 * a);
			d = 1 / u;
			p = (a - d) / (u - d);
		}

		double q = 1 - p;
		double p_ = p * r_inv;
		double q_ = q * r_inv;

		//memorization
		double* V = new double[2 * no_of_divisions + 1];
		double* S = new double[2 * no_of_divisions + 1];
		S[no_of_divisions] = initial_stock_price;

		// terminal stock value
		for (int i = 1; i <= no_of_divisions; i++) {
			S[no_of_divisions - i] = S[no_of_divisions - i + 1] * d;
			S[no_of_divisions + i] = S[no_of_divisions + i - 1] * u;
		}

		// terminal option value
		if ((option_type == "AC") || (option_type == "EC")) {
			for (int i = 0; i <= 2 * no_of_divisions; i += 2) {
				V[i] = max((S[i] - strike_price), 0);
			}
		}
		else if ((option_type == "AP") || (option_type == "EP")) {
			for (int i = 0; i <= 2 * no_of_divisions; i += 2) {
				V[i] = max((strike_price - S[i]), 0);
			}
		}
		else {
			cout << "Unknown/Unsupport option type" << endl;
		}

		//backward pricing
		if (option_type == "AC") {

			//American call option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					if (i == no_of_divisions - 1) {
						V[j + no_of_divisions] = max(S[j + no_of_divisions] - strike_price,
							option_price_call_black_scholes(S[j + no_of_divisions], strike_price, risk_free_rate,
								dividend_rate, volatility, delta_t));
					}
					else {
						V[j + no_of_divisions] = max(p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1], S[j + no_of_divisions] - strike_price);
					}
				}
			}

			cout << "BBS Binomial Price of an American Call Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else if (option_type == "AP") {

			//American put option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					if (i == no_of_divisions - 1) {
						V[j + no_of_divisions] = max(strike_price - S[j + no_of_divisions],
							option_price_put_black_scholes(S[j + no_of_divisions], strike_price, risk_free_rate,
								dividend_rate, volatility, delta_t));
					}
					else {
						V[j + no_of_divisions] = max(p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1], strike_price - S[j + no_of_divisions]);
					}
					V[j + no_of_divisions] = max(p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1], strike_price - S[j + no_of_divisions]);
				}
			}

			cout << "BBS Binomial Price of an American Put Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else if (option_type == "EC") {

			// European call option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					if (i == no_of_divisions - 1) {
						V[j + no_of_divisions] = option_price_call_black_scholes(S[j + no_of_divisions], strike_price, risk_free_rate,
								dividend_rate, volatility, delta_t);
					}
					else {
						V[j + no_of_divisions] = p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1];
					}					
				}
			}

			cout << "BBS Binomial Price of an European Call Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else if (option_type == "EP") {

			//  European put option
			for (int i = no_of_divisions - 1; i >= 0; i--) {
				for (int j = -i; j <= i; j += 2) {
					if (i == no_of_divisions - 1) {
						V[j + no_of_divisions] = option_price_put_black_scholes(S[j + no_of_divisions], strike_price, risk_free_rate,
							dividend_rate, volatility, delta_t);
					}
					else {
						V[j + no_of_divisions] = p_ * V[j + no_of_divisions + 1] + q_ * V[j + no_of_divisions - 1];
					}
				}
			}

			cout << "BBS Binomial Price of an European Put Option = " << V[no_of_divisions] << endl;
			clock_t end = clock();
			cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
			cout << "--------------------------------------" << endl;
		}
		else {
			cout << "Unknown/Unsupport option type" << endl;
		}
		double option_price = V[no_of_divisions];
		delete []V;
		delete []S;
		return option_price;
	}

	double parameter() {
		cout << "BBS Binomial Option Pricing(CRR)" << endl;
		cout << "Expiration Time (Years) = " << expiration_time << endl;
		cout << "Number of Divisions = " << no_of_divisions << endl;
		cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
		cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
		cout << "Initial Stock Price = " << initial_stock_price << endl;
		cout << "Strike Price = " << strike_price << endl;
		cout << "Up Factor (CRR) = " << u << endl;
		cout << "Down Factor (CRR) = " << d << endl;
		cout << "Uptick Probability (CRR) = " << p << endl;
		cout << "Dividend Rate = " << dividend_rate << endl;
		cout << "--------------------------------------" << endl;
		return 0;
	}
};

class backward_binomial_option_pricing_BBSR {
public:

	string option_type;
	double dividend_rate, risk_free_rate, strike_price;
	double initial_stock_price, expiration_time, volatility;
	int no_of_divisions;
	bool is_simplified_u;
	// we are using BBS method in this class instead of write it twice
	double option_pricing()
	{
		//no_of_divisions
		backward_binomial_option_pricing_BBS option0;
		option0.expiration_time = expiration_time;
		option0.no_of_divisions = no_of_divisions;
		option0.risk_free_rate = risk_free_rate;
		option0.volatility = volatility;
		option0.initial_stock_price = initial_stock_price;
		option0.strike_price = strike_price;
		option0.dividend_rate = dividend_rate;
		option0.option_type = option_type;
		option0.is_simplified_u = is_simplified_u;

		// 1/2 no_of_divisions
		backward_binomial_option_pricing_BBS option1;
		option1.expiration_time = expiration_time;
		option1.no_of_divisions = floor(0.5 * no_of_divisions);
		option1.risk_free_rate = risk_free_rate;
		option1.volatility = volatility;
		option1.initial_stock_price = initial_stock_price;
		option1.strike_price = strike_price;
		option1.dividend_rate = dividend_rate;
		option1.option_type = option_type;
		option1.is_simplified_u = is_simplified_u;
		
		double BBSR_price = 2 * option0.option_pricing() - option1.option_pricing();
		if (option_type == "AC") {
			cout << "BBSR Binomial Price of an American Call Option = " << BBSR_price << endl;
		}
		else if (option_type == "AP") {
			cout << "BBSR Binomial Price of an American Put Option = " << BBSR_price << endl;
		}
		else if (option_type == "EC") {
			cout << "BBSR Binomial Price of an European Call Option = " << BBSR_price << endl;
		}
		else if (option_type == "EP") {
			cout << "BBSR Binomial Price of an European Put Option = " << BBSR_price << endl;
		}
		return 2 * option0.option_pricing() - option1.option_pricing();
	}


};


/*int main(int argc, char* argv[])
{

	backward_binomial_option_pricing_optimized option2;

	option2.expiration_time = 0.5;
	option2.no_of_divisions = 10000;
	option2.risk_free_rate = 0.05;
	option2.volatility = 0.2;
	option2.initial_stock_price = 100;
	option2.strike_price = 100;
	option2.dividend_rate = 0.02;
	option2.option_type = "AC";
	option2.is_simplified_u = true;

	clock_t start = clock();
	for (int i = 1; i <= 100; i++) {
		option2.option_pricing();
	}

	clock_t end = clock();
	cout << "It took-- " << (long double)(end - start) / CLOCKS_PER_SEC * 0.01 << " seconds for average" << endl;
	
	system("pause");
}
*/
