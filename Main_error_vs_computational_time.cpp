#include "Binomial Option Pricing.cpp"
#include "Trinomial Option Pricing.cpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;


int main(int argc, char* argv[])
{
	//define parameter sets for short to long term, and ATM/OTM/ITM option
	double parameter_sets[7][9];
	for (int i = 0; i <= 2; i++)
	{
		parameter_sets[0][3*i] = 0.25; //expiration_time-three month
		parameter_sets[0][3 * i + 1] = 0.5; //expiration_time-half year
		parameter_sets[0][3 * i + 2] = 1; //expiration_time-one year

		parameter_sets[1][3 * i] = 0.012675; //risk_free_rate-three month
		parameter_sets[1][3 * i + 1] = 0.015713; //risk_free_rate-half year
		parameter_sets[1][3 * i + 2] = 0.026138; //risk_free_rate-one year
		//https://www.global-rates.com/en/interest-rates/libor/american-dollar/american-dollar.aspx

		parameter_sets[2][3 * i] = 0.137488582; //volatility-three month
		parameter_sets[2][3 * i + 1] = 0.137488582; //volatility-half year
		parameter_sets[2][3 * i + 2] = 0.137488582; //volatility-one year
		//https://finance.yahoo.com/quote/%5EGSPC/history?period1=1601337600&period2=1632960000&interval=1d&filter=history&frequency=1d&includeAdjustedClose=true

		parameter_sets[3][3 * i] = 100; //initial_stock_price-three month
		parameter_sets[3][3 * i + 1] = 100; //initial_stock_price-half year
		parameter_sets[3][3 * i + 2] = 100; //initial_stock_price-one year

		parameter_sets[5][3 * i] = 0.0138; //dividend_rate-three month
		parameter_sets[5][3 * i + 1] = 0.0138; //dividend_rate-half year
		parameter_sets[5][3 * i + 2] = 0.0138; //dividend_rate-one year
		//https://ycharts.com/indicators/sp_500_dividend_yield

	}
	for (int i = 0; i <= 8; i++)
	{
		if (i < 3) {
			parameter_sets[4][i] = 90; //strike_price-in the money
		}
		else if (i < 6){
			parameter_sets[4][i] = 100; //strike_price-at the money
		}
		else if (i < 9) {
			parameter_sets[4][i] = 110; //strike_price-out of money
		}
	}

	parameter_sets[6][0] = 250; //three month-in the money
	parameter_sets[6][1] = 650; //half year-in the money
	parameter_sets[6][2] = 1250; //one year-in the money
	parameter_sets[6][3] = 2100; //three month-at the money
	parameter_sets[6][4] = 6100; //half year-at the money
	parameter_sets[6][5] = 12100; //one year-at the money
	parameter_sets[6][6] = 2200; //three month-out of money
	parameter_sets[6][7] = 6200; //half year-out of money
	parameter_sets[6][8] = 12200; //one year-out of money

	//test each parameter set
	for (int option_status = 0; option_status < 9; option_status++) {

		//initialize each pricing method
		backward_binomial_option_pricing_optimized option0;
		option0.expiration_time = parameter_sets[0][option_status];
		option0.risk_free_rate = parameter_sets[1][option_status];
		option0.volatility = parameter_sets[2][option_status];
		option0.initial_stock_price = parameter_sets[3][option_status];
		option0.strike_price = parameter_sets[4][option_status];
		option0.dividend_rate = parameter_sets[5][option_status];
		option0.option_type = "AC";

		backward_binomial_option_pricing_BBS option1;
		option1.expiration_time = parameter_sets[0][option_status];
		option1.risk_free_rate = parameter_sets[1][option_status];
		option1.volatility = parameter_sets[2][option_status];
		option1.initial_stock_price = parameter_sets[3][option_status];
		option1.strike_price = parameter_sets[4][option_status];
		option1.dividend_rate = parameter_sets[5][option_status];
		option1.option_type = "AC";
		option1.is_simplified_u = true;

		backward_binomial_option_pricing_BBSR option2;
		option2.expiration_time = parameter_sets[0][option_status];
		option2.risk_free_rate = parameter_sets[1][option_status];
		option2.volatility = parameter_sets[2][option_status];
		option2.initial_stock_price = parameter_sets[3][option_status];
		option2.strike_price = parameter_sets[4][option_status];
		option2.dividend_rate = parameter_sets[5][option_status];
		option2.option_type = "AC";
		option2.is_simplified_u = true;

		backward_Trinomial_option_pricing option3;
		option3.expiration_time = parameter_sets[0][option_status];
		option3.risk_free_rate = parameter_sets[1][option_status];
		option3.volatility = parameter_sets[2][option_status];
		option3.initial_stock_price = parameter_sets[3][option_status];
		option3.strike_price = parameter_sets[4][option_status];
		option3.dividend_rate = parameter_sets[5][option_status];

		//error compare
		//calculate true value using BBSR
		option2.no_of_divisions = 10000;
		double BBSR_true_price = option2.option_pricing();

		double simplified_price, complex_price, BBS_price, BBSR_price, Tri_price;

		string file_name;
		ofstream oFile_error;
		if (parameter_sets[6][option_status] == 250) {
			file_name = "Error Compare/Error Compare-three month-in the money.csv";
		}
		else if (parameter_sets[6][option_status] == 650) {
			file_name = "Error Compare/Error Compare-half year-in the money.csv";
		}
		else if (parameter_sets[6][option_status] == 1250) {
			file_name = "Error Compare/Error Compare-one year-in the money.csv";
		}
		else if (parameter_sets[6][option_status] == 2100) {
			file_name = "Error Compare/Error Compare-three month-at the money.csv";
		}
		else if (parameter_sets[6][option_status] == 6100) {
			file_name = "Error Compare/Error Compare-half year-at the money.csv";
		}
		else if (parameter_sets[6][option_status] == 12100) {
			file_name = "Error Compare/Error Compare-one year-at the money.csv";
		}
		else if (parameter_sets[6][option_status] == 2200) {
			file_name = "Error Compare/Error Compare-three month-out of money.csv";
		}
		else if (parameter_sets[6][option_status] == 6200) {
			file_name = "Error Compare/Error Compare-half year-out of money.csv";
		}
		else if (parameter_sets[6][option_status] == 12200) {
			file_name = "Error Compare/Error Compare-one year-out of money.csv";
		}

		oFile_error.open(file_name, ios::out | ios::ate);
		oFile_error << "Error compare for lattice methods" << endl;
		oFile_error << "expiration_time" << "," << parameter_sets[0][option_status] << ","
			<< "risk_free_rate" << "," << parameter_sets[1][option_status] << ","
			<< "volatility" << "," << parameter_sets[2][option_status] << ","
			<< "initial_stock_price" << "," << parameter_sets[3][option_status] << ","
			<< "strike_price" << "," << parameter_sets[4][option_status] << ","
			<< "dividend_rate" << "," << parameter_sets[5][option_status] << ","
			<< endl;

		oFile_error << "N" << ","
			<< "Appendix B (Simplified U)" << "," << "Appendix B (Simplified U) Run time" << ","
			<< "Appendix B (Complex U)" << "," << "Appendix B (Complex U) Run time" << ","
			<< "BBS(Simplified U)" << "," << "BBS(Simplified U) Run time" << ","
			<< "BBSR(Simplified U)" << "," << "BBSR(Simplified U) Run time" << ","
			<< "Trinomial" << "," << "Trinomial Run time" << ","
			<< "True Value*" << endl;
		
		//even number N, step size is 50
		for (int i = 50; i <= 1000; i += 50)
		{
			//Appendix B Simplified U
			clock_t start0 = clock();

			option0.is_simplified_u = true;
			option0.no_of_divisions = i;
			for (int j = 1; j <= 100; j++) {
				simplified_price = option0.option_pricing();
			}
			clock_t end0 = clock();
			//cout << "It took-- " << (long double)(end0 - start0) / CLOCKS_PER_SEC * 0.01 << " seconds for average(Appendix B simplified u)" << endl;

			//Appendix B Complex U
			clock_t start1 = clock();
			option0.is_simplified_u = false;
			option0.no_of_divisions = i;
			for (int j = 1; j <= 100; j++) {
				complex_price = option0.option_pricing();
			}
			clock_t end1 = clock();
			//cout << "It took-- " << (long double)(end1 - start1) / CLOCKS_PER_SEC * 0.01 << " seconds for average(Appendix B complex u)" << endl;

			//BBS Simplified U
			clock_t start2 = clock();
			option1.no_of_divisions = i;
			for (int j = 1; j <= 100; j++) {
				BBS_price = option1.option_pricing();
			}
			clock_t end2 = clock();
			//cout << "It took-- " << (long double)(end2 - start2) / CLOCKS_PER_SEC * 0.01 << " seconds for average(Appendix B complex u)" << endl;

			//BBSR Simplified U
			clock_t start3 = clock();
			option2.no_of_divisions = i;
			for (int j = 1; j <= 100; j++) {
				BBSR_price = option2.option_pricing();
			}
			clock_t end3 = clock();
			//cout << "It took-- " << (long double)(end3 - start3) / CLOCKS_PER_SEC * 0.01 << " seconds for average(Appendix B complex u)" << endl;

			//Trinomial Simplified U
			clock_t start4 = clock();
			option3.no_of_divisions = i;
			for (int j = 1; j <= 100; j++) {
				Tri_price = option3.option_pricing();
			}
			clock_t end4 = clock();
			//cout << "It took-- " << (long double)(end3 - start3) / CLOCKS_PER_SEC * 0.01 << " seconds for average(Appendix B complex u)" << endl;

			oFile_error << i << ","
				<< simplified_price << "," << (long double)(end0 - start0) / CLOCKS_PER_SEC * 0.01 << ","
				<< complex_price << "," << (long double)(end1 - start1) / CLOCKS_PER_SEC * 0.01 << ","
				<< BBS_price << "," << (long double)(end2 - start2) / CLOCKS_PER_SEC * 0.01 << ","
				<< BBSR_price << "," << (long double)(end3 - start3) / CLOCKS_PER_SEC * 0.01 << ","
				<< Tri_price << "," << (long double)(end4 - start4) / CLOCKS_PER_SEC * 0.01 << ","
				<< BBSR_true_price << endl;
			cout << i << endl;

		}
		oFile_error.close();
	}
	system("pause");
}
