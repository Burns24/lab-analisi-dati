#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "asset.h"

using namespace std;

// Metodi Classe Asset
Asset :: Asset(){asset_price = 0;}

Asset :: ~Asset(){}

double Asset :: GetAssetPrice(){
	return asset_price;
};

void Asset :: SetAssetPrice(double my_asset_price){
	asset_price = my_asset_price;
};

// Metodi Classe GBM
GBM :: GBM(double my_drift, double my_volatility){
	drift = my_drift;
	volatility = my_volatility;
};

GBM :: ~GBM(){}

void GBM :: SetGBMGaussianVar(double my_gaussian_var){
	gaussian_var = my_gaussian_var;
};

double GBM :: GetDrift(){
	return drift;
};

double GBM :: GetVolatility(){
	return volatility;
};

void GBM :: UpdateAssetPrice(double initial_date, double expire_date){
	double old_asset_price = this->GetAssetPrice();
	double new_asset_price = old_asset_price*exp( (drift - 0.5*pow(volatility,2.) )*(expire_date - initial_date) + volatility*gaussian_var*sqrt(expire_date - initial_date) );
	this->SetAssetPrice(new_asset_price);
};
