#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "option.h"

using namespace std;

double Option :: GetCallOptionProfit(){
	return call_option_profit;
};

void Option :: SetCallOptionProfit(double my_call_option_profit){
    call_option_profit = my_call_option_profit;
}

double Option :: GetPutOptionProfit(){
	return put_option_profit;
};

void Option :: SetPutOptionProfit(double my_put_option_profit){
    put_option_profit = my_put_option_profit;
}

European :: European(){
	asset_price = 0;
    strike_price = 0;
};

European :: ~European(){};

void European :: SetStrikePrice(double my_strike_price){
    strike_price = my_strike_price;
}

void European :: SetAssetPrice(double my_asset_price){
    asset_price = my_asset_price;
}

void European :: UpdateCallOptionProfit(){
    double diff = asset_price - strike_price;
    if (diff <=0 ){
        this->SetCallOptionProfit(0);
    } else {
        this->SetCallOptionProfit(diff);
    };
};

void European :: UpdatePutOptionProfit(){
    double diff = strike_price - asset_price;
    if (diff <=0 ){
        this->SetPutOptionProfit(0);
    } else {
        this->SetPutOptionProfit(diff);
    }
};
