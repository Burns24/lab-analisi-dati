#ifndef ASSET_H
#define ASSET_H

class Asset{

private:
	double asset_price;
protected:

public:
	Asset();
	~Asset();

	// Metodi
	virtual void UpdateAssetPrice(double initial_date, double expire_date){};
	void SetAssetPrice(double my_asset_value);
	double GetAssetPrice();
};

class GBM :public Asset{
	/*
	Classe cheimplementa il moto geometrico browniano con media=drift e std.dev=volatility
	per il metodo UpdateAssetValue(t).
	*/
private:
	double drift, volatility;
	double gaussian_var;
protected:

public:
	GBM(double my_drift=0, double my_volatility=1);
	~GBM();

	// Metodi
	void UpdateAssetPrice(double initial_date, double expire_date);
	void SetGBMGaussianVar(double my_gaussian_var);
	double GetDrift();
	double GetVolatility();
};

#endif 
