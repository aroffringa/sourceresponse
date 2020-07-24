#ifndef POLYNOMIAL_FITTER_H
#define POLYNOMIAL_FITTER_H

#include <aocommon/uvector.h>

#include <array>

class PolynomialFitter
{
public:
	void Clear() { _dataPoints.clear(); }
	
	void AddDataPoint(double x, double y, double w)
	{
		_dataPoints.emplace_back(std::array<double,3>{{x, y, w}});
	}
	
	void Fit(aocommon::UVector<double>& terms, size_t nTerms);
	
	static double Evaluate(double x, const aocommon::UVector<double>& terms)
	{
		double val = terms[0];
		double f = 1.0;
		for(size_t i=1; i!=terms.size(); ++i)
		{
			f *= x;
			val += f * terms[i];
		}
		return val;
	}
	
	size_t size() const { return _dataPoints.size(); }
	
private:
	aocommon::UVector<std::array<double,3>> _dataPoints;
};

#endif
