#include "Limiter.h"



Limiter::Limiter()
{
}


Limiter::~Limiter()
{
}

StateVector Limiter::limiter(const StateVector & r)
{
	StateVector y;
	y << 1, 1, 1, 1, 1;
	return limiter(r, y);
}


//Safer form of the minmod limiter
StateVector LimiterMinmod::limiter(const StateVector &x, const StateVector &y)
{
	StateVector phis;
	for (int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		double yi = y[i];
		double phi = 0;

		if (xi*yi > 0 && std::isfinite(xi) && std::isfinite(yi))
		{
			if (std::abs(xi) < std::abs(yi))
				phi = xi / yi;
			else
				phi = 1;
		}
		phis[i] = phi;// *2 / (1 + r);
	}
	return phis;
}


//Monotone centered limiter
StateVector LimiterMonotoneCentered::limiter(const StateVector &x, const StateVector &y)
{
	StateVector phis;
	for (int i = 0; i < x.size(); i++)
	{
		double r = x[i]/y[i];
		double m = std::min(2 * r, 2.0);
		double phi = std::min(m, 0.5*(1 + r));
		phi = std::max(phi, 0.0);
		phis[i] = phi;
		if (!std::isfinite(r))
			phis[i] = 0;
	}
	return phis;
}

StateVector LimiterVanAlbada::limiter(const StateVector &x, const StateVector &y)
{
	StateVector phis;
	for (int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		double yi = y[i];
		double phi = 0;

		if (xi*yi > 0 && std::isfinite(xi) && std::isfinite(yi))
		{
			phi = 2 / (xi / yi + yi / xi);
		}
		phis[i] = phi;// *2 / (1 + r);
	}
	return phis;
}