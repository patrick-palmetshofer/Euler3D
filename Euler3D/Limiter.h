#pragma once
#include "GlobalTypes.h"
class Limiter
{
public:
	Limiter();
	virtual ~Limiter();

	//Limiter function using the ratio. Can lead to problems with 0/0 inputs
	StateVector2D limiter(const StateVector2D &r);

	//Limiter function using explicit values. Returns Limiter function in terms of r
	virtual StateVector2D limiter(const StateVector2D & x, const StateVector2D & y) = 0;
};

class LimiterMinmod :
	public Limiter
{
public:
	StateVector2D limiter(const StateVector2D & x, const StateVector2D & y);
};

class LimiterVanAlbada :
	public Limiter
{
	StateVector2D limiter(const StateVector2D & x, const StateVector2D & y);
};

class LimiterMonotoneCentered :
	public Limiter
{
	StateVector2D limiter(const StateVector2D & x, const StateVector2D & y);
};

class NoLimiter :
	public Limiter
{
	StateVector2D limiter(const StateVector2D & x, const StateVector2D & y) { StateVector2D ret = { 1,1,1,1 }; return ret; };
};