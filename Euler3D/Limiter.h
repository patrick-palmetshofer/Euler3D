#pragma once
#include "GlobalTypes.h"
class Limiter
{
public:
	Limiter();
	virtual ~Limiter();

	//Limiter function using the ratio. Can lead to problems with 0/0 inputs
	StateVector limiter(const StateVector &r);

	//Limiter function using explicit values. Returns Limiter function in terms of r
	virtual StateVector limiter(const StateVector & x, const StateVector & y) = 0;
};

class LimiterMinmod :
	public Limiter
{
public:
	StateVector limiter(const StateVector & x, const StateVector & y);
};

class LimiterVanAlbada :
	public Limiter
{
	StateVector limiter(const StateVector & x, const StateVector & y);
};

class LimiterMonotoneCentered :
	public Limiter
{
	StateVector limiter(const StateVector & x, const StateVector & y);
};

class NoLimiter :
	public Limiter
{
	StateVector limiter(const StateVector & x, const StateVector & y) { StateVector ret; ret << 1, 1, 1, 1, 1; return ret; };
};