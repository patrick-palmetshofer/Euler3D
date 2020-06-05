#include "FluxRoe.h"



FluxRoe::FluxRoe(int new_dim) : Flux(new_dim)
{
}


FluxRoe::~FluxRoe()
{
}

//Roe solver
StateVector FluxRoe::calcDissip(std::pair<StateVector, StateVector> leftrightstates, DirVector &n)
{
	prim_left = fluid->cons2prim(leftrightstates.first);
	prim_right = fluid->cons2prim(leftrightstates.second);

	double sqrtrho_right = std::sqrt(prim_right[0]);
	double sqrtrho_left = std::sqrt(prim_left[0]);

	double rho = sqrtrho_right * sqrtrho_left;

	// Calculate Roe averaged variables
	Eigen::Vector3d uvec = (sqrtrho_right*prim_right.segment(1,3) + sqrtrho_left * prim_left.segment(1, 3)) / (sqrtrho_right + sqrtrho_left);
	Eigen::Vector3d duvec = (prim_right - prim_left).segment(1, 3);

	double pleft = fluid->calcPprim(prim_left);
	double pright = fluid->calcPprim(prim_right);

	double Hleft = prim_left[4] + pleft / prim_left[0];
	double Hright = prim_right[4] + pright / prim_right[0];
	//double Hright = prim_right[3] + (gamma - 1)*(prim_right[3] - 0.5*(prim_right[1] * prim_right[1] + prim_right[2] * prim_right[2]));
	//double Hleft = prim_left[3] + (gamma - 1)*(prim_left[3] - 0.5*(prim_left[1] * prim_left[1] + prim_left[2] * prim_left[2]));

	double h = (sqrtrho_right*Hright + sqrtrho_left * Hleft) / (sqrtrho_right + sqrtrho_left);

	double unorm = uvec.dot(n);
	Eigen::Vector3d upar = uvec.cross(n); //Check this
	Eigen::Vector3d uzero = upar.cross(n);

	double c = std::sqrt((fluid->getGamma() - 1)*(h - 0.5*uvec.squaredNorm()));

	double drho = prim_right[0] - prim_left[0];
	double dp = pright - pleft;
	//Check this
	double dunorm = duvec.dot(n);

	Eigen::Vector3d dupar = duvec.cross(n);
	Eigen::Vector3d duzero = dupar.cross(n);

	//Calculate transformation matrices
	roevectors <<
		1,
		uvec[0] - c * n[0],
		uvec[1] - c * n[1],
		uvec[2] - c * n[2],
		h - unorm * c,

		1,
		uvec[0],
		uvec[1],
		uvec[2],
		0.5*uvec.squaredNorm(),

		0,
		dupar[0],
		dupar[1],
		dupar[2],
		upar.dot(dupar),

		0,
		duzero[0],
		duzero[1],
		duzero[2],
		uzero.dot(duzero),

		1,
		uvec[0] + c * n[0],
		uvec[1] + c * n[1],
		uvec[2] + c * n[2],
		h + unorm * c;


	roefactors[0] = (dp - rho * c*dunorm) / (2 * c*c);
	roefactors[1] = -(dp - c * c * drho) / (c*c);
	roefactors.segment(2,2) = rho;
	roefactors[4] = (dp + rho * c*dunorm) / (2 * c*c);

	//Entropy fix
	double unorm_L = prim_left.segment(1,3).matrix().dot(n);
	double unorm_R = prim_right.segment(1, 3).matrix().dot(n);

	eigenvals[0] = unorm - c;// std::min(, unorm_L - c);
	eigenvals[1] = unorm;
	eigenvals[2] = unorm;
	eigenvals[3] = unorm;
	eigenvals[4] = unorm + c;

	dissip = roevectors.transpose() * (eigenvals.abs()*roefactors).matrix();

	return dissip;
}
