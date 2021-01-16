/*
 * vortexparticles.h
 *
 *  Created on: Oct 7, 2008
 *      Author: bzhao
 */

#ifndef VORTEXPARTICLES_H_
#define VORTEXPARTICLES_H_

#include <list>
using std::list;

#include "voxel.h"
#include "particle.h"


class VortexParticles{

public:
	VortexParticles(u_int n, float s, float t, int X, int Y, int Z)
	: N(n), support(s), initVorticity(t), DimX(X), DimY(Y), DimZ(Z){

	}
	~VortexParticles(){}
	void AddVortexParticle(const VortexParticle &p){
		Particles.push_back(p);
	}
	void Seed(const Voxel &voxel, int i, int j, int k);
	
	void VortexForce(const Voxel &voxel, 
						float *u0, float *v0, float *w0,
						float *phi_u, float *phi_v, float *phi_w,
						float *phi_u_obj, float *phi_v_obj, float *phi_w_obj);
	list<VortexParticle> Particles;
	
private:
	u_int N;
	float support;
	float initVorticity;	
	int DimX;
	int DimY;
	int DimZ;
};


#endif /* VORTEXPARTICLES_H_ */
