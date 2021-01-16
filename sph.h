/*
 * sph.h
 *
 *  Created on: Sep 5, 2008
 *      Author: bzhao
 */

#ifndef SPH_H_
#define SPH_H_

#ifdef WIN32
// VC++  < v7.1 in header "hash_map" in namespace std.
// VC++ >= v7.1 in header "hash_map" in namespace stdext (std version deprecated in VC++).
#include <hash_map>
using stdext::hash_map;
#else
#include <tr1/unordered_map>
using namespace std::tr1;
#endif

#include "geometry.h"
#include "voxel.h"
#include "timemanager.h"

class VectorField3D;

class SPH{
public:
	SPH( float _h, TimeManager *_tm, float _g, 
	     float ap, int nx, int ny, int nz, int i, int j, int k)
	  : h(_h), tm(_tm), g(_g), 
	    DimX(nx), DimY(ny), DimZ(nz), I(i), J(j), K(k), 
	    averagePeriod(ap), influenceRadius(1.5f*_h),
	    searchRadius(2.f*_h){		
		
		radiusSquared = searchRadius * searchRadius;
		cWeights = new float[nx*ny*nz];
		uWeights = new float[nx*ny*nz];
		vWeights = new float[nx*ny*nz];
		wWeights = new float[nx*ny*nz];
//		u = new float[nx*ny*nz];
//		v = new float[nx*ny*nz];
//		w = new float[nx*ny*nz];
		u0 = new float[nx*ny*nz];
		v0 = new float[nx*ny*nz];
		w0 = new float[nx*ny*nz];
		SetZero(cWeights);
		SetZero(uWeights);
		SetZero(vWeights);
		SetZero(wWeights);
//		SetZero(u);
//		SetZero(v);
//		SetZero(w);
		SetZero(u0);
		SetZero(v0);
		SetZero(w0);
	}
	~SPH(){
		delete [] cWeights;
		delete [] uWeights;
		delete [] vWeights;
		delete [] wWeights;
//		delete [] u;
//		delete [] v;
//		delete [] w;
		delete [] u0;
		delete [] v0;
		delete [] w0;
	}
	void AddWaterParticles(const WaterParticle &p);
	void TimeStepping(const Voxel *voxel, const VectorField3D *vel);
	void ApplyForce();
	void ComputeWeights(const Voxel *voxel);
	void ComputeTargetDensity(const Voxel *voxel, int n);
	void ApplyParticleSlip(const Voxel *voxel);
	void VelocityChangeFLIP(const Voxel *voxel, const Point &pos, 
					float &uvel, float &vvel, float &wvel) const;
	void FloodFill(int ii, int jj, int kk, char *color) const;
	void Projection(const Voxel *voxel, const float *upls, const float *vpls, const float *wpls);
	void LinearSolver( const Voxel *voxel, float *x,
#ifdef WIN32
					   hash_map<u_int, u_int> &a, hash_map<u_int, float> &x0,
#else
					   unordered_map<u_int, u_int> &a, unordered_map<u_int, float> &x0,
#endif
					   u_int n );
	void OutputWaterParticlesBinary(const Voxel *voxel, const char *name);
	list<WaterParticle> waterParticles;
private:
	void SetZero(float *data){
		FOR_EACH_CELL
			data[INDEX(i,j,k)] = 0.f;
		END_FOR
	}
	void SetZero(char *data){
		FOR_EACH_CELL
			data[INDEX(i,j,k)] = 0;
		END_FOR
	}
	float h;
	TimeManager *tm;
	float g;
	int DimX;
	int DimY;
	int DimZ;
	int I;
	int J;
	int K;
	float averagePeriod;
	float incompressibilityNumberDensity;
	float influenceRadius;	
	float searchRadius;
	float radiusSquared;
	float *cWeights,*uWeights,*vWeights,*wWeights;
	
//	float *u, *v, *w, *u0, *v0, *w0;
	float *u0, *v0, *w0;

};

#endif /* SPH_H_ */
