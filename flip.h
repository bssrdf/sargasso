/*
 * flip.h
 *
 *  Created on: Oct 29, 2010
 *      Author: bzhao
 */

#ifndef FLIP_H_
#define FLIP_H_
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

class ScalarField3D;
class VectorField3D;

class FLIP{
public:
	FLIP(float _h, const TimeManager *_tm, float _g,
	     float pf,  const ScalarField3D *obj,
	     int nx, int ny, int nz, int i, int j, int k)
	  : h(_h), mTimeManager(_tm), g(_g),
	    DimX(nx), DimY(ny), DimZ(nz), I(i), J(j), K(k),
	    mPICFLIPRatio(pf), mObject(obj){

		invh = 1.f / _h;
		mWeights = new float[nx*ny*nz];
		u 		 = new float[nx*ny*nz];
		v 		 = new float[nx*ny*nz];
		w 		 = new float[nx*ny*nz];
		du 		 = new float[nx*ny*nz];
		dv 		 = new float[nx*ny*nz];
		dw 	     = new float[nx*ny*nz];
		mMarker  = new  char[nx*ny*nz];
		phi      = new float[nx*ny*nz];
		tmp      = new float[nx*ny*nz];
		SetZero(mWeights);
		SetZero(u);
		SetZero(v);
		SetZero(w);
		SetZero(du);
		SetZero(dv);
		SetZero(dw);
		SetZero(mMarker);
		SetZero(phi);
		SetZero(tmp);

	}
	~FLIP(){
		delete [] mWeights;
		delete [] u;
		delete [] v;
		delete [] w;
		delete [] du;
		delete [] dv;
		delete [] dw;
		delete [] mMarker;
		delete [] phi;
		delete [] tmp;
	}


	void TimeStepping(const Voxel *voxel, const VectorField3D *vel);
	void MoveParticlesInGrid(const Voxel *vox);
	void TransferParticleValuesToGrid(const Voxel *voxel);
	void SaveGridVelocities();
	void ApplyForce(float dt);
	void UpdateGridVelocityIncrement();
	void UpdateParticleVelocityFromGridValues();
	void ComputeGridLevelset();
	void ExtrapolateGridVelocity();
	void ApplyBoundaryConditions(const Voxel *voxel, const VectorField3D *vel);
	void Projection(const Voxel *voxel, float dt);
	void AddWaterParticles(const WaterParticle &p);
	void OutputWaterParticlesBinary(const Voxel *voxel, const char *name);


//	void ComputeWeights(const Voxel *voxel);
//	void ComputeTargetDensity(const Voxel *voxel, int n);
//	void ApplyParticleSlip(const Voxel *voxel);
//	void VelocityChangeFLIP(const Voxel *voxel, const Point &pos,
//					float &uvel, float &vvel, float &wvel) const;

//	void Projection(const Voxel *voxel, const float *upls, const float *vpls, const float *wpls);
//	void LinearSolver( const Voxel *voxel, float *x,
//#ifdef WIN32
//					   hash_map<u_int, u_int> &a, hash_map<u_int, float> &x0,
//#else
//					   unordered_map<u_int, u_int> &a, unordered_map<u_int, float> &x0,
//#endif
//					   u_int n );



	WaterParticleList mWaterParticles;

private:
	void AdvectRK3Subcycling(const Voxel *voxel, float dts);

	float TriInterp(const Voxel *voxel, const Point &p, float *phi) const;
	void  TriInterp(float &ui, float &vi, float &wi,
			const Point &p,	float *uvel, float *vvel, float *wvel) const;
	Vector ObjectNormalAt(const Voxel *voxel, const Point &p, float *phi_c_obj) const;

	void AccumulateValues(const Point &pos, float vel, int index, float *wh, float *vc);

	void InitializePhi();
	void SweepPhi();
	void SweepVelocity();
	void SweepU(int i0, int i1, int j0, int j1, int k0, int k1);
	void SweepV(int i0, int i1, int j0, int j1, int k0, int k1);
	void SweepW(int i0, int i1, int j0, int j1, int k0, int k1);
	void CheckVelocity() const;
	void CheckNewmannBCSurroundedRegion(const Voxel *voxel);
	void FloodFill(int ii, int jj, int kk, char *color) const;

	void LinearSolver(const Voxel *voxel, float *x,
#ifdef WIN32
					  hash_map<u_int, u_int> &a, hash_map<u_int, float> &x0,
#else
					  unordered_map<u_int, u_int> &a, unordered_map<u_int, float> &x0,
#endif
		  			 u_int n );

	void PrintOneCellVel(int i, int j, int k, const char *str) const{
		u_int pos = INDEX(i,j,k);
		printf(" %s: at (%d, %d, %d) vel = (%f,%f,%f)  vel1 = (%f,%f,%f) \n",
				str, i, j, k, u[pos], v[pos], w[pos],
				u[INDEX(i-1,j,k)], v[INDEX(i,j-1,k)], w[INDEX(i,j,k-1)]);
	}
	void PrintOneCellPhi(int i, int j, int k, const char *str) const{
		u_int pos = INDEX(i,j,k);
		printf(" %s: at (%d, %d, %d) phi = %f  phi(i-1)= %f phi(i+1) = %f"
				" phi(j-1) = %f phi(j+1) = %f phi(k-1) = %f phi(k+1) = %f \n",
				str, i, j, k, phi[pos], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
				phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
				phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
	}


	void SetZero(float *data){
		FOR_EACH_CELL
			data[INDEX(i,j,k)] = 0.f;
		END_FOR
	}
	void SetZero(float *data, u_int N){
		for(u_int i = 0; i < N; ++i)
			data[i] = 0.f;

	}
	void SetZero(char *data){
		FOR_EACH_CELL
			data[INDEX(i,j,k)] = 0;
		END_FOR
	}

	void SetEqual(float *x, float *x0){
		for(u_int i = 0; i < DimX*DimY*DimZ; ++i)
			x0[i]  = x[i];
	}

	void SetScalar(float rhs, float *val){
		for(u_int i = 0; i < DimX*DimY*DimZ; ++i)
			val[i] = rhs;
	}

	bool IndexOutofBounds(int i, int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ)
			return true;
		else
			return false;
	}

	float h;
	float invh;
	const TimeManager *mTimeManager;
	const ScalarField3D *mObject;
	float g;
	int DimX;
	int DimY;
	int DimZ;
	int I;
	int J;
	int K;
	float mPICFLIPRatio; // regularization parameter tuning PIC/FLIP

	float *mWeights;
	float *u, *v, *w, *du, *dv, *dw;
	float *phi;
	float *tmp;
	char  *mMarker;

};

#endif /* FLIP_H_ */
