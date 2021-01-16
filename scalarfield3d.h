#ifndef SCALARFIELD3D_H_
#define SCALARFIELD3D_H_

#include <string.h>
#include "field3d.h"
#include "particle.h"
#include "voxel.h"
#include "vectorfield3d.h"
#include "minheap.h"
#include "geometry.h"
#include "blendparticles.h"
#include "CIsoSurface.h"
#include "wavegenerator.h"

#ifdef COUPLED_SPH
#include "sph.h"
#endif


#define MAX_CURVATURE  0.5f

#define RADIUS_MAX 0.5f*delta
#define RADIUS_MIN 0.1f*delta

#define SEED_PARTICLE_ITERATIONS 15

class ScalarField3D : public Field3D{

public:
	ScalarField3D(PhysicsWorld *pw, int nx, int ny, int nz,
				  int i, int j, int k)
		: Field3D(nx,ny,nz,pw),
		  I(i), J(j), K(k){
		phi = new float[nx*ny*nz];
//		memset(phi, 0, nx*ny*nz);
		SetInfinity(phi);
		phi0 = new float[nx*ny*nz];
//		memset(phi0, 0, nx*ny*nz);
		SetInfinity(phi0);
		phiPos = new float[nx*ny*nz];
//		memset(phiPos, 0, nx*ny*nz);
		SetInfinity(phiPos);
		phiNeg = new float[nx*ny*nz];
//		memset(phiNeg, 0, nx*ny*nz);
		SetInfinity(phiNeg);
		source = NULL;
		LiquidLastStep = 0;
//		phiPos = new float[(nx+1)*(ny+1)*(nz+1)];
//		memset(phiPos, 0, (nx+1)*(ny+1)*(nz+1));
//		phiNeg = new float[(nx+1)*(ny+1)*(nz+1)];
//		memset(phiNeg, 0, (nx+1)*(ny+1)*(nz+1));

#ifdef CUDA
		obj     = new char[nx*ny*nz];
		mov_obj = new char[nx*ny*nz];
		memset(obj,     0, nx*ny*nz*sizeof(char));
		memset(mov_obj, 0, nx*ny*nz*sizeof(char));
		dims[0] = nx;
		dims[1] = ny;
		dims[2] = nz;
#endif
	}
	~ScalarField3D(){
		if(phi){
			delete [] phi;
			phi = NULL;
		}
		if(phi0){
			delete [] phi0;
			phi0 = NULL;
		}
		if(phiPos){
			delete [] phiPos;
			phiPos = NULL;
		}
		if(phiNeg){
			delete [] phiNeg;
			phiNeg = NULL;
		}
		if(source)
			delete source;
#ifdef CUDA
		if(obj){
			delete [] obj;
			obj = NULL;
		}
		if(mov_obj){
			delete [] mov_obj;
			mov_obj = NULL;
		}
#endif
	}
	void Print() const;
	void OutputGridData(const Voxel &voxel, char *) const;
	void OutputBinaryGridData(const Voxel &voxel, char *filename) const;
	void ReadBinaryGridData(char *filename);
	void ReadinGridData(char *);
	void WriteRestart(FILE* fp) const;
	bool ReadRestart(FILE* fp);
	void EvaluatePhi(const vector<Particle> &particles,
					map<u_int, KnownPoint> &band,
					const Voxel &v);
	void EvaluatePhi(const Voxel &voxel, const Point &pos, float radius);
	void EvaluatePhi(const Voxel &voxel,
					int start_x0, int start_y0, int start_z0,
					int end_x0,  int end_y0,   int end_z0);
	void EvaluatePhiAirWater(const Voxel &voxel,
					int start_x0, int start_y0, int start_z0,
					int end_x0,  int end_y0,   int end_z0);
	void EvaluatePhiAirWater(const Voxel &voxel,
					int start_x0, int start_y0, int start_z0,
					int end_x0,  int end_y0,   int end_z0,
					int start_x1, int start_y1, int start_z1,
					int end_x1,  int end_y1,   int end_z1);
	void EvaluatePhiAirWater(const Voxel &voxel,float d, float H);
	void FastMarching(map<u_int, KnownPoint> &posband,
 				      map<u_int, KnownPoint> &negband,
					  Voxel &v, bool inital);
	void FastMarchingAirWater(map<u_int, KnownPoint> &posband,
						 	map<u_int, KnownPoint> &negband,
						    Voxel &v, bool initial, float *phi_tmp);

	float UpdateTentativeValue(const KnownPoint &p, u_int index,
						 	   map<u_int, KnownPoint> &band,
						 	   const Voxel &v, float *ph);
	float ProcessOneTerm(int axis, float phi, const Voxel &v);

	float ProcessTwoTerms(int axis1, int axis2,
							  float phi1, float phi2,
							  const Voxel &v);

  	float ProcessThreeTerms(int axis1, int axis2, int axis3,
						    float phi1, float phi2, float phi3,
						    const Voxel &v);

  	void InitializeFromWaveGen(const WaveGenerator *wave, const Voxel *voxel){
  		wave->SetInitialWaveLevelset(voxel, DimX, DimY, DimZ,
  				WALL_THICKNESS+1, WALL_THICKNESS+1, WALL_THICKNESS+1, phi);
  	}

#ifdef SPMD
  	bool NearInterface(int i, int j, int k, float e, float limit, float *tmp) const;
	void ReInitialize(Voxel &voxel, const ScalarField3D &solid, bool ini);
	void ReInitialize(Voxel &voxel, const ScalarField3D &solid, bool ini, float *value, float *value0);
#if 0
	void ReInitializeForRendering(const Voxel &voxel, float *value, float *value0);
#endif
	void FindSmoothFactor(const Voxel &voxel, char *needed, float *value0, float *S);
	void FindSolidNeighbors(const Voxel &voxel,
			const char *needed, const float *x0,
			char *axis_x, char *axis_y, char *axis_z) const;
	void EulerStep(const Voxel &voxel, float dtau, const char *need, const float *S,
			const char *axis_x, const char *axis_y, const char *axis_z,
			float *value, const float *value0);
#else
	void ReInitialize(Voxel &voxel,
					  map<u_int, KnownPoint> &pos_band,
					  map<u_int, KnownPoint> &neg_band,
					  const ScalarField3D &solid);
#endif
	void Initialize(Voxel &voxel,
					map<u_int, KnownPoint> &posband,
					map<u_int, KnownPoint> &negband,
					const ScalarField3D &solid);
	void InitializeObject( Voxel &voxel,
					map<u_int, KnownPoint> &pos_band,
					map<u_int, KnownPoint> &neg_band);
	void InitializeInterface(Voxel &voxel,
					map<u_int, KnownPoint> &pos_band,
					map<u_int, KnownPoint> &neg_band,
					const ScalarField3D &solid,
					float *phi_tmp);
	float LossasoDistance(const Voxel &voxel, const Point &p, float phimin) const;

	void Recompute(const Voxel &voxel, float *phi_tmp) const{
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
//		 	TentativeGridPoint tp(0.f, i, j, k);
			if(voxel.InSolid(i,j,k)){
				phi_tmp[pos] = max(phi[pos], phi_tmp[pos]);
//				if(phi[INDEX(i,j,k)] > phi_tmp[INDEX(i,j,k)])
//					phi_tmp[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
			}
			// this modification will make the boundary between solid and liquid
		    // look smooth, only for visualizing purpose
//			else if(voxel.CloseToSolid(tp) && phi_tmp[INDEX(i,j,k)] <= 0.f)
//			else if(voxel.CloseToSolid(tp))
			else if(voxel.CloseToSolidTwoRings(i,j,k))
				phi_tmp[pos] = max(phi[pos], phi_tmp[pos]);
		END_FOR
	}

	void RecomputeMovingObj(const Voxel &voxel, float *phi_tmp) const{
		FOR_EACH_CELL
		    // if a cell is occupied by a "moving" solid object, we
		    // reset its levelset to phi_c_obj (a positive quantity)
			if(voxel.InMovingSolid(i,j,k)){
				if(phi[INDEX(i,j,k)] > phi_tmp[INDEX(i,j,k)])
					phi_tmp[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
			}
		END_FOR
	}
	void ExtrapolateIntoObject(const Voxel &voxel, const ScalarField3D &obj, bool limit){
//		u_int NM = 0;
//		FOR_EACH_CELL
//			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//				if(phi[INDEX(i,j,k)] <= 0.f)
//					++NM;
//			}
//		END_FOR
//		printf("\nbefore extrapolate into object there are %ld liquid points \n\n", NM);
#ifdef SPMD
		SetEqual(phi, phi0);
		obj.Extrapolate(voxel, phi, phi0, limit);
#else
		obj.Extrapolate(voxel, phi, limit);
#endif
//		NM = 0;
//		FOR_EACH_CELL
//			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//				if(phi[INDEX(i,j,k)] <= 0.f)
//					++NM;
//			}
//		END_FOR
//		printf("\nafter extrapolate into object there are %ld liquid points \n\n", NM);
		FOR_EACH_CELL
			if(i == I && j == J && k == K){
			    printf("after extrapolate into object: at i = %d, j = %d, k = %d,  phi = %20.16f, "
			   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
				 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
				 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
				 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
			}
		END_FOR
		obj.RecomputeMovingObj(voxel, phi);
	}
	void ExtrapolateIntoObject1(const Voxel &voxel, const ScalarField3D &obj, bool limit){
//		u_int NM = 0;
//		FOR_EACH_CELL
//			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//				if(phi[INDEX(i,j,k)] <= 0.f)
//					++NM;
//			}
//		END_FOR
//		printf("\nbefore extrapolate into object there are %ld liquid points \n\n", NM);
#ifdef SPMD
		SetEqual(phi, phi0);
		obj.Extrapolate(voxel, phi, phi0, limit);
#else
		obj.Extrapolate(voxel, phi, limit);
#endif
//		NM = 0;
//		FOR_EACH_CELL
//			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//				if(phi[INDEX(i,j,k)] <= 0.f)
//					++NM;
//			}
//		END_FOR
//		printf("\nafter extrapolate into object there are %ld liquid points \n\n", NM);
		FOR_EACH_CELL
			if(i == I && j == J && k == K){
				printf("after extrapolate into object: at i = %d, j = %d, k = %d,  phi = %20.16f, "
					"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
						I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
						phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
						phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
			}
		END_FOR
	}

	void ResetObject(const Voxel &voxel, const ScalarField3D &obj){
		obj.Recompute(voxel, phi);
//#ifdef SPMD
//		SetEqual(phi, phi0);
//		ReInitializeForRendering(voxel, phi, phi0);
//#endif
	}

	void CopyFrom(float *s){
		FOR_EACH_CELL
			phi[INDEX(i,j,k)] = s[INDEX(i,j,k)];
		END_FOR
	}

	void CopyTo(ScalarField3D &s){
		s.CopyFrom(phi);
	}


	void BuilSignedDistanceMap(const Voxel *voxel){
		mPhysWorld->BuildSignedDistanceMap(voxel, 0, phi);
	}

#ifdef SPMD
	void Extrapolate(const Voxel &voxel, float *phi_tmp, float *phi_tmp0, bool limit) const;
#else
	void Extrapolate(const Voxel &voxel, float *phi_tmp, bool limit) const;
#endif

	void print_grad_phi(int I, int J, int K) const;

	float TriInterp(const Voxel &v, const Point &p) const;
	float TriInterpDebug(const Voxel &voxel, const Point &p) const;

	inline u_int Index(const TentativeGridPoint &p){
		return p.kk*DimX*DimY+p.jj*DimX+p.ii;
	}

	inline u_int Index(int i, int j, int k){
		return k*DimX*DimY+j*DimX+i;
	}

	bool IndexOutofBounds(int i, int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ)
			return true;
		else
			return false;
	}

	void setScalar(int i, int j, int k, float rhs){
		phi[k*DimX*DimY+j*DimX+i] = rhs;
	}

	float getScalar(int i, int j, int k) const{
		return phi[k*DimX*DimY+j*DimX+i];
	}

	float* getScalar() const{
		return phi;
	}

	void SetEqual(float *x, float *x0){
		for(int i = 0; i < DimX*DimY*DimZ; ++i)
			x0[i]  = x[i];
	}


	void RecomputeObjectLevelset(const Voxel &voxel, VectorField3D &vel){
		vel.RecomputeObjectLevelset(voxel, phi);
	}

#ifdef SPMD
	bool CFLExceedsOne(const Voxel &voxel, VectorField3D &vel, float dts, float *cfl) const{
		return vel.CFLExceedsOne(voxel, dts, phi, cfl);
	}
	float FindTimeStep(const Voxel &voxel, VectorField3D &vel) const{
		return vel.FindTimeStep(voxel, phi);
	}
#endif

	void Advect(const Voxel &voxel, VectorField3D &vel){
		SWAP(phi0, phi);
		u_int N = 0;
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi0[INDEX(i,j,k)] <= 0.f)
				++N;
		END_FOR
		printf("Before advection, there are %u liquid points \n", N);
		u_int N2 = 0;
		FOR_EACH_CELL
			if(voxel.InSolid(i,j,k) && phi0[INDEX(i,j,k)] < 0.f)
				++N2;
		END_FOR
		printf("\n before advection there are %u solid cells having < 0 phi\n\n", N2);
#ifdef SPMD
//		vel.Advect(voxel, phi, phi0);
//#else
		vel.AdvectSemiLagrangian(voxel, phi, phi0);
#endif
		N = 0;
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f){
//				printf("in ScalarField::Advect() at (%d, %d, %d) phi = %16.13f \n", i, j, k, phi[INDEX(i,j,k)]);
				++N;
			}
		END_FOR
		printf("After advection, there are %u liquid points \n", N);
	}
#ifdef SPMD
	void AdvectSubcycling(const Voxel &voxel, float dts, VectorField3D &vel){
		SWAP(phi0, phi);
		u_int N = 0;
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi0[INDEX(i,j,k)] <= 0.f)
				++N;
		END_FOR
		printf("Before advection, there are %u liquid points \n", N);
		u_int N2 = 0;
		FOR_EACH_CELL
			if(voxel.InSolid(i,j,k) && phi0[INDEX(i,j,k)] < 0.f)
				++N2;
		END_FOR
		printf("\n before advection there are %u solid cells having < 0 phi\n\n", N2);

		vel.AdvectSubcycling(voxel, dts, phi, phi0);
		N = 0;
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f){
//				printf("in ScalarField::Advect() at (%d, %d, %d) phi = %16.13f \n", i, j, k, phi[INDEX(i,j,k)]);
				++N;
			}
		END_FOR
		printf("After advection, there are %u liquid points \n", N);
	}
#endif

	void AdvectEscapedParticles(const Voxel &voxel,
			          list<Particle> &posParticles,
			          list<Particle> &negParticles,
			          float dt) const;

	void CorrectNeighborLevelSet(const Voxel &voxel, const Point &pos,
		                         const Particle &ps, float delta, char sign);
	void ErrorCorrection(const Voxel &voxel,
//						list<Particle> &particles
						ParticleList &particles
#ifdef POSITIVE_PARTICLES
//						,list<Particle> &posParticles
						,ParticleList &posParticles
#endif
								);

	float Curvature(const Voxel &voxel, int i, int j, int k) const;

	void UpdateVoxel(Voxel &voxel);

	void AttractParticles(const Voxel &voxel, const ScalarField3D &object,
							list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
							list<Particle> &posParticles,
#endif
							int nPerCell) const;
	void SeedParticles(const Voxel &voxel,
			ParticleList &negParticles,
#ifdef POSITIVE_PARTICLES
		    ParticleList &posParticles,
#endif
			int nPerCell) const;

	void SeedParticles(const Voxel &voxel, const ScalarField3D &object,
						list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
						list<Particle> &posParticles,
#endif
						int nPerCell) const;

	void SeedParticlesAtSource(const Voxel &voxel,
							list<Particle> &negParticles,
							int start_x0, int start_y0,  int start_z0,
							int end_x0,   int end_y0,   int  end_z0,
							//list<Particle> &posParticles,
							int nPerCell) const;
	void NumOfNonEscapedParticles(const Voxel &voxel, int *cellParticles, int s,
				list<Particle> &particles) const;
	void ReSeedSourceParticles(const Voxel &voxel, const ScalarField3D &object,
				list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
				list<Particle> &posParticles,
#endif
				int nPerCell) const;
	void ReSeedParticles(const Voxel &voxel,
				const ScalarField3D &object,
				ParticleList &negParticles,
#ifdef POSITIVE_PARTICLES
				ParticleList &posParticles,
#endif
				int nPerCell) const;
	void ReSeedParticles2(const Voxel &voxel, const ScalarField3D &object,
					list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
					list<Particle> &posParticles,
#endif
					int nPerCell) const;
	void ReSeedParticlesNew(const Voxel &voxel, const ScalarField3D &object,
					list<Particle> &negParticles,
	#ifdef POSITIVE_PARTICLES
					list<Particle> &posParticles,
	#endif
					int nPerCell) const;

	void AdjustParticleRadii(const Voxel &voxel, const VectorField3D &vel,
				int n, int iters,
//				list<Particle> &negParticles
				ParticleList &negParticles
#ifdef POSITIVE_PARTICLES
//				,list<Particle> &posParticles
				,ParticleList &posParticles
#endif
#ifdef COUPLED_SPH
				, SPH *water
#else
#ifdef COUPLED_FLIP
				, WaterParticleList &wpl
#else
				,BlendParticles *water
#endif
#endif
//				,list<WaterParticle> &absorb
				,WaterParticleList &absorb
				) const;
	void AdjustWaterParticleRadii(const Voxel &voxel,
							list<WaterParticle> &currentWaterParticles,
							BlendParticles *water) const;

	void GradPhi(const Voxel &voxel, const Point &p,  Vector &g) const;
	void NormalAt(const Voxel &voxel, const Point &p,  Vector &g) const;

	void SourceTerm(const Voxel &voxel,
					int start_x0, int start_y0,  int start_z0,
					int end_x0,   int end_y0,   int  end_z0,
					int nx, int ny, int nz, float speed) {

		source = new SourceRegion( start_x0, start_y0,  start_z0,
				 					end_x0,  end_y0,   end_z0,
				 					nx, ny, nz,
				 					speed);
		int i, j, k;
		k = start_z0;
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i)
				source->AddSource(i,j,k);
		k = end_z0;
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i)
				source->AddSource(i,j,k);
//		k = start_z0;
//		i = start_x0-1;
//		for(j = start_y0; j <= end_y0; ++j){
//			source->AddSource(i,j,k);
//		}
//		i = end_x0+1;
//		for(j = start_y0; j <= end_y0; ++j){
//			source->AddSource(i,j,k);
//		}
//		j = start_y0-1;
//		for(i = start_x0; i <= end_x0; ++i){
//			source->AddSource(i,j,k);
//		}
//		j = end_y0+1;
//		for(i = start_x0; i <= end_x0; ++i){
//			source->AddSource(i,j,k);
//		}

	}

	void ApplySourceTerm(const Voxel &voxel, bool initial){
		for(int n = 0; n < sources.size();++n){
			sources[n]->Merge(voxel, phi, initial);
		}
		FOR_EACH_CELL
			if(i == I && j == J && k == K){
				printf("after apply source before advection: at i = %d, j = %d, k = %d,  phi = %20.16f, "
					"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
						I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
						phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
						phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
			}
		END_FOR
	}

	void SetZero(float *data){
		FOR_EACH_CELL
			data[INDEX(i,j,k)] = 0.f;
		END_FOR
	}
	void SetInfinity(float *data){
		FOR_EACH_CELL
			data[INDEX(i,j,k)] = INFINITY;
		END_FOR
	}
	void SetZero(char *data){
//		FOR_EACH_CELL
//			data[INDEX(i,j,k)] = 0;
//		END_FOR
		memset(data, 0, DimX*DimY*DimZ*sizeof(char));
	}

	bool CheckIndex(int i, int j, int k) const {
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
			printf("Error in (ScalarField3D::CheckIndex())! index out of range i=%d, j=%d, k=%d \n",
					i,j,k);
			exit(1);
//			return false;
		}
		else
			return true;
	}

	bool CheckIndexNoExit(int i, int j, int k) const {
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ ){
			printf("Error in (ScalarField3D::CheckIndexNoExit())! index out of range i=%d, j=%d, k=%d \n",
					i,j,k);
			exit(1);
			return false;
		}
		else
			return true;
	}

	void ConstructMesh(const Voxel &voxel){
		float h = voxel.VoxelDelta();
		mc.GenerateSurface(phi, 0, DimX-1, DimY-1, DimZ-1, h, h, h);
	}
	void DestroyMesh(){
		mc.DeleteSurface();
	}
	void OutputGridDataRib(char *fname){
		mc.RIBDump(fname);
	}

	void MergeWithWaterParticles(const BlendParticles *waterSpray){
		FOR_EACH_CELL
			if(waterSpray->phi[INDEX(i,j,k)] < 0.f)
				phi[INDEX(i,j,k)] = phi[INDEX(i,j,k)] < waterSpray->phi[INDEX(i,j,k)] ?
									phi[INDEX(i,j,k)] : waterSpray->phi[INDEX(i,j,k)];
		END_FOR
	}
	void MergeWith(const Voxel &voxel, const ScalarField3D &s){
		s.AndOperator(voxel, phi);
	}
	void AndOperator(const Voxel &voxel, float *phi_tmp) const{
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k))
				phi_tmp[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)] < phi[INDEX(i,j,k)] ?
									phi_tmp[INDEX(i,j,k)] : phi[INDEX(i,j,k)];
		END_FOR
	}
	void CheckLocalMinima(const Voxel &voxel, int m, int n, int t) const{
		float delta = voxel.VoxelDelta();
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
				if( phi[INDEX(i,j,k)] > 0.f && phi[INDEX(i,j,k)]+E_EPSIL < FASTMARCH_LIMIT ){
					if( phi[INDEX(i,j,k)] < phi[INDEX(i+1,j,k)] &&
						phi[INDEX(i,j,k)] < phi[INDEX(i-1,j,k)] &&
						phi[INDEX(i,j,k)] < phi[INDEX(i,j+1,k)] &&
						phi[INDEX(i,j,k)] < phi[INDEX(i,j-1,k)] &&
						phi[INDEX(i,j,k)] < phi[INDEX(i,j,k+1)] &&
						phi[INDEX(i,j,k)] < phi[INDEX(i,j,k-1)] ){
						printf("at m = %d, n = %d, t = %d, (%d, %d, %d) found local minima level set = %f \n",
								m, n, t, i, j, k, phi[INDEX(i,j,k)]);
						printf("at i = %d, j = %d, k = %d,  phi = %f, "
							"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
								i, j, k, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
								phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
								phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);

						exit(1);
					}
				}
			}
		END_FOR
	}

	void UpdateStatus(const Voxel &voxel){
		u_int N = 0;
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f)
					++N;
		END_FOR
		printf("\nat this step, there are %u liquid points while last step with %u difference is %d \n\n",
					N, LiquidLastStep, N-LiquidLastStep);
		LiquidLastStep = N;
	}

private:
//	int DimX;
//	int DimY;
//	int DimZ;
	float *phi, *phi0;
	float *phiPos, *phiNeg;
	int I, J, K;
	SourceRegion *source;
	CIsoSurface<float> mc;
	mutable u_int LiquidLastStep;

#ifdef CUDA
	char *obj;
	char *mov_obj;
	int dims[3];
#endif

	ScalarField3D(const ScalarField3D &a)
	{
		assert(false);
	}

	ScalarField3D& operator=(const ScalarField3D &a)
	{ assert(false); return *this; }

};


#endif /*SCALARFIELD3D_H_*/
