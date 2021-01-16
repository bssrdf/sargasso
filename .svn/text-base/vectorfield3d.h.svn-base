#ifndef VECTORFIELD3D_H_
#define VECTORFIELD3D_H_

#include <map>
#include <stack>
#include <algorithm>
#ifdef WIN32
// VC++  < v7.1 in header "hash_map" in namespace std.
// VC++ >= v7.1 in header "hash_map" in namespace stdext (std version deprecated in VC++).
#include <hash_map>
using stdext::hash_map;
#else
#include <tr1/unordered_map>
using namespace std::tr1;
#endif
//using namespace std;




#include "geometry.h"
#include "field3d.h"
#include "voxel.h"
#include "minheap.h"
#include "movingobj.h"
#include "sourceregion.h"
#include "vortexparticles.h"
#include "CIsoSurface.h"
#include "wavegenerator.h"

#ifdef COUPLED_SPH
#include "sph.h"
#endif

class Container;

class VectorField3D : public Field3D{
public:
	VectorField3D(PhysicsWorld *pw, int nx, int ny, int nz,
//			       float _dt, float _g, float _visc,
			       TimeManager *tm, float _g, float _visc,
			       bool mov, Container *c,
			       int i, int j, int k,
			       const WaveGenerator *wg=NULL)
		: Field3D(nx, ny, nz, pw),
		  g(_g), visc(_visc),
		  tmg(tm), hasMovingObj(mov), container(c),
		  I(i), J(j), K(k), mWaveGenerator(wg){
		if(_visc < 0.f){
			printf("Wrong! viscosity (%f) must be >= 0 \n", _visc);
			exit(1);
		}

		u = new float[nx*ny*nz];
		v = new float[nx*ny*nz];
		w = new float[nx*ny*nz];
		u0 = new float[nx*ny*nz];
		v0 = new float[nx*ny*nz];
		w0 = new float[nx*ny*nz];
		u1 = new float[nx*ny*nz];
		v1 = new float[nx*ny*nz];
		w1 = new float[nx*ny*nz];
		phi_u = new float[nx*ny*nz];
		phi_v = new float[nx*ny*nz];
		phi_w = new float[nx*ny*nz];
		if(hasMovingObj){
			phi_c_fixed_obj = new float[nx*ny*nz];
			phi_u_fixed_obj = new float[nx*ny*nz];
			phi_v_fixed_obj = new float[nx*ny*nz];
			phi_w_fixed_obj = new float[nx*ny*nz];
		}
		temp = new float[nx*ny*nz];
		needRecomputePhi= new char[nx*ny*nz];
#ifdef AIR_DIV_FREE
//		phiu_star = new float[nx*ny*nz];
//		phiv_star = new float[nx*ny*nz];
//		phiw_star = new float[nx*ny*nz];
//		u_star = new float[nx*ny*nz];
//		v_star = new float[nx*ny*nz];
//		w_star = new float[nx*ny*nz];
#endif
		SetZero(u);
		SetZero(v);
		SetZero(w);
		SetZero(u0);
		SetZero(v0);
		SetZero(w0);
		SetZero(u1);
		SetZero(v1);
		SetZero(w1);
		SetZero(phi_u);
		SetZero(phi_v);
		SetZero(phi_w);
		if(hasMovingObj){
			SetZero(phi_c_fixed_obj);
			SetZero(phi_u_fixed_obj);
			SetZero(phi_v_fixed_obj);
			SetZero(phi_w_fixed_obj);
		}
		SetZero(temp);
		SetZero(needRecomputePhi);
#ifdef AIR_DIV_FREE
//		memset(u_star, 0, nx*ny*nz);
//		memset(v_star, 0, nx*ny*nz);
//		memset(w_star, 0, nx*ny*nz);
//		memset(phiu_star, 0, nx*ny*nz);
//		memset(phiv_star, 0, nx*ny*nz);
//		memset(phiw_star, 0, nx*ny*nz);
//		SetZero(u_star);
//		SetZero(v_star);
//		SetZero(w_star);
//		SetZero(phiu_star);
//		SetZero(phiv_star);
//		SetZero(phiw_star);
#endif
		source = NULL;
		u_max = 0.f;
		v_max = 0.f;
		w_max = 0.f;
		u_max_step = 0;
		v_max_step = 0;
		w_max_step = 0;
#ifdef AIR_DIV_FREE
		residual_projection_air = 0.0;
#endif
		mc = new CIsoSurface<float>();
	}
	~VectorField3D(){
		delete [] u;
		delete [] v;
		delete [] w;
		delete [] u0;
		delete [] v0;
		delete [] w0;
		delete [] u1;
		delete [] v1;
		delete [] w1;
		delete [] phi_u;
		delete [] phi_v;
		delete [] phi_w;
//		delete [] phi_u_obj;
//		delete [] phi_v_obj;
//		delete [] phi_w_obj;
		if(hasMovingObj){
			delete [] phi_c_fixed_obj;
			delete [] phi_u_fixed_obj;
			delete [] phi_v_fixed_obj;
			delete [] phi_w_fixed_obj;
		}
		delete [] temp;
		delete [] needRecomputePhi;
#ifdef AIR_DIV_FREE
//		delete [] phiu_star;
//		delete [] phiv_star;
//		delete [] phiw_star;
//		delete [] u_star;
//		delete [] v_star;
//		delete [] w_star;
#endif
		if(source)
			delete source;
		delete mc;
	}
	void LinearCombine(const VectorField3D &rhs1,
					   const VectorField3D &rhs2,
					   float t1, float t2)	{
		for(int i=0;i<DimX*DimY*DimZ;i++){
			u[i] = t1*rhs1.u[i] + t2*rhs2.u[i];
			v[i] = t1*rhs1.v[i] + t2*rhs2.v[i];
			w[i] = t1*rhs1.w[i] + t2*rhs2.w[i];
		}
	}

	void TimeStepping(u_int m, u_int restart, Voxel &voxel, float *phi,
			                ParticleList &negParticles
#ifdef POSITIVE_PARTICLES
							,ParticleList &posParticles
#endif
							){
//		AssignPhi(phi);
//		FixedVelocity(voxel);
//		SetEqual(u, u1);
//		SetEqual(v, v1);
//		SetEqual(w, w1);
//		printf("Velocity updated \n");
//		if(hasMovingObj)
//			UpdateObjectBoundary(voxel);
//		printf(" at the beginning of TimeStepping u = %f, v = %f, w = %f \n",
//						u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)]);
//		SetEqual(u, u1);
//		SetEqual(v, v1);
//		SetEqual(w, w1);
		AssignPhi(voxel, phi);
#ifdef SPMD
		// according to Tariq Aslam, we can just use
		// phi_u, phi_v and phi_w interpolated from phi_c and there is no need
		// to reinitialize them
//		ReInitialize(voxel, 1, phi_u, u1);
//		ReInitialize(voxel, 2, phi_v, v1);
//		ReInitialize(voxel, 3, phi_w, w1);
#else
//		ReInitialize(voxel);
#endif
//		ErrorCorrection(voxel, negParticles
//#ifdef POSITIVE_PARTICLES
//				,posParticles
//#endif
//				);
		UpdateVoxel(voxel);

//		if(m == restart){
//			Extrapolate(voxel);
//		}
#ifdef AIR_DIV_FREE
		//ResetLquidVelocity();
#endif

//		printf(" after add gravity force u = %f, v = %f, w = %f \n",
//							u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)]);
//		Extrapolate(voxel);

//		printf(" after diffuse u = %f, v = %f, w = %f \n",
//					u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)]);
//		ApplySourceTerm(voxel);

		//ideally, advection of velocity is better done within a divergence-free
		//background velocity field.

		// we are not doing projection here because wrong (extrapolated) velocities
		// from last time step may be used to compute the divergence.
		// instead, the advection of velocity is performed after force is applied,
		// such that the projection uses valid fluid velocities (advected and forced).
		// at last, projection is done to make the final velocity divergence free.
//		Projection(voxel);

//		if(source)
//			ApplySourceTerm();
//		printf("before diffuse vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//						I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//						phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//						            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//						            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);

		/* we first compute velocity diffusion here */
		Diffuse(voxel);

		// now u1, v1 and w1 contain the intermediate velocity after diffusion
//		printf("before advect vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//								I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//								phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//								            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//								            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
		// we then do advection
		Advect(voxel);
//		AdvectMacCormack(voxel);

//		printf("before addforce vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//								I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//								phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//								            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//								            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
		/* we apply gravity and other body forces here */
		GravityForce(voxel);
		AddSource(voxel);


//		printf(" before diffuse u = %f, v = %f, w = %f \n",
//					u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)]);




		// see Bridson's book p.94 for why UpdateObjectBoundary() is called here
		UpdateObjectBoundary(voxel, false);

	}

	void TimeStepping1(u_int m, u_int restart, Voxel &voxel, float *phi){

		AssignPhi(voxel, phi);
		UpdateVoxel(voxel);

		ApplySourceTerm(voxel);

		Projection(voxel);

#ifdef AIR_DIV_FREE
//		SetEqual(phi_u, phiu_star);
//		SetEqual(phi_v, phiv_star);
//		SetEqual(phi_w, phiw_star);
#endif

//		K = 9;
		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
				I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
				phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
				            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
				            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
//		K = 10;
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//				I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//				phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//				            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//				            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
//		K = 12;
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//						I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//						phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//						            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//						            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
//		K = 14;
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//						I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//						phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//						            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//						            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
//		K = 16;
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//						I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//						phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//						            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//						            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
//		K = 18;
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//						I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//						phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//						            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//						            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
//		K = 20;
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//						I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
//						phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//						            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//						            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);

		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
				I, J-1, K, u[INDEX(I,J-1,K)], v[INDEX(I,J-1,K)], w[INDEX(I,J-1,K)],
				phi_c[INDEX(I,J-1,K)], u[INDEX(I,J-1,K)]-u[INDEX(I-1,J-1,K)]
				            		 +v[INDEX(I,J-1,K)]-v[INDEX(I,J-2,K)]
				            	     +w[INDEX(I,J-1,K)]-w[INDEX(I,J-1,K-1)]);

		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
						I, J+1, K, u[INDEX(I,J+1,K)], v[INDEX(I,J+1,K)], w[INDEX(I,J+1,K)],
						phi_c[INDEX(I,J+1,K)],u[INDEX(I,J+1,K)]-u[INDEX(I-1,J+1,K)]
		   				            		 +v[INDEX(I,J+1,K)]-v[INDEX(I,J,K)]
		   				            	     +w[INDEX(I,J+1,K)]-w[INDEX(I,J+1,K-1)]);

		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
					I+1, J, K, u[INDEX(I+1,J,K)], v[INDEX(I+1,J,K)], w[INDEX(I+1,J,K)],
					phi_c[INDEX(I+1,J,K)],u[INDEX(I+1,J,K)]-u[INDEX(I,J,K)]
	   				            		 +v[INDEX(I+1,J,K)]-v[INDEX(I+1,J-1,K)]
	   				            	     +w[INDEX(I+1,J,K)]-w[INDEX(I+1,J,K-1)]);

		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
					I-1, J, K, u[INDEX(I-1,J,K)], v[INDEX(I-1,J,K)], w[INDEX(I-1,J,K)],
					phi_c[INDEX(I-1,J,K)],u[INDEX(I-1,J,K)]-u[INDEX(I-2,J,K)]
	   				            		 +v[INDEX(I-1,J,K)]-v[INDEX(I-1,J-1,K)]
	   				            	     +w[INDEX(I-1,J,K)]-w[INDEX(I-1,J,K-1)]);

		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
					I, J, K+1, u[INDEX(I,J,K+1)], v[INDEX(I,J,K+1)], w[INDEX(I,J,K+1)],
					phi_c[INDEX(I,J,K+1)],u[INDEX(I,J,K+1)]-u[INDEX(I-1,J,K+1)]
	   				            		 +v[INDEX(I,J,K+1)]-v[INDEX(I,J-1,K+1)]
	   				            	     +w[INDEX(I,J,K+1)]-w[INDEX(I,J,K)]);

		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
					I, J, K-1, u[INDEX(I,J,K-1)], v[INDEX(I,J,K-1)], w[INDEX(I,J,K-1)],
					phi_c[INDEX(I,J,K-1)],u[INDEX(I,J,K-1)]-u[INDEX(I-1,J,K-1)]
	   				            		 +v[INDEX(I,J,K-1)]-v[INDEX(I,J-1,K-1)]
	   				            	     +w[INDEX(I,J,K-1)]-w[INDEX(I,J,K-2)]);

//		for(int i=I; i<35; ++i)
//		for(int k=K; k>20; --k)
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, phi_obj = %f \n",
//						I, J, k, u[INDEX(I,J,k)], v[INDEX(I,J,k)], w[INDEX(I,J,k)],
//						phi_c[INDEX(I,J,k)],phi_c_obj[INDEX(I,J,k)]);

//		printf("index = %ld, index1 = %ld \n", INDEX(I,J,K), INDEX(I,J,K-1));

//		FOR_EACH_CELL
//			//if(voxel.InLiquid(i,j,k) || voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] < 5.f){
//			if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
////			if( k == K && voxel.InAir(i,j,k) && voxel.InLiquid(i-1,j,k) && voxel.InLiquid(i,j-1,k)
////			&& voxel.InLiquid(i+1,j,k) && voxel.InLiquid(i,j+1,k) && voxel.InLiquid(i,j,k+1) ){
//				float div = u[INDEX(i,j,k)]-u[INDEX(i-1,j,k)]
//		            		 +v[INDEX(i,j,k)]-v[INDEX(i,j-1,k)]
//		            	     +w[INDEX(i,j,k)]-w[INDEX(i,j,k-1)];
//				if(fabs(div) > 1.e-4){
////					printf("at (%d, %d, %d) divergence = %f \n",
////							i,j,k, div);
//					printf("(%d,%d,%d) phi = %f, u = %f, v = %f, w = %f, div = %f\n",
//										i, j, k, phi_c[INDEX(i,j,k)],
//										u[INDEX(i,j,k)], v[INDEX(i,j,k)], w[INDEX(i,j,k)], div);
//				}
//			}
////			if( !voxel.InSolid(i,j,k) && w[INDEX(i,j,k)] > 0.f ){
////				printf("(%d,%d,%d) phi = %f, phi(k-1) = %f, u = %f, v = %f, w = %f \n",
////						i, j, k, phi_c[INDEX(i,j,k)],  phi_c[INDEX(i,j,k-1)],
////						u[INDEX(i,j,k)], v[INDEX(i,j,k)], w[INDEX(i,j,k)]);
//////				exit(1);
////			}
//		END_FOR
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//				I+source->end_x0-source->start_x0, J, K, u[INDEX(I+source->end_x0-source->start_x0,J,K)], v[INDEX(I+source->end_x0-source->start_x0,J,K)], w[INDEX(I+source->end_x0-source->start_x0,J,K)],
//				phi_c[INDEX(I+source->end_x0-source->start_x0,J,K)], u[INDEX(I+source->end_x0-source->start_x0,J,K)]-u[INDEX(I+source->end_x0-source->start_x0-1,J,K)]
//				            		 +v[INDEX(I+source->end_x0-source->start_x0,J,K)]-v[INDEX(I+source->end_x0-source->start_x0,J-1,K)]
//				            	     +w[INDEX(I+source->end_x0-source->start_x0,J,K)]-w[INDEX(I+source->end_x0-source->start_x0,J,K-1)]);
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//				I+source->end_x0-source->start_x0, J+source->end_y0-source->start_y0, K,
//				u[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K)],
//				v[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K)],
//				w[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K)],
//				phi_c[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K)],
//				u[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K)]-u[INDEX(I+source->end_x0-source->start_x0-1,J+source->end_y0-source->start_y0,K)]
//                +v[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K)]-v[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0-1,K)]
//    	        +w[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K)]-w[INDEX(I+source->end_x0-source->start_x0,J+source->end_y0-source->start_y0,K-1)]);
//		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
//				I, J+source->end_y0-source->start_y0, K,
//				u[INDEX(I,J+source->end_y0-source->start_y0,K)],
//				v[INDEX(I,J+source->end_y0-source->start_y0,K)],
//				w[INDEX(I,J+source->end_y0-source->start_y0,K)],
//				phi_c[INDEX(I,J+source->end_y0-source->start_y0,K)],
//				u[INDEX(I,J+source->end_y0-source->start_y0,K)]-u[INDEX(I-1,J+source->end_y0-source->start_y0,K)]
//                +v[INDEX(I,J+source->end_y0-source->start_y0,K)]-v[INDEX(I,J+source->end_y0-source->start_y0-1,K)]
//    	        +w[INDEX(I,J+source->end_y0-source->start_y0,K)]-w[INDEX(I,J+source->end_y0-source->start_y0,K-1)]);
//		 printf("divergence = %16.13f \n",
//				 u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
//		                   +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
//	                                 +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
//		 if(voxel.InLiquid(I,J,K))
//			 printf(" (%d, %d, %d) in liquid \n", I, J, K);
//		 if(voxel.InAir(I,J,K))
//			printf(" (%d, %d, %d) in air \n", I, J, K);
//		 if(voxel.InSurface(I,J,K))
//			printf(" (%d, %d, %d) in surface \n", I, J, K);
//		 if(voxel.InSolid(I,J,K))
//	 		printf(" (%d, %d, %d) in solid \n", I, J, K);
//		CheckDiv(voxel);
		CheckVelocity(voxel, m);

	}


	void InitializeFromWaveGen(const WaveGenerator *wave, const Voxel *voxel){
		wave->SetInitialWaveVelocity(voxel, DimX, DimY, DimZ,
				WALL_THICKNESS+1, WALL_THICKNESS+1, WALL_THICKNESS+1, u, w);
		SetEqual(u, u1);
		SetEqual(w, w1);
	}

	void FixedVelocity(const Voxel &voxel){
		FOR_EACH_CELL
			Point up = voxel.VelPosition(1,i,j,k);
			u[INDEX(i,j,k)] = (M_PI/314)*(50-up[1]);
			Point vp = voxel.VelPosition(2,i,j,k);
			v[INDEX(i,j,k)] = (M_PI/314)*(vp[0]-50);
			w[INDEX(i,j,k)] = 0.f;
		END_FOR
	}

	void SetEqual(const float *x, float *x0){
		for(u_int i = 0; i < DimX*DimY*DimZ; ++i)
			x0[i]  = x[i];
	}

	void AssignPhi(const Voxel &voxel, float *phi){

		phi_c = phi;
//		u_int N = 0;
//		FOR_EACH_CELL
//			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi_c[INDEX(i,j,k)] < 0.f)
//				++N;
//		END_FOR
//		printf("after AssignPhi(), there are %u liquid points \n", N);
		FOR_EACH_CELL
			if( i < DimX-1 )
				phi_u[INDEX(i,j,k)] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i+1,j,k)]);
			else
				phi_u[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
			if( j < DimY-1 )
				phi_v[INDEX(i,j,k)] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j+1,k)]);
			else
				phi_v[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
			if( k < DimZ-1 )
				phi_w[INDEX(i,j,k)] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j,k+1)]);
			else
				phi_w[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
			if(i==I && j==J && k==K){
		    	printf(" \n before reinit cell phi = (%f, %f, %f) \n\n",
		    			phi_u[INDEX(i,j,k)], phi_v[INDEX(i,j,k)], phi_w[INDEX(i,j,k)]);
			 }
		END_FOR

	}

	void AssignPhiObject(float *phi, int axis){
		if(axis == 0)
			phi_c_obj = phi;
		else if(axis == 1)
			phi_u_obj = phi;
		else if(axis == 2)
			phi_v_obj = phi;
		else if(axis == 3)
			phi_w_obj = phi;


	}

	void RecomputePhiObject(float *phi){

		FOR_EACH_CELL
			if( i < DimX-1 )
				phi_u_obj[INDEX(i,j,k)] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i+1,j,k)]);
			else
				phi_u_obj[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
			if( j < DimY-1 )
				phi_v_obj[INDEX(i,j,k)] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j+1,k)]);
			else
				phi_v_obj[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
			if( k < DimZ-1 )
				phi_w_obj[INDEX(i,j,k)] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j,k+1)]);
			else
				phi_w_obj[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
//			if(i == I && j == J && k == K){
//				printf("phi_c = %f \n", phi_c_obj[INDEX(i,j,k)]);
//				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
//						"phi[k-1] = %f, phi[k+1] = %f \n\n",
//   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
//   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
//   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
//			}

		END_FOR
	}

#ifdef AIR_DIV_FREE
//	void ResetLquidVelocity(){
//		FOR_EACH_CELL
//			if(phi_u[INDEX(i,j,k)] < 0.f && phiu_star[INDEX(i,j,k)] > 0.f)
//				u[INDEX(i,j,k)]= u_star[INDEX(i,j,k)];
//			if(phi_v[INDEX(i,j,k)] < 0.f && phiv_star[INDEX(i,j,k)] > 0.f)
//				v[INDEX(i,j,k)]= v_star[INDEX(i,j,k)];
//			if(phi_w[INDEX(i,j,k)] < 0.f && phiw_star[INDEX(i,j,k)] > 0.f)
//				w[INDEX(i,j,k)]= w_star[INDEX(i,j,k)];
//		END_FOR
//	}
#endif

	void UpdateVoxel(Voxel &voxel) const;
	void CorrectLevelSet(const Voxel &voxel, const Point &pos,
				        float *phiNeg, float* phiPos,
		                const Particle &ps, char sign,
		                unsigned char vel_index) const;
	void ErrorCorrection(const Voxel &voxel,
							list<Particle> &particles
#ifdef POSITIVE_PARTICLES
							,list<Particle> &posParticles
#endif
									);
#ifdef SPMD
	void FindSmoothFactor(const Voxel &voxel, char *needed, char *valid,
			float *value0, float *S);
	void EulerStep(const Voxel &voxel, char *needed, char *valid,
			float *S, float *value, float *value0);
	void ReInitialize(const Voxel &voxel, int vel_index, float *value, float *value0);
	void EulerStepNeeded(const Voxel &voxel, const char *needed, float dtau,
			float *value, float *value0, const float *S);
	void ReInitializeNeeded(Voxel &voxel, const char *needed, float *value, float *value0);
#else
	void ReInitialize(const Voxel &voxel);
#endif
	void ReInitializeObject(const Voxel &voxel);
	void InitializeInterface(const Voxel &voxel, char *mask, float *phi, int index);
	float ClosestDistance(const Voxel &voxel, float *phi,
			const char *valid, int i, int j, int k, int vel_index) const;
	float RefinedClosestDistance(const Voxel &voxel, float *phi,
				int i, int j, int k, int vel_index) const;
	void ProcessNonSolidPoint(const Voxel &voxel, float *phi_tmp,
			                  float *phi, const char *valid,
			                  int i, int j, int k, int vel_index) const;
	void SetValidPoint(const Voxel &voxel, char *valid, int vel_index) const;
	void InitializeInterface( const Voxel &voxel,
						map<u_int, KnownPoint> &pos_band,
						map<u_int, KnownPoint> &neg_band,
						float *phi, const float *phi_obj,
						char *valid, int vel_index);
	void FastMarching(char *mask, const Voxel &voxel, float *phi);

	void FastMarching(map<u_int, KnownPoint> &posband,
						 map<u_int, KnownPoint> &negband,
						 const Voxel &voxel, float *phi,
						 const float *phi_obj, const char *valid);
	float UpdateTentativeValue(const KnownPoint &p, u_int index,
							const char *mask, char MASK,
							 const Voxel &v,
							 float *phi);
	float UpdateTentativeValue(const KnownPoint &p, u_int index,
							 map<u_int, KnownPoint> &band,
							 const Voxel &v,
							 float *phi);
	float ProcessOneTerm(int axis, float phi, const Voxel &v);
	float ProcessTwoTerms(int axis1, int axis2,
							  float phi1, float phi2,
							  const Voxel &v);
	float ProcessThreeTerms(int axis1, int axis2, int axis3,
						    float phi1, float phi2, float phi3,
						    const Voxel &v);

	void Diffuse(const Voxel &voxel);
	void ExplicitDiffusion(const Voxel &voxel, float *vel, float *vel0, float *phi_tmp, int vel_index);
	void ImplicitDiffusion(const Voxel &voxel, float *vel, float *vel0, float *phi_tmp, int vel_index);
	void AddSource(const Voxel &voxel);
	float EstimateFluidVolume(int ii, int jj, int kk, float delta) const;
	void FloodFillProjection(int ii, int jj, int kk, char *color) const;
	void Projection(const Voxel &voxel);
#ifdef AIR_DIV_FREE
	void ProjectionInAirCell(const Voxel &voxel);
#endif
	Vector ObjectNormalAt(float delta, int i, int j, int k) const;
	Vector ObjectNormalAt(const Voxel &voxel, const Point &p) const;
	Vector ObjectNormalAt(const Voxel &voxel, int i, int j, int k) const;
	Vector GeometricNormalAt(const Voxel &voxel, int i, int j, int k) const;
	void RecomputeObjectLevelset(const Voxel &voxel, float *phi);
	void SetBoundary( const Voxel &voxel, int b, float *x );
	void SetSurfaceBoundary(const Voxel &voxel);

#ifdef COUPLED_FLIP
	void SetGridVelocity(int i, int j, int k, int index, float *vel) const;
#endif
	void SetSolidBoundary(const Voxel &voxel, int b, float *x);
	void SetSolidBoundaryHouston2003(const Voxel &voxel, char *valid);
	void SetSolidBoundary(const Voxel &voxel, char *c);
	void SetMovingSolidBoundary(const Voxel &voxel);
	void SetMovingSolidBoundary2(const Voxel &voxel);
	void SetSolidBoundaryForAdvection(const Voxel &voxel,
			map<u_int, KnownPoint> &knownPointU,
			map<u_int, KnownPoint> &knownPointV,
			map<u_int, KnownPoint> &knownPointW);
	void SetSolidBoundaryForAdvectionHouston2003(const Voxel &voxel, const char *valid);
	void SetSolidBoundaryForAdvection(const Voxel &voxel);
	void SetSolidBoundaryNoExtrapolate(const Voxel &voxel, int b, float *x);
	void SetSolidBoundaryNoExtrapolate(const Voxel &voxel);
	void SetAirBoundary(const Voxel &voxel, int b, float *x);
	void GravityForce(const Voxel &voxel){
		SetZero(u0);
		SetZero(v0);
		SetZero(w0);
		FOR_EACH_CELL
		//if(!voxel.InSolid(i,j,k)){
		if(phi_w_obj[INDEX(i,j,k)] < 0.f){
			if(phi_w[INDEX(i,j,k)] <= 0.f)
				w0[INDEX(i,j,k)] = -g;
			else if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
				w0[INDEX(i,j,k)] = -g;
			else if(voxel.InAir(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
				w0[INDEX(i,j,k)] = -g;
			else if(voxel.InSurface(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
				w0[INDEX(i,j,k)] = -g;
		}
		END_FOR

	}

	void VortexForce(const Voxel &voxel, VortexParticles *VP);

	void ClearData(float *data);
	void SetZero(float *data){
		for(int i = 0; i < DimX*DimY*DimZ; ++i)
			data[i]  = 0.f;
	}
	void SetZero(float *data, u_int N){
		for(u_int n=0; n < N; ++n)
			data[n] = 0.f;
	}
	void SetZero(char *data){
		for(int i = 0; i < DimX*DimY*DimZ; ++i)
			data[i]  = 0;
	}

	float TriInterp(const Voxel &voxel, const Point &p, float *phi) const;
	float TriInterp(const Voxel &voxel, const Point &p, float *phi,
					unsigned char index) const;

	void AssignEightVelocities(int i0, int j0, int k0,
			int i1, int j1, int k1, const float *vel, float f[]) const;
	void HandleEightRaysForVelocityInterpolation(const Voxel &voxel,
				int i0, int j0, int k0,
				int i1, int j1, int k1,
				int vel_index, const Point &pi,
			    float f[]) const;
	void TriInterp(const Voxel &voxel, const Point &pi,
				float &ui, float &vi, float &wi,
		 		float inv_delta, float *uvel, float *vvel, float *wvel) const;

	void GetVelocityAtPos(const Voxel &voxel, const Point &p, float delta,
			float &uu, float &vv, float &ww) const;

#if defined(COUPLED_SPH) || defined(COUPLED_FLIP)
	void MomentumConservation(const Voxel &voxel,
//			list<WaterParticle> &absorbed
			WaterParticleList &absorbed
			);
#endif
#ifdef COUPLED_FLIP
	void AccumulateValues(const Point &pos, float invh,
					float vel, int index,
					float *mWeights, float *vc);
	void ConvertParticlesToLevelSet(const Voxel &voxel,
									WaterParticleList &wpl,
									ParticleList &pl,
									int nPerCell);
#endif
	void LinearSolver( const Voxel &voxel, int b, float * x, float * x0, float a, float c );
	void LinearSolver( const Voxel &voxel,
						int b, float * x, float * x0, u_int n );
	void LinearSolver( const Voxel &voxel, float * x,
#ifdef WIN32
					 hash_map<u_int, u_int> &a,
#else
					 unordered_map<u_int, u_int> &a,
#endif
					 const float *x0,
					 const char *border, u_int n );
	void LinearSolverForDiffusion( const Voxel &voxel,
				 float * x, map<u_int, u_int> &a,
				 map<u_int, float> &x0 , u_int n, int vel_index );
#ifdef AIR_DIV_FREE
	void LinearSolverAtAirCell( const Voxel &voxel,
								float * x, map<u_int, u_int> &a,
								 map<u_int, float> &b, char *c,
								 u_int n );
	void SetBoundaryForAirCell(const Voxel &voxel,
		                       float *x);
	void FloodFill(const Voxel &voxel, char *color,
								   int ii, int jj, int kk,
			                       float delta) const;
#endif
	//friend void ScalarField3D::ExtractVelocity(const Voxel &voxel, VectorField3D &vel);
	 void Extrapolate(const Voxel &voxel, const char *valid = NULL);
	 void ExtrapolateOneVelocity(const Voxel &voxel,
			 float *phi_tmp, float *phi_obj, float *vel, int b, float delta);



#ifdef SPMD
	 void ExtrapolateOneVelocityIntoObject(const Voxel &voxel,
	 		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta);
#else
	 void UpdateCells(MinHeap<TentativeGridPoint, float> &heap, int vel_index,
	 			 float *vel, float *phi_tmp);
	 void ExtrapolateOneVelocityIntoObject(const Voxel &voxel,
	  		float *phi_tmp, float *vel, int vel_index, float delta);
#endif

#ifdef SPMD
	 void ExtrapolateOneAirVelocityIntoObject(const Voxel &voxel,
	 		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta);
#else
	 void ExtrapolateOneAirVelocityIntoObject(const Voxel &voxel,
	  		float *phi_tmp, float *vel, int vel_index, float delta);
#endif


#ifdef SPMD
	 void ExtrapolateOneVelocityIntoLiquid(const Voxel &voxel,
	 		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta);
#else
	 void ExtrapolateOneVelocityIntoLiquid(const Voxel &voxel,
	 	 	float *phi_tmp, float *vel, int vel_index, float delta);
#endif
#ifdef SPMD
	 void ExtrapolateOneVelocityIntoNonSolid(const Voxel &voxel,
	 		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta);
#else
	 void ExtrapolateOneVelocityIntoNonSolid(const Voxel &voxel,
	 	 	float *phi_tmp, float *vel, int vel_index, float delta);
#endif


#ifdef SPMD
	 void ExtrapolateOneVelocityIntoObjectForAdvection(const Voxel &voxel,
	 		//map<u_int, KnownPoint> &knownPoint,
	 		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta);
	 void ExtrapolateOneVelocityIntoObjectForAdvectionHouston2003(const Voxel &voxel,
	 		//map<u_int, KnownPoint> &knownPoint,
	 		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta);
	 void ExtrapolateOneVelocityIntoObjectForAdvection1(const Voxel &voxel,
	 	 		//map<u_int, KnownPoint> &knownPoint,
	 	 	float *phi_tmp, const char *valid, float *vel, float *vel0, int vel_index, float delta);
#else
	 void ExtrapolateOneVelocityIntoObjectForAdvection(const Voxel &voxel,
	 	 	//map<u_int, KnownPoint> &knownPoint,
	 	 	float *phi_tmp, float *vel, int vel_index, float delta);
#endif

	 void ExtrapolateOneForProjection(const Voxel &voxel,
                     map<u_int, KnownPoint> &knownPoint,
                    float *phi_tmp, float *phi_obj, float *vel, int vel_index);
#ifdef SPMD
	 void ExtrapolateOneForProjection(const Voxel &voxel,
	 		         map<u_int, KnownPoint> &knownPoint,
	 		         float *phi_tmp, float *phi_obj, float *vel, float *vel0, int vel_index);
	 void ExtrapolateOneForProjection(const Voxel &voxel, float *phi_tmp,
	 						float *vel, float *vel0, int vel_index);
	 float ProcessTrialPoints(const Voxel &voxel, const Point &p,
	 				const float *phi_tmp, char *needed, char *valid, float*vel, int vel_index);
	 void ProcessFarAwayPoints(const Voxel &voxel, float *phi_tmp,
	 		                   char *valid, char *needed, int vel_index,
							 char *mask, float *vel);
	 bool TestVelPointIsLiquid(const Voxel &voxel, int i, int j, int k, int vel_index,
	 		                   const char *valid) const;
	 void DetectIsolatedPoints(const Voxel &voxel,
			 const char *air, const float *phi_tmp, float *vel, int vel_index);

	 void ExtrapolateOneUsingFM(const Voxel &voxel,  const char *air,
	 		                        float *phi_tmp, float *vel, int vel_index);
#endif
	 void ExtrapolateOneForProjection(const Voxel &voxel,
	 		                        float *phi_tmp, float *vel, int vel_index);

#ifdef SPMD
	 u_int HasNonUpdatedPoints(const Voxel &voxel, char *valid, float *phi_tmp, int vel_index,
#ifdef WIN32
							hash_map<u_int, TentativeGridPoint> &heap) const;
#else
							unordered_map<u_int, TentativeGridPoint> &heap) const;
#endif

#endif
	 u_int HasNonUpdatedPoints(const Voxel &voxel, char *valid, float *phi_tmp, int vel_index,
	 			 MinHeap<TentativeGridPoint, float> &heap) const;

	 void StretchVortex(const Voxel &voxel, VortexParticles *vortexParticles);
	 void AdvectRK3(const Voxel &voxel, list<Particle> &particles, int sign) const;
	 void AdvectRK3NoSubcycling(const Voxel &voxel, list<Particle> &particles, int sign) const;
	 bool ClipSemiLagrangianRay(const Point &pold, Point &pos) const;
	 void AdvectRK3Subcycling(const Voxel &voxel, float dts,
//			 list<Particle> &particles,
			 ParticleList &particles,
			 int s) const;
	 void AdvectRK2(const Voxel &voxel, list<Particle> &particles, int sign) const;
	 void AdvectForwardEuler(const Voxel &voxel, list<Particle> &particles, int sign) const;
	 void AdvectForwardEulerSubcycling(const Voxel &voxel, float dts, list<Particle> &particles, int sign) const;
	 void AdjustOneVelocity(const Voxel &voxel, int vel_index,
	 					list<WaterParticle> &waterParticles, float *vel);
	 void AdjustVelocity(const Voxel &voxel, list<WaterParticle> &waterParticles);
	 void Advect(const Voxel &voxel, list<WaterParticle> &particles) const;
#ifdef COUPLED_SPH
	 void AdvectSPHParticles(const Voxel *voxel, float dts, list<WaterParticle> &waterParticles) const;
#endif
	 void AdvectSubcycling(const Voxel &voxel, float dts, WaterParticleList &waterParticles) const;
	 void Advect(const Voxel &voxel, VortexParticles *VP ) const;
#ifdef SPMD
	 void ExtrapolatePhiIntoObject(const Voxel &voxel,
	 							   float *phi_tmp,
	 							    float *phi_tmp0);
	 bool CFLExceedsOne(const Voxel &voxel, float dts, const float *x0, float *cfl);
	 float FindTimeStep(const Voxel &voxel, const float *x0);
	 void SemiLagrangian(const Voxel &voxel, int i, int j, int k,
	 				float dt0, float *d, float *d0);
	 // 3rd order RK and 5th order WENO for levelset
	 void Euler(const Voxel &voxel, float dts, char *needed, float *x, float *x0, char *cfl);
	 void Advect_RK3_WENO(const Voxel &voxel, char *needed, float *x, float *x0);
	 void Advect_FE_WENO(const Voxel &voxel, char *needed, float *x, float *x0);
	 void Advect(const Voxel &voxel, float *x, float *x0);
	 void AdvectSubcycling(const Voxel &voxel, float dts, float *x, float *x0) ;
	 void Advect_RK3_WENO_Subcycling(const Voxel &voxel, float dts, char *needed, float *x, float *x0);
//#else


//	  semi-lagrangian for levelset
	 void AdvectSemiLagrangian(const Voxel &voxel, float *d, float *d0);
#endif
	 void Advect(const Voxel &voxel);
	 void AdvectHouston2003(const Voxel &voxel);
	 bool HasNewExtrema(const Voxel &voxel, int i, int j, int k, int b, float d, float *d0) const;
	 void AdvectMacCormack(const Voxel &voxel);
	 void AdvectOneVelocity(const Voxel &voxel, int b, float *d, const float *d0);



	 void Print(int I, int J, int K){
		 printf("at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
		 			I,J,K, u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)]);
	 }

	 inline u_int  Index(const TentativeGridPoint &p){
	 		return p.kk*DimX*DimY+p.jj*DimX+p.ii;
	 }

	 bool CheckIndex(int i, int j, int k) const;

	 void CheckVelocity(const Voxel &voxel, u_int m) {

		 float delta = voxel.VoxelDelta();

		 float max_u = 0.f, max_v = 0.f, max_w = 0.f;
		 int max_w_index = 0, w_ii = 0, w_jj = 0, w_kk = 0;
		 int u_ii = 0, u_jj = 0, u_kk = 0;
		 int v_ii = 0, v_jj = 0, v_kk = 0;
		 int max_u_index = 0, max_v_index = 0;

		 int mag_ii = 0, mag_jj = 0, mag_kk = 0;
		 float max_vel = 0.f;
		 float max_ud = 0.f, max_vd = 0.f, max_wd = 0.f;
 		 int max_ud_i = 1, max_ud_j = 1, max_ud_k = 1;
 		 int max_vd_i = 1, max_vd_j = 1, max_vd_k = 1;
 		 int max_wd_i = 1, max_wd_j = 1, max_wd_k = 1;
 		 FOR_EACH_CELL
 		 	if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
 		 		 voxel.InSource(i,j,k) ){
 		 		float ud = u1[INDEX(i,j,k)] - u1[INDEX(i-1,j,k)];
 		 		float vd = v1[INDEX(i,j,k)] - v1[INDEX(i,j-1,k)];
 		 		float wd = w1[INDEX(i,j,k)] - w1[INDEX(i,j,k-1)];
 		 		if(fabsf(max_ud) < fabsf(ud)){
 		 			max_ud = ud;
 		 			max_ud_i = i;
 		 			max_ud_j = j;
 		 			max_ud_k = k;
 		 		}
 		 		if(fabsf(max_vd) < fabsf(vd)){
 		 			max_vd = vd;
 		 			max_vd_i = i;
 		 			max_vd_j = j;
 		 			max_vd_k = k;
 		 		}
 		 		if(fabsf(max_wd) < fabsf(wd)){
 		 			max_wd = wd;
 		 			max_wd_i = i;
 		 			max_wd_j = j;
 		 			max_wd_k = k;
 		 		}
 		 	}
 		 END_FOR

 		printf("\n so far max u diff at (%d,%d,%d) u = %f u1 = %f diff = %f phic = %f\n",
				 max_ud_i, max_ud_j, max_ud_k,
				 u1[INDEX(max_ud_i, max_ud_j, max_ud_k)],
				 u1[INDEX(max_ud_i-1, max_ud_j, max_ud_k)],
				 max_ud, phi_c[INDEX(max_ud_i, max_ud_j, max_ud_k)]);
		 printf(" u1[i+1] = %f phiu[i+1] = %f, u[j-1] = %f phiu[j-1] = %f, u[j+1] = %f phiu[j+1] = %f \n",
 				  u1[INDEX(max_ud_i+1, max_ud_j, max_ud_k)], phi_u[INDEX(max_ud_i+1, max_ud_j, max_ud_k)],
 				  u1[INDEX(max_ud_i, max_ud_j-1, max_ud_k)], phi_u[INDEX(max_ud_i, max_ud_j-1, max_ud_k)],
 				  u1[INDEX(max_ud_i, max_ud_j+1, max_ud_k)], phi_u[INDEX(max_ud_i, max_ud_j+1, max_ud_k)]);
		 printf(" u[k+1] = %f phiu[k+1] = %f, u[k-1] = %f phiu[k-1] = %f\n",
  				  u1[INDEX(max_ud_i, max_ud_j, max_ud_k+1)], phi_u[INDEX(max_ud_i, max_ud_j, max_ud_k+1)],
  				  u1[INDEX(max_ud_i, max_ud_j, max_ud_k-1)], phi_u[INDEX(max_ud_i, max_ud_j, max_ud_k-1)]);
		 printf(" so far max v diff at (%d,%d,%d) v = %f v1 = %f diff = %f phic = %f \n",
 				max_vd_i, max_vd_j, max_vd_k,
 				 v1[INDEX(max_vd_i, max_vd_j, max_vd_k)],
 				 v1[INDEX(max_vd_i, max_vd_j-1, max_vd_k)],
 				 max_vd, phi_c[INDEX(max_vd_i, max_vd_j, max_vd_k)]);
		 printf(" v[i+1] = %f phiv[i+1] = %f, v[i-1] = %f phiv[i-1] = %f, v[j+1] = %f phiv[j+1] = %f \n",
  				  v1[INDEX(max_vd_i+1, max_vd_j, max_vd_k)], phi_v[INDEX(max_vd_i+1, max_vd_j, max_vd_k)],
  				  v1[INDEX(max_vd_i-1, max_vd_j, max_vd_k)], phi_v[INDEX(max_vd_i-1, max_vd_j, max_vd_k)],
  				  v1[INDEX(max_vd_i, max_vd_j+1, max_vd_k)], phi_v[INDEX(max_vd_i, max_vd_j+1, max_vd_k)]);
 		 printf(" v[k+1] = %f phiv[k+1] = %f, v[k-1] = %f phiv[k-1] = %f\n",
   				  v1[INDEX(max_vd_i, max_vd_j, max_vd_k+1)], phi_v[INDEX(max_vd_i, max_vd_j, max_vd_k+1)],
   				  v1[INDEX(max_vd_i, max_vd_j, max_vd_k-1)], phi_v[INDEX(max_vd_i, max_vd_j, max_vd_k-1)]);
		 printf(" so far max w diff at (%d,%d,%d) w = %f w1 = %f diff = %f phic = %f\n",
 				 max_wd_i, max_wd_j, max_wd_k,
 				 w1[INDEX(max_wd_i, max_wd_j, max_wd_k)],
 				 w1[INDEX(max_wd_i, max_wd_j, max_wd_k-1)],
 				 max_wd, phi_c[INDEX(max_wd_i, max_wd_j, max_wd_k)]);
		 printf(" w[i+1] = %f phiw[i+1] = %f, w[i-1] = %f phiw[i-1] = %f, w[j+1] = %f phiw[j+1] = %f \n",
  				  w1[INDEX(max_wd_i+1, max_wd_j, max_wd_k)], phi_w[INDEX(max_wd_i+1, max_wd_j, max_wd_k)],
  				  w1[INDEX(max_wd_i-1, max_wd_j, max_wd_k)], phi_w[INDEX(max_wd_i-1, max_wd_j, max_wd_k)],
  				  w[INDEX(max_wd_i, max_wd_j+1, max_wd_k)], phi_w[INDEX(max_wd_i, max_wd_j+1, max_wd_k)]);
 		 printf(" w[k+1] = %f phiw[k+1] = %f, w[j-1] = %f phiw[j-1] = %f\n\n",
   				  w1[INDEX(max_wd_i, max_wd_j, max_wd_k+1)], phi_w[INDEX(max_wd_i, max_wd_j, max_wd_k+1)],
   				  w1[INDEX(max_wd_i, max_wd_j-1, max_wd_k)], phi_w[INDEX(max_wd_i, max_wd_j-1, max_wd_k)]);

		 FOR_EACH_CELL
		 	if(!voxel.InSolid(i,j,k)){
		 		u0[INDEX(i,j,k)] = u1[INDEX(i,j,k)];
		 		v0[INDEX(i,j,k)] = v1[INDEX(i,j,k)];
		 		w0[INDEX(i,j,k)] = w1[INDEX(i,j,k)];
		 	}
		 END_FOR

		 FOR_EACH_CELL
		 	if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
		 		(voxel.InSurface(i,j,k)) ||
		 		(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] < EXTRAPOLATE_VEL_LIMIT) ){
		 		float tmp = sqrt( u0[INDEX(i,j,k)]*u0[INDEX(i,j,k)] +
		 						  v0[INDEX(i,j,k)]*v0[INDEX(i,j,k)] +
		 						  w0[INDEX(i,j,k)]*w0[INDEX(i,j,k)]);
		 		if(max_vel < tmp){
		 			max_vel = tmp;
		 			mag_ii = i;
		 			mag_jj = j;
		 			mag_kk = k;
		 		}
		 	}
		 END_FOR
		 printf("\n max vel magitude at (%d,%d,%d) phic = %f, phi_u = %f, phi_v = %f, phi_w = %f \n",
		  		mag_ii, mag_jj, mag_kk, phi_c[INDEX(mag_ii, mag_jj, mag_kk)],
		  		phi_u[INDEX(mag_ii, mag_jj, mag_kk)],
		  		phi_v[INDEX(mag_ii, mag_jj, mag_kk)],
		  		phi_w[INDEX(mag_ii, mag_jj, mag_kk)]);
		 printf(" max vel magitude = %f at (%d,%d,%d)  u = %f, v = %f, w = %f \n\n",
				 max_vel, mag_ii, mag_jj, mag_kk, u0[INDEX(mag_ii, mag_jj, mag_kk)],
		 		  v0[INDEX(mag_ii, mag_jj, mag_kk)],
		 		  w0[INDEX(mag_ii, mag_jj, mag_kk)]);
		 FOR_EACH_CELL
		 //for(int n=0; n<DimX*DimY*DimZ; ++n){
		   u_int n = INDEX(i,j,k);
		 	if(phi_c[INDEX(i,j,k)] <= 0.f && phi_u_obj[INDEX(i,j,k)] <= 0.f && max_u < fabs(u[n])){
		 		max_u = fabs(u[n]);
		 		max_u_index = n;
		 		 u_ii = i;
		 	     u_jj = j;
		 	     u_kk = k;
		 	}
		 	if(phi_c[INDEX(i,j,k)] <= 0.f && phi_v_obj[INDEX(i,j,k)] <= 0.f && max_v < fabs(v[n])){
		 		 max_v = fabs(v[n]);
		 		 max_v_index = n;
		 		 v_ii = i;
		 	     v_jj = j;
		 	     v_kk = k;
		 	}
		 	if(phi_c[INDEX(i,j,k)] <= 0.f && phi_w_obj[INDEX(i,j,k)] <= 0.f && max_w < fabs(w[n])){
		 		max_w = fabs(w[n]);
		 	    max_w_index = n;
		 	    w_ii = i;
		 	    w_jj = j;
		 	    w_kk = k;
		 	}
		 	if( isinf(phi_c[INDEX(i,j,k)]) || isnan(phi_c[INDEX(i,j,k)]) ||
	 			isinf(phi_u[INDEX(i,j,k)]) || isnan(phi_u[INDEX(i,j,k)]) ||
				isinf(phi_v[INDEX(i,j,k)]) || isnan(phi_v[INDEX(i,j,k)]) ||
                isinf(phi_w[INDEX(i,j,k)]) || isnan(phi_w[INDEX(i,j,k)]) ||
                isinf(u[INDEX(i,j,k)]) || isnan(u[INDEX(i,j,k)])      ||
                isinf(v[INDEX(i,j,k)]) || isnan(v[INDEX(i,j,k)])      ||
                isinf(w[INDEX(i,j,k)]) || isnan(w[INDEX(i,j,k)])
                ){

//	 		if( fabs(phi_c[INDEX(i,j,k)]) > 1.e10 || fabs(phi_c[INDEX(i,j,k)]) > 1.e10 ||
//	 			fabs(phi_u[INDEX(i,j,k)]) > 1.e10 || fabs(phi_u[INDEX(i,j,k)]) > 1.e10 ||
//				fabs(phi_v[INDEX(i,j,k)]) > 1.e10 || fabs(phi_v[INDEX(i,j,k)]) > 1.e10 ||
//				fabs(phi_w[INDEX(i,j,k)]) > 1.e10 || fabs(phi_w[INDEX(i,j,k)]) > 1.e10 ){
                	 printf("\n nan or inf phi at (%d,%d,%d) phi = %f and phiu = %f, phiv = %f, phiw = %f \n",
                			 i,j,k, phi_c[INDEX(i,j,k)], phi_u[INDEX(i,j,k)],
                			 phi_v[INDEX(i,j,k)], phi_w[INDEX(i,j,k)]);
                	 printf(" u = %f, v = %f, w = %f \n",
                			 u[INDEX(i,j,k)], v[INDEX(i,j,k)], w[INDEX(i,j,k)] );
                	 I = i; J = j; K = k;
                	 printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f, divergence = %16.13f \n",
        	 				I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
        	 				phi_c[INDEX(I,J,K)], u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
        	 				            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
        	 				            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);

        	 		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f \n",
        	 						I, J-1, K, u[INDEX(I,J-1,K)], v[INDEX(I,J-1,K)], w[INDEX(I,J-1,K)],
        	 						phi_c[INDEX(I,J-1,K)]);

        	 		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f \n",
        	 						I, J+1, K, u[INDEX(I,J+1,K)], v[INDEX(I,J+1,K)], w[INDEX(I,J+1,K)],
        	 						phi_c[INDEX(I,J+1,K)]);

        	 		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f \n",
        	 					I+1, J, K, u[INDEX(I+1,J,K)], v[INDEX(I+1,J,K)], w[INDEX(I+1,J,K)],
        	 					phi_c[INDEX(I+1,J,K)]);

        	 		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f \n",
        	 					I-1, J, K, u[INDEX(I-1,J,K)], v[INDEX(I-1,J,K)], w[INDEX(I-1,J,K)],
        	 					phi_c[INDEX(I-1,J,K)]);

        	 		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f \n",
        	 					I, J, K+1, u[INDEX(I,J,K+1)], v[INDEX(I,J,K+1)], w[INDEX(I,J,K+1)],
        	 					phi_c[INDEX(I,J,K+1)]);

        	 		printf("vel at (%d, %d, %d) is (%f, %f, %f) with phi = %f \n\n",
        	 					I, J, K-1, u[INDEX(I,J,K-1)], v[INDEX(I,J,K-1)], w[INDEX(I,J,K-1)],
        	 					phi_c[INDEX(I,J,K-1)]);
                	 exit(1);
                }
//		 	 printf("\n i = %d, j = %d, k = %d, u = %f, v = %f, w = %f \n",
//		 			i, j, k, u[n], v[n], w[n] );
		 END_FOR

//		 printf("max_u at %d, max_v at %d ,max_w at %d \n ",
//				 max_u_index, max_v_index, max_w_index);

		 u_max_step = u_max < max_u ? m : u_max_step;
		 u_max = u_max < max_u ? max_u : u_max;
		 v_max_step = v_max < max_v ? m : v_max_step;
		 v_max = v_max < max_v ? max_v : v_max;
		 w_max_step = w_max < max_w ? m : w_max_step;
		 w_max = w_max < max_w ? max_w : w_max;


		 printf("\n max_u = %f, max_v = %f, max_w = %f \n",
			 u[max_u_index], v[max_v_index], w[max_w_index] );

		 int type;
		 if(phi_c_obj[INDEX(u_ii, u_jj, u_kk)] >= 0.f)
			 type = 0;
		 else{
			 if(phi_c[INDEX(u_ii, u_jj, u_kk)] > 0.f)
				 type = 2;
			 else
				 type = 1;
		 }
		 printf("max_u at (%d,%d,%d) type = %d, phic = %f, phi_u = %f, phi_v = %f, phi_w = %f \n\n",
 				 u_ii, u_jj, u_kk, type, phi_c[max_u_index], phi_u[max_u_index],phi_v[max_u_index],phi_w[max_u_index]);
		 if(phi_c_obj[INDEX(v_ii, v_jj, v_kk)] >= 0.f)
 			 type = 0;
 		 else{
 			 if(phi_c[INDEX(v_ii, v_jj, v_kk)] > 0.f)
 				 type = 2;
 			 else
 				 type = 1;
 		 }
		 printf("max_v at (%d,%d,%d) type = %d, phic = %f, phi_u = %f, phi_v = %f, phi_w = %f \n\n",
				 v_ii, v_jj, v_kk, type, phi_c[max_v_index], phi_u[max_v_index],phi_v[max_v_index],phi_w[max_v_index]);
		 if(phi_c_obj[INDEX(w_ii, w_jj, w_kk)] >= 0.f)
 			 type = 0;
 		 else{
 			 if(phi_c[INDEX(w_ii, w_jj, w_kk)] > 0.f)
 				 type = 2;
 			 else
 				 type = 1;
 		 }
		 printf("max_w at (%d,%d,%d) type = %d, phic = %f, phi_u = %f, phi_v = %f, phi_w = %f \n\n",
				 w_ii, w_jj, w_kk, type, phi_c[max_w_index], phi_u[max_w_index],phi_v[max_w_index],phi_w[max_w_index]);
		printf("so far max_u = %f at step = %d, max_v = %f at step = %d, max_w = %f at step = %d \n ",
				u_max, u_max_step, v_max, v_max_step, w_max, w_max_step);

	 }

	 void CheckDiv(const Voxel &voxel) const{
		 float north = 0.f, south = 0.f,  east = 0.f,
		       west = 0.f, top = 0.f, bottom = 0.f;
		 float top1 = 0.f;
		 FOR_EACH_CELL
		 	if(voxel.InLiquid(i,j,k) && k == K){
//		 		if(v[INDEX(i,j,k)] >= 0.f)
		 			north += v[INDEX(i,j,k)];
//		 		else
		 			south += v[INDEX(i,j-1,k)];
//		 		if(u[INDEX(i,j,k)] >= 0.f)
		 			east += u[INDEX(i,j,k)];
//		 		else
		 			west += u[INDEX(i-1,j,k)];
		 		top += w[INDEX(i,j,k)];
		 		bottom += w[INDEX(i,j,k-1)];
		 	}
		 	if(k == K && voxel.InLiquid(i,j,k) && phi_w[INDEX(i,j,k)] <= 0.f)
		 			top1 += w[INDEX(i,j,k)];
		 END_FOR
		 printf("north = %f, south = %f, east = %f, west = %f, top = %f, bottom =%f \n",
				 north, south, east, west, top, bottom);
		 printf("horizontal = %f, vertical = %f top1 = %f\n",
				 north-south+east-west, top, top1);
	 }
	 void OutputGridData(char *filename, int vel_index, int kk) const{
	 	FILE *fp=fopen(filename, "w");
	 	if(!fp){
	       printf("Couldn't open file to write \n");
	       exit(-1);
    	}

    	if(vel_index == 1){
//    		for(int k=0;k<DimZ;k++){
    			for(int j=0;j<DimY;j++){
		 			for(int i=0;i<DimX;i++){
		 				for(int k=0;k<DimZ;k++)
		 					fprintf(fp,"%f ", u[INDEX(i,j,k)]);
		 			fprintf(fp,"\n");
		    	}
    		}
    	}
    	else if(vel_index == 2){
//    		for(int k=0;k<DimZ;k++){
		    	for(int j=0;j<DimY;j++){
		 			for(int i=0;i<DimX;i++)	 {
		 				for(int k=0;k<DimZ;k++)
		 					fprintf(fp,"%f ", v[INDEX(i,j,k)]);
		 			fprintf(fp,"\n");
		 		}
    		}
    	}
    	else if(vel_index == 3){
//    		for(int k=0;k<DimZ;k++){
	    		for(int j=0;j<DimY;j++){
		 			for(int i=0;i<DimX;i++){
		 				for(int k=0;k<DimZ;k++)
		 					fprintf(fp,"%f ", w[INDEX(i,j,k)]);
		 			fprintf(fp,"\n");
		 		}
    		}
    	}
    	else if(vel_index == 4){
//    		for(int k=0;k<DimZ;k++){
		    	for(int j=0;j<DimY;j++){
		 			for(int i=0;i<DimX;i++)	{
		 				for(int k=0;k<DimZ;k++)
		 					fprintf(fp,"%f ", phi_u[INDEX(i,j,k)]);
		 			fprintf(fp,"\n");
		 		}
    		}
    	}
    	else if(vel_index == 5){
//    		for(int k=0;k<DimZ;k++){
		    	for(int j=0;j<DimY;j++){
		 			for(int i=0;i<DimX;i++)	{
		 				for(int k=0;k<DimZ;k++)
		 					fprintf(fp,"%f ", phi_v[INDEX(i,j,k)]);
		 			fprintf(fp,"\n");
		 		}
    		}
    	}
    	else if(vel_index == 6){
//    		for(int k=0;k<DimZ;k++){
		    	for(int j=0;j<DimY;j++){
		 			for(int i=0;i<DimX;i++)	{
		 				for(int k=0;k<DimZ;k++)
		 					fprintf(fp,"%f ", phi_w[INDEX(i,j,k)]);
		 			fprintf(fp,"\n");
		 		}
    		}
    	}
    	else if(vel_index == 7){
//    		for(int k=0;k<DimZ;k++){
		    	for(int j=0;j<DimY;j++){
		 			for(int i=0;i<DimX;i++){
		 				for(int k=0;k<DimZ;k++)
		 					fprintf(fp,"%f ", phi_c[INDEX(i,j,k)]);
		 			fprintf(fp,"\n");
		 		}
    		}
    	}
    	else{
    		printf("Error! vel_index must be in [1...7]\n");
    		exit(1);
    	}
	 	fclose(fp);

	 }
	 void OutputBindaryData(const Voxel &voxel, char *filename) const;

	 void WriteRestart(FILE* fp) const;

	 bool ReadRestart(FILE* fp);

	 void SourceTerm(const Voxel &voxel,
	 				int start_x0, int start_y0,  int start_z0,
 					int end_x0,   int end_y0,   int  end_z0,
 					int nx, int ny, int nz,
 					float speed) {

		 source = new SourceRegion( start_x0, start_y0,  start_z0,
				 					end_x0,  end_y0,   end_z0,
				 					nx, ny, nz,
				 					speed);

 		int i, j, k;
// 		for(k = start_z0; k <= end_z0; ++k)
 			k= start_z0;
 			for(j = start_y0; j <= end_y0; ++j)
 				for(i = start_x0; i <= end_x0; ++i)
 					source->AddSource(i,j,k);
   	    k = end_z0;
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i)
 				source->AddSourceNeighbor(i,j,k);

// 		printf("(SourceTerm) vel at (%d, %d, %d) is (%f, %f, %f) with divergence = %16.13f \n",
// 						I, J, K, u[INDEX(I,J,K)], v[INDEX(I,J,K)], w[INDEX(I,J,K)],
// 						u[INDEX(I,J,K)]-u[INDEX(I-1,J,K)]
// 						            		 +v[INDEX(I,J,K)]-v[INDEX(I,J-1,K)]
// 						            	     +w[INDEX(I,J,K)]-w[INDEX(I,J,K-1)]);
 	}

	void ApplySourceTerm(const Voxel &voxel){
		for(int n = 0; n < sources.size();++n){
			sources[n]->SetVelocity(voxel, u, v, w);
		}

	}

	void AddMovingObject(MovingObject *obj){
		movingObjects.push_back(obj);
	}

	bool HasMovingObject(){
		bool stillMoving = false;
		if(hasMovingObj){
			for(int i = 0; i < movingObjects.size(); ++i){
				MovingObject *obj = movingObjects[i];
				if(!obj->StopNow()){
					stillMoving = true;
					break;
				}
			}
		}
		return stillMoving;
	}

#define EXTENT 5

	bool NearInterface(int ii, int jj, int kk, float e, float limit, float *tmp) const{

		for(int i = ii-EXTENT; i <= ii+EXTENT; ++i)
			for(int j = jj-EXTENT; j <= jj+EXTENT; ++j)
				for(int k = kk-EXTENT; k <= kk+EXTENT; ++k){
					if(fabsf(tmp[INDEX(i,j,k)])+e < limit)
						return true;
				}
		return false;

	}

	void UpdateObjectBoundary(Voxel &voxel, bool initial){
		float delta = voxel.VoxelDelta();
		if(hasMovingObj || initial){
			bool stillMoving = false;
			for(int i = 0; i < movingObjects.size(); ++i){
				MovingObject *obj = movingObjects[i];
				if(!obj->StopNow()){
					stillMoving = true;
					break;
				}
			}
			if(initial || stillMoving){
			float *obj_b = temp;
			SetZero(needRecomputePhi);
			SetEqual(phi_c_obj, obj_b);
			if(hasMovingObj){
			// at every time step, reset phi_c_obj to the static object levelset
				SetEqual(phi_c_fixed_obj, phi_c_obj);
				SetEqual(phi_u_fixed_obj, phi_u_obj);
				SetEqual(phi_v_fixed_obj, phi_v_obj);
				SetEqual(phi_w_fixed_obj, phi_w_obj);
			}
	//		FOR_EACH_CELL
	//		if(i == I && j == J && k == K){
	//				printf("phi_c = %f \n", phi_c_obj[INDEX(i,j,k)]);
	//				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
	//						"phi[k-1] = %f, phi[k+1] = %f \n\n",
	//   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
	//   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
	//   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
	//			}
	//		END_FOR
			for(int i = 0; i < movingObjects.size(); ++i){
				MovingObject *obj = movingObjects[i];
				obj->Update(voxel, phi_c_obj, phi_u_obj, phi_v_obj, phi_w_obj);
			}


//			for(int i = 0; i < movingObjects.size(); ++i){
//				MovingObject *obj = movingObjects[i];
//				obj->Intersection(phi_c_obj);
//			}

			FOR_EACH_CELL
			if(i == I && j == J && k == K){
					printf("phi_c = %f \n", phi_c_obj[INDEX(i,j,k)]);
					printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
							"phi[k-1] = %f, phi[k+1] = %f \n\n",
	   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
	   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
	   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
					printf("phi_u_obj = %f, phi_v_obj = %f,phi_w_obj = %f\n ",
							phi_u_obj[INDEX(i,j,k)], phi_v_obj[INDEX(i,j,k)],phi_w_obj[INDEX(i,j,k)]);
				}
			END_FOR
			u_int NM = 0;
			FOR_EACH_CELL
				if(voxel.InSolid(i,j,k)){
					++NM;
				}
//				if(phi_c != NULL){
//					if(!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)] <= 0.f)
//						u0[INDEX(i,j,k)] = -999.f;
//				}
//				else{
//					if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
//						u0[INDEX(i,j,k)] = -999.f;
//				}
			END_FOR
			printf("\nbefore moving object update there are %ld solid points \n\n", NM);
			u_int n1 = 0, n2 = 0, n3 = 0, n4 = 0;
			FOR_EACH_CELL
				if(phi_c_obj[INDEX(i,j,k)] > 0.f){  // within moving object boundary
					if(!voxel.InSolid(i,j,k))
						voxel.UpdateCellType(i,j,k,SET,SOLID);
					if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
						++n1;
					}
					if(voxel.InLiquid(i,j,k)){
						voxel.UpdateCellType(i,j,k,CLEAR,LIQUID);
					}
					if(voxel.InAir(i,j,k)){
						voxel.UpdateCellType(i,j,k,CLEAR,EMPTY);
						++n4;
					}
					if(voxel.InSurface(i,j,k)){
						voxel.UpdateCellType(i,j,k,CLEAR,SURFACE);
						++n4;
					}

				}
				else{
					if(voxel.InSolid(i,j,k)){
						voxel.UpdateCellType(i,j,k,CLEAR,SOLID);
//						if(phi_c[INDEX(i,j,k)] <= 0.f){
//							voxel.UpdateCellType(i,j,k,SET,LIQUID);
//							++n2;
//						}
//						else{
						voxel.UpdateCellType(i,j,k,SET,EMPTY);
						phi_c[INDEX(i,j,k)] = phi_c[INDEX(i,j,k)] > phi_c_obj[INDEX(i,j,k)] ?
									phi_c[INDEX(i,j,k)] : phi_c_obj[INDEX(i,j,k)] ;
						++n3;
						needRecomputePhi[INDEX(i,j,k)] = 1;
//						}
					}
				}
			END_FOR
			printf("\n after moving object %u liquid -> solid, and %u solid -> liquid, and %u solid -> air, %u air -> solid\n\n", n1, n2, n3, n4);
	//		AssignPhiObject(phi_c_obj);
			NM = 0;
			int NP = 0, NN = 0;
			FOR_EACH_CELL
				if(voxel.InSolid(i,j,k)){
						++NM;
				}
				if(voxel.InSolid(i,j,k) && u0[INDEX(i,j,k)] == -999.f){
					//NP++;
//					printf(" liquid cell engalved by the moving object at (%d, %d, %d) \n", i,j,k);
				}
				if(obj_b[INDEX(i,j,k)] < 0.f && phi_c_obj[INDEX(i,j,k)] > 0.f)
					NP++;
				if(obj_b[INDEX(i,j,k)] > 0.f && phi_c_obj[INDEX(i,j,k)] < 0.f)
					NN++;
				if(i == I && j == J && k == K){
					if(voxel.InSolid(i,j,k))
						printf(" grid point (%d, %d, %d)in solid \n ",i,j,k);
					else{
						printf(" grid point (%d, %d, %d) not in solid \n ",i,j,k);
						if(voxel.InLiquid(I,J,K))
						  printf(" \n (%d, %d, %d) in liquid \n\n", I, J, K);
					    if(voxel.InAir(I,J,K))
						  printf(" \n (%d, %d, %d) in air \n\n", I, J, K);
					    if(voxel.InSurface(I,J,K))
						  printf("\n  (%d, %d, %d) in surface \n\n", I, J, K);
					}
				}
			END_FOR
			printf("\nafter moving object update there are %ld solid points \n\n", NM);
			printf("\nafter moving object update there are %d points nonsolid->solid, %d points solid->nonsolid \n\n",
					NP, NN);
			//ReInitializeObjectPhic(voxel);
//			RecomputePhiObject(phi_c_obj);
//			ReInitializeObject(voxel);
//			delete [] obj_b;
#ifdef SPMD
			if(!initial){
				NN = 0;
				FOR_EACH_CELL
					if(voxel.InMovingSolid(i,j,k)){
//						phi_c[INDEX(i,j,k)] = phi_c[INDEX(i,j,k)] > phi_c_obj[INDEX(i,j,k)] ?
//								phi_c[INDEX(i,j,k)] : phi_c_obj[INDEX(i,j,k)] ;
						phi_c[INDEX(i,j,k)] = phi_c_obj[INDEX(i,j,k)];
					}
//					if(!voxel.InSolid(i,j,k) && fabsf(phi_c[INDEX(i,j,k)])+E_EPSIL < FASTMARCH_LIMIT){
					if( !voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)
							&& NearInterface(i, j, k, E_EPSIL, FASTMARCH_LIMIT, phi_c) ){
						needRecomputePhi[INDEX(i,j,k)] = 1;
						++NN;
					}
				END_FOR
				SetEqual(phi_c, w0);
				SetEqual(phi_c, temp);
				printf("\n Reintializing %u points including those (%u) which are just occupied by moving solids\n\n", NN, n3);
				ReInitializeNeeded(voxel, needRecomputePhi, w0, temp);
				FOR_EACH_CELL
					if(!voxel.InSolid(i,j,k) && fabs(w0[INDEX(i,j,k)]) > FASTMARCH_LIMIT){
						w0[INDEX(i,j,k)] = w0[INDEX(i,j,k)] > 0.f ? FASTMARCH_LIMIT : -FASTMARCH_LIMIT;
//						++N1;
					}
					else if(!voxel.InSolid(i,j,k) && FASTMARCH_LIMIT - fabs(w0[INDEX(i,j,k)]) < delta){
						w0[INDEX(i,j,k)] = w0[INDEX(i,j,k)] > 0.f ? FASTMARCH_LIMIT : -FASTMARCH_LIMIT;
//						++N1;
					}
				END_FOR
				FOR_EACH_CELL
//					if(needRecomputePhi[INDEX(i,j,k)])
						phi_c[INDEX(i,j,k)] = w0[INDEX(i,j,k)];
				END_FOR
			}
#endif

			}
		}
	}

	void PrintMaxVelocity(){
		printf("MAX u = %f, v = %f, w = %f \n",
				u_max, v_max, w_max);
	}

	// this function only gets called once at the beginning of the simulation
	// subsequently, phi_c_obj is reset every time step based on
	// values set here
	void FixedObjectBoundary(){
		if(hasMovingObj){
			SetEqual(phi_c_obj, phi_c_fixed_obj);
			SetEqual(phi_u_obj, phi_u_fixed_obj);
			SetEqual(phi_v_obj, phi_v_fixed_obj);
			SetEqual(phi_w_obj, phi_w_fixed_obj);
		}
	}

	void SetConstraintVelocity(const Voxel &voxel, Vector &obj_vel, int i, int j, int k) const{
		for(int n = 0; n < movingObjects.size(); ++n){
			MovingObject *obj = movingObjects[n];
			obj->SetConstraintVelocity(voxel, obj_vel, i, j, k);
		}
	}
	void SetConstraintVelocity(const Voxel &voxel, Vector &obj_vel, int i, int j, int k, int vel_index) const{
		for(int n = 0; n < movingObjects.size(); ++n){
			MovingObject *obj = movingObjects[n];
			obj->SetConstraintVelocity(voxel, obj_vel, i, j, k, vel_index);
		}
	}
	void ConstraintFriction(const Voxel &voxel,  int vel_index, int i, int j, int k, float *f) const{
		for(int n = 0; n < movingObjects.size(); ++n){
			MovingObject *obj = movingObjects[n];
			obj->Friction(voxel, vel_index, i, j, k, f);
		}
	}

	bool IsDynamicBorderCell(const Voxel &voxel, int i, int j, int k) const{
		bool is = false;
		for(int n = 0; n < movingObjects.size(); ++n){
			MovingObject *obj = movingObjects[n];
			is = obj->IsDynamicBorderCell(voxel, i, j, k);
		}
		return is;
	}


	void NumofCornerInSolid(float *f, int &n, int *s) const{
		for(int l=0; l < 4; ++l){
			if( f[l] - E_EPSIL > 0.f){
				++n;
				s[l] = l;
			}
		}
	}

	float OneCornerInSolid(float *f, int *s) const{
		int one = 0;
		for(int l=0; l < 4; ++l){
			if(s[l] != -99){
				one = l;
				break;
			}
		}
		float f1, f2;
		if(one == 0){
			f1 = f[0]/(f[0]-f[1]);
			f2 = f[0]/(f[0]-f[3]);
		}
		else if(one == 1 || one == 2){
			f1 = f[one]/(f[one]-f[one-1]);
			f2 = f[one]/(f[one]-f[one+1]);
		}
		else if(one == 3){
			f1 = f[3]/(f[3]-f[2]);
			f2 = f[3]/(f[3]-f[0]);
		}
		return 1.f - 0.5f*f1*f2;
	}

	float TwoCornersInSolid(float *f, int *s) const{
		float f1, f2;
		if(s[0] != -99 && s[1] != -99){
			f1 = f[0]/(f[0]-f[3]);
			f2 = f[1]/(f[1]-f[2]);
		}
		else if(s[1] != -99 && s[2] != -99){
			f1 = f[1]/(f[1]-f[0]);
			f2 = f[2]/(f[2]-f[3]);
		}
		else if(s[2] != -99 && s[3] != -99){
			f1 = f[2]/(f[2]-f[1]);
			f2 = f[3]/(f[3]-f[0]);
		}
		else if(s[3] != -99 && s[0] != -99){
			f1 = f[3]/(f[3]-f[2]);
			f2 = f[0]/(f[0]-f[1]);
		}
		return 1.f - 0.5f*(f1 + f2);
	}

	float ThreeCornersInSolid(float *f, int *s) const{
		int one = 0;
		for(int l=0; l < 4; ++l){
			if(s[l] == -99){
				one = l;
				break;
			}
		}
		float f1, f2;
		if(one == 0){
			f1 = f[0]/(f[0]-f[1]);
			f2 = f[0]/(f[0]-f[3]);
		}
		else if(one == 1 || one == 2){
			f1 = f[one]/(f[one]-f[one-1]);
			f2 = f[one]/(f[one]-f[one+1]);
		}
		else if(one == 3){
			f1 = f[3]/(f[3]-f[2]);
			f2 = f[3]/(f[3]-f[0]);
		}
		return 0.5f*f1*f2;
	}

	void CutCellFaceArea(const Voxel &voxel, int i, int j, int k, int m, float *c) const;
	float MarchingCubeArea(const Voxel &voxel, int i, int j, int k) const;

	void UpdateNonWaterVelocity(const Voxel &voxel){
		char *maskAirVelSet = new char[DimX*DimY*DimZ];
		SetZero(maskAirVelSet);
#ifdef SPMD
//		ExtrapolateOneForProjection(voxel, knownPointUAir, phi_u, phi_u_obj, u, u0, 1);
//		ExtrapolateOneForProjection(voxel, knownPointVAir, phi_v, phi_v_obj, v, v0, 2);
//		ExtrapolateOneForProjection(voxel, knownPointWAir, phi_w, phi_w_obj, w, w0, 3);
//		Extrapolate(voxel, valid);
//#else
		ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_u, u, 1);
		ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_v, v, 2);
		ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_w, w, 3);
//		ExtrapolateOneForProjection(voxel, knownPointUAir, phi_u, phi_u_obj, u, 1);
//		ExtrapolateOneForProjection(voxel, knownPointVAir, phi_v, phi_v_obj, v, 2);
//		ExtrapolateOneForProjection(voxel, knownPointWAir, phi_w, phi_w_obj, w, 3);
#endif
		delete [] maskAirVelSet;
//#ifdef AIR_DIV_FREE
//	ApplySourceTerm(voxel);
//	SetEqual(u, u1);
//	SetEqual(v, v1);
//	SetEqual(w, w1);
//	FOR_EACH_CELL
//		if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
//			u1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,1);
//			}
//			if(obj_vel.x != 0.f)
//				u1[INDEX(i,j,k)] = obj_vel.x;
//		}
//		else if(i < DimX-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i+1,j,k) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i+1,j,k)) )
//			u1[INDEX(i,j,k)] = 0.f;
//
//		if(phi_v_obj[INDEX(i,j,k)] >= 0.f){
//			v1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel = Vector(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,2);
//			}
//			if(obj_vel.y != 0.f)
//				v1[INDEX(i,j,k)] = obj_vel.y;
//		}
//		else if(j < DimY-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i,j+1,k) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i,j+1,k)) )
//			v1[INDEX(i,j,k)] = 0.f;
//
//		if(phi_w_obj[INDEX(i,j,k)] >= 0.f){
//			w1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel = Vector(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,3);
//			}
//			if(obj_vel.z != 0.f)
//				w1[INDEX(i,j,k)] = obj_vel.z;
//		}
//		else if( k < DimZ-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i,j,k+1) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i,j,k+1)) )
//			w1[INDEX(i,j,k)] = 0.f;
//	END_FOR
//
//	ProjectionInAirCell(voxel);
//
//#endif
//		SetSolidBoundaryForAdvection(voxel, knownPointU, knownPointV, knownPointW);
		SetSolidBoundaryForAdvection(voxel);
		for(int n = 0; n < sources.size();++n)
			sources[n]->SetVelocity(voxel, u1, v1, w1);

		printf("after UpdateNonWaterVelocity point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
					I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );
		printf("after UpdateNonWaterVelocity point (%d, %d, %d)  u = %f, v = %f, w = %f \n",
				I, J, K, u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)] );

	}

	void UpdateTimeStep() {
		dt = tmg->GetDt();
	}

#ifdef COUPLED_SPH
	void ProjectionSPH(const Voxel *voxel, SPH *sph) const{
		sph->Projection(voxel, u, v, w);
	}
#endif

private:
//	int DimX;
//	int DimY;
//	int DimZ;
	float *u, *v, *w;
	float *u0, *v0, *w0;
	float *u1, *v1, *w1;  // used in ErrorCorrection
	float *phi_u, *phi_v, *phi_w;
	float *phi_c;
	float *phi_c_obj, *phi_c_fixed_obj;
	float *phi_u_obj, *phi_v_obj, *phi_w_obj;
	float *phi_u_fixed_obj, *phi_v_fixed_obj, *phi_w_fixed_obj;
	char *needRecomputePhi;
	float *temp;
#ifdef AIR_DIV_FREE
//	float *phiu_star, *phiv_star, *phiw_star;
//	float *u_star, *v_star, *w_star;
#endif
	float dt;
	TimeManager *tmg;
	float g;
	float visc;
	bool hasMovingObj;
	Container *container;
	int I, J, K;
	float u_max, v_max, w_max;
	int   u_max_step, v_max_step, w_max_step;
#ifdef AIR_DIV_FREE
	double residual_projection_air;
#endif
	SourceRegion *source;
	vector<MovingObject *> movingObjects;
	CIsoSurface<float> *mc;
	const WaveGenerator *mWaveGenerator;

	VectorField3D(const VectorField3D &a)
	{
		assert(false);
	}

	VectorField3D& operator=(const VectorField3D &a)
	{ assert(false); return *this; }

};

#endif /*VECTORFIELD3D_H_*/

