/*
 * flip.cpp
 *
 *  Created on: Oct 29, 2010
 *      Author: bzhao
 */

#include <queue>
using std::queue;

#include "flip.h"
#include "mtrand.h"
#include "pcg_solver.h"
#include "scalarfield3d.h"


#define TARGET_CELL if(I==i && J==j && K==k)

#ifdef COUPLED_FLIP

void FLIP::AddWaterParticles(const WaterParticle &p){
	mWaterParticles.push_back(p);
}

void FLIP::TimeStepping(const Voxel *voxel, const VectorField3D *vel){
	float dt = mTimeManager->GetDt();
	printf(" There are %u FLIP particles \n", mWaterParticles.size());
	if(mWaterParticles.size() == 0)
		return;
	else{
		MoveParticlesInGrid(voxel);
		TransferParticleValuesToGrid(voxel);
		SaveGridVelocities();
		ApplyForce(dt);
		ComputeGridLevelset();
		ExtrapolateGridVelocity();
		ApplyBoundaryConditions(voxel, vel);
		Projection(voxel, dt);
		ExtrapolateGridVelocity();
		UpdateGridVelocityIncrement();
		UpdateParticleVelocityFromGridValues();
	}
}

void FLIP::Projection(const Voxel *voxel, float dt){

#ifdef WIN32
	hash_map<u_int, float> rhs;
	hash_map<u_int, u_int> A;
	hash_map<u_int, u_int>::iterator found;
#else
	unordered_map<u_int, float> rhs;
	unordered_map<u_int, u_int> A;
	unordered_map<u_int, u_int>::iterator found;
#endif


	char *color = new char[DimX*DimY*DimZ];
	SetZero(color);
	FOR_EACH_CELL
		if(voxel->InSolid(i,j,k))
			color[INDEX(i,j,k)] = 2;
		if(voxel->InSource(i,j,k))
			color[INDEX(i,j,k)] = 2;
		if(voxel->InLiquid(i,j,k))
			color[INDEX(i,j,k)] = 2;
	END_FOR
	TentativeGridPoint floodSeed(0.f, DimX/2, DimY/2, DimZ-(WALL_THICKNESS+1)-1);
	FloodFill(floodSeed.ii, floodSeed.jj, floodSeed.kk, color);
	u_int air_pockets = 0, air_points = 0, others = 0;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(color[pos] == 0){
			++air_pockets;
			if(mMarker[pos] == 1)
				mMarker[pos] = 0;
			printf("at (%d, %d, %d) grid point is in air pocket\n", i,j,k);

		}
		if(color[pos] == 1)
			++air_points;
		if(color[pos] == 2)
			++others;
	END_FOR
	printf("\nthere are %u points in air pocket and %u points in air \n\n", air_pockets, air_points);
	delete [] color;
//	CheckNewmannBCSurroundedRegion(voxel);

	u_int N = 0;
	SetZero(tmp);
	float dxdt = h / dt;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mMarker[pos] == 1){
			tmp[pos] = -(u[pos] - u[INDEX(i-1,j,k)] +
						 v[pos] - v[INDEX(i,j-1,k)] +
						 w[pos] - w[INDEX(i,j,k-1)])*dxdt;
			rhs.insert( make_pair(pos, tmp[pos]) );
			A.insert( make_pair(pos, N) );
			++N;
//			TARGET_CELL
//			{
//				printf("Before Projection: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f,"
//						" u_1 = %f, v_1 = %f, w_1 = %f, div = %f\n",
//							i,j,k, u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)],
//							u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)], tmp[INDEX(i,j,k)]);
//				printf("Before Projection: marker = %d, "
//						" %d, %d,  %d, %d, %d, %d \n",
//								mMarker[INDEX(i,j,k)], mMarker[INDEX(i+1,j,k)],mMarker[INDEX(i,j+1,k)],
//								mMarker[INDEX(i,j,k+1)], mMarker[INDEX(i-1,j,k)],mMarker[INDEX(i,j-1,k)],
//								mMarker[INDEX(i,j,k-1)]);
//			}
		}
		TARGET_CELL
			PrintOneCellVel(i,j,k, "before projection");
	END_FOR
	if(N == 0)
		return;
	float *p = new float[N];
	SetZero(p, N);
	printf("FLIP particles occupying %u grid cells \n", N);
	LinearSolver(voxel, p, A, rhs, N);
	SetZero(tmp);
	found = A.begin();
	for(; found != A.end(); ++found){
		tmp[found->first] = p[found->second];
	}
	delete [] p;

	float dtdx = dt * invh;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//		TARGET_CELL
//		printf(" p[i,j,k] = %f, p[i-1,j,k]= %f, p[i,j-1,k]= %f, p[i,j,k-1]= %f"
//				"               p[i+1,j,k]= %f, p[i,j+1,k]= %f, p[i,j,k+1]= %f \n",
//				tmp[INDEX(i,j,k)], tmp[INDEX(i-1,j,k)],tmp[INDEX(i,j-1,k)],tmp[INDEX(i,j,k-1)],
//				tmp[INDEX(i+1,j,k)],tmp[INDEX(i,j+1,k)],tmp[INDEX(i,j,k+1)]);
		if(i < DimX-1){
			if(mMarker[pos] == 1 || mMarker[INDEX(i+1,j,k)] == 1)
				if(!voxel->InLiquid(i,j,k)   && !voxel->InSource(i,j,k)   && !voxel->InSolid(i,j,k) &&
				   !voxel->InLiquid(i+1,j,k) && !voxel->InSource(i+1,j,k) && !voxel->InSolid(i+1,j,k)){
					u[pos] -= dtdx * (tmp[INDEX(i+1,j,k)] - tmp[pos]);
//					TARGET_CELL
//						printf("A. p1 = %f, p = %f \n",tmp[INDEX(i+1,j,k)], tmp[pos]);
				}
		}
		if(j < DimY-1){
			if(mMarker[pos] == 1 || mMarker[INDEX(i,j+1,k)] == 1)
				if(!voxel->InLiquid(i,j,k)   && !voxel->InSource(i,j,k)   && !voxel->InSolid(i,j,k) &&
				   !voxel->InLiquid(i,j+1,k) && !voxel->InSource(i,j+1,k) && !voxel->InSolid(i,j+1,k)){
					v[pos] -= dtdx * (tmp[INDEX(i,j+1,k)] - tmp[pos]);
//					TARGET_CELL
//						printf("B. p1 = %f, p = %f \n",tmp[INDEX(i,j+1,k)], tmp[pos]);
				}
		}
		if(k < DimZ-1){
			if(mMarker[pos] == 1 || mMarker[INDEX(i,j,k+1)] == 1)
				if(!voxel->InLiquid(i,j,k)   && !voxel->InSource(i,j,k)   && !voxel->InSolid(i,j,k) &&
				   !voxel->InLiquid(i,j,k+1) && !voxel->InSource(i,j,k+1) && !voxel->InSolid(i,j,k+1)){
					w[pos] -= dtdx * (tmp[INDEX(i,j,k+1)] - tmp[pos]);
//					TARGET_CELL
//						printf("C. p1 = %f, p = %f \n",tmp[INDEX(i,j,k+1)], tmp[pos]);
				}
		}
		TARGET_CELL
			PrintOneCellVel(i,j,k, "after projection");
	END_FOR

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mMarker[pos]){
				float div =  u[pos]-u[INDEX(i-1,j,k)]
		  				    +v[pos]-v[INDEX(i,j-1,k)]
						    +w[pos]-w[INDEX(i,j,k-1)];
				if(fabs(div) > 1.e-3){
					printf("(%d,%d,%d) phi = %f, u = %f, v = %f, w = %f,"
							" u_1 = %f, v_1 = %f, w_1 = %f div = %f\n",
								i, j, k, phi[INDEX(i,j,k)],
						u[INDEX(i,j,k)], v[INDEX(i,j,k)], w[INDEX(i,j,k)],
						u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)], div);
					exit(1);
				}
		}
	END_FOR

	CheckVelocity();

}


void FLIP::FloodFill(int ii, int jj, int kk, char *color) const {
	queue<TentativeGridPoint> Q;
	TentativeGridPoint it(0.f, ii, jj, kk);
	u_int index = INDEX(ii,jj,kk);
//	printf("current position (%d, %d, %d) color = %d \n", ii, jj, kk,color[index]);
	if(color[index] == 2){
		printf("seed point in solid\n");
		return;
	}
//	else if(found != visitedPoints.end())
	else if(color[index] == 1){
		printf("seed point already set \n");
		return;
	}
	else{
		// put the seed point into queue to start the loop
		Q.push(it);
		color[index] = 1;
		do{
			// get the current point
			TentativeGridPoint tp = Q.front();
			vector<TentativeGridPoint> neighbors;
			tp.Neigbor(neighbors, DimX, DimY, DimZ);
			for(int m=0; m<neighbors.size(); m++){
				// if a neighbor point has not been set to 1, set it to 1 and put it into queue
				// these neighbors are immediately set to 1 here such that they will not be
				// duplicated in queue later on
				int i = neighbors[m].ii;
				int j = neighbors[m].jj;
				int k = neighbors[m].kk;
				if(color[INDEX(i, j, k)] == 0){
				   color[INDEX(i, j, k)] = 1;
				   TentativeGridPoint tmp(0.f, i, j, k);
				   Q.push(tmp);
				}
			}
			Q.pop();
//			printf("there are %lu iterms in queues \n", Q.size());
		}while(!Q.empty());
	}
}

void FLIP::CheckNewmannBCSurroundedRegion(const Voxel *voxel){
	FOR_EACH_CELL
	   u_int index = INDEX(i,j,k);
	   if( mMarker[index] == 1 ){
		 TentativeGridPoint tp(0.f,i,j,k);
		 vector<TentativeGridPoint> neighbors;
		 tp.Neigbor(neighbors, DimX, DimY, DimZ);
		 int solidNeighbors=0;
		 for(u_int m=0; m<neighbors.size(); m++){
			 if(voxel->InSolid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ||
				voxel->InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
				++solidNeighbors;
			 }
			 else if(voxel->InSource(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
				 ++solidNeighbors;
			 }
		 }
		 if(solidNeighbors == 6){
			printf("at (%d, %d, %d) Newmann b.c. for all 6 directions \n", i,j,k);
//			I = i; J = j; K = k;
			mMarker[index] = 0;
		 }
	 }
	END_FOR
}

void FLIP::CheckVelocity() const{

		 int mag_ii = 0, mag_jj = 0, mag_kk = 0;
		 float max_vel = 0.f;

		 FOR_EACH_CELL
			 u_int pos = INDEX(i,j,k);
//			 if(mMarker[pos]){
		  		float tmp = sqrt( u[pos]*u[pos] +
		 						  v[pos]*v[pos] +
		 						  w[pos]*w[pos]);
		 		if(max_vel < tmp){
		 			max_vel = tmp;
		 			mag_ii = i;
		 			mag_jj = j;
		 			mag_kk = k;
		 		}
//			 }
		 END_FOR

		 printf("\n max vel magitude at (%d,%d,%d) \n",
		  		mag_ii, mag_jj, mag_kk );
		 printf(" max vel magitude = %f at (%d,%d,%d)  u = %f, v = %f, w = %f \n\n",
				 max_vel, mag_ii, mag_jj, mag_kk, u[INDEX(mag_ii, mag_jj, mag_kk)],
		 		  v[INDEX(mag_ii, mag_jj, mag_kk)],
		 		  w[INDEX(mag_ii, mag_jj, mag_kk)]);
		 if(mag_ii > 0 && mag_jj > 0 && mag_kk > 0)
			 printf("u0 = %f, v0 = %f, w0 = %f \n\n",
					  u[INDEX(mag_ii-1, mag_jj, mag_kk)],
					  v[INDEX(mag_ii, mag_jj-1, mag_kk)],
					  w[INDEX(mag_ii, mag_jj, mag_kk-1)]);


	 }

void FLIP::LinearSolver(const Voxel *voxel, float *x,
#ifdef WIN32
					 hash_map<u_int, u_int> &a, hash_map<u_int, float> &x0,
#else
					 unordered_map<u_int, u_int> &a, unordered_map<u_int, float> &x0,
#endif
					 u_int n ){

#ifdef WIN32
	hash_map<u_int, u_int>::iterator found;
	hash_map<u_int, float>::iterator found1;
#else
	unordered_map<u_int, u_int>::iterator found;
	unordered_map<u_int, float>::iterator found1;
#endif
	PCGSolver<double> solver;
	int cgIters = 100;
	if(n > 500000){
		cgIters *= int(round((float)n/500000.0f));
	}
	solver.set_solver_parameters(1e-13, cgIters, 0.97, 0.25);

	double residual;
	int iterations;
	SparseMatrixd A(n);
	vector<double> rhs;
	vector<double> p;
	for(u_int i=0; i<n; ++i)
		p.push_back((double)x[i]);
	u_int diag;

	FOR_EACH_CELL
		 u_int index = INDEX(i,j,k);
		 if( mMarker[index] == 1 ){
			 found = a.find(index);
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back((double)found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 int solidNeighbors=0;
			 for(u_int m=0; m<neighbors.size(); m++){
				 if(voxel->InSolid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ||
					voxel->InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
					++solidNeighbors;
				 }
				 else if(voxel->InSource(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
					 ++solidNeighbors;
				 }
				 else if(mMarker[INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)] == 1){
					 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
					 A.set_element(diag,found->second,-1);
				 }
			 }
			 if(solidNeighbors == 6){
				printf("at (%d, %d, %d) Newmann b.c. for all directions \n", i,j,k);
				exit(1);
			 }

			 A.set_element(diag,diag,(float)(neighbors.size()-solidNeighbors));
//			 TARGET_CELL
//				 printf("non solid neighbors = %d at(%d, %d,%d) \n",
//						 (neighbors.size()-solidNeighbors), i,j,k);
		 }
	END_FOR


	A.check_symmetry();

	u_int M = 0;
	if(!A.check_diagnal(M)){
		FOR_EACH_CELL
			u_int index = INDEX(i,j,k);
			found = a.find(index);
			if(found != a.end()){
				diag = found->second;
				if( diag == M ){
					printf(" at (%d, %d, %d) matrix diagnal is too small or negative\n", i, j, k);
					exit(1);
				}
			}
		END_FOR
	}

	solver.solve(A, rhs, p, residual, iterations);

	for(u_int i=0; i<n; ++i)
		x[i] = (float)p[i];

	printf("\nPCG iterations = %d, and residual = %16.13f \n\n ",
			iterations, residual);

}

void FLIP::ApplyBoundaryConditions(const Voxel *voxel, const VectorField3D *vel){
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//		if(mMarker[pos] == 1){  // fluid cell
			if(i > 0){
				if(voxel->InLiquid(i-1,j,k) || voxel->InSource(i-1,j,k) || voxel->InSolid(i-1,j,k))
					vel->SetGridVelocity(i-1, j, k, 1, u);
			}
			if(i < DimX-1){
				if(voxel->InLiquid(i+1,j,k) || voxel->InSource(i+1,j,k) || voxel->InSolid(i+1,j,k))
					vel->SetGridVelocity(i, j, k, 1, u);
			}
			if(j > 0){
				if(voxel->InLiquid(i,j-1,k) || voxel->InSource(i,j-1,k) || voxel->InSolid(i,j-1,k))
					vel->SetGridVelocity(i, j-1, k, 2, v);
			}
			if(j < DimY-1){
				if(voxel->InLiquid(i,j+1,k) || voxel->InSource(i,j+1,k) || voxel->InSolid(i,j+1,k))
					vel->SetGridVelocity(i, j, k, 2, v);
			}
			if(k > 0){
				if(voxel->InLiquid(i,j,k-1) || voxel->InSource(i,j,k-1) || voxel->InSolid(i,j,k-1))
					vel->SetGridVelocity(i, j, k-1, 3, w);
			}
			if(k < DimZ-1){
				if(voxel->InLiquid(i,j,k+1) || voxel->InSource(i,j,k+1) || voxel->InSolid(i,j,k+1))
					vel->SetGridVelocity(i, j, k, 3, w);
			}
			if(i > 0){
				if(voxel->InSolid(i-1,j,k))
					u[INDEX(i-1,j,k)] = 0.f;
			}
			if(i < DimX-1){
				if(voxel->InSolid(i+1,j,k))
					u[pos] = 0.f;
			}
			if(j > 0){
				if(voxel->InSolid(i,j-1,k))
					v[INDEX(i,j-1,k)] = 0.f;
			}
			if(j < DimY-1){
				if(voxel->InSolid(i,j+1,k))
					v[pos] = 0.f;
			}
			if(k > 0){
				if(voxel->InSolid(i,j,k-1))
					w[INDEX(i,j,k-1)] = 0.f;
			}
			if(k < DimZ-1){
				if(voxel->InSolid(i,j,k+1))
					w[pos] = 0.f;
			}

//		}
	END_FOR
}

void FLIP::ExtrapolateGridVelocity(){
//	PrintOneCellVel(I,J,K,"before extrapolation");
	for(int i=0; i<8; ++i){
		SweepVelocity();
//		char s[18] = "Sweep ";
//		char f[6];
//		sprintf(f, "f%d",i);
//		strcat(s,f);
//		PrintOneCellVel(I,J,K,s);
//		PrintOneCellPhi(I,J,K,s);
	}
}

void FLIP::SweepVelocity(){
	SweepU(1, DimX-1, 1, DimY-1, 1, DimZ-1);
	SweepU(1, DimX-1, DimY-2, 0, 1, DimZ-1);
	SweepU(1, DimX-1, 1, DimY-1, DimZ-2, 0);
	SweepU(1, DimX-1, DimY-2, 0, DimZ-2, 0);
	SweepU(DimX-2, 0, 1, DimY-1, 1, DimZ-1);
	SweepU(DimX-2, 0, DimY-2, 0, 1, DimZ-1);
	SweepU(DimX-2, 0, 1, DimY-1, DimZ-2, 0);
	SweepU(DimX-2, 0, DimY-2, 0, DimZ-2, 0);

	SweepV(1, DimX-1, 1, DimY-1, 1, DimZ-1);
	SweepV(1, DimX-1, DimY-2, 0, 1, DimZ-1);
	SweepV(1, DimX-1, 1, DimY-1, DimZ-2, 0);
	SweepV(1, DimX-1, DimY-2, 0, DimZ-2, 0);
	SweepV(DimX-2, 0, 1, DimY-1, 1, DimZ-1);
	SweepV(DimX-2, 0, DimY-2, 0, 1, DimZ-1);
	SweepV(DimX-2, 0, 1, DimY-1, DimZ-2, 0);
	SweepV(DimX-2, 0, DimY-2, 0, DimZ-2, 0);

	SweepW(1, DimX-1, 1, DimY-1, 1, DimZ-1);
	SweepW(1, DimX-1, DimY-2, 0, 1, DimZ-1);
	SweepW(1, DimX-1, 1, DimY-1, DimZ-2, 0);
	SweepW(1, DimX-1, DimY-2, 0, DimZ-2, 0);
	SweepW(DimX-2, 0, 1, DimY-1, 1, DimZ-1);
	SweepW(DimX-2, 0, DimY-2, 0, 1, DimZ-1);
	SweepW(DimX-2, 0, 1, DimY-1, DimZ-2, 0);
	SweepW(DimX-2, 0, DimY-2, 0, DimZ-2, 0);
}

void FLIP::SweepU(int i0, int i1, int j0, int j1, int k0, int k1)
{
   int di=(i0<i1) ? 1 : -1;
   int dj=(j0<j1) ? 1 : -1;
   int dk=(k0<k1) ? 1 : -1;
   float dp, dq, dr, alpha, beta;
   for(int k=k0; k!=k1; k+=dk)
	   for(int j=j0; j!=j1; j+=dj)
		   for(int i=i0; i!=i1; i+=di){
			  if(mMarker[INDEX(i,j,k)] != 1 && mMarker[INDEX(i+1,j,k)] != 1){
				 dp=di*(phi[INDEX(i+1,j,k)]-phi[INDEX(i,j,k)]);
				 if(dp<0) continue; // not useful on this sweep direction
				 dq=0.5*(phi[INDEX(i+1,j,k)]+phi[INDEX(i,j,k)]-phi[INDEX(i+1,j-dj,k)]-phi[INDEX(i,j-dj,k)]);
				 if(dq<0) continue; // not useful on this sweep direction
				 dr=0.5*(phi[INDEX(i+1,j,k)]+phi[INDEX(i,j,k)]-phi[INDEX(i+1,j,k-dk)]-phi[INDEX(i,j,k-dk)]);
				 if(dr<0) continue; // not useful on this sweep direction
				 float t=dp+dq+dr;
				 if(t==0) {
					 alpha = beta = 1.f/3;
				 }
				 else {
					 alpha = dp/t;
					 beta  = dq/t;
				 }
				 u[INDEX(i,j,k)]=alpha*u[INDEX(i-di,j,k)]+
								 beta *u[INDEX(i,j-dj,k)]+
								 (1-alpha-beta)*u[INDEX(i,j,k-dk)];
			  }
		   }
}

void FLIP::SweepV(int i0, int i1, int j0, int j1, int k0, int k1)
{
   int di=(i0<i1) ? 1 : -1;
   int dj=(j0<j1) ? 1 : -1;
   int dk=(k0<k1) ? 1 : -1;
   float dp, dq, dr, alpha, beta;
   for(int k=k0; k!=k1; k+=dk)
	   for(int j=j0; j!=j1; j+=dj)
		   for(int i=i0; i!=i1; i+=di){
			  if(mMarker[INDEX(i,j,k)] != 1 && mMarker[INDEX(i,j+1,k)] != 1){
				 dp= 0.5*(phi[INDEX(i,j+1,k)]+phi[INDEX(i,j,k)]-phi[INDEX(i-di,j+1,k)]-phi[INDEX(i-di,j,k)]);
				 if(dp<0) continue; // not useful on this sweep direction
				 dq=di*(phi[INDEX(i,j+1,k)]-phi[INDEX(i,j,k)]);
				 if(dq<0) continue; // not useful on this sweep direction
				 dr=0.5*(phi[INDEX(i,j+1,k)]+phi[INDEX(i,j,k)]-phi[INDEX(i,j+1,k-dk)]-phi[INDEX(i,j,k-dk)]);
				 if(dr<0) continue; // not useful on this sweep direction
				 float t=dp+dq+dr;
				 if(t==0) {
					 alpha = beta = 1.f/3;
				 }
				 else {
					 alpha = dp/t;
					 beta  = dq/t;
				 }
				 v[INDEX(i,j,k)]=alpha*v[INDEX(i-di,j,k)]+
								 beta *v[INDEX(i,j-dj,k)]+
								 (1-alpha-beta)*v[INDEX(i,j,k-dk)];
			  }
		   }
}

void FLIP::SweepW(int i0, int i1, int j0, int j1, int k0, int k1)
{
   int di=(i0<i1) ? 1 : -1;
   int dj=(j0<j1) ? 1 : -1;
   int dk=(k0<k1) ? 1 : -1;
   float dp, dq, dr, alpha, beta;
   for(int k=k0; k!=k1; k+=dk)
	   for(int j=j0; j!=j1; j+=dj)
		   for(int i=i0; i!=i1; i+=di){
			  if(mMarker[INDEX(i,j,k)] != 1 && mMarker[INDEX(i,j,k+1)] != 1){
				 dp= 0.5*(phi[INDEX(i,j,k+1)]+phi[INDEX(i,j,k)]-phi[INDEX(i-di,j,k+1)]-phi[INDEX(i-di,j,k)]);
				 if(dp<0) continue; // not useful on this sweep direction
				 dq=0.5*(phi[INDEX(i,j,k+1)]+phi[INDEX(i,j,k)]-phi[INDEX(i,j-dj,k+1)]-phi[INDEX(i,j-dj,k)]);
				 if(dq<0) continue; // not useful on this sweep direction
				 dr=di*(phi[INDEX(i,j,k+1)]-phi[INDEX(i,j,k)]);
				 if(dr<0) continue; // not useful on this sweep direction
				 float t=dp+dq+dr;
				 if(t==0) {
					 alpha = beta = 1.f/3;
				 }
				 else {
					 alpha = dp/t;
					 beta  = dq/t;
				 }
				 w[INDEX(i,j,k)]=alpha*w[INDEX(i-di,j,k)]+
								 beta *w[INDEX(i,j-dj,k)]+
								 (1-alpha-beta)*w[INDEX(i,j,k-dk)];
			  }
		   }
}

void FLIP::ComputeGridLevelset(){
	InitializePhi();
	// in 3D, we need 2^3(=8) GS sweeps to have converged distance values
	for(int i=0; i<8; ++i)
	    SweepPhi();
}

static float SolveDistance(float o, float p, float q, float s){
	vector<float> a;
	a.push_back(o);
	a.push_back(p);
	a.push_back(q);
	std::sort(a.begin(), a.end());
	float d = a[0]+1;
    if(d > a[1]){
    	 float sq = a[0]-a[1];
    	 d = (a[0]+a[1]+sqrt(2-sq*sq))/2;
    	 if(d > a[2]){
    		 sq = a[1]-a[2];
    		 d = (a[1]+a[2]+sqrt(2-sq*sq))/2;
    	 }
    }
    if(d < s)
    	return d;
    else
    	return s;

}

void FLIP::SweepPhi(){
	// fast sweeping outside the fluid in all 8 sweep directions

	for(int k=1;k<DimZ;k++) {
		for(int j=1;j<DimY;j++) {
			for(int i=1;i<DimX;i++){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i-1,j,k)], phi[INDEX(i,j-1,k)],
													  phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k)]);

			}
		}
	}
	for(int k=1;k<DimZ;k++) {
		for(int j=DimY-2;j>=0;j--) {
			for(int i=1;i<DimX;i++){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i-1,j,k)], phi[INDEX(i,j+1,k)],
													  phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k)]);

			}
		}
	}
	for(int k=1;k<DimZ;k++) {
		for(int j=1;j<DimY;j++) {
			for(int i=DimX-2;i>=0;i--){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i+1,j,k)], phi[INDEX(i,j-1,k)],
													  phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k)]);

			}
		}
	}
	for(int k=1;k<DimZ;k++) {
		for(int j=DimY-2;j>=0;j--) {
			for(int i=DimX-2;i>=0;i--){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i+1,j,k)], phi[INDEX(i,j+1,k)],
													  phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k)]);

			}
		}
	}

	for(int k=DimZ-2;k>=0;k--) {
		for(int j=1;j<DimY;j++) {
			for(int i=1;i<DimX;i++){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i-1,j,k)], phi[INDEX(i,j-1,k)],
													  phi[INDEX(i,j,k+1)], phi[INDEX(i,j,k)]);

			}
		}
	}
	for(int k=DimZ-2;k>=0;k--) {
		for(int j=DimY-2;j>=0;j--) {
			for(int i=1;i<DimX;i++){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i-1,j,k)], phi[INDEX(i,j+1,k)],
													  phi[INDEX(i,j,k+1)], phi[INDEX(i,j,k)]);

			}
		}
	}
	for(int k=DimZ-2;k>=0;k--) {
		for(int j=1;j<DimY;j++) {
			for(int i=DimX-2;i>=0;i--){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i+1,j,k)], phi[INDEX(i,j-1,k)],
													  phi[INDEX(i,j,k+1)], phi[INDEX(i,j,k)]);

			}
		}
	}
	for(int k=DimZ-2;k>=0;k--) {
		for(int j=DimY-2;j>=0;j--) {
			for(int i=DimX-2;i>=0;i--){
				if(mMarker[INDEX(i,j,k)] != 1)  // non_fluid cell
					phi[INDEX(i,j,k)] = SolveDistance(phi[INDEX(i+1,j,k)], phi[INDEX(i,j+1,k)],
													  phi[INDEX(i,j,k+1)], phi[INDEX(i,j,k)]);

			}
		}
	}

}

void FLIP::InitializePhi(){
	float large_distance=DimX+DimY+DimZ+3;
	SetScalar(large_distance, phi);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mMarker[pos] == 1)  // fluid cell
			phi[pos] = -0.5f;
	END_FOR

}

void FLIP::UpdateParticleVelocityFromGridValues(){
	float up, vp, wp, dup, dvp, dwp;
	float old_up, old_vp, old_wp;

//	list<WaterParticle>::iterator particle;
	WaterParticleList::iterator particle;

	for(particle  = mWaterParticles.begin();
		particle != mWaterParticles.end();){
		WaterParticle &ps = *particle;
		Point pos = ps.Position();
		old_up = ps.UVelocity();
		old_vp = ps.VVelocity();
		old_wp = ps.WVelocity();
		TriInterp(up, vp, wp, pos, u, v, w);
		TriInterp(dup, dvp, dwp, pos, du, dv, dw);
		up = mPICFLIPRatio * up + (1.f - mPICFLIPRatio) * (old_up + dup);
		vp = mPICFLIPRatio * vp + (1.f - mPICFLIPRatio) * (old_vp + dvp);
		wp = mPICFLIPRatio * wp + (1.f - mPICFLIPRatio) * (old_wp + dwp);
		ps.SetVelocity(up, vp, wp);
		++particle;
	}
}

void FLIP::UpdateGridVelocityIncrement(){
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		du[pos] = u[pos] - du[pos];
		dv[pos] = v[pos] - dv[pos];
		dw[pos] = w[pos] - dw[pos];
	END_FOR
}

void FLIP::ApplyForce(float dt){
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		w[pos] -= dt*g;
	END_FOR
}



void FLIP::AccumulateValues(const Point &pos, float vel, int index,
		float *wh, float *vc){

	int i0, j0, k0, i1, j1, k1;
	int x0, y0, z0;
	float xp0, yp0, zp0;

	float s0, t0, s1, t1;
	float r0, r1;
	float weight;

	float xp = pos.x * invh;
	float yp = pos.y * invh;
	float zp = pos.z * invh;

	xp = xp < 0.5f      ? 0.5f      : xp;
	xp = xp > DimX-1.5f ? DimX-1.5f : xp;
	yp = yp < 0.5f      ? 0.5f      : yp;
	yp = yp > DimY-1.5f ? DimY-1.5f : yp;
	zp = zp < 0.5f      ? 0.5f      : zp;
	zp = zp > DimZ-1.5f ? DimZ-1.5f : zp;

	if(index == 1){

		i0=floor(xp)-1; if(i0 < 0) i0 = 0;
		i1 = i0 + 1;
		s1 = xp - 1 - i0;
		s0 = 1 - s1;
		j0 = floor(yp); y0 = j0;
		if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		j1=j0+1;
		yp0 = yp-y0;
		t1 = yp0 < 0.5f ? 0.5f+yp0 : yp0-0.5f;
		t0 = 1-t1;
		k0 = floor(zp); z0 = k0;
		if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
		k1=k0+1;
		zp0 = zp-z0;
		r1 = zp0 < 0.5f ? 0.5f+zp0 : zp0-0.5f;
		r0 = 1-r1;

	}

	if(index == 2){

		j0=floor(yp)-1; if(j0 < 0) j0 = 0;
		j1 = j0 + 1;
		t1 = yp - 1 - j0;
		t0 = 1 - t1;
		i0 = floor(xp); x0 = i0;
		if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		i1=i0+1;
		xp0 = xp - x0;
		s1 = xp0 < 0.5f ? 0.5f+xp0 : xp0-0.5f;
		s0 = 1-s1;
		k0 = floor(zp); z0 = k0;
		if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
		k1=k0+1;
		zp0 = zp-z0;
		r1 = zp0 < 0.5f ? 0.5f+zp0 : zp0-0.5f;
		r0 = 1-r1;

	}


	if(index == 3){

		k0=floor(zp)-1; if(k0 < 0) k0 = 0;
		k1 = k0+1;
		r1 = zp - 1 - k0;
		r0 = 1 - r1;
		i0 = floor(xp); x0 = i0;
		if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		i1=i0+1;
		xp0 = xp - x0;
		s1 = xp0 < 0.5f ? 0.5f+xp0 : xp0-0.5f;
		s0 = 1-s1;
		j0 = floor(yp); y0 = j0;
		if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		j1=j0+1;
		yp0 = yp-y0;
		t1 = yp0 < 0.5f ? 0.5f+yp0 : yp0-0.5f;
		t0 = 1-t1;

	}


	weight = r0*s0*t0;
	vc[INDEX(i0,j0,k0)] += weight * vel;
	wh[INDEX(i0,j0,k0)] += weight;
	weight = r0*s0*t1;
	vc[INDEX(i0,j1,k0)] += weight * vel;
	wh[INDEX(i0,j1,k0)] += weight;
	weight = r0*s1*t0;
	vc[INDEX(i1,j0,k0)] += weight * vel;
	wh[INDEX(i1,j0,k0)] += weight;
	weight = r0*s1*t1;
	vc[INDEX(i1,j1,k0)] += weight * vel;
	wh[INDEX(i1,j1,k0)] += weight;
	weight = r1*s0*t0;
	vc[INDEX(i0,j0,k1)] += weight * vel;
	wh[INDEX(i0,j0,k1)] += weight;
	weight = r1*s0*t1;
	vc[INDEX(i0,j1,k1)] += weight * vel;
	wh[INDEX(i0,j1,k1)] += weight;
	weight = r1*s1*t0;
	vc[INDEX(i1,j0,k1)] += weight * vel;
	wh[INDEX(i1,j0,k1)] += weight;
	weight = r1*s1*t1;
	vc[INDEX(i1,j1,k1)] += weight * vel;
	wh[INDEX(i1,j1,k1)] += weight;

}

void FLIP::TransferParticleValuesToGrid(const Voxel *voxel){
//	list<WaterParticle>::iterator particle;
	WaterParticleList::iterator particle;

	SetZero(u);
	SetZero(v);
	SetZero(w);

	SetZero(mWeights);
	for(particle  = mWaterParticles.begin();
	    particle != mWaterParticles.end();){
		WaterParticle &ps = *particle;
		Point pos = ps.Position();
		Vector vel = ps.Velocity();
		AccumulateValues(pos, vel.x, 1, mWeights, u);
		++particle;
	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mWeights[pos] != 0.f)
			u[pos] /= mWeights[pos];
	END_FOR

	SetZero(mWeights);
	for(particle  = mWaterParticles.begin();
		particle != mWaterParticles.end();){
		WaterParticle &ps = *particle;
		Point pos = ps.Position();
		Vector vel = ps.Velocity();
		AccumulateValues(pos, vel.y, 2, mWeights, v);
		++particle;
	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mWeights[pos] != 0.f)
			v[pos] /= mWeights[pos];
	END_FOR

	SetZero(mWeights);
	for(particle  = mWaterParticles.begin();
		particle != mWaterParticles.end();){
		WaterParticle &ps = *particle;
		Point pos = ps.Position();
		Vector vel = ps.Velocity();
		AccumulateValues(pos, vel.z, 3, mWeights, w);
		++particle;
	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//		TARGET_CELL
//			printf("mWeights = %f, w = %f \n", mWeights[pos],w[pos]);
		if(mWeights[pos] != 0.f)
			w[pos] /= mWeights[pos];
	END_FOR



	PrintOneCellVel(I,J,K, "after particle transfer");


	SetZero(mMarker);
	for(particle  = mWaterParticles.begin();
		particle != mWaterParticles.end();){
		WaterParticle &ps = *particle;
		Point pos = ps.Position();
		Vector vel = ps.Velocity();
		int i = (int)(pos.x * invh);
		int j = (int)(pos.y * invh);
		int k = (int)(pos.z * invh);
		if(!voxel->InLiquid(i,j,k) && !voxel->InSource(i,j,k) && !voxel->InSolid(i,j,k))
			mMarker[INDEX(i,j,k)] = 1;
		TARGET_CELL
			printf("in cell (%d,%d,%d), marker=%d, particle at (%f,%f,%f) with vel(%f,%f,%f) \n",
					i,j,k, mMarker[INDEX(i,j,k)],
					pos.x, pos.y, pos.z, vel.x, vel.y, vel.z);
		++particle;
	}

}


void FLIP::SaveGridVelocities(){

	SetEqual(u, du);
	SetEqual(v, dv);
	SetEqual(w, dw);

}


void FLIP::MoveParticlesInGrid(const Voxel *voxel){
//	float dts = mTimeManager->GetSubcyclingDt();
//	int iterations = mTimeManager->GetSubcyclingIterations();

	float dt = mTimeManager->GetDt();
//	list<WaterParticle>::iterator particle;
	WaterParticleList::iterator particle;

	float velmax = 0.f;
	Vector pvel(0.f);
	for(particle=mWaterParticles.begin(); particle != mWaterParticles.end(); ){
		WaterParticle &ps = *particle;
		Vector vel = ps.Velocity();
		float velmag = vel.Length();
		if(velmag > velmax){
			velmax = velmag;
			pvel = vel;
		}
		++particle;
	}
	printf(" max particle vel = (%f, %f, %f) mag = %f \n",
			pvel.x, pvel.y, pvel.z, velmax);
	int iterations;
	float dts;
	if(velmax > 0.f)
		dts = h / velmax;
	else
		dts = dt;
	iterations = floor(dt/dts);
	printf("\n Advection step with dts = %f, iters = %d \n\n", dts, iterations);
	for(int n = 0; n < iterations; ++n){
		AdvectRK3Subcycling(voxel, dts);
	}
	if(dts*iterations < dt){
		dts = dt-dts*iterations;
		AdvectRK3Subcycling(voxel, dts);
	}
}

void FLIP::AdvectRK3Subcycling(const Voxel *voxel, float dts){

	Vector pvel;
	Point pold;

	int N = 0;

	float delta = h;
//	list<WaterParticle>::iterator particle;
	WaterParticleList::iterator particle;

	float *phi_c_obj = mObject->getScalar();

	for(particle=mWaterParticles.begin(); particle != mWaterParticles.end(); ){
		WaterParticle &ps = *particle;
		Point pos = ps.Position();
		float colldist = ps.CollisionDistance();
		pold  = pos;
		TentativeGridPoint o_tp = voxel->ContainsPoint(pos);
		if(IndexOutofBounds(o_tp.ii, o_tp.jj, o_tp.kk))
			particle = mWaterParticles.erase(particle);
		else{
			TriInterp(pvel.x, pvel.y, pvel.z, pos, u, v, w);
			pos = pold + 0.5f*dts*pvel;
			float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
			if( phi_i_new >= 0.f ||	fabsf(phi_i_new) < colldist ){
				Vector o_normal = ObjectNormalAt(voxel, pos, phi_c_obj);
				Point projP = pos + (colldist + phi_i_new) * o_normal;
				pos = projP;
			}
			TriInterp(pvel.x, pvel.y, pvel.z, pos, u, v, w);
			pos = pold + dts*pvel;
			TentativeGridPoint tp = voxel->ContainsPoint(pos);
			if(IndexOutofBounds(tp.ii, tp.jj, tp.kk)){
//				printf(" particle at (%f, %f, %f) with vel (%f,%f,%f) out of bounds \n",
//					pos.x, pos.y, pos.z, up, vp, wp);
//				printf(" particle originally at (%f,%f,%f) object levelset = %f \n ",
//						o_xp, o_yp, o_zp, phi_i_new);
//				printf(" up1 = %f, vp1 = %f, wp1 = %f \n",
//						up1, vp1, wp1);
//				printf(" particle original vel = (%f, %f, %f) dts = %f\n",
//						old_up, old_vp, old_zp, dts);
				particle = mWaterParticles.erase(particle);
			}
			else{
				phi_i_new = TriInterp(voxel, pos, phi_c_obj);
				if( phi_i_new >= 0.f ||	fabsf(phi_i_new) < colldist ){
					Vector o_normal = ObjectNormalAt(voxel, pos, phi_c_obj);
					Point projP = pos + (colldist + phi_i_new) * o_normal;
					pos = projP;
				}
				ps.SetParticle(pos);
				if(tp.ii == I && tp.jj == J && tp.kk == K)
					++N;
				++particle;
			}
		}
	}

}

Vector FLIP::ObjectNormalAt(const Voxel *voxel, const Point &p, float *phi_c_obj) const{
	float gx, gy, gz;

	float delta = voxel->VoxelDelta();

	Point px0(p.x-0.5f*delta, p.y, p.z);
	Point px1(p.x+0.5f*delta, p.y, p.z);
	float phix0 = TriInterp(voxel, px0, phi_c_obj);
	float phix1 = TriInterp(voxel, px1, phi_c_obj);
	gx = (phix1 - phix0) / delta;

	Point py0(p.x, p.y-0.5f*delta, p.z);
	Point py1(p.x, p.y+0.5f*delta, p.z);
	float phiy0 = TriInterp(voxel, py0, phi_c_obj);
	float phiy1 = TriInterp(voxel, py1, phi_c_obj);
	gy = (phiy1 - phiy0) / delta;

	Point pz0(p.x, p.y, p.z-0.5f*delta);
	Point pz1(p.x, p.y, p.z+0.5f*delta);
	float phiz0 = TriInterp(voxel, pz0, phi_c_obj);
	float phiz1 = TriInterp(voxel, pz1, phi_c_obj);
	gz = (phiz1 - phiz0) / delta;

	Vector N(gx, gy, gz);
	if(N.Length() == 0.f){
		MTRand mt;
		double rn1 = mt();
		double rn2 = mt();
		double rn3 = mt();
		N.x = (float)rn1;
		N.y = (float)rn2;
		N.z = (float)rn3;
	}
	N = Normalize(N);
	return -1*N;  // multiply by -1 because phi at object grid points is positive

}

float FLIP::TriInterp(const Voxel *voxel, const Point &p, float *phi) const{

	float delta = voxel->VoxelDelta();

	float x = p.x/delta - 0.5f;
    float y = p.y/delta - 0.5f;
    float z = p.z/delta - 0.5f;
    if(x < 0) x = 0.f;
    if(y < 0) y = 0.f;
    if(z < 0) z = 0.f;
	int xi = floor(x);
	int yi = floor(y);
	int zi = floor(z);
//	printf(" xi = %d, yi = %d, zi = %d \n", xi, yi, zi);
	//Point point = voxel.VoxelCenterPosition(xi,yi,zi);
	Point point = voxel->VoxelCornerPosition(xi,yi,zi,1);
	float r = x - point.x/delta;
    float s = y - point.y/delta;
    float t = z - point.z/delta;
    int xj, yj, zj;
    if( xi >= DimX-1 ){
    	xi = DimX-1;
    	xj = xi;
    	r = 1.f;
    }
    else
    	xj = xi + 1;
    if( yi >= DimY-1 ){
    	yi = DimY-1;
        yj = yi;
        s = 1.f;
    }
    else
    	yj = yi + 1;
    if( zi >= DimZ-1 ){
    	zi = DimZ-1;
        zj = zi;
        t = 1.f;
    }
    else
    	zj = zi + 1;

//	printf(" r = %f,  s = %f,  t = %f \n",  r, s, t);
	return trilerp(phi[INDEX(xi, yi, zi)], phi[INDEX(xj, yi, zi)],
				   phi[INDEX(xi, yi, zj)], phi[INDEX(xj, yi, zj)],
				   phi[INDEX(xi, yj, zi)], phi[INDEX(xj, yj, zi)],
				   phi[INDEX(xi, yj, zj)], phi[INDEX(xj, yj, zj)],
				   r, s, t);
}


void FLIP::TriInterp(float &ui, float &vi, float &wi,
		const Point &p,	float *uvel, float *vvel, float *wvel) const {

	int i0, j0, k0, i1, j1, k1;
	float s0, t0, s1, t1;
	int x0, y0, z0;
	float xp0, yp0, zp0;
	float r0, r1;

	float xp = p.x * invh;
	float yp = p.y * invh;
	float zp = p.z * invh;

	xp = xp < 0.5f      ? 0.5f      : xp;
	xp = xp > DimX-1.5f ? DimX-1.5f : xp;
	yp = yp < 0.5f      ? 0.5f      : yp;
	yp = yp > DimY-1.5f ? DimY-1.5f : yp;
	zp = zp < 0.5f      ? 0.5f      : zp;
	zp = zp > DimZ-1.5f ? DimZ-1.5f : zp;

	i0=floor(xp)-1; if(i0 < 0) i0 = 0;
	i1 = i0 + 1;
	s1 = xp - 1 - i0;
	s0 = 1 - s1;
	j0 = floor(yp); y0 = j0;
	if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
	j1=j0+1;
	yp0 = yp-y0;
	t1 = yp0 < 0.5f ? 0.5f+yp0 : yp0-0.5f;
	t0 = 1-t1;
	k0 = floor(zp); z0 = k0;
	if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
	k1=k0+1;
	zp0 = zp-z0;
	r1 = zp0 < 0.5f ? 0.5f+zp0 : zp0-0.5f;
	r0 = 1-r1;



//	CheckIndex(i0,j0,k0);
//	CheckIndex(i0,j1,k0);
//	CheckIndex(i1,j0,k0);
//	CheckIndex(i1,j1,k0);
//	CheckIndex(i0,j0,k1);
//	CheckIndex(i0,j1,k1);
//	CheckIndex(i1,j0,k1);
//	CheckIndex(i1,j1,k1);

	ui = r0*(s0*(t0*uvel[INDEX(i0,j0,k0)]+t1*uvel[INDEX(i0,j1,k0)]) +
		     s1*(t0*uvel[INDEX(i1,j0,k0)]+t1*uvel[INDEX(i1,j1,k0)]))+
 		 r1*(s0*(t0*uvel[INDEX(i0,j0,k1)]+t1*uvel[INDEX(i0,j1,k1)]) +
			 s1*(t0*uvel[INDEX(i1,j0,k1)]+t1*uvel[INDEX(i1,j1,k1)]));


	j0=floor(yp)-1; if(j0 < 0) j0 = 0;
	j1 = j0 + 1;
	t1 = yp - 1 - j0;
	t0 = 1 - t1;
	i0 = floor(xp); x0 = i0;
	if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
	i1=i0+1;
	xp0 = xp - x0;
	s1 = xp0 < 0.5f ? 0.5f+xp0 : xp0-0.5f;
	s0 = 1-s1;
	k0 = floor(zp); z0 = k0;
	if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
	k1=k0+1;
	zp0 = zp-z0;
	r1 = zp0 < 0.5f ? 0.5f+zp0 : zp0-0.5f;
	r0 = 1-r1;


//	CheckIndex(i0,j0,k0);
//	CheckIndex(i0,j1,k0);
//	CheckIndex(i1,j0,k0);
//	CheckIndex(i1,j1,k0);
//	CheckIndex(i0,j0,k1);
//	CheckIndex(i0,j1,k1);
//	CheckIndex(i1,j0,k1);
//	CheckIndex(i1,j1,k1);

    vi = r0*(s0*(t0*vvel[INDEX(i0,j0,k0)]+t1*vvel[INDEX(i0,j1,k0)])+
		     s1*(t0*vvel[INDEX(i1,j0,k0)]+t1*vvel[INDEX(i1,j1,k0)])) +
		 r1*(s0*(t0*vvel[INDEX(i0,j0,k1)]+t1*vvel[INDEX(i0,j1,k1)])+
		     s1*(t0*vvel[INDEX(i1,j0,k1)]+t1*vvel[INDEX(i1,j1,k1)]));

    k0=floor(zp)-1; if(k0 < 0) k0 = 0;
	k1 = k0+1;
	r1 = zp - 1 - k0;
	r0 = 1 - r1;
	i0 = floor(xp); x0 = i0;
	if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
	i1=i0+1;
	xp0 = xp - x0;
	s1 = xp0 < 0.5f ? 0.5f+xp0 : xp0-0.5f;
	s0 = 1-s1;
	j0 = floor(yp); y0 = j0;
	if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
	j1=j0+1;
	yp0 = yp-y0;
	t1 = yp0 < 0.5f ? 0.5f+yp0 : yp0-0.5f;
	t0 = 1-t1;


//	CheckIndex(i0,j0,k0);
//	CheckIndex(i0,j1,k0);
//	CheckIndex(i1,j0,k0);
//	CheckIndex(i1,j1,k0);
//	CheckIndex(i0,j0,k1);
//	CheckIndex(i0,j1,k1);
//	CheckIndex(i1,j0,k1);
//	CheckIndex(i1,j1,k1);

    wi = r0*(s0*(t0*wvel[INDEX(i0,j0,k0)]+t1*wvel[INDEX(i0,j1,k0)])+
	     	 s1*(t0*wvel[INDEX(i1,j0,k0)]+t1*wvel[INDEX(i1,j1,k0)])) +
	     r1*(s0*(t0*wvel[INDEX(i0,j0,k1)]+t1*wvel[INDEX(i0,j1,k1)])+
			 s1*(t0*wvel[INDEX(i1,j0,k1)]+t1*wvel[INDEX(i1,j1,k1)]));
//    if(i0 == I && j0 == J && k0 == K){
//		printf("particle at (%f, %f, %f) with w = %f \n ",
//			   	xp, yp, zp, wi);
//		printf("t0 = %f, t1 = %f, s0 = %f, s1 = %f, r0 = %f, r1 = %f \n",
//				t0, t1, s0, s1, r0, r1);
//		printf("w(i0,j0,k0) = %f,w(i0,j1,k0) = %f \n",
//				wvel[INDEX(i0,j0,k0)], wvel[INDEX(i0,j1,k0)]);
//		printf("w(i1,j0,k0) = %f,w(i1,j1,k0) = %f \n",
//				wvel[INDEX(i1,j0,k0)], wvel[INDEX(i1,j1,k0)]);
//		printf("w(i0,j0,k1) = %f,w(i0,j1,k1) = %f \n",
//				wvel[INDEX(i0,j0,k1)], wvel[INDEX(i0,j1,k1)]);
//		printf("w(i1,j0,k1) = %f,w(i1,j1,k1) = %f \n",
//				wvel[INDEX(i1,j0,k1)], wvel[INDEX(i1,j1,k1)]);
//
//	}
}

void FLIP::OutputWaterParticlesBinary(const Voxel *voxel, const char *name){

	float delta = voxel->VoxelDelta();

//	list<WaterParticle>::iterator iter_particle;
	WaterParticleList::iterator iter_particle;


	FILE *fp = fopen(name, "wb");
	if(!fp){
      printf("Couldn't open file [ %s ] to write \n", name);
      exit(-1);
   	}
	u_int num = mWaterParticles.size();
	Point p0 = voxel->GetP0();
	//   	fprintf(fp,"%f %f %f \n",p0.x, p0.y, p0.z);
   	fwrite(&(p0.x), sizeof(float), 1, fp);
   	fwrite(&(p0.y), sizeof(float), 1, fp);
   	fwrite(&(p0.z), sizeof(float), 1, fp);
   	Point p1 = voxel->GetP1();
//    fprintf(fp,"%f %f %f \n",p1.x, p1.y, p1.z);
    fwrite(&(p1.x), sizeof(float), 1, fp);
    fwrite(&(p1.y), sizeof(float), 1, fp);
    fwrite(&(p1.z), sizeof(float), 1, fp);

    fwrite(&delta, sizeof(float), 1, fp);
	fwrite(&num, sizeof(u_int), 1, fp);
	for(iter_particle  = mWaterParticles.begin();
		iter_particle != mWaterParticles.end();
		){
		WaterParticle &wps = *iter_particle;
		float r = wps.Radius();
		Point pos = wps.Position();
		Vector vel = wps.Velocity();
		fwrite(&(pos.x), sizeof(float), 1, fp);
		fwrite(&(pos.y), sizeof(float), 1, fp);
		fwrite(&(pos.z), sizeof(float), 1, fp);
		fwrite(&(vel.x), sizeof(float), 1, fp);
		fwrite(&(vel.y), sizeof(float), 1, fp);
		fwrite(&(vel.z), sizeof(float), 1, fp);
		fwrite(&r, sizeof(float), 1, fp);
		++iter_particle;
	}
	fclose(fp);
}
#endif
