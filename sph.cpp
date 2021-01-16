/*
 * sph.cpp
 *
 *  Created on: Sep 5, 2008
 *      Author: bzhao
 */
#include <queue>
using std::queue;

#include "sph.h"
#include "kdtree.h"
#include "pcg_solver.h"
#include "vectorfield3d.h"

#ifdef COUPLED_SPH

struct FloatPair{
	FloatPair(){
		nFoundParticles = 0;
		smallestR = INFINITY;
		value = INFINITY;
	}
	mutable vector<WaterParticle> foundParticles;
	mutable u_int nFoundParticles;
	mutable float smallestR, value;
	void operator()(const WaterParticle &p,
					 float distSquare,
					 float &maxDistSquared) const{
		foundParticles.push_back(p);
		nFoundParticles++;
		float dist = sqrt(distSquare);
		if( dist <= smallestR){
			smallestR = dist;
			value = smallestR;
		}
		return;
	}
	bool Empty() const{
		return foundParticles.empty();
	}
	u_int Size() const {
		return (unsigned int)foundParticles.size();
	}
	void Clear(){
		nFoundParticles = 0;
		smallestR = INFINITY;
		value = INFINITY;
		foundParticles.clear();
	}
};

void SPH::AddWaterParticles(const WaterParticle &p){
	waterParticles.push_back(p);
}

void SPH::TimeStepping(const Voxel *voxel, const VectorField3D *vel){
	float dt = tm->GetDt();
	printf(" There are %u SPH particles \n", waterParticles.size());
	if(waterParticles.size() == 0)
		return;
	else{
		ApplyForce();
		ComputeWeights(voxel);
		vel->ProjectionSPH(voxel, this);
//		ApplyParticleSlip(voxel);
		vel->AdvectSPHParticles(voxel, dt, waterParticles);
	}
}


void SPH::FloodFill(int ii, int jj, int kk, char *color) const {
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
				if(color[INDEX(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk)] == 0){
					color[INDEX(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk)] = 1;
					TentativeGridPoint tmp(0.f, neighbors[m].ii, neighbors[m].jj, neighbors[m].kk);
					Q.push(tmp);
				}
			}
			Q.pop();
//			printf("there are %lu iterms in queues \n", Q.size());
		}while(!Q.empty());
	}
}

void SPH::Projection(const Voxel *voxel, const float *upls, const float *vpls, const float *wpls){

	float delta = voxel->VoxelDelta();
	float delta2 = delta * delta;
	float deltar = 1.f / delta;
	float delta3 = delta * delta2;
	float dt = tm->GetDt();
	float apr = 1.f / averagePeriod;
	float dtdx = dt / delta;
#ifdef WIN32
	hash_map<u_int, float> rhs;
	hash_map<u_int, u_int> A;
	hash_map<u_int, u_int>::iterator found;
#else
	unordered_map<u_int, float> rhs;
	unordered_map<u_int, u_int> A;
	unordered_map<u_int, u_int>::iterator found;
#endif
	float *div = new float[DimX*DimY*DimZ];
	SetZero(div);

	char *color = new char[DimX*DimY*DimZ];
	SetZero(color);
	FOR_EACH_CELL
		if(voxel->InSolid(i,j,k))
			color[INDEX(i,j,k)] = 2;
		if(voxel->InSource(i,j,k))
			color[INDEX(i,j,k)] = 2;
		if(voxel->InLiquid(i,j,k) && !voxel->InSurface(i,j,k))
			color[INDEX(i,j,k)] = 2;
	END_FOR
	TentativeGridPoint floodSeed(0.f, 35, 35, 285);
	FloodFill(floodSeed.ii, floodSeed.jj, floodSeed.kk, color);
	u_int air_pockets = 0, air_points = 0, others = 0;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(color[pos] == 0){
			++air_pockets;
//			printf("at (%d, %d, %d) grid point is in air pocket\n", i,j,k);
		}
		if(color[pos] == 1)
			++air_points;
		if(color[pos] == 2)
			++others;
	END_FOR
	printf("\nthere are %u points in air pocket and %u points in air \n\n", air_pockets, air_points);

	u_int N = 0;
	float max_rhs = 0.f;
	int max_rhs_i = 0, max_rhs_j = 0,max_rhs_k = 0;
	float dxdt  = delta / dt;
	float dx2dt = delta * delta / dt;
	// compute the r.h.s. of the pressure equation
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if( (voxel->InAir(i,j,k) || voxel->InSurface(i,j,k)) && color[pos] != 0 &&
			cWeights[pos] > 0.15f * incompressibilityNumberDensity ){
			float uvel1, uvel0, vvel1, vvel0, wvel1, wvel0;
			if( (voxel->InLiquid(i+1,j,k) && !voxel->InSurface(i+1,j,k)) ||
				 voxel->InSolid (i+1,j,k) ||  voxel->InSource (i+1,j,k) )
				uvel1 = upls[INDEX(i,j,k)];
			else
				uvel1 = u0[INDEX(i,j,k)];
			if( (voxel->InLiquid(i-1,j,k) && !voxel->InSurface(i-1,j,k)) ||
				 voxel->InSolid (i-1,j,k) ||  voxel->InSource (i-1,j,k) )
				uvel0 = upls[INDEX(i-1,j,k)];
			else
				uvel0 = u0[INDEX(i-1,j,k)];
			if( (voxel->InLiquid(i,j+1,k) && !voxel->InSurface(i,j+1,k)) ||
				 voxel->InSolid (i,j+1,k) ||  voxel->InSource (i,j+1,k) )
				vvel1 = vpls[INDEX(i,j,k)];
			else
				vvel1 = v0[INDEX(i,j,k)];
			if( (voxel->InLiquid(i,j-1,k) && !voxel->InSurface(i,j-1,k)) ||
				 voxel->InSolid (i,j-1,k) ||  voxel->InSource (i,j-1,k) )
				vvel0 = vpls[INDEX(i,j-1,k)];
			else
				vvel0 = v0[INDEX(i,j-1,k)];
			if( (voxel->InLiquid(i,j,k+1) && !voxel->InSurface(i,j,k+1)) ||
				 voxel->InSolid (i,j,k+1) ||  voxel->InSource (i,j,k+1) )
				wvel1 = wpls[INDEX(i,j,k)];
			else
				wvel1 = w0[INDEX(i,j,k)];
			if( (voxel->InLiquid(i,j,k-1) && !voxel->InSurface(i,j,k-1)) ||
				 voxel->InSolid (i,j,k-1) ||  voxel->InSource (i,j,k-1) )
				wvel0 = wpls[INDEX(i,j,k-1)];
			else
				wvel0 = w0[INDEX(i,j,k-1)];

			div[pos] = -(uvel1 - uvel0 + vvel1 - vvel0 + wvel1 - wvel0)*dxdt;
			div[pos] -= apr * dx2dt * logf(cWeights[pos] / incompressibilityNumberDensity);
		   if(fabs(max_rhs) < fabs(div[pos])){
			   max_rhs = div[pos];
			   max_rhs_i = i; max_rhs_j = j; max_rhs_k = k;
		   }
		   rhs.insert( make_pair(pos, div[pos]) );
		   A.insert( make_pair(pos, N) );
		   ++N;
		}

	END_FOR
	printf("Before SPH Projection: there are %u points \n", N);

	if(N == 0){
		delete [] div;
		delete [] color;
		return;
	}

	printf("\nmaximum rhs (%f) occurs at (%d, %d, %d) \n", max_rhs,
				max_rhs_i,max_rhs_j,max_rhs_k);

	float *p = new float[N];
//	memset(p, 0, N);
	for(int n=0; n<N; ++n)
		p[n] = 0.f;
    printf("pressure solver stage 1, N = \n", N);
	LinearSolver(voxel, p, A, rhs, N);
	SetZero(div);
	for(found = A.begin(); found != A.end(); ++found)
		div[found->first] = p[found->second];
	delete [] p;

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if( (voxel->InAir(i,j,k) || voxel->InSurface(i,j,k)) &&
			cWeights[pos] > 0.15f * incompressibilityNumberDensity ){
			if( (voxel->InLiquid(i+1,j,k) && !voxel->InSurface(i+1,j,k)) ||
				 voxel->InSolid (i+1,j,k) ||  voxel->InSource (i+1,j,k) )
				u0[pos] = upls[pos] - u0[pos];
			else
				u0[pos] =  dtdx * (div[INDEX(i+1,j,k)] - div[pos]);
			if( (voxel->InLiquid(i-1,j,k) && !voxel->InSurface(i-1,j,k)) ||
				 voxel->InSolid (i-1,j,k) ||  voxel->InSource (i-1,j,k) )
				u0[INDEX(i-1,j,k)] = upls[INDEX(i-1,j,k)] - u0[INDEX(i-1,j,k)];
			else
				u0[INDEX(i-1,j,k)] = dtdx * (div[pos] - div[INDEX(i-1,j,k)]);
			if( (voxel->InLiquid(i,j+1,k) && !voxel->InSurface(i,j+1,k)) ||
				 voxel->InSolid (i,j+1,k) ||  voxel->InSource (i,j+1,k) )
				v0[pos] = vpls[pos] - v0[pos];
			else
				v0[pos] =  dtdx * (div[INDEX(i,j+1,k)] - div[pos]);
			if( (voxel->InLiquid(i,j-1,k) && !voxel->InSurface(i,j-1,k)) ||
				 voxel->InSolid (i,j-1,k) ||  voxel->InSource (i,j-1,k) )
				v0[INDEX(i,j-1,k)] = vpls[INDEX(i,j-1,k)] - v0[INDEX(i,j-1,k)];
			else
				v0[INDEX(i,j-1,k)] = dtdx * (div[pos] - div[INDEX(i,j-1,k)]);
			if( (voxel->InLiquid(i,j,k+1) && !voxel->InSurface(i,j,k+1)) ||
				 voxel->InSolid (i,j,k+1) ||  voxel->InSource (i,j,k+1) )
				w0[pos] = wpls[pos] - w0[pos];
			else
				w0[pos] =  dtdx * (div[INDEX(i,j,k+1)] - div[pos]);
			if( (voxel->InLiquid(i,j,k-1) && !voxel->InSurface(i,j,k-1)) ||
				 voxel->InSolid (i,j,k-1) ||  voxel->InSource (i,j,k-1) )
				w0[INDEX(i,j,k-1)] = wpls[INDEX(i,j,k-1)] - w0[INDEX(i,j,k-1)];
			else
				w0[INDEX(i,j,k-1)] = dtdx * (div[pos] - div[INDEX(i,j,k-1)]);
		}
	END_FOR

	delete [] div;
	delete [] color;

	ApplyParticleSlip(voxel);

}

void SPH::LinearSolver( const Voxel *voxel, float *x,
#ifdef WIN32
					 hash_map<u_int, u_int> &a, hash_map<u_int, float> &x0,
#else
					 unordered_map<u_int, u_int> &a, unordered_map<u_int, float> &x0,
#endif
					 u_int n ){

	float delta = voxel->VoxelDelta();
	float delta2 = delta * delta;
	float delta3 = delta * delta2;
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
	printf("pressure solver stage 2, n = %u \n", n);
	double residual;
	int iterations;
	SparseMatrixd A(n);
	vector<double> rhs;
	vector<double> p;
	for(u_int i=0; i<n; ++i)
		p.push_back((double)x[i]);
	u_int diag;
	printf("pressure solver stage 3 \n");
	FOR_EACH_CELL
		 u_int index = INDEX(i,j,k);
		 found = a.find(index);
		 if(found != a.end()){
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back((double)found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 int solidNeighbors=0;
			 for(u_int m=0; m<neighbors.size(); m++){
					 if(voxel->InSolid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
						 ++solidNeighbors;
					 }
					 else if(voxel->InSource(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
						 ++solidNeighbors;
					 }

					 else if(voxel->InAir(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ||
							 voxel->InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
						 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
						 if(found != a.end())
							 A.set_element(diag,found->second,-1);
					 }
					 else if(voxel->InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) &&
							!voxel->InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)
							){
						 ++solidNeighbors;
					 }
			 }

			 if(solidNeighbors == 6){
				printf("at (%d, %d, %d) Newmann b.c. for all directions \n", i,j,k);
				exit(1);
			 }
			 A.set_element(diag,diag,(double)(neighbors.size()-solidNeighbors));
		}
	END_FOR
	printf("pressure solver stage 4 \n");
	//A.write_matlab(cout,"A");
//	printf("\nmatrix A is ready \n\n ");
	A.check_symmetry();
	printf("pressure solver stage 5 \n");
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
	printf("pressure solver stage 6 \n");
	solver.solve(A, rhs, p, residual, iterations);
	printf("pressure solver stage 7 \n");
	for(u_int i=0; i<n; ++i)
		x[i] = (float)p[i];
	printf("pressure solver stage 8 \n");
	printf("\nPCG iterations = %d, and residual = %16.13f \n\n ",
			iterations, residual);


}


void SPH::ComputeWeights(const Voxel *voxel){

	list<WaterParticle>::iterator iter_particle;
	vector<WaterParticle> wp;
	list<WaterParticle> &sphParticles = waterParticles;
	for(iter_particle = sphParticles.begin();
		iter_particle != sphParticles.end();
		){
			WaterParticle &wps = *iter_particle;
			wp.push_back(wps);
			++iter_particle;
	}
	KdTree<WaterParticle, FloatPair> dist(wp);
	FloatPair fp;
	SetZero(cWeights);
	SetZero(uWeights);
	SetZero(vWeights);
	SetZero(wWeights);
	SetZero(u0);
	SetZero(v0);
	SetZero(w0);
	FOR_EACH_CELL
		if( voxel->InAir(i,j,k) || voxel->InSurface(i,j,k) ){
			Point p = voxel->VoxelCenterPosition(i,j,k);
			dist.Lookup(p, fp, radiusSquared);
		    float totalWeight = 0.f;
		    u_int N = fp.nFoundParticles;
		    if(N > 0){
		    	for(int n=0; n < N; ++n){
					WaterParticle &wps = fp.foundParticles[n];
					totalWeight += wps.Influence(p, influenceRadius);
				}
		    }
		    cWeights[INDEX(i,j,k)] = totalWeight;
		    fp.Clear();
		    if(i == I && j == J && k == K)
		    	printf("\nat (%d, %d, %d), cWeights = %f \n", i,j,k,cWeights[INDEX(i,j,k)]);

		    p = voxel->VelPosition(1,i,j,k);
			dist.Lookup(p, fp, radiusSquared);
			totalWeight = 0.f;
			float totalVel = 0.f;
			N = fp.nFoundParticles;
			if(N > 0){
				for(int n=0; n < N; ++n){
					WaterParticle &wps = fp.foundParticles[n];
					float uvel = wps.UVelocity();
					float weight = wps.Influence(p, influenceRadius);
					totalWeight += weight;
					totalVel += weight * uvel;
				}
			}
			uWeights[INDEX(i,j,k)] = totalWeight;
			if(totalWeight > 0.f)
				u0[INDEX(i,j,k)] = totalVel / totalWeight;
			else
				u0[INDEX(i,j,k)] = 0.f;
			fp.Clear();

			p = voxel->VelPosition(2,i,j,k);
			dist.Lookup(p, fp, radiusSquared);
			totalWeight = 0.f;
			totalVel = 0.f;
			N = fp.nFoundParticles;
			if(N > 0){
				for(int n=0; n < N; ++n){
					WaterParticle &wps = fp.foundParticles[n];
					float vvel = wps.VVelocity();
					float weight = wps.Influence(p, influenceRadius);
					totalWeight += weight;
					totalVel += weight * vvel;
				}
			}
			vWeights[INDEX(i,j,k)] = totalWeight;
			if(totalWeight > 0.f)
				v0[INDEX(i,j,k)] = totalVel / totalWeight;
			else
				v0[INDEX(i,j,k)] = 0.f;
			fp.Clear();

			p = voxel->VelPosition(3,i,j,k);
			dist.Lookup(p, fp, radiusSquared);
			totalWeight = 0.f;
			totalVel = 0.f;
			N = fp.nFoundParticles;
			if(N > 0){
				for(int n=0; n < N; ++n){
					WaterParticle &wps = fp.foundParticles[n];
					float wvel = wps.WVelocity();
					float weight = wps.Influence(p, influenceRadius);
					totalWeight += weight;
					totalVel += weight * wvel;
				}
			}
			wWeights[INDEX(i,j,k)] = totalWeight;
			if(totalWeight > 0.f)
				w0[INDEX(i,j,k)] = totalVel / totalWeight;
			else
				w0[INDEX(i,j,k)] = 0.f;
			fp.Clear();
		}
	END_FOR
}

void SPH::ComputeTargetDensity(const Voxel *voxel, int nPerVoxel) {
	list<Particle> testParticles;
	list<Particle>::iterator iter_particle;
	for(int m=0; m<nPerVoxel; ++m){
		Particle ps;
		voxel->PlaceParticles(I, J, K, ps, 0);
		Point pos = ps.Position();
		printf(" seeded particle %d at (%f, %f, %f) \n ", m, pos.x, pos.y, pos.z);
		testParticles.push_back(ps);
	}

	Point p = voxel->VoxelCenterPosition(I,J,K);
    float totalWeight = 0.f;
    printf(" Point p at (%f, %f, %f) \n", p.x, p.y, p.z);
    for(iter_particle  = testParticles.begin();
    	iter_particle != testParticles.end();
    	){
    	Particle &wps = *iter_particle;
		totalWeight += wps.Influence(p, influenceRadius);
//		printf(" totalWeight = %f \n", totalWeight);
		++iter_particle;
    }
    incompressibilityNumberDensity = totalWeight;
    printf("\n SPH incompressibilityNumberDensity = %f \n\n", incompressibilityNumberDensity);


}

void SPH::ApplyForce(){
	Vector G(0, 0, -g);
	float dt = tm->GetDt();
	list<WaterParticle>::iterator iter_particle;
	for(iter_particle = waterParticles.begin();
		iter_particle != waterParticles.end();
		){
		WaterParticle &wps = *iter_particle;
		Vector velocity = wps.Velocity();
	    velocity += G * dt;
		wps.SetVelocity(velocity);
		++iter_particle;
	}
}

void SPH::VelocityChangeFLIP(const Voxel *voxel, const Point &pos,
				float &uvel, float &vvel, float &wvel) const{

	TentativeGridPoint tp = voxel->ContainsPoint(pos);
	float f = (uWeights[INDEX(tp.ii,tp.jj,tp.kk)] + uWeights[INDEX(tp.ii-1,tp.jj,tp.kk)]);
	if( f > 0.f )
		uvel = (uWeights[INDEX(tp.ii,  tp.jj,tp.kk)]*u0[INDEX(tp.ii,  tp.jj,tp.kk)] +
			    uWeights[INDEX(tp.ii-1,tp.jj,tp.kk)]*u0[INDEX(tp.ii-1,tp.jj,tp.kk)]) / f;
	else
		uvel = 0.f;
	f = (vWeights[INDEX(tp.ii,tp.jj,tp.kk)] + vWeights[INDEX(tp.ii,tp.jj-1,tp.kk)]);
	if( f > 0.f )
		vvel =  (vWeights[INDEX(tp.ii,tp.jj,  tp.kk)]*v0[INDEX(tp.ii,tp.jj,  tp.kk)] +
				 vWeights[INDEX(tp.ii,tp.jj-1,tp.kk)]*v0[INDEX(tp.ii,tp.jj-1,tp.kk)]) / f;
	else
		vvel = 0.f;
	f = (wWeights[INDEX(tp.ii,tp.jj,tp.kk)] + wWeights[INDEX(tp.ii,tp.jj,tp.kk-1)]);
	if( f > 0.f )
		wvel = (wWeights[INDEX(tp.ii,tp.jj,tp.kk  )]*w0[INDEX(tp.ii,tp.jj,tp.kk  )] +
				wWeights[INDEX(tp.ii,tp.jj,tp.kk-1)]*w0[INDEX(tp.ii,tp.jj,tp.kk-1)]) / f;
	else
		wvel = 0.f;

}


void SPH::ApplyParticleSlip(const Voxel *voxel){
	list<WaterParticle>::iterator iter_particle;
	vector<WaterParticle> wp;
	list<WaterParticle> &sphParticles = waterParticles;
	for(iter_particle = sphParticles.begin();
		iter_particle != sphParticles.end();
		){
			WaterParticle &wps = *iter_particle;
			wp.push_back(wps);
			++iter_particle;
	}
	KdTree<WaterParticle, FloatPair> dist(wp);
	FloatPair fp;

	for(iter_particle = waterParticles.begin();
		iter_particle != waterParticles.end();
		){
		WaterParticle &wps = *iter_particle;
		Vector velocity = wps.Velocity();
	    Point pos = wps.Position();
	    dist.Lookup(pos, fp, radiusSquared);
	    u_int N = fp.nFoundParticles;
	    float totalWeight = 0.f;
		if(N > 0){
			for(int n=0; n < N; ++n){
				WaterParticle &wps = fp.foundParticles[n];
				if(wps.DistanceTo(pos) > 1.e-6f) // ignore the particel itself
					totalWeight += wps.Influence(pos, influenceRadius);
			}
		}
		float slip = totalWeight / incompressibilityNumberDensity;
		float uc = 0.f, vc = 0.f, wc = 0.f;
		VelocityChangeFLIP(voxel, pos, uc, vc, wc);
		velocity += slip * Vector(uc, vc, wc);
		wps.SetVelocity(velocity);
		++iter_particle;
	}
}

void SPH::OutputWaterParticlesBinary(const Voxel *voxel, const char *name){

	float delta = voxel->VoxelDelta();

	list<WaterParticle>::iterator iter_particle;


	FILE *fp = fopen(name, "wb");
	if(!fp){
      printf("Couldn't open file [ %s ] to write \n", name);
      exit(-1);
   	}
	u_int num = waterParticles.size();
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
	for(iter_particle = waterParticles.begin();
		iter_particle != waterParticles.end();
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
