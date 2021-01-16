#include "vortexparticles.h"
#include "kdtree.h"
#include "mtrand.h"

struct FloatPair{
	FloatPair(){
		nFoundParticles = 0;
		smallestR = INFINITY;
		value = INFINITY;
	}
	mutable vector<VortexParticle> foundParticles;
	mutable u_int nFoundParticles;
	mutable float smallestR, value;
	void operator()(const VortexParticle &p,
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

static float ComputeVortexForce(const Point &p, int vel_index, float radiusSquared, float s,
		FloatPair &fp, const KdTree<VortexParticle, FloatPair> &dist){
		
	Vector force(0.f, 0.f, 0.f);
	dist.Lookup(p, fp, radiusSquared);      
    u_int N = 0;
    if(fp.nFoundParticles > 0){
    	N = fp.Size();    	
		for(int n=0; n < N; ++n){
			VortexParticle &ps = fp.foundParticles[n];
			Point pos = ps.Position();
			Vector vor = ps.Vorticity();
			Vector Np = p - pos;
			Np = Normalize(Np);
			Vector omega = ps.Influence(p, s) * vor;
			force += Cross(Np, omega);
		}
	}	
    fp.Clear();
	if(vel_index == 1)
		return force.x;
	else if(vel_index == 2)
		return force.y;
	else if(vel_index == 3)
		return force.z;
}

void VortexParticles::VortexForce(const Voxel &voxel, 
								float *u0, float *v0, float *w0,
								float *phi_u, float *phi_v, float *phi_w,
								float *phi_u_obj, float *phi_v_obj, float *phi_w_obj) {
	
	if(Particles.size() == 0)
		return;
	list<VortexParticle>::iterator iter_particle;
	vector<VortexParticle> wp;
	
	for(iter_particle  = Particles.begin();
		iter_particle != Particles.end();
		){
		VortexParticle &wps = *iter_particle;
		wp.push_back(wps);
		++iter_particle;		
	}
	
	KdTree<VortexParticle, FloatPair> dist(wp);
//	printf("got here (1.3)\n");
	FloatPair fp;	
	float radiusSquared = support * support;
	
	FOR_EACH_CELL
		//if(!voxel.InSolid(i,j,k)){
		if(phi_u_obj[INDEX(i,j,k)] < 0.f){
			Point pos = voxel.VelPosition(1,i,j,k);
			if(phi_u[INDEX(i,j,k)] <= 0.f)
				u0[INDEX(i,j,k)] = ComputeVortexForce(pos, 1, radiusSquared, support, fp, dist);
			else if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
				u0[INDEX(i,j,k)] = ComputeVortexForce(pos, 1, radiusSquared, support, fp, dist);
			else if(voxel.InAir(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
				u0[INDEX(i,j,k)] = ComputeVortexForce(pos, 1, radiusSquared, support, fp, dist);
			else if(voxel.InSurface(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
				u0[INDEX(i,j,k)] = ComputeVortexForce(pos, 1, radiusSquared, support, fp, dist);
		}
		if(phi_v_obj[INDEX(i,j,k)] < 0.f){
			Point pos = voxel.VelPosition(2,i,j,k);
			if(phi_v[INDEX(i,j,k)] <= 0.f)
				v0[INDEX(i,j,k)] = ComputeVortexForce(pos, 2, radiusSquared, support, fp, dist);
			else if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
				v0[INDEX(i,j,k)] = ComputeVortexForce(pos, 2, radiusSquared, support, fp, dist);
			else if(voxel.InAir(i,j,k) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
				v0[INDEX(i,j,k)] = ComputeVortexForce(pos, 2, radiusSquared, support, fp, dist);
			else if(voxel.InSurface(i,j,k) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
				v0[INDEX(i,j,k)] = ComputeVortexForce(pos, 2, radiusSquared, support, fp, dist);
		}
		if(phi_w_obj[INDEX(i,j,k)] < 0.f){
			Point pos = voxel.VelPosition(3,i,j,k);
			if(phi_w[INDEX(i,j,k)] <= 0.f)
				w0[INDEX(i,j,k)] = ComputeVortexForce(pos, 3, radiusSquared, support, fp, dist);
			else if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
				w0[INDEX(i,j,k)] = ComputeVortexForce(pos, 3, radiusSquared, support, fp, dist);
			else if(voxel.InAir(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
				w0[INDEX(i,j,k)] = ComputeVortexForce(pos, 3, radiusSquared, support, fp, dist);
			else if(voxel.InSurface(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
				w0[INDEX(i,j,k)] = ComputeVortexForce(pos, 3, radiusSquared, support, fp, dist);
		}
		END_FOR

	
}

void VortexParticles::Seed(const Voxel &voxel, int i, int j, int k){
	
	float delta = voxel.VoxelDelta();
	
	for(int n=0; n<N; ++n){
		MTRand mt;
		double rn = mt();
		double rn1 = mt();
		double rn2 = mt();
		double rn3 = mt();
		Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
		voxel.PlaceParticles(i, j, k, ps, 0);
		Point pos = ps.Position();
		float cl = ps.CollisionDistance();
		Vector dir((float)rn1, (float)rn2, (float)rn);
		dir = Normalize(dir);
		VortexParticle wps(pos, 0, -1, cl, initVorticity, dir.x, dir.y, dir.z);
		AddVortexParticle(wps);
	}
	
}