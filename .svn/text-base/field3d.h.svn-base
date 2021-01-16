#ifndef FIELD3D_H_
#define FIELD3D_H_

#define FASTMARCH_LIMIT	10.f*delta
#define LS_ADV_LIMIT    5.f*delta
#define EXTRAPOLATE_VEL_LIMIT  10.f*delta
#define PARTICLE_SEED_LIMIT 3*delta

#ifdef AIR_DIV_FREE
#define AIR_DIV_FREE_LIMIT  5.f*delta
#endif


#define MONTE_CARLO_ITERATIONS 2048

#ifdef SPMD
#ifdef TBB
#define THREADS 2
#define CHUNK_SIZE 10
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
//using namespace tbb;
#endif
#endif

#include<vector>
//using namespace std;
#include "geometry.h"
#include "physicsworld.h"
#include "source.h"

class Voxel;

class Field3D{
public:
	Field3D(){
		DimX = 0;
		DimY = 0;
		DimZ = 0;
		mPhysWorld = NULL;
	}
	Field3D(int nx, int ny, int nz, PhysicsWorld *pw)
	: DimX(nx), DimY(ny), DimZ(nz), mPhysWorld(pw) {
	}

	virtual ~Field3D(){

	}
	void AddSourceObject(Source *obj){
		sources.push_back(obj);
	}
	Field3D(const Field3D &a){
		DimX = a.DimX;
		DimY = a.DimY;
		DimZ = a.DimZ;
		sources = a.sources;
		mPhysWorld = a.mPhysWorld;
	}

protected:
	vector<Source *> sources;
	PhysicsWorld *mPhysWorld;
	int ThreeRingNeighborVisible(const Voxel &voxel, const float *d0, int i0, int j0, int k0, int index,
				const Point &pi, float &fa) const;
	int TwoRingNeighborVisible(const Voxel &voxel, const float *d0, int i0, int j0, int k0, int index,
				const Point &pi, float &fa) const;
	int OneRingNeighborVisible(const Voxel &voxel, const float *d0, int i0, int j0, int k0, int index,
				const Point &pi, float &fa) const;
	void ProcessInvisibleNeightborForLevelsetAdvection(const Voxel &voxel,
				const float *d0, int i0, int j0, int k0, int index, double o,
				const Point &pc, const Point &pi, float &f) const;

	int DimX;
	int DimY;
	int DimZ;

private:

};



#endif /*FIELD3D_H_*/
