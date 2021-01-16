#ifndef BLENDPARTICLES_H_
#define BLENDPARTICLES_H_
#include <vector>
using std::vector;
#include "geometry.h"
#include "voxel.h"
#include "kdtree.h"
#include "particle.h"


class BlendParticles{
public:
	BlendParticles(int nx, int ny, int nz, float th)
		: DimX(nx), DimY(ny), DimZ(nz), threshold(th){
		if(th  <= 0.f){
			printf("Error! threshold must be > 0\n");
			exit(1);
		}
		phi = new float[nx*ny*nz];
	}
	~BlendParticles(){
		waterParticles.clear();
		delete [] phi;
	}
	void AddWaterParticles(const WaterParticle &p);
	void EvaluatePhi(const Voxel &voxel);
//	float Kernel(const Point &p0, const Point &p1, float r);
	void Print() const{
		printf("\n there are %u water particles \n\n", waterParticles.size());
	}
	u_int Size() const{
		return waterParticles.size();
	}
	void OutputBinaryGridData(const Voxel &voxel, char *filename) const;
	void OutputWaterParticles(const Voxel &voxel, char *filename);
	void OutputWaterParticlesBinary(const Voxel &voxel, char *name);
	void ReInitialize(int steps, const Voxel &voxel, float *value);
	void MergeWith(const Voxel &voxel, float *phi0) const;
//	list<WaterParticle> waterParticles;
	WaterParticleList waterParticles;
	float *phi;
private:
	int DimX;
	int DimY;
	int DimZ;
	float threshold;  // must be negative
};
#endif /*BLENDPARTICLES_H_*/
