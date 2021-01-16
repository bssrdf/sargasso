#ifndef CUBESOURCE_H_
#define CUBESOURCE_H_
#include "geometry.h"
#include "source.h"

class CubeSource : public Source{
public:
	CubeSource(const Voxel &voxel, int _start_x0=0, int _start_y0=0, int _start_z0=0,
	int _end_x0=0, int _end_y0=0, int _end_z0=0,
	int X = 0, int Y = 0, int Z = 0, int i=0, int j=0, int k=0,
	float _speed = 0.f, float _h = 0.f,
	char _axis = 0, char _side = -1)
	: Source(i,j,k), start_x0(_start_x0), start_y0(_start_y0), start_z0(_start_z0),
	  end_x0(_end_x0), end_y0(_end_y0), end_z0(_end_z0),
	  DimX(X), DimY(Y), DimZ(Z),
	  speed(_speed), h(_h),
	  axis(_axis), side(_side){
		stop = false;
		dimx = end_x0 - start_x0 + 1;
		dimy = end_y0 - start_y0 + 1;
		dimz = end_z0 - start_z0 + 1;
		phi = new float[dimx*dimy*dimz];
		Initialize(voxel, phi);
	}
	~CubeSource(){
		delete [] phi;
	}

	void Merge(const Voxel &voxel, float *phi0, bool initial) const;
	void Initialize(const Voxel &voxel, float *f);
	void UpdateVoxel(Voxel &voxel) const;
	void RestoreVoxel(Voxel &voxel) const;
	void SetVelocity(const Voxel &voxel, float *u, float *v, float *w) const;
	void Destroy(Voxel &voxel);

private:
	int start_x0, start_y0, start_z0;
	int end_x0, end_y0, end_z0;
	float *phi;
	int DimX, DimY, DimZ;
	float speed;
	float h;
	int dimx, dimy, dimz;
	char axis, side; // axis = 1  ==>  x
					 // axis = 2  ==>  y
					 // axis = 3  ==>  z
	 	             // side = -1  ==> emits water at start_x0, start_y0 or start_z0
     				 // side = 1   ==> emits water at end_x0, end_y0 or end_z0
	bool stop;
};
#endif /*CUBESOURCE_H_*/
