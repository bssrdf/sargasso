#ifndef SOURCE_H_
#define SOURCE_H_

#include "voxel.h"
#include "transform.h"

class Source{
public:
	Source(int i, int j, int k):
	I(i), J(j), K(k){}
	virtual ~Source(){}
	virtual void Merge(const Voxel &voxel, float *phi, bool initial) const = 0;
//	virtual void Initialize(const Voxel &voxel, float *phi) = 0;
	virtual void UpdateVoxel(Voxel &voxel) const = 0;
	virtual void RestoreVoxel(Voxel &voxel) const = 0;
	virtual void SetVelocity(const Voxel &voxel, float *u, float *v, float *w) const = 0;
	virtual void Destroy(Voxel &voxel) = 0;
protected:
	int I, J, K;
};


#endif /*SOURCE_H_*/
