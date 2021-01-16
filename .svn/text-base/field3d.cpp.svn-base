/*
 * field3d.cpp
 *
 *  Created on: Nov 10, 2010
 *      Author: bzhao
 */

#include "field3d.h"
#include "voxel.h"

int Field3D::
	ThreeRingNeighborVisible(const Voxel &voxel, const float *d0, int i0, int j0, int k0, int index,
			const Point &pi, float &fa) const{
	int nvis = 0;
	Point pc(0.f);
	double o, q;
	if(index == 0){
	    pc = voxel.VoxelCenterPosition(i0+1,j0+1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0+1,k0+1)];
		}
	}
	else if(index == 1){
		pc = voxel.VoxelCenterPosition(i0+1,j0-1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0-1,k0+1)];
		}
	}
	else if(index == 2){
		pc = voxel.VoxelCenterPosition(i0-1,j0+1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0+1,k0+1)];
		}
	}
	else if(index == 3){
		pc = voxel.VoxelCenterPosition(i0-1,j0-1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0-1,k0+1)];
		}
	}
	else if(index == 4){
		pc = voxel.VoxelCenterPosition(i0+1,j0+1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0+1,k0-1)];
		}
	}
	else if(index == 5){
		pc = voxel.VoxelCenterPosition(i0+1,j0-1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0-1,k0-1)];
		}
	}
	else if(index == 6){
		pc = voxel.VoxelCenterPosition(i0-1,j0+1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0+1,k0-1)];
		}
	}
	else if(index == 7){
		pc = voxel.VoxelCenterPosition(i0-1,j0-1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0-1,k0-1)];
		}
	}

	return nvis;
}

int Field3D::
	TwoRingNeighborVisible(const Voxel &voxel, const float *d0, int i0, int j0, int k0, int index,
			const Point &pi, float &fa) const{
	int nvis = 0;
	Point pc(0.f);
	double o, q;
	if(index == 0){
	    pc = voxel.VoxelCenterPosition(i0+1,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0+1,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0+1,k0)];
		}
	}
	else if(index == 1){
		pc = voxel.VoxelCenterPosition(i0+1,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0+1,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0-1,k0)];
		}
	}
	else if(index == 2){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0-1,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0+1,k0)];
		}
	}
	else if(index == 3){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0+1)];
		}
		pc = voxel.VoxelCenterPosition(i0-1,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0-1,k0)];
		}
	}
	else if(index == 4){
		pc = voxel.VoxelCenterPosition(i0+1,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0+1,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0+1,k0)];
		}
	}
	else if(index == 5){
		pc = voxel.VoxelCenterPosition(i0+1,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0+1,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0-1,k0)];
		}
	}
	else if(index == 6){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0-1,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0+1,k0)];
		}
	}
	else if(index == 7){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0-1)];
		}
		pc = voxel.VoxelCenterPosition(i0-1,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0-1,k0)];
		}
	}

	return nvis;
}

int Field3D::
	OneRingNeighborVisible(const Voxel &voxel, const float *d0, int i0, int j0, int k0, int index,
			const Point &pi, float &fa) const{
	int nvis = 0;
	Point pc(0.f);
	double o, q;
	if(index == 0){
	    pc = voxel.VoxelCenterPosition(i0+1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0+1)];
		}
	}
	else if(index == 1){
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0+1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0+1)];
		}
	}
	else if(index == 2){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0+1)];
		}
	}
	else if(index == 3){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0+1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0+1)];
		}
	}
	else if(index == 4){
		pc = voxel.VoxelCenterPosition(i0+1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0-1)];
		}
	}
	else if(index == 5){
		pc = voxel.VoxelCenterPosition(i0+1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0+1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0-1)];
		}
	}
	else if(index == 6){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0+1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0+1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0-1)];
		}
	}
	else if(index == 7){
		pc = voxel.VoxelCenterPosition(i0-1,j0,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0-1,j0,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0-1,k0);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0-1,k0)];
		}
		pc = voxel.VoxelCenterPosition(i0,j0,k0-1);
		if(!mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
			++nvis;
			fa += d0[INDEX(i0,j0,k0-1)];
		}
	}

	return nvis;
}

void Field3D::
	ProcessInvisibleNeightborForLevelsetAdvection(const Voxel &voxel, const float *d0,
			int i0, int j0, int k0, int index, double o,
			const Point &pc, const Point &pi, float &f) const{
	float faccum = 0.f;
	int nvis = 0;
	nvis = OneRingNeighborVisible(voxel, d0, i0, j0, k0, index, pi, faccum);
	if(nvis > 0)
		f = faccum / nvis;
	else{
		nvis = TwoRingNeighborVisible(voxel, d0, i0, j0, k0, index, pi, faccum);
		if(nvis > 0)
			f = faccum / nvis;
		else{
			nvis = ThreeRingNeighborVisible(voxel, d0, i0, j0, k0, index, pi, faccum);
			if(nvis > 0)
				f = faccum / nvis;
			else{ // use positive distance to object along this line segment
				Vector pci = o * (pc-pi);
				f = pci.Length();
			}
		}
	}
}
