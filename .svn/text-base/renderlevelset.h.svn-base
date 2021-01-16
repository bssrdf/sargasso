/*
 * renderlevelset.h
 *
 *  Created on: Nov 11, 2010
 *      Author: bzhao
 */

#ifndef RENDERLEVELSET_H_
#define RENDERLEVELSET_H_

#include "field3d.h"
#include "scalarfield3d.h"
#include "voxel.h"


class RenderLevelSet : Field3D{

public:
	RenderLevelSet(PhysicsWorld *pw, int nx, int ny, int nz,
				  int i, int j, int k)
		: Field3D(nx,ny,nz,pw),
		  I(i), J(j), K(k){

		phi = new float[nx*ny*nz];
		SetZero(phi);

	}
	~RenderLevelSet(){
		if(phi){
			delete [] phi;
			phi = NULL;
		}
	}

	void CopyFrom(float *s){
		FOR_EACH_CELL
			u_int index = INDEX(i,j,k);
			phi[index] = s[index];
		END_FOR
	}

	void ResetObject(const Voxel &voxel, const ScalarField3D &obj){
		obj.Recompute(voxel, phi);
	}

	void OutputBinaryGridData(const Voxel &voxel, const char *filename) const;

private:

	void SetZero(float *data){
		FOR_EACH_CELL
			data[INDEX(i,j,k)] = 0.f;
		END_FOR
	}

	float *phi;
	int I, J, K;
};
#endif /* RENDERLEVELSET_H_ */
