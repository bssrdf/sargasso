/*
 * signeddistancecontainer.h
 *
 *  Created on: Feb 11, 2009
 *      Author: bzhao
 */

#ifndef SIGNEDDISTANCECONTAINER_H_
#define SIGNEDDISTANCECONTAINER_H_

#include "container.h"

class SignedDistanceContainer : public Container{

public:
	SignedDistanceContainer(const char *filename, Container *s, int X, int Y, int Z, float f)
	: Container(X, Y, Z, f), surrogate(s)
	{
		LoadPhi(filename);
	}
	~SignedDistanceContainer(){}

	bool IsSignedDistance() const{
		return true;
	}
	bool IsBorderCell(const Voxel &voxel, int i, int j, int k) const{
		return surrogate->IsBorderCell(voxel, i, j, k);
	}
	void SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const{
		surrogate->SetBorderCellPhi(voxel, i, j, k, n, f);
	}

	void SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const{
		surrogate->SetNormal(voxel, N, i, j, k);
	}
private:
//	void EvaluatePhi(char *file);
	Container *surrogate;
	int start_x0, start_y0,  start_z0;
	int end_x0,   end_y0,    end_z0;

};

#endif /* SIGNEDDISTANCECONTAINER_H_ */
