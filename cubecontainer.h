/*
 * cubecontainer.h
 *
 *  Created on: Sep 25, 2008
 *      Author: bzhao
 */

#ifndef CUBECONTAINER_H_
#define CUBECONTAINER_H_

#include "container.h"

class CubeContainer : public Container{

public:
	CubeContainer(const Voxel &voxel, int X, int Y, int Z, float f,
			int s_x0, int s_y0, int s_z0,
			int e_x0, int e_y0, int e_z0,
			const char *filename = NULL)
	: Container(X, Y, Z, f),
	start_x0(s_x0), start_y0(s_y0),  start_z0(s_z0),
	end_x0(e_x0),  end_y0(e_y0), end_z0(e_z0){

		h          = voxel.VoxelDelta();
		leftWall   = voxel.FacePostion(1, false, start_x0, start_y0, start_z0);
		rightWall  = voxel.FacePostion(1, true, end_x0, start_y0, start_z0);
		backWall   = voxel.FacePostion(2, false, start_x0, start_y0, start_z0);
		frontWall  = voxel.FacePostion(2, true, start_x0, end_y0, start_z0);
		bottomWall = voxel.FacePostion(3, false, start_x0, start_y0, start_z0);
		topWall    = voxel.FacePostion(3, true, start_x0, start_y0, end_z0);

		if(filename){
			LoadPhi(filename);
			mHasSignedDistanceMap = true;
		}
		else{
			EvaluatePhi(voxel, 0, phi);
			EvaluatePhi(voxel, 1, phiu);
			EvaluatePhi(voxel, 2, phiv);
			EvaluatePhi(voxel, 3, phiw);
			mHasSignedDistanceMap = false;
		}
	}
	~CubeContainer(){}

	bool IsSignedDistance() const{
		return mHasSignedDistanceMap;
	}

	bool IsBorderCell(const Voxel &voxel, int i, int j, int k) const;
	void SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const;
	void SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const;

private:
	void EvaluatePhi(const Voxel &voxel, int axis, float *phitmp);
	float EvaluatePhi(const Point &p) const;
	int start_x0, start_y0,  start_z0;
	int end_x0,   end_y0,    end_z0;
	float h;
	float leftWall;
	float rightWall;
	float backWall;
	float frontWall;
	float bottomWall;
	float topWall;

};

#endif /* CUBECONTAINER_H_ */
