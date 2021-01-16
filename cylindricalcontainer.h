/*
 * cylindricalcontainer.h
 *
 *  Created on: Sep 25, 2008
 *      Author: bzhao
 */

#ifndef CYLINDRICALCONTAINER_H_
#define CYLINDRICALCONTAINER_H_

#include "container.h"

class CylinderContainer : public Container{

public:
	CylinderContainer(const Voxel &voxel, int X, int Y, int Z,
			float f,
			const Point &p, const Vector &dir,
        	 float r, float _h,
        	 const char *filename = NULL)
	  :  Container(X, Y, Z, f),
	     pos(p), ydir(dir), radius(r), h(_h){

		CoordinateSystem(ydir, &xdir, &zdir);
		trans = Translate(Vector(-pos.x, -pos.y, -pos.z));
		Matrix4x4 m(xdir.x, xdir.y, xdir.z, 0,
					ydir.x, ydir.y, ydir.z, 0,
					zdir.x, zdir.y, zdir.z, 0,
					0, 0, 0, 1);
		trans = Transform(m) * trans;
		trans_inv = trans.GetInverse();

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
	~CylinderContainer(){}

	bool IsSignedDistance() const{
		return mHasSignedDistanceMap;
	}
	void SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const;
	bool IsBorderCell(const Voxel &voxel, int i, int j, int k) const;
	void SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const;
	void SetNormal1(const Voxel &voxel, Vector &N, int i, int j, int k) const;

private:
	void EvaluatePhi(const Voxel &voxel, int axis, float *phitmp);
	float EvaluatePhi(const Point &p) const;
	float TriInterp(const Voxel &voxel, const Point &p) const;
	Point pos;
	Vector xdir, ydir, zdir;
	float radius; // radius of cylinder
	float h;
	Transform trans, trans_inv;

};

#endif /* CYLINDRICALCONTAINER_H_ */
