#ifndef CYLINDRICALSOURCE_H_
#define CYLINDRICALSOURCE_H_

#include "geometry.h"
#include "source.h"

class CylinderSource : public Source{
public:
	CylinderSource(const Voxel &voxel, const Point &p, const Vector &dir,
	int X, int Y, int Z, float r, float _h, int i, int j, int k,
	float _speed = 0.f) :
	Source(i, j, k), pos(p), ydir(dir), radius(r), h(_h), velocity(_speed), DimX(X), DimY(Y), DimZ(Z){
		stop = false;
//		if(fabs(ydir.x) <= fabs(ydir.y) && fabs(ydir.x) <= fabs(ydir.z))
//			xdir = Vector(0, -ydir.z, ydir.y);
//		if(fabs(ydir.y) <= fabs(ydir.x) && fabs(ydir.y) <= fabs(ydir.z))
//			xdir = Vector(-ydir.z, 0, ydir.x);
//		if(fabs(ydir.z) <= fabs(ydir.x) && fabs(ydir.z) <= fabs(ydir.y))
//			xdir = Vector(-ydir.y, ydir.x, 0);
//		xdir = Normalize(xdir);
//		zdir = Cross(xdir, ydir);
//		zdir = Normalize(zdir);
		CoordinateSystem(ydir, &xdir, &zdir);
		trans = Translate(Vector(-pos.x, -pos.y, -pos.z));
		Matrix4x4 m(xdir.x, xdir.y, xdir.z, 0,
					ydir.x, ydir.y, ydir.z, 0,
					zdir.x, zdir.y, zdir.z, 0,
		            0, 0, 0, 1);
		trans = Transform(m) * trans;
		trans_inv = trans.GetInverse();
	}

	~CylinderSource(){
	}

	void Merge(const Voxel &voxel, float *phi0, bool initial) const;
	void UpdateVoxel(Voxel &voxel) const;
	void RestoreVoxel(Voxel &voxel) const;
	void SetVelocity(const Voxel &voxel, float *u, float *v, float *w) const;
	void Destroy(Voxel &voxel);

private:
	Point pos;
	Vector xdir, ydir, zdir;
	float radius; // radius of cylinder
	float h;
	float velocity; // along -dir
	Transform trans, trans_inv;
	int DimX, DimY, DimZ;
	bool stop;
};

#endif /*CYLINDRICALSOURCE_H_*/
