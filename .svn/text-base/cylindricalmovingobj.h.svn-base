/*
 * cylindricalmovingobj.h
 *
 *  Created on: Sep 25, 2008
 *      Author: bzhao
 */

#ifndef CYLINDRICALMOVINGOBJ_H_
#define CYLINDRICALMOVINGOBJ_H_

#include "geometry.h"
#include "movingobj.h"

class CylinderMovingObject : public MovingObject{

public:
	CylinderMovingObject(const Voxel &voxel,
			           int nx, int ny, int nz, float f,
			           const Point &p, const Vector &dir,
			           float ye,
			           float u, float v, float w,
			           float r, TimeManager *tg,
			           float _h,
			           int i, int j, int k)
	           : MovingObject(nx, ny, nz, _h, f, tg, i, j, k),
	           pos(p), ydir(dir), y_extent(ye), radius(r), us(u), vs(v), ws(w)
	           {

		if(u != 0.f || v != 0.f || w != 0.f){
			printf("This cylinder moving object serves as a container, \n so set its velocity to zero \n");
			exit(1);
		}
		CoordinateSystem(ydir, &xdir, &zdir);
		trans = Translate(Vector(-pos.x, -pos.y, -pos.z));
		Matrix4x4 m(xdir.x, xdir.y, xdir.z, 0,
					ydir.x, ydir.y, ydir.z, 0,
					zdir.x, zdir.y, zdir.z, 0,
					0, 0, 0, 1);
		trans = Transform(m) * trans;
		trans_inv = trans.GetInverse();
		phi = new float[nx*ny*nz];
		memset(phi, 0, (nx*ny*nz)*sizeof(float));
	}

	~CylinderMovingObject(){
		if(phi)
			delete [] phi;
	}

	void Update(Voxel &voxel, float *phic, float *phiu, float *phiv, float *phiw);
	void SetConstraintVelocity(Vector &vel, int i, int j, int k) const;
	void SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k) const;
	void SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k, int vel_index) const;
	void SetObjectLevelSet(const Voxel &voxel, int ax);
	void SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const;
	void Friction(const Voxel &voxel, int vel_index, int i, int j, int k, float *f) const;
	bool StopNow() const{
		if(us == 0.f && vs == 0.f && ws == 0.f)
			return true;
		else
			return false;
	}
private:
	Point pos;
	Vector xdir, ydir, zdir;
	float y_extent;
	float radius; // radius of cylinder
	float us, vs, ws;
//	float dt;
	Transform trans_inv;
};

#endif /* CYLINDRICALMOVINGOBJ_H_ */
