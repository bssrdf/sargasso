#ifndef SPHEREMOVINGOBJ_H_
#define SPHEREMOVINGOBJ_H_

#include "geometry.h"
#include "movingobj.h"

class SphereMovingObject : public MovingObject{

public:
	SphereMovingObject(const Voxel &voxel,
			           int nx, int ny, int nz,
			           float x, float y, float z,
			           float u, float v, float w,
			           float r,  TimeManager *tg,
			           float _h, float f,
			           int i, int j, int k)
	           : MovingObject(nx, ny, nz, _h, f, tg, i, j, k),
	           px(x), py(y), pz(z),
	           us(u), vs(v), ws(w),
	           rad(r){
		phi = new float[nx*ny*nz];
		memset(phi, 0, nx*ny*nz);
		Point p(x, y, z);
		SetObjectLevelSet(voxel, p, r, 0, phi);
	}
	~SphereMovingObject(){
		if(phi)
			delete [] phi;
	}

	void Update(Voxel &voxel, float *phic, float *phiu, float *phiv, float *phiw);
	void SetConstraintVelocity(Vector &vel, int i, int j, int k) const;
	void SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k) const;
	void SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k, int vel_index) const;
	void SetObjectLevelSet(const Voxel &voxel, const Point &p, float r, char ax, float *phi);
	void SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const;
	void Friction(const Voxel &voxel, int vel_index, int i, int j, int k, float *f) const;
	bool IsBorderCell(const Voxel &voxel, int i, int j, int k) const;
	bool IsDynamicBorderCell(const Voxel &voxel, int i, int j, int k) const;
	void SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const;
	bool StopNow() const{
		if(fabsf(us) < E_EPSIL && fabsf(vs) < E_EPSIL && fabsf(ws) < E_EPSIL)
			return true;
		else
			return false;
	}
private:
	float EvaluatePhi(const Point &q) const;
	float px, py, pz;
	float us, vs, ws;
	float rad;
//	float dt;
};

#endif /*SPHEREMOVINGOBJ_H_*/
