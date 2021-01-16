#ifndef MOVINGOBJ_H_
#define MOVINGOBJ_H_

#include <vector>
#include <map>
#include <algorithm>
//using namespace std;
using std::map;
#include "minheap.h"
#include "geometry.h"
#include "voxel.h"
#include "transform.h"
#include "timemanager.h"

#define POS(x) (x).kk*(DimX*DimY)+(x).jj*DimX+(x).ii

// abstract class for all moving objects
class MovingObject{
public:
	MovingObject(int nx, int ny, int nz, float _h, float f, TimeManager *tg, int i, int j, int k)
	: DimX(nx), DimY(ny), DimZ(nz),
	  h(_h), friction(f), tmg(tg), I(i), J(j), K(k)
	{
	}

	virtual ~MovingObject(){

	}
	virtual void Update(Voxel &voxel, float *phic, float *phiu, float *phiv, float *phiw) = 0;
	virtual void SetConstraintVelocity(Vector &vel, int i, int j, int k) const = 0;
	virtual void SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k) const = 0;
	virtual void SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k, int vel_index) const = 0;
	virtual bool StopNow() const = 0;
	virtual void SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const = 0;
	void Intersection(float *phi0) const{
		FOR_EACH_CELL
//			if(i == 49 && j == 49 && k == 49 )
//						printf(" phi0 = %f, phi = %f \n",
//								phi0[INDEX(i,j,k)],phi[INDEX(i,j,k)]);

			phi0[INDEX(i,j,k)] = phi0[INDEX(i,j,k)] < phi[INDEX(i,j,k)] ?
									phi[INDEX(i,j,k)] : phi0[INDEX(i,j,k)];
//			if(i == 49 && j == 49 && k == 49 )
//					printf(" phi0 = %f, phi = %f \n",
//							phi0[INDEX(i,j,k)],phi[INDEX(i,j,k)]);


		END_FOR
	}
	virtual void Friction(const Voxel &voxel, int vel_index, int i, int j, int k, float *f) const = 0;
	virtual bool IsBorderCell(const Voxel &voxel, int i, int j, int k) const = 0;
	virtual bool IsDynamicBorderCell(const Voxel &voxel, int i, int j, int k) const = 0;
	virtual void SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const = 0;
protected:
	float *phi;
	int DimX, DimY, DimZ;
	float h;
	float friction;
	TimeManager *tmg;
	int I, J, K;
	Transform trans;
};
#endif /*MOVINGOBJ_H_*/
