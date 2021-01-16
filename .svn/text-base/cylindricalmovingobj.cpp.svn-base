/*
 * cylindricalmovingobj.cpp
 *
 *  Created on: Sep 25, 2008
 *      Author: bzhao
 */

#include "cylindricalmovingobj.h"


void CylinderMovingObject::Update(Voxel &voxel, float *phic, float *phiu, float *phiv, float *phiw){

	SetObjectLevelSet(voxel, 1);
	Intersection(phiu);
	SetObjectLevelSet(voxel, 2);
	Intersection(phiv);
	SetObjectLevelSet(voxel, 3);
	Intersection(phiw);

	SetObjectLevelSet(voxel, 0);
	FOR_EACH_CELL
		if(phi[INDEX(i,j,k)] > 0.f){
			if(!voxel.InMovingSolid(i,j,k))
				voxel.UpdateCellType(i,j,k,SET,MOVSOLID);
		}
		else{
			if(voxel.InMovingSolid(i,j,k))
				voxel.UpdateCellType(i,j,k,CLEAR,MOVSOLID);
		}
	END_FOR
	// if moving solid objects stop, they become part of the static solid boundary
	// then the reinitialize routine will treat them accordingly
	if(StopNow()){
		FOR_EACH_CELL
			if(voxel.InMovingSolid(i,j,k))
				voxel.UpdateCellType(i,j,k,CLEAR,MOVSOLID);
		END_FOR
	}
	Intersection(phic);
	int NM = 0;
	FOR_EACH_CELL
		if(phi[INDEX(i,j,k)] > 0.f)
			NM++;
	END_FOR
	printf("moving object velocity = (%f, %f, %f)\n", us, vs, ws);
	printf("moving object occupies %d cells \n", NM);


//	FOR_EACH_CELL
//		if(phi[INDEX(i,j,k)] >= 0.f){  // in moving object
//			if(!voxel.InSolid(i,j,k))
//				voxel.UpdateCellType(i,j,k,SET,SOLID);
//			if(voxel.InLiquid(i,j,k))
//				voxel.UpdateCellType(i,j,k,CLEAR,LIQUID);
//			if(voxel.InAir(i,j,k))
//				voxel.UpdateCellType(i,j,k,CLEAR,EMPTY);
//			if(voxel.InSurface(i,j,k))
//				voxel.UpdateCellType(i,j,k,CLEAR,SURFACE);
//		}
//		else{
//			if(voxel.InSolid(i,j,k))
//				voxel.UpdateCellType(i,j,k,CLEAR,SOLID);
//		}
//	END_FOR
	return;
}

void CylinderMovingObject::SetObjectLevelSet(const Voxel &voxel, int axis){
	FOR_EACH_CELL
		Point q;
		if(axis == 0)
			q = voxel.VoxelCenterPosition(i,j,k);
		else
			q = voxel.VelPosition(axis,i,j,k);
		Point newp = trans(q);
		float d2 = newp.x * newp.x + newp.z * newp.z;
		float d = sqrt(d2);
		float dist, dist1, dist2;
		if(newp.y <= y_extent && newp.y >= -y_extent && d <= radius){
			dist1 = radius - d;
			dist2 = min(y_extent-newp.y, newp.y+y_extent);
			dist = min(dist1, dist2);
			phi[INDEX(i,j,k)] = -dist;  // < 0 : inside cylinder
		}
		else{
			dist1 = d - radius;
			dist2 = min(fabs(newp.y-y_extent), fabs(-y_extent-newp.y));
			if(d <= radius)
				dist = dist2;
			else if(newp.y <= y_extent && newp.y >= -y_extent)
				dist = dist1;
			else
				dist = min(dist1, dist2);
			phi[INDEX(i,j,k)] = dist;  // > 0 : outside of cylinder
		}
		if(i == I && j == J && k == K )
			printf(" CylinderMovingObject::SetObjectLevelSet  point at (%f, %f, %f) transformed to (%f, %f, %f) phi = %f \n",
					q.x, q.y, q.z, newp.x, newp.y, newp.z, phi[INDEX(i,j,k)]);
	END_FOR
	return;
}

void CylinderMovingObject::SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const {

	if(phi[INDEX(i,j,k)] >= 0.f){
		Point q = voxel.VoxelCenterPosition(i,j,k);
		Point newp = trans(q);
		float d2 = newp.x * newp.x + newp.z * newp.z;
		float d = sqrt(d2);
		float dist, dist1, dist2, dist3;
		dist1 = y_extent - newp.y;
		dist2 = newp.y + y_extent;
		dist = radius - d;
		dist3 = min(dist1, dist2);
		if(dist <= dist3)
			N = Vector(newp.x, 0, newp.z);
		else{
			N = dist1 <= dist2 ? Vector(0, 1, 0) : Vector(0, -1, 0);
		}
		N = trans_inv(N);
		if(N.Length() != 0.f)
			N = Normalize(N);
	}
	return;
}

void CylinderMovingObject::SetConstraintVelocity(Vector &vel, int i, int j, int k) const{

	if(phi[INDEX(i,j,k)] >= 0.f){
		vel.x = us;
		vel.y = vs;
		vel.z = ws;
	}
	return;
}

void CylinderMovingObject::SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k) const{
	float dt = tmg->GetDt();
	if(phi[INDEX(i,j,k)] >= 0.f){
		Point q = voxel.VoxelCenterPosition(i,j,k);
		Point old_pos = trans(q);
		Point new_pos = trans(q);
		Vector v = new_pos - old_pos;
		vel = v / dt;
	}

}

void CylinderMovingObject::SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k, int vel_index) const{
		float dt = tmg->GetDt();
		Point old_pos = voxel.VelPosition(vel_index,i,j,k);
		Point new_pos = trans(old_pos);
		Vector v = new_pos - old_pos;
		vel = v / dt;
}

void CylinderMovingObject::Friction(const Voxel &voxel, int vel_index, int i, int j, int k, float *f) const{
	Point q = voxel.VelPosition(vel_index,i,j,k);
	Point newp = trans(q);
	float d2 = newp.x * newp.x + newp.z * newp.z;
	float d = sqrt(d2);
	float dist, dist1, dist2;
	if(newp.y < y_extent && newp.y > -y_extent && d < radius);
	else
		*f = friction;

}

