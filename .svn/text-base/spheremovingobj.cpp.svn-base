#include "spheremovingobj.h"


void SphereMovingObject::Update(Voxel &voxel, float *phic, float *phiu, float *phiv, float *phiw){

	float dt = tmg->GetDt();
	px += dt * us;
	py += dt * vs;
	pz += dt * ws;

	if(pz - (WALL_THICKNESS+1) * h <= rad + 2 * h){
		us = 0.f;
		vs = 0.f;
		ws = 0.f;
	}
	trans = Translate(dt*Vector(us, vs, ws));
	Point p(px, py, pz);

	SetObjectLevelSet(voxel, p, rad, 1, phi);
	Intersection(phiu);
	SetObjectLevelSet(voxel, p, rad, 2, phi);
	Intersection(phiv);
	SetObjectLevelSet(voxel, p, rad, 3, phi);
	Intersection(phiw);

	SetObjectLevelSet(voxel, p, rad, 0, phi);
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

void SphereMovingObject::SetObjectLevelSet(const Voxel &voxel, const Point &p, float r, char axis, float *_phi){
	FOR_EACH_CELL
		Point q;
		if(axis == 0)
			q = voxel.VoxelCenterPosition(i,j,k);
		else
			q = voxel.VelPosition(axis,i,j,k);
		float d = Distance(p, q);
		_phi[INDEX(i,j,k)] = r - d;  // > 0 : inside object
		                             // < 0 : outside object
 		if(axis == 0 && i == I && j == J && k == K){
			printf("p at (%f, %f, %f), q at (%f, %f, %f) \n",
					p.x, p.y, p.z, q.x, q.y, q.z);
			printf(" d = %f, r = %f, phi = %f \n",
					d, r, _phi[INDEX(i,j,k)]);
//			Vector N = q - p;
//			N = Normalize(N);
//			printf("object normal = (%f, %f, %f) \n", N.x, N.y, N.z);
		}
	END_FOR
//	int minz = 79;
//	int mini, minj;
//	FOR_EACH_CELL
//		if(phi[INDEX(i,j,k)] > 0.f)
//			if(k < minz){
//				minz = k;
//				minj = j;
//				mini = i;
//			}
//	END_FOR
//	printf("lowest sphere point at (%d, %d, %d) \n", mini, minj, minz);
}

float SphereMovingObject::EvaluatePhi(const Point &q) const {
	Point p(px,py,pz);
	float f = 0;
	float d = Distance(p, q);
	f = rad - d;  // > 0 : inside object
	              // < 0 : outside object
	return f;
}


void SphereMovingObject::SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const {
	if(phi[INDEX(i,j,k)]+E_EPSIL > 0.f){
//	if(IsBorderCell(voxel,i,j,k)){
		Point q = voxel.VoxelCenterPosition(i,j,k);
		Point p(px, py, pz);
		N = q - p;
		if(N.Length() != 0.f)
			N = Normalize(N);
	}
}
void SphereMovingObject::SetConstraintVelocity(Vector &vel, int i, int j, int k) const{

	if(phi[INDEX(i,j,k)] > 0.f){
		vel.x = us;
		vel.y = vs;
		vel.z = ws;
	}

}

void SphereMovingObject::SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k) const{
	float dt = tmg->GetDt();
	if(phi[INDEX(i,j,k)]+E_EPSIL > 0.f){
//	if(IsBorderCell(voxel,i,j,k)){
		Point old_pos = voxel.VoxelCenterPosition(i,j,k);
		Point new_pos = trans(old_pos);
		Vector v = new_pos - old_pos;
		vel = v / dt;
	}

}


bool SphereMovingObject::IsBorderCell(const Voxel &voxel, int i, int j, int k) const{
//	Point p(px,py,pz);
	if(phi[INDEX(i,j,k)] > E_EPSIL){
		for(int n=0; n<6; ++n){
			for(int m=0; m<4; ++m){
				Point q = voxel.SubFacePosition(n,m,i,j,k);
//				float d = Distance(p, q);
				if( EvaluatePhi(q) + E_EPSIL < 0.f )
					return true;
			}
		}
		return false;
	}
	else{
		for(int n=0; n<6; ++n){
			for(int m=0; m<4; ++m){
				Point q = voxel.SubFacePosition(n,m,i,j,k);
//				float d = Distance(p, q);
//				if(i == 39 && j == 39 && k == 47){
//					printf("n = %d, m = %d, rad = %f, d = %f, val = %f \n", n, m, rad, d, EvaluatePhi(q));
//					printf("p at (%f, %f, %f), q at (%f, %f, %f) \n", p.x, p.y, p.z, q.x, q.y, q.z);
//				}
				if(  EvaluatePhi(q) + E_EPSIL > 0.f )
					return true;
			}
		}
		return false;
	}
}

bool SphereMovingObject::IsDynamicBorderCell(const Voxel &voxel, int i, int j, int k) const{
	if(IsBorderCell(voxel,i,j,k) && !StopNow())
		return true;
	else
		return false;
}

void SphereMovingObject::SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const {
	Point p;
	p = voxel.SubFacePosition(n,0,i,j,k);
	f[0] = EvaluatePhi(p);
	p = voxel.SubFacePosition(n,1,i,j,k);
	f[1] = EvaluatePhi(p);
	p = voxel.SubFacePosition(n,2,i,j,k);
	f[2] = EvaluatePhi(p);
	p = voxel.SubFacePosition(n,3,i,j,k);
	f[3] = EvaluatePhi(p);
//	if(i == 39 && j == 40 && k == 47){
//		printf("n = %d, f[0] = %f  f[1] = %f f[2] = %f f[3] = %f \n", n, f[0], f[1], f[2], f[3]);
//	}

}


void SphereMovingObject::SetConstraintVelocity(const Voxel &voxel, Vector &vel, int i, int j, int k, int vel_index) const{
		float dt = tmg->GetDt();
		Point old_pos = voxel.VelPosition(vel_index,i,j,k);
		float d = Distance(Point(px, py, pz), old_pos);
		float dist = rad - d;  // > 0 : inside object
				               // < 0 : outside object
		if(dist+E_EPSIL > 0.f){
			Point new_pos = trans(old_pos);
			Vector v = new_pos - old_pos;
			vel = v / dt;
		}
}
void SphereMovingObject::Friction(const Voxel &voxel, int vel_index, int i, int j, int k, float *f) const{
	Point old_pos = voxel.VelPosition(vel_index,i,j,k);
	float d = Distance(Point(px, py, pz), old_pos);
	float dist = rad - d;  // > 0 : inside object
			               // < 0 : outside object
	if(dist+E_EPSIL > 0.f){
		*f = friction;
	}
}


