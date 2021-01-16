/*
 * cylindricalcontainer.cpp
 *
 *  Created on: Sep 25, 2008
 *      Author: bzhao
 */

#include "cylindricalcontainer.h"
#include "mtrand.h"

void CylinderContainer::EvaluatePhi(const Voxel &voxel, int axis, float *phitmp){
		FOR_EACH_CELL
			Point p;
			if( axis == 0 )
				p = voxel.VoxelCenterPosition(i,j,k);
			else
				p = voxel.VelPosition(axis,i,j,k);
			float dist, dist1, dist2;
			Point newp = trans(p);
			float d2 = newp.x * newp.x + newp.z * newp.z;
			if(newp.y <= h && newp.y >= -h && d2 <= radius*radius){
				dist1 = radius - sqrt(d2);
				dist2 = min(h-newp.y, newp.y+h);
				dist = min(dist1, dist2);
				phitmp[INDEX(i,j,k)] =  -dist;
//					if(i == 33 && j == 20 && k ==91 )
//					if( k == 91 )
//						printf(" at (%d, %d, %d) p at (%f, %f, %f) newp at (%f, %f, %f) phi is source = %f \n",
//								i,j,k, p.x, p.y, p.z, newp.x, newp.y, newp.z, -dist);
			}
			else{
				float d  = sqrt(d2);
				dist1 = d - radius;
				dist2 = min(fabs(newp.y-h), fabs(-h-newp.y));
				if(d <= radius)
					dist = dist2;
				else if(newp.y <= h && newp.y >= -h)
					dist = dist1;
				else
					dist = max(dist1, dist2);
				phitmp[INDEX(i,j,k)] =  dist;  // > 0 : outside of cylinder
			}
		END_FOR

}
float CylinderContainer::EvaluatePhi(const Point &p) const {
		
		float dist, dist1, dist2;
		Point newp = trans(p);
		float d2 = newp.x * newp.x + newp.z * newp.z;
		if(newp.y <= h && newp.y >= -h && d2 <= radius*radius){
			dist1 = radius - sqrt(d2);
			dist2 = min(h-newp.y, newp.y+h);
			dist = min(dist1, dist2);
			return -dist;
		}
		else{
			float d  = sqrt(d2);
			dist1 = d - radius;
			dist2 = min(fabs(newp.y-h), fabs(-h-newp.y));
			if(d <= radius)
				dist = dist2;
			else if(newp.y <= h && newp.y >= -h)
				dist = dist1;
			else
				dist = max(dist1, dist2);
			return dist;  // > 0 : outside of cylinder
		}
}

bool CylinderContainer::IsBorderCell(const Voxel &voxel, int i, int j, int k) const{
	for(int n=0; n<8; ++n){
		Point q = voxel.VoxelCornerPosition(i,j,k,n+1);
		float f = EvaluatePhi(q);
		if( f > E_EPSIL )
			return true;
	}
	return false;
}

void CylinderContainer::SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const{

	Point p;
	p = voxel.SubFacePosition(n,0,i,j,k);
	f[0] = EvaluatePhi(p);
	p = voxel.SubFacePosition(n,1,i,j,k);
	f[1] = EvaluatePhi(p);
	p = voxel.SubFacePosition(n,2,i,j,k);
	f[2] = EvaluatePhi(p);
	p = voxel.SubFacePosition(n,3,i,j,k);
	f[3] = EvaluatePhi(p);


}

float CylinderContainer::TriInterp(const Voxel &voxel, const Point &p) const{

	float delta = voxel.VoxelDelta();

	float x = p.x/delta - 0.5f;
    float y = p.y/delta - 0.5f;
    float z = p.z/delta - 0.5f;
    if(x < 0) x = 0.f;
    if(y < 0) y = 0.f;
    if(z < 0) z = 0.f;
	int xi = floor(x);
	int yi = floor(y);
	int zi = floor(z);
//	printf(" xi = %d, yi = %d, zi = %d \n", xi, yi, zi);
	//Point point = voxel.VoxelCenterPosition(xi,yi,zi);
	Point point = voxel.VoxelCornerPosition(xi,yi,zi,1);
	float r = x - point.x/delta;
    float s = y - point.y/delta;
    float t = z - point.z/delta;
    int xj, yj, zj;
    if( xi >= DimX-1 ){
    	xi = DimX-1;
    	xj = xi;
    	r = 1.f;
    }
    else
    	xj = xi + 1;
    if( yi >= DimY-1 ){
    	yi = DimY-1;
        yj = yi;
        s = 1.f;
    }
    else
    	yj = yi + 1;
    if( zi >= DimZ-1 ){
    	zi = DimZ-1;
        zj = zi;
        t = 1.f;
    }
    else
    	zj = zi + 1;

//	printf(" r = %f,  s = %f,  t = %f \n",  r, s, t);
	return trilerp(phi[INDEX(xi, yi, zi)], phi[INDEX(xj, yi, zi)],
				   phi[INDEX(xi, yi, zj)], phi[INDEX(xj, yi, zj)],
				   phi[INDEX(xi, yj, zi)], phi[INDEX(xj, yj, zi)],
				   phi[INDEX(xi, yj, zj)], phi[INDEX(xj, yj, zj)],
				   r, s, t);
}

void CylinderContainer::SetNormal1(const Voxel &voxel, Vector &N, int i, int j, int k) const{
	
		float gx, gy, gz;
		float delta = voxel.VoxelDelta();
		Point p = voxel.VoxelCenterPosition(i,j,k);
		Point px0(p.x-0.5f*delta, p.y, p.z);
		Point px1(p.x+0.5f*delta, p.y, p.z);
		float phix0 = TriInterp(voxel, px0);
		float phix1 = TriInterp(voxel, px1);
		gx = (phix1 - phix0) / delta;

		Point py0(p.x, p.y-0.5f*delta, p.z);
		Point py1(p.x, p.y+0.5f*delta, p.z);
		float phiy0 = TriInterp(voxel, py0);
		float phiy1 = TriInterp(voxel, py1);
		gy = (phiy1 - phiy0) / delta;

		Point pz0(p.x, p.y, p.z-0.5f*delta);
		Point pz1(p.x, p.y, p.z+0.5f*delta);
		float phiz0 = TriInterp(voxel, pz0);
		float phiz1 = TriInterp(voxel, pz1);
		gz = (phiz1 - phiz0) / delta;

		N = Vector(gx, gy, gz);
		if(N.Length() == 0.f){
			printf("at (%d, %d, %d) object normal undefined!\n", i,j,k);
//			printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
//			printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
//					"phi[k-1] = %f, phi[k+1] = %f \n\n",
//			   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
//				 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
//				 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			MTRand mt;
			double rn1 = mt();
			double rn2 = mt();
			double rn3 = mt();
			N.x = (float)rn1;
			N.y = (float)rn2;
			N.z = (float)rn3;
		}
		N = -1 * Normalize(N); // multiply by -1 because phi at object grid points is positive

}


void CylinderContainer::SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const{
	Point p = voxel.VoxelCenterPosition(i,j,k);
	float dist, dist1, dist2;
	Point newp = trans(p);
	float d2 = newp.x * newp.x + newp.z * newp.z;
	if(newp.y <= h && newp.y >= -h && d2 <= radius*radius){
		dist1 = radius - sqrt(d2);
		dist2 = min(h-newp.y, newp.y+h);
//		if(dist1 < dist2){
		if(IsBorderCell(voxel,i,j,k)){
			Vector N1(newp.x,0,newp.z);
			N1 = Normalize(N1);
			N = trans_inv(N1);
			N = -1*Normalize(N);
		}
		else{
			if(newp.y < 0.f){
				Vector N1(0,1,0);
				N = trans_inv(N1);
			}
			else{
				Vector N1(0,-1,0);
				N = trans_inv(N1);
			}
			N = Normalize(N);
		}
	}
	else{
		float d  = sqrt(d2);
		dist1 = d - radius;
		dist2 = min(fabs(newp.y-h), fabs(-h-newp.y));
		if(d <= radius){
			if(newp.y > h){
				Vector N1(0,-1,0);
				N = trans_inv(N1);
			}
			else if(newp.y < -h){
				Vector N1(0,1,0);
				N = trans_inv(N1);
			}
			N = Normalize(N);
		}
		else if(newp.y <= h && newp.y >= -h){
			Vector N1(newp.x,0,newp.z);
			N1 = Normalize(N1);
			N = trans_inv(N1);
			N = -1*Normalize(N);
		}
		else
			SetNormal1(voxel, N, i,j,k);
	}
}
