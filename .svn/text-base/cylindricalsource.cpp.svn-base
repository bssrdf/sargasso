#include "cylindricalsource.h"


void CylinderSource::UpdateVoxel(Voxel &voxel) const{
	FOR_EACH_CELL
		Point p = voxel.VoxelCenterPosition(i,j,k);
		Point newp = trans(p);
		float d2 = newp.x * newp.x + newp.z * newp.z;
		if(newp.y <= h && newp.y >= -h && d2 <= radius *radius){
			if(!voxel.InSource(i,j,k)){
				voxel.UpdateCellType(i,j,k,SET,SOURCE);
				if(voxel.InLiquid(i,j,k))
					voxel.UpdateCellType(i,j,k,CLEAR,LIQUID);
				if(voxel.InAir(i,j,k))
					voxel.UpdateCellType(i,j,k,CLEAR,EMPTY);
				if(voxel.InSurface(i,j,k))
					voxel.UpdateCellType(i,j,k,CLEAR,SURFACE);
				if(voxel.InSolid(i,j,k))
					voxel.UpdateCellType(i,j,k,CLEAR,SOLID);
				if(i == I && j == J && k == K )
					printf(" (%d, %d, %d) is included in source \n", i, j, k);
			}
		}
	END_FOR
}

void CylinderSource::Merge(const Voxel &voxel, float *phi0, bool initial) const{
	int i,j,k;
	if(!stop){
			FOR_EACH_CELL
				Point p = voxel.VoxelCenterPosition(i,j,k);
				float dist, dist1, dist2;
				Point newp = trans(p);
				float d2 = newp.x * newp.x + newp.z * newp.z;
				if(newp.y <= h && newp.y >= -h && d2 <= radius*radius){
					dist1 = radius - sqrt(d2);
					dist2 = min(h-newp.y, newp.y+h);
					dist = min(dist1, dist2);
					phi0[INDEX(i,j,k)] =  -dist;
					if(i == I && j == J && k == K )
//					if( k == 91 )
						printf(" at (%d, %d, %d) p at (%f, %f, %f) newp at (%f, %f, %f) phi is source = %f \n",
								i,j,k, p.x, p.y, p.z, newp.x, newp.y, newp.z, -dist);
				}
				else if(initial){
					float d = sqrt(d2);
					dist1 = d - radius;
					dist2 = min(fabs(newp.y-h), fabs(-h-newp.y));
					if(d <= radius)
						dist = dist2;
					else if(newp.y <= h && newp.y >= -h)
						dist = dist1;
					else
						dist = max(dist1, dist2);
					float old_phi = phi0[INDEX(i,j,k)];
					phi0[INDEX(i,j,k)] = phi0[INDEX(i,j,k)] <= dist ?
							phi0[INDEX(i,j,k)] : dist;  // > 0 : outside of cylinder
					if(i == I && j == J && k == K ){
						printf(" at (%d, %d, %d) p at (%f, %f, %f) newp at (%f, %f, %f) phi is out of source = %f \n",
							i,j,k, p.x, p.y, p.z, newp.x, newp.y, newp.z, dist);
						printf("old_phi = %f, dist1 = %f, dist2 = %f d = %f radius = %f, h = %f\n",
								old_phi, dist1, dist2, d, radius, h);
					}
				}
			END_FOR

	}
//	else{
//		FOR_EACH_CELL
//			Point p = voxel.VoxelCenterPosition(i,j,k);
//			float dist, dist1, dist2;
//			Point newp = trans(p);
//			float d2 = newp.x * newp.x + newp.z * newp.z;
//			if(newp.y <= h && newp.y >= -h && d2 <= radius *radius){
//				dist1 = fabs(sqrt(d2) - radius);
//				dist2 = min(h-newp.y, newp.y+h);
//				dist = min(dist1, dist2);
//				phi0[INDEX(i,j,k)] = phi0[INDEX(i,j,k)] >= dist ?
//							phi0[INDEX(i,j,k)] : dist;
//			}
//		END_FOR
//
//	}

}

void  CylinderSource::SetVelocity(const Voxel &voxel, float *u, float *v, float *w) const {
	if(!stop){
		Vector vel(0, velocity, 0);
		Vector new_vel = trans_inv(vel);
		FOR_EACH_CELL
//			Point p = voxel.VoxelCenterPosition(i,j,k);
//			Point newp = trans(p);
//			float d2 = newp.x * newp.x + newp.z * newp.z;
//			if(newp.y <= h && newp.y >= -h && d2 <= radius*radius){
			if(voxel.InSource(i,j,k)){
				u[INDEX(i,j,k)] = new_vel.x;
				v[INDEX(i,j,k)] = new_vel.y;
				w[INDEX(i,j,k)] = new_vel.z;
				if(!voxel.InSource(i-1,j,k))
					u[INDEX(i-1,j,k)] = new_vel.x;
				if(!voxel.InSource(i,j-1,k))
					v[INDEX(i,j-1,k)] = new_vel.y;
				if(!voxel.InSource(i,j,k-1))
					w[INDEX(i,j,k-1)] = new_vel.z;
			}
		END_FOR
	}
}

void CylinderSource::Destroy(Voxel &voxel){
	//int i,j,k;
	if(!stop){
		stop = true;
		RestoreVoxel(voxel);
	}

}

void CylinderSource::RestoreVoxel(Voxel &voxel) const{
	FOR_EACH_CELL
//		Point p = voxel.VoxelCenterPosition(i,j,k);
//		Point newp = trans(p);
//		float d2 = newp.x * newp.x + newp.z * newp.z;
//		if(newp.y <= h && newp.y >= -h && d2 <= radius*radius){
				if(voxel.InSource(i,j,k)){
					voxel.UpdateCellType(i,j,k,CLEAR,SOURCE);
					if(!voxel.InLiquid(i,j,k))
						voxel.UpdateCellType(i,j,k, SET,LIQUID);
					if(voxel.InAir(i,j,k))
						voxel.UpdateCellType(i,j,k,CLEAR,EMPTY);
					if(voxel.InSurface(i,j,k))
						voxel.UpdateCellType(i,j,k,CLEAR,SURFACE);
					if(voxel.InSolid(i,j,k))
						voxel.UpdateCellType(i,j,k,CLEAR,SOLID);
				}
//			}
	END_FOR
}
