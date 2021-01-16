#include "cubesource.h"

void CubeSource::Initialize(const Voxel &voxel, float *f) {

	float leftWall = voxel.FacePostion(1, false, start_x0, start_y0, start_z0);
    float rightWall = voxel.FacePostion(1, true, end_x0, start_y0, start_z0);
    float backWall = voxel.FacePostion(2, false, start_x0, start_y0, start_z0);
    float frontWall = voxel.FacePostion(2, true, start_x0, end_y0, start_z0);
    float bottomWall =  voxel.FacePostion(3, false, start_x0, start_y0, start_z0);
    float topWall =  voxel.FacePostion(3, true, start_x0, start_y0, end_z0);
    int i,j,k;
	for(k = start_z0; k <= end_z0; ++k)
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i){
				Point p = voxel.VoxelCenterPosition(i,j,k);

    		if( bottomWall <= p[2] && p[2] <= topWall
    			&& leftWall <= p[0] && p[0] <= rightWall
    			&& backWall <= p[1] && p[1] <= frontWall) {
    	       float dist[6];
    	       dist[0] = p[0] - leftWall;
    	       dist[1] = rightWall - p[0];
    	       dist[2] = p[1] - backWall;
    	       dist[3] = frontWall - p[1];
    	       dist[4] = p[2] - bottomWall;
    	       dist[5] = topWall - p[2];
    	       float d = FindMinimum(dist, 6);
    	       int ii = i - start_x0;
    	       int jj = j - start_y0;
    	       int kk = k - start_z0;
    	       f[kk*(dimx*dimy)+jj*dimx+ii] = -d;
    	    }
    		else{
    			printf("Error in CubeSource::Initialize()! Point outside of source\n ");
    			exit(1);
    		}

//    		if( i == I && j == J && k == K)
//    			printf("p at (%f,%f,%f) with phi = %f \n", p[0], p[1], p[2], phi[INDEX(i,j,k)]);
    	}

}

void CubeSource::UpdateVoxel(Voxel &voxel) const{
	int i,j,k;
	for(k = start_z0; k <= end_z0; ++k)
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i){
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
				}
			}
}

void CubeSource::RestoreVoxel(Voxel &voxel) const{
	int i,j,k;
	for(k = start_z0; k <= end_z0; ++k)
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i){
				if(voxel.InSource(i,j,k)){
					voxel.UpdateCellType(i,j,k,CLEAR,SOURCE);
					if(!voxel.InLiquid(i,j,k))
						voxel.UpdateCellType(i,j,k,SET,LIQUID);
					if(voxel.InAir(i,j,k))
						voxel.UpdateCellType(i,j,k,CLEAR,EMPTY);
					if(voxel.InSurface(i,j,k))
						voxel.UpdateCellType(i,j,k,CLEAR,SURFACE);
					if(voxel.InSolid(i,j,k))
						voxel.UpdateCellType(i,j,k,CLEAR,SOLID);
				}
			}
}


void CubeSource::Merge(const Voxel &voxel, float *phi0, bool initial) const{
	int i,j,k;
	if(!stop){
		for(k = start_z0; k <= end_z0; ++k)
			for(j = start_y0; j <= end_y0; ++j)
				for(i = start_x0; i <= end_x0; ++i){
					int ii = i - start_x0;
					int jj = j - start_y0;
					int kk = k - start_z0;
					phi0[INDEX(i,j,k)] = phi0[INDEX(i,j,k)] <= phi[kk*(dimx*dimy)+jj*dimx+ii] ?
							phi0[INDEX(i,j,k)] : phi[kk*(dimx*dimy)+jj*dimx+ii];
				}
	}
//	else{
//		for(k = start_z0; k <= end_z0; ++k)
//			for(j = start_y0; j <= end_y0; ++j)
//				for(i = start_x0; i <= end_x0; ++i){
//					int ii = i - start_x0;
//					int jj = j - start_y0;
//					int kk = k - start_z0;
//					phi0[INDEX(i,j,k)] = phi0[INDEX(i,j,k)] >= phi[kk*(dimx*dimy)+jj*dimx+ii] ?
//							phi0[INDEX(i,j,k)] : phi[kk*(dimx*dimy)+jj*dimx+ii];
//				}
//
//	}

}

void  CubeSource::SetVelocity(const Voxel &voxel, float *u, float *v, float *w) const {
	int i,j,k;

//	for(k = start_z0; k <= end_z0; ++k)
//		for(j = start_y0; j <= end_y0; ++j)
//			for(i = start_x0; i <= end_x0; ++i){
//				u[INDEX(i,j,k)] = speed;
//				u[INDEX(i-1,j,k)] = speed;
//				v[INDEX(i,j,k)] = 0.f;
//				v[INDEX(i,j-1,k)] = 0.f;
//				w[INDEX(i,j,k)] = 0.f;
//				w[INDEX(i,j,k-1)] = 0.f;
//			}
//	i = start_x0;
//	for(k = start_z0; k <= end_z0; ++k)
//			for(j = start_y0; j <= end_y0; ++j){
//				u[INDEX(i,j,k)] = 0.f;
//				u[INDEX(i-1,j,k)] = 0.f;
//			}

	// a source with flow out of the bottom face
//	for(k = start_z0; k <= end_z0; ++k)
//		for(j = start_y0; j <= end_y0; ++j)
//			for(i = start_x0; i <= end_x0; ++i){
//				u[INDEX(i,j,k)] = 0.f;
//				u[INDEX(i-1,j,k)] = 0.f;
//				v[INDEX(i,j,k)] = 0.f;
//				v[INDEX(i,j-1,k)] = 0.f;
//				w[INDEX(i,j,k)] = speed;
//				w[INDEX(i,j,k-1)] = speed;
//			}
//	k = end_z0;
//	for(j = start_y0; j <= end_y0; ++j)
//		for(i = start_x0; i <= end_x0; ++i){
//			w[INDEX(i,j,k)] = 0.f;
//			w[INDEX(i,j,k-1)] = 0.f;
//		}
	if(!stop){
		for(k = start_z0; k <= end_z0; ++k)
			for(j = start_y0; j <= end_y0; ++j)
				for(i = start_x0; i <= end_x0; ++i){
					if(axis == 1){
						u[INDEX(i,j,k)] = speed;
						u[INDEX(i-1,j,k)] = speed;
						v[INDEX(i,j,k)] = 0.f;
						v[INDEX(i,j-1,k)] = 0.f;
						w[INDEX(i,j,k)] = 0.f;
						w[INDEX(i,j,k-1)] = 0.f;

					}
					else if(axis == 2){
						u[INDEX(i,j,k)] = 0.f;
						u[INDEX(i-1,j,k)] = 0.f;
						v[INDEX(i,j,k)] = speed;
						v[INDEX(i,j-1,k)] = speed;
						w[INDEX(i,j,k)] = 0.f;
						w[INDEX(i,j,k-1)] = 0.f;
					}
					else if(axis == 3){
						u[INDEX(i,j,k)] = 0.f;
						u[INDEX(i-1,j,k)] = 0.f;
						v[INDEX(i,j,k)] = 0.f;
						v[INDEX(i,j-1,k)] = 0.f;
						w[INDEX(i,j,k)] = speed;
						w[INDEX(i,j,k-1)] = speed;
					}

				}

		if(axis == 1){
			if(side == -1)
				i = end_x0;
			else
				i = start_x0;
			for(k = start_z0; k <= end_z0; ++k)
				for(j = start_y0; j <= end_y0; ++j){
					u[INDEX(i,j,k)] = 0.f;
					u[INDEX(i-1,j,k)] = 0.f;
				}

		}
		else if(axis == 2){
			if(side == -1)
				j = end_y0;
			else
				j = start_y0;
			for(k = start_z0; k <= end_z0; ++k)
				for(i = start_x0; i <= end_x0; ++i){
					v[INDEX(i,j,k)] = 0.f;
					v[INDEX(i,j-1,k)] = 0.f;
				}

		}
		else if(axis == 3){
			if(side == -1)
				k = end_z0;
			else
				k = start_z0;
			for(j = start_y0; j <= end_y0; ++j)
				for(i = start_x0; i <= end_x0; ++i){
					w[INDEX(i,j,k)] = 0.f;
					w[INDEX(i,j,k-1)] = 0.f;
				}
		}
	}

}

void CubeSource::Destroy(Voxel &voxel){
	int i,j,k;
	if(!stop){
		stop = true;
		RestoreVoxel(voxel);
	}
//	for(k = start_z0; k <= end_z0; ++k)
//		for(j = start_y0; j <= end_y0; ++j)
//			for(i = start_x0; i <= end_x0; ++i){
//				int ii = i - start_x0;
//				int jj = j - start_y0;
//				int kk = k - start_z0;
//				phi[kk*(dimx*dimy)+jj*dimx+ii] = -1 * phi[kk*(dimx*dimy)+jj*dimx+ii];
//			}

}















