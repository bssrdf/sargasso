/*
 * cubecontainer.cpp
 *
 *  Created on: Sep 25, 2008
 *      Author: bzhao
 */

#include "cubecontainer.h"


void CubeContainer::EvaluatePhi(const Voxel &voxel, int axis, float *phi_tmp){

	    printf("left wall at %f, right wall at %f \n",leftWall, rightWall);
	    printf("back wall at %f, front wall at %f \n",backWall, frontWall);
	    printf("bottom wall at %f, top wall at %f \n",bottomWall, topWall);

		FOR_EACH_CELL
			Point p;
			if( axis == 0 )
				p = voxel.VoxelCenterPosition(i,j,k);
			else
				p = voxel.VelPosition(axis,i,j,k);

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
		       phi_tmp[INDEX(i,j,k)] = -d;
		    }
			else if(p[0] < leftWall)
				phi_tmp[INDEX(i,j,k)] = leftWall - p[0];
			else if(p[0] > rightWall)
				phi_tmp[INDEX(i,j,k)] = p[0] - rightWall;
			else if(p[1] < backWall)
				phi_tmp[INDEX(i,j,k)] = backWall - p[1];
			else if(p[1] > frontWall)
				phi_tmp[INDEX(i,j,k)] = p[1] - frontWall;
			else if(p[2] < bottomWall)
				phi_tmp[INDEX(i,j,k)] = bottomWall - p[2];
			else if(p[2] > topWall)
				phi_tmp[INDEX(i,j,k)] = p[2] - topWall;
		END_FOR
}

float CubeContainer::EvaluatePhi(const Point &p) const{
	float f;
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
       f = -d;
    }
	else if(p[0] < leftWall)
		f = leftWall - p[0];
	else if(p[0] > rightWall)
		f = p[0] - rightWall;
	else if(p[1] < backWall)
		f = backWall - p[1];
	else if(p[1] > frontWall)
		f = p[1] - frontWall;
	else if(p[2] < bottomWall)
		f = bottomWall - p[2];
	else if(p[2] > topWall)
		f = p[2] - topWall;
	return f;
}

bool CubeContainer::IsBorderCell(const Voxel &voxel, int i, int j, int k) const{

//	for(int n=0; n<8; ++n){
//		Point p = voxel.VoxelCornerPosition(i,j,k,n+1);
//		float distance;
//		if( bottomWall <= p[2] && p[2] <= topWall
//			&& leftWall <= p[0] && p[0] <= rightWall
//			&& backWall <= p[1] && p[1] <= frontWall) {
//		   float dist[6];
//		   dist[0] = p[0] - leftWall;
//		   dist[1] = rightWall - p[0];
//		   dist[2] = p[1] - backWall;
//		   dist[3] = frontWall - p[1];
//		   dist[4] = p[2] - bottomWall;
//		   dist[5] = topWall - p[2];
//		   float d = FindMinimum(dist, 6);
//		   distance = -d;
//		}
//		else if(p[0] < leftWall)
//			distance = leftWall - p[0];
//		else if(p[0] > rightWall)
//			distance = p[0] - rightWall;
//		else if(p[1] < backWall)
//			distance = backWall - p[1];
//		else if(p[1] > frontWall)
//			distance = p[1] - frontWall;
//		else if(p[2] < bottomWall)
//			distance = bottomWall - p[2];
//		else if(p[2] > topWall)
//			distance = p[2] - topWall;
//		if( distance > E_EPSIL )
//			return true;
//	}
	return false;

}


void CubeContainer::SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const{

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

void CubeContainer::SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const{
	TentativeGridPoint p(0.f, i,j,k);
	vector<TentativeGridPoint> neighbors;
	p.Neigbor(neighbors, DimX, DimY, DimZ);
	for(int m=0; m<neighbors.size(); ++m){
		if(!voxel.InSolid(neighbors[m])){
			if(p.LeftNeighbor(neighbors[m])){
				N = Vector(-1.f, 0.f, 0.f);
			}
			if(p.RightNeighbor(neighbors[m])){
				N = Vector(1.f, 0.f, 0.f);
			}
			if(p.FrontNeighbor(neighbors[m])){
				N = Vector(0.f, 1.f, 0.f);
			}
			if(p.BackNeighbor(neighbors[m])){
				N = Vector(0.f, -1.f, 0.f);
			}
			if(p.TopNeighbor(neighbors[m])){
				N = Vector(0.f, 0.f, 1.f);
			}
			if(p.BottomNeighbor(neighbors[m])){
				N = Vector(0.f, 0.f, -1.f);
			}
		}
	}
}