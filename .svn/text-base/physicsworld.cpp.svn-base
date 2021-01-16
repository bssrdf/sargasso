/*
 * physicsworld.cpp
 *
 *  Created on: Nov 5, 2010
 *      Author: bzhao
 */

#include <queue>
using std::queue;

#include "physicsworld.h"
#include "voxel.h"
#include "mtrand.h"
#include "configfile.h"
#include "aliaswavefrontobj.h"
#include "vec.h"
#include "util.h"
#include "bounding_box.h"

void PhysicsWorld::BuildPhysicsWorld(const char *descfile, const Voxel *voxel){
	ConfigFile cf(descfile);
	unsigned int n = 0;
	while(cf.HasNextObject()){
		Transform o2w = cf.ReadTransform(n);
		string objfile = cf.ReadObjectFileName(n);
		++n;
		if(cf.IsAWObj(objfile)){
			AliasWavefrontObj *object = new AliasWavefrontObj;
			object->LoadObject(objfile.c_str());
			object->TransformO2W(o2w);
			mStaticTriangleMeshObjects.push_back(object);
			int num_verts = object->NumVertices();
			const double *verts = object->Vertices();
			int num_tris  = object->NumTriangles();
			const int *tris     = object->TriangleIndices();
			//LabelCellsCloseToMesh(voxel, num_verts, verts, num_tris, tris);
			MeshObject *meshObj = construct_mesh_object(num_verts, verts, num_tris, tris);
			mStaticMeshObjects.push_back(meshObj);
		}
	}
}

inline bool PhysicsWorld::IndexInBounds(const TentativeGridPoint &tp) const{
	if( tp.ii < 0 || tp.ii >= DimX || tp.jj < 0 || tp.jj >= DimY  || tp.kk < 0 || tp.kk >= DimZ )
		return false;
	else
		return true;
}

void PhysicsWorld::FloodFillSolid(int ii, int jj, int kk, char *color) const {
	queue<TentativeGridPoint> Q;
	TentativeGridPoint it(0.f, ii, jj, kk);
	u_int index = INDEX(ii,jj,kk);
//	printf("current position (%d, %d, %d) color = %d \n", ii, jj, kk,color[index]);
	if(color[index] == 2){
		printf("seed point in solid\n");
		return;
	}
//	else if(found != visitedPoints.end())
	else if(color[index] == 1){
		printf("seed point already set \n");
		return;
	}
	else{
		// put the seed point into queue to start the loop
		Q.push(it);
		color[index] = 1;
		do{
			// get the current point
			TentativeGridPoint tp = Q.front();
			vector<TentativeGridPoint> neighbors;
			tp.Neigbor(neighbors, DimX, DimY, DimZ);
			for(int m=0; m<neighbors.size(); m++){
				// if a neighbor point has not been set to 1, set it to 1 and put it into queue
				// these neighbors are immediately set to 1 here such that they will not be
				// duplicated in queue later on
				int i = neighbors[m].ii;
				int j = neighbors[m].jj;
				int k = neighbors[m].kk;
				if(color[INDEX(i, j, k)] == 0){
				   color[INDEX(i, j, k)] = 1;
				   TentativeGridPoint tmp(0.f, i, j, k);
				   Q.push(tmp);
				}
			}
			Q.pop();
//			printf("there are %lu iterms in queues \n", Q.size());
		}while(!Q.empty());
	}
}

void PhysicsWorld::BuildSignedDistanceMap(const Voxel *voxel, int index, float *phi) const{
	float delta = voxel->VoxelDelta();
	Point p0;
	double t;
	double P0[3], P1[3];
	float dist[6];
	bool  inside[6];
	float big = 1.e16f;
//	(DimX+DimY+DimZ)*delta;

	char *outsideflag = new char[DimX*DimY*DimZ];
	memset(outsideflag, 0, DimX*DimY*DimZ*sizeof(char));
	int I = 77;
	int J = 97;
	int K = 76;
 	FOR_EACH_CELL
//	int i = I;
//	int j = J;
//	int k = K;
			if(index == 0){
				p0 = voxel->VoxelCenterPosition(i,j,k);
			}
			else{
				p0 = voxel->VelPosition(index,i,j,k);
			}
			P0[0] = p0.x;P0[1] = p0.y;P0[2] = p0.z;
			for(unsigned int n=0; n<6; ++n){
				dist[n] = big;
				inside[n] = false;
			}
			bool pInside = false;
			for(unsigned int n=0; n < mStaticMeshObjects.size(); ++n) {
				double a,b,c;
				int tri;
				float d = big;
				P1[0] = 1.0;P1[1] = 0.0;P1[2] = 0.0;
				if(ray_intersects_mesh(P0, P1, mStaticMeshObjects[n], &tri, &t, &a, &b, &c, inside[0]))
					d = t;
				if(d < dist[0]) dist[0] = d;

				d = big;
				P1[0] = -1.0; P1[1] = 0.0;P1[2] = 0.0;
				if(ray_intersects_mesh(P0, P1, mStaticMeshObjects[n], &tri, &t, &a, &b, &c, inside[1]))
					d = t;
				if(d < dist[1]) dist[1] = d;

				d = big;
				P1[0] = 0.0;P1[1] = 1.0;P1[2] = 0.0;
				if(ray_intersects_mesh(P0, P1, mStaticMeshObjects[n], &tri, &t, &a, &b, &c,inside[2]))
					d = t;
				if(d < dist[2]) dist[2] = d;

				d = big;
				P1[0] = 0.0;P1[1] = -1.0;P1[2] = 0.0;
				if(ray_intersects_mesh(P0, P1, mStaticMeshObjects[n], &tri, &t, &a, &b, &c,inside[3]))
					d = t;
				if(d < dist[3]) dist[3] = d;

				d = big;
				P1[0] = 0.0;P1[1] = 0.0;P1[2] = 1.0;
				if(ray_intersects_mesh(P0, P1, mStaticMeshObjects[n], &tri, &t, &a, &b, &c, inside[4]))
					d = t;
				if(d < dist[4]) dist[4] = d;

				d = big;
				P1[0] = 0.0;P1[1] = 0.0;P1[2] = -1.0;
				if(ray_intersects_mesh(P0, P1, mStaticMeshObjects[n], &tri, &t, &a, &b, &c, inside[5]))
					d = t;
				if(d < dist[5]) dist[5] = d;

				int nin = 0, nout = 0;
				for(unsigned int m=0; m<6; ++m){
					if(inside[m])
						++nin;
					else
						++nout;
				}
				if(nin < nout)
					pInside = false;
				else if(nin > nout)
					pInside = true;
				else{
					Vector N;
					MTRand mt;
					double rn1 = mt();
					double rn2 = mt();
					double rn3 = mt();
					N.x = (float)rn1;
					N.y = (float)rn2;
					N.z = (float)rn3;
					N = Normalize(N);
					P1[0] = N.x;P1[1] = N.y;P1[2] = N.z;
					if(!ray_intersects_mesh(P0, P1, mStaticMeshObjects[n], &tri, &t, &a, &b, &c, pInside))
						pInside = false;
				}
//				if(point_inside_mesh(P0, mStaticMeshObjects[n]))
//					pInside = true;
			}
//			float df = min(dist[0],dist[1],dist[2],dist[3],dist[4],dist[5]);
			float df = big;
			for(unsigned int m=0; m<6; ++m){
				if(pInside){
					if(inside[m])
						if(dist[m] < df)
							df = dist[m];
				}
				else{
					if(!inside[m])
						if(dist[m] < df)
							df = dist[m];
				}
				if(index == 0 && i == I && j == J && k == K){
					printf(" dist[%d] = %f \n", m, dist[m]);
					if(inside[m])
						printf(" inside[%d] is true \n", m);
					else
						printf(" inside[%d] is false \n", m);
				}
			}
			if(!pInside){
				df = -1 * df;
				outsideflag[INDEX(i,j,k)] = 1;
			}
			phi[INDEX(i,j,k)] = df;

//			if(index == 0 && pInside){
			if(index == 0 && i == I && j == J && k == K){
				printf(" at cell (%d, %d, %d) with pos (%f, %f, %f) phi = %f\n ",
					i,j,k, p0.x, p0.y, p0.z, phi[INDEX(i,j,k)]);
			}
	END_FOR

	unsigned int NM = 0;
	FOR_EACH_CELL
			if(!voxel->InSolid(i,j,k)){
				if(phi[INDEX(i,j,k)] > 0.f){
//					printf("CELL (%d, %d, %d) inside the mesh \n", i, j, k);
					++NM;
				}
			}
			if(index == 0 && i == I && j == J && k == K){
				printf("at (%d,%d,%d) phi = %f \n", i,j,k, phi[INDEX(i,j,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i-1,j,k, phi[INDEX(i-1,j,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i+1,j,k, phi[INDEX(i+1,j,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j-1,k, phi[INDEX(i,j-1,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j+1,k, phi[INDEX(i,j+1,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j,k-1, phi[INDEX(i,j,k-1)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j,k+1, phi[INDEX(i,j,k+1)]);
			}
	END_FOR
//	if(index == 0){
//		for(int k=0; k<DimZ; ++k)
//			printf(" at (%d, %d, %d) phi = %f \n", I, J, k, phi[INDEX(I,J,k)]);
//	}
	printf("\n Before floodfill there are %u cells inside the mesh for index %d \n", NM, index);


	char *color = new char[DimX*DimY*DimZ];
	memset(color, 0, DimX*DimY*DimZ*sizeof(char));
	FOR_EACH_CELL
		if(phi[INDEX(i,j,k)] >= 0.f)
			color[INDEX(i,j,k)] = 2;
		if(voxel->InSolid(i,j,k))
			color[INDEX(i,j,k)] = 2;
	END_FOR
	//TentativeGridPoint floodSeed(0.f, 159, 190, 94);
	TentativeGridPoint floodSeed(0.f, DimX/2, 190, DimZ-(WALL_THICKNESS+1)-1);
	FloodFillSolid(floodSeed.ii, floodSeed.jj, floodSeed.kk, color);
	u_int numEnslavePoints = 0;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(outsideflag[pos] == 1 && color[pos] == 0){
			phi[pos] *= -1;
			++numEnslavePoints;
		}
	END_FOR
	delete [] color;
	delete [] outsideflag;
	printf("\n After floodfill found %u cells enslaved by mesh for index %d \n",
			numEnslavePoints, index);

	NM = 0;
	FOR_EACH_CELL
			if(!voxel->InSolid(i,j,k)){
				if(phi[INDEX(i,j,k)] > 0.f){
//					printf("CELL (%d, %d, %d) inside the mesh \n", i, j, k);
					++NM;
				}
			}
			if(index == 0 && i == I && j == J && k == K){
				printf("at (%d,%d,%d) phi = %f \n", i,j,k, phi[INDEX(i,j,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i-1,j,k, phi[INDEX(i-1,j,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i+1,j,k, phi[INDEX(i+1,j,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j-1,k, phi[INDEX(i,j-1,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j+1,k, phi[INDEX(i,j+1,k)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j,k-1, phi[INDEX(i,j,k-1)]);
				printf("at (%d,%d,%d) phi = %f \n", i,j,k+1, phi[INDEX(i,j,k+1)]);
			}
	END_FOR
//	if(index == 0){
//		for(int k=0; k<DimZ; ++k)
//			printf(" at (%d, %d, %d) phi = %f \n", I, J, k, phi[INDEX(I,J,k)]);
//	}
	printf("\n After floodfill there are %u cells inside the mesh for index %d \n", NM, index);
}

void PhysicsWorld::LabelCellsCloseToMesh(const Voxel *voxel, int num_verts, const double *verts,
		int num_tris, const int *tris){
#ifdef WIN32
	hash_map<u_int, char>::iterator iter_cell;
#else
	std::tr1::unordered_map<u_int, char>::iterator iter_cell;
#endif
	if(num_tris == 0) return;
	const Vec3d *x = (const Vec3d *)verts;
	const Vec3i *tri = (const Vec3i *)tris;
	std::vector<BoundingBox> box(num_tris);
	for(int t=0; t<num_tris; ++t){
		int a, b, c; assign(tri[t], a, b, c);
		box[t].build_from_points(x[a], x[b], x[c]);
		box[t].enlarge(1.e-5);
		Point p0(box[t].xmin[0], box[t].xmin[1], box[t].xmin[2]);
		Point p1(box[t].xmax[0], box[t].xmax[1], box[t].xmax[2]);
		TentativeGridPoint tp0 = voxel->ContainsPoint(p0);
		TentativeGridPoint tp1 = voxel->ContainsPoint(p1);
		if(IndexInBounds(tp0) && IndexInBounds(tp1)){
			for(int k=tp0.kk; k<=tp1.kk; ++k)
				for(int j=tp0.jj; j<=tp1.jj; ++j)
					for(int i=tp0.ii; i<=tp1.ii; ++i){
						unsigned int index = k*DimX*DimY+j*DimX+i;
						iter_cell = mCloseToMeshCells.find(index);
						if(iter_cell == mCloseToMeshCells.end())
							mCloseToMeshCells.insert(make_pair(index, (char)1));
					}
		}
	}
	printf("\n there are %u lablled cells \n", mCloseToMeshCells.size());

}

bool PhysicsWorld::IsCellCloseToMesh(int i, int j, int k){
#ifdef WIN32
	hash_map<u_int, char>::iterator iter_cell;
#else
	std::tr1::unordered_map<u_int, char>::iterator iter_cell;
#endif
	unsigned int index = k*DimX*DimY+j*DimX+i;
	iter_cell = mCloseToMeshCells.find(index);
	if(iter_cell != mCloseToMeshCells.end())
		return true;
	else
		return false;
}

bool PhysicsWorld::SegmentIntersectsMesh(const Point &p0, const Point &p1, double *s, double *t) const{
	for(unsigned int n=0; n < mStaticMeshObjects.size(); ++n) {
		double P0[3], P1[3];
		double a,b,c;
		int tri;
		P0[0] = p0.x;P0[1] = p0.y;P0[2] = p0.z;
		P1[0] = p1.x;P1[1] = p1.y;P1[2] = p1.z;
		if(segment_intersects_mesh(P0, P1, mStaticMeshObjects[n],
									&tri, s, t, &a, &b, &c))
			return true;
	}
	return false;
}

bool PhysicsWorld::PointInsideMesh(const Point &p) const{
	for(unsigned int n=0; n < mStaticMeshObjects.size(); ++n){
		double P[3];
		P[0] = p.x;P[1] = p.y;P[2] = p.z;
		if(point_inside_mesh(P, mStaticMeshObjects[n]))
			return true;
	}
	return false;
}

bool PhysicsWorld::PointDistToMesh(const Point &p, float &rad, Vector &objNormal) const{
	bool hasCloser = false;
	for(unsigned int n=0; n < mStaticMeshObjects.size(); ++n){
		double P[3];
		P[0] = p.x;P[1] = p.y;P[2] = p.z;
		double N[3];
		double dist = rad;
      	if(point_to_mesh(P, dist, N, mStaticMeshObjects[n])){
		   if(dist < rad){
			rad = dist;
			objNormal.x = N[0];
			objNormal.y = N[1];
			objNormal.z = N[2];
		   }
		   hasCloser = true;
      	}
	}
	return hasCloser;
}

void PhysicsWorld::DestroyPhysicsWorld(){
	for(unsigned int n=0; n < mStaticMeshObjects.size(); ++n)
		destroy_mesh_object(mStaticMeshObjects[n]);
	for(unsigned int n=0; n < mStaticTriangleMeshObjects.size(); ++n)
		delete mStaticTriangleMeshObjects[n];
}
