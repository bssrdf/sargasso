/*
 * physicsworld.h
 *
 *  Created on: Nov 4, 2010
 *      Author: bzhao
 */

#ifndef PHYSICSWORLD_H_
#define PHYSICSWORLD_H_

#include <vector>
#ifdef WIN32
#include <hash_map>
using namespace stdext;
#else
#include <tr1/unordered_map>
//using namespace std::tr1;
#endif
#include "mesh_query.h"
#include "geometry.h"


class TriangleMeshObject;
class Voxel;

class PhysicsWorld{
public:
	PhysicsWorld(const Point &p0, const Point &p1, int nx, int ny, int nz)
	: mPmin(p0), mPmax(p1), DimX(nx), DimY(ny), DimZ(nz){

	}

	~PhysicsWorld(){
		DestroyPhysicsWorld();
	}

	void BuildPhysicsWorld(const char *descfile, const Voxel *voxel);


	bool SegmentIntersectsMesh(const Point &p0, const Point &p1, double *s, double *t) const;

	bool PointInsideMesh(const Point &p) const;

	bool PointDistToMesh(const Point &p, float &rad, Vector &objNormal) const;

	void DestroyPhysicsWorld();

	void BuildSignedDistanceMap(const Voxel *voxel, int ind, float *phi) const;

	bool IsCellCloseToMesh(int i, int j, int k);
	void LabelCellsCloseToMesh(const Voxel *voxel, int num_verts, const double *verts,
			int num_tris, const int *tris);
private:
	inline bool IndexInBounds(const TentativeGridPoint &tp) const;
	void FloodFillSolid(int ii, int jj, int kk, char *color) const;
	Point mPmin, mPmax;
	int DimX;
	int DimY;
	int DimZ;
	std::vector<MeshObject*> mStaticMeshObjects;
	std::vector<TriangleMeshObject*> mStaticTriangleMeshObjects;
#ifdef WIN32
	hash_map<u_int, char> mCloseToMeshCells;
#else
	std::tr1::unordered_map<u_int, char> mCloseToMeshCells;
#endif

};

#endif /* PHYSICSWORLD_H_ */
