/*
 * aliaswavefrontobj.cpp
 *
 *  Created on: Nov 5, 2010
 *      Author: bzhao
 */

#include <iostream>
#include "aliaswavefrontobj.h"
#include "mesh.h"
#include "geometry.h"
#include "transform.h"

void  AliasWavefrontObj::LoadObject(const char *obj_file){
	mesh *meshObj = new mesh;
	meshObj->LoadMesh(string(obj_file));
	mNumVerts = meshObj->vTotal;
	mNumTris = meshObj->fTotal;
	mVerts = new double[3*mNumVerts];
	mTris  = new int[3*mNumTris];
	for(unsigned int n=0; n<mNumVerts; ++n){
		mVerts[n*3]   = meshObj->vList[n][0];
		mVerts[n*3+1] = meshObj->vList[n][1];
		mVerts[n*3+2] = meshObj->vList[n][2];
	}
	for(unsigned int n=0; n<mNumTris; ++n){
		mTris[n*3]   = meshObj->faceList[n][0].v;
		mTris[n*3+1] = meshObj->faceList[n][1].v;
		mTris[n*3+2] = meshObj->faceList[n][2].v;
	}
	delete meshObj;
}

void AliasWavefrontObj::TransformO2W(const Transform &o2w){
	BBox bobj;
	for(unsigned int n=0; n<mNumVerts; ++n){
		Point pb((float)mVerts[n*3], (float)mVerts[n*3+1], (float)mVerts[n*3+2]);
		Point pt = o2w(pb);
		bobj = Union(bobj, pt);
		mVerts[n*3]   = pt.x;
		mVerts[n*3+1] = pt.y;
		mVerts[n*3+2] = pt.z;
	}

	printf("world bound box  at (%f,%f,%f) -> (%f,%f,%f) \n",
    bobj.pMin.x, bobj.pMin.y, bobj.pMin.z,
    bobj.pMax.x, bobj.pMax.y, bobj.pMax.z);
    printf("object to world transform is: \n");
    std::cout << o2w << std::endl;

}

int  AliasWavefrontObj::NumVertices() const{
	return mNumVerts;
}

const double* AliasWavefrontObj::Vertices()  const{
	return mVerts;
}

int AliasWavefrontObj::NumTriangles()    const{
	return mNumTris;
}

const int*  AliasWavefrontObj::TriangleIndices() const{
	return mTris;
}
