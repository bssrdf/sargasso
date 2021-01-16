/*
 * aliaswavefrontobj.h
 *
 *  Created on: Nov 5, 2010
 *      Author: bzhao
 */

#ifndef ALIASWAVEFRONTOBJ_H_
#define ALIASWAVEFRONTOBJ_H_

#include "trianglemeshobject.h"

class Transform;

class AliasWavefrontObj : public TriangleMeshObject{
public:
	AliasWavefrontObj(){

	}
	virtual ~AliasWavefrontObj(){
		if(mVerts)
			delete mVerts;
		if(mTris)
			delete mTris;
	}

	void          LoadObject(const char *s);
	void          TransformO2W(const Transform &o2w);
	int           NumVertices()     const;
	const double* Vertices()        const;
	int           NumTriangles()    const;
	const int*    TriangleIndices() const;

	double *mVerts;
	int    *mTris;

private:
	int mNumVerts;
	int mNumTris;

};
#endif /* ALIASWAVEFRONTOBJ_H_ */
