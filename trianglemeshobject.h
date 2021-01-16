/*
 * trianglemeshobject.h
 *
 *  Created on: Nov 5, 2010
 *      Author: bzhao
 */

#ifndef TRIANGLEMESHOBJECT_H_
#define TRIANGLEMESHOBJECT_H_

class Transform;

class TriangleMeshObject{
public:
	TriangleMeshObject(){
	}
	virtual ~TriangleMeshObject(){
	}
	virtual void          LoadObject(const char *s)   = 0;
	virtual void          TransformO2W(const Transform &o2w) = 0;
	virtual int           NumVertices()     const = 0;
	virtual const double* Vertices()        const = 0;
	virtual int           NumTriangles()    const = 0;
	virtual const int*    TriangleIndices() const = 0;

};

#endif /* TRIANGLEMESHOBJECT_H_ */
