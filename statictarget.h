/*
 * statictarget.h
 *
 *  Created on: Mar 6, 2009
 *      Author: bzhao
 */

#ifndef STATICTARGET_H_
#define STATICTARGET_H_

#include <map>
using std::map;
using std::make_pair;
#include "targetshape.h"

class StaticTargetShape : public TargetShape{
public:
	StaticTargetShape(int X, int Y, int Z, float a, float b, float c, int g, const char *filename)
	: TargetShape(X, Y, Z, a, b, c, g){
		ReadinPhi(filename);
	}
	void VelocityFeedbackForce(const Voxel &voxel, const float *u, const float *v, const float *w,
				  float *u0, float *v0, float *w0) const;
	void ShapeFeedbackForce(const Voxel &voxel, const float *phic, float *u0, float *v0, float *w0) const;
	void GeometricPotentialForce(const Voxel &voxel, float *u0, float *v0, float *w0) const;

private:
	void ReadinPhi(const char *file);
	void LinearSolver( const Voxel &voxel,
		  			 float * x, map<u_int, u_int> &a,
					 map<u_int, float> &x0, u_int n ) const;
};


#endif /* STATICTARGET_H_ */
