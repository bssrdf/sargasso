/*
 * targetshape.h
 *
 *  Created on: Mar 6, 2009
 *      Author: bzhao
 */

#ifndef TARGETSHAPE_H_
#define TARGETSHAPE_H_

#include "voxel.h"

class TargetShape{
public:
	TargetShape(int X, int Y, int Z, float a, float b, float c, int g)
	: DimX(X), DimY(Y), DimZ(Z),
	  alpha(a), beta(b), C(c), gamma(g){
		  phi = new float[X*Y*Z];
		  memset(phi, 0, X*Y*Z*sizeof(float));
	}

	virtual ~TargetShape(){
		  delete [] phi;
	}

	void ControlForce(const Voxel &voxel, const float *u, const float *v, const float *w, const float *phic,
				float *u0, float *v0, float *w0) const{
		VelocityFeedbackForce(voxel, u, v, w, u0, v0, w0);
		ShapeFeedbackForce(voxel, phic, u0, v0, w0);
		GeometricPotentialForce(voxel, u0, v0, w0);
	 }

	virtual void VelocityFeedbackForce(const Voxel &voxel, const float *u, const float *v, const float *w,
			  float *u0, float *v0, float *w0) const = 0;
	virtual void ShapeFeedbackForce(const Voxel &voxel, const float *phic, float *u0, float *v0, float *w0) const = 0;
	virtual void GeometricPotentialForce(const Voxel &voxel, float *u0, float *v0, float *w0) const = 0;

protected:
	float *phi;
	int DimX, DimY, DimZ;
	float alpha, beta, C;
	int gamma;
};


#endif /* TARGETSHAPE_H_ */
