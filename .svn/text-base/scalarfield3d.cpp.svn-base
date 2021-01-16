#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>
//using namespace std;
using std::map;

#include "scalarfield3d.h"
#include "kdtree.h"
#include "mtrand.h"
#include "wallclocktime.h"

#ifdef CUDA
#include "reinit_cuda.h"
#endif

#define POS(x) (x).kk*(DimX*DimY)+(x).jj*DimX+(x).ii

struct FloatPair{
	FloatPair(){
		nFoundParticles = 0;
		smallestR = INFINITY;
		value = INFINITY;
	}
	mutable vector<Particle> foundParticles;
	mutable u_int nFoundParticles;
	mutable float smallestR, value;
	void operator()(const Particle &p,
					 float distSquare,
					 float &maxDistSquared) const;
	bool Empty() const{
		return foundParticles.empty();
	}
	u_int Size() const {
		return (unsigned int)foundParticles.size();
	}
	void Clear(){
		nFoundParticles = 0;
		smallestR = INFINITY;
		value = INFINITY;
		foundParticles.clear();
	}
};

void FloatPair::operator()(const Particle &p,
							float distSquare,
							float &maxDistSquared) const {
	foundParticles.push_back(p);
	nFoundParticles++;
	float dist = sqrt(distSquare);
	if( dist <= smallestR){
		smallestR = dist;
		value = p.EvaluatePhiP(smallestR);
	}
	return;
}

void ScalarField3D::Print() const{
	for(int i=0; i<DimX*DimY*DimZ; i++)
			cout << phi[i] << endl;
}

float ScalarField3D::TriInterp(const Voxel &voxel, const Point &p) const{

	float delta = voxel.VoxelDelta();
	float f[8];
	float x = p.x/delta - 0.5f;
    float y = p.y/delta - 0.5f;
    float z = p.z/delta - 0.5f;
    if(x < 0) x = 0;
    if(y < 0) y = 0;
    if(z < 0) z = 0;
	int xi = floor(x);
	int yi = floor(y);
	int zi = floor(z);
	//Point point = voxel.VoxelCenterPosition(xi,yi,zi);
	Point point = voxel.VoxelCornerPosition(xi,yi,zi,1);
	float r = x - point.x/delta;
    float s = y - point.y/delta;
    float t = z - point.z/delta;
    int xj, yj, zj;
    if( xi >= DimX-1 ){
    	xi = DimX-1;
    	xj = xi;
    	r = 1.f;
    }
    else
    	xj = xi + 1;
    if( yi >= DimY-1 ){
    	yi = DimY-1;
        yj = yi;
        s = 1.f;
    }
    else
    	yj = yi + 1;
    if( zi >= DimZ-1 ){
    	zi = DimZ-1;
        zj = zi;
        t = 1.f;
    }
    else
    	zj = zi + 1;

#ifdef USE_MESH
    double o, q;

    f[0] =  phi[INDEX(xi,yi,zi)];
	f[1] =  phi[INDEX(xi,yj,zi)];
	f[2] =  phi[INDEX(xj,yi,zi)];
	f[3] =  phi[INDEX(xj,yj,zi)];
	f[4] =  phi[INDEX(xi,yi,zj)];
	f[5] =  phi[INDEX(xi,yj,zj)];
	f[6] =  phi[INDEX(xj,yi,zj)];
	f[7] =  phi[INDEX(xj,yj,zj)];

//	TentativeGridPoint tp = voxel.ContainsPoint(p);
//	if(mPhysWorld->IsCellCloseToMesh(tp.ii, tp.jj, tp.kk)){

	Point pc = voxel.VoxelCenterPosition(xi, yi, zi);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xi, yi, zi, 0, o, pc, p, f[0]);
	}
	pc = voxel.VoxelCenterPosition(xi, yj, zi);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xi, yj, zi, 1, o, pc, p, f[1]);
	}
	pc = voxel.VoxelCenterPosition(xj, yi, zi);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xj, yi, zi, 2, o, pc, p, f[2]);
	}
	pc = voxel.VoxelCenterPosition(xj, yj, zi);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xj, yj, zi, 3, o, pc, p, f[3]);
	}
	pc = voxel.VoxelCenterPosition(xi, yi, zj);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xi, yi, zj, 4, o, pc, p, f[4]);
	}
	pc = voxel.VoxelCenterPosition(xi, yj, zj);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xi, yj, zj, 5, o, pc, p, f[5]);
	}
	pc = voxel.VoxelCenterPosition(xj, yi, zj);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xj, yi, zj, 6, o, pc, p, f[6]);
	}
	pc = voxel.VoxelCenterPosition(xj, yj, zj);
	if(mPhysWorld->SegmentIntersectsMesh(p, pc, &o, &q)){
		ProcessInvisibleNeightborForLevelsetAdvection(voxel, phi,
				xj, yj, zj, 7, o, pc, p, f[7]);
	}
//	}
//    printf("\npoint at (%f,%f,%f) \n", p.x, p.y, p.z);
//    printf(" u = %f, v = %f, w = %f\n", u, v, w);
//	printf("xi = %d, yi = %d, zi = %d \n\n", xi, yi, zi);
//	return trilerp(phi[zi*(DimX*DimY)+yi*DimX+xi], phi[zi*(DimX*DimY)+yi*DimX+xi+1],
//				   phi[(zi+1)*(DimX*DimY)+yi*DimX+xi], phi[(zi+1)*(DimX*DimY)+yi*DimX+xi+1],
//				   phi[zi*(DimX*DimY)+(yi+1)*DimX+xi], phi[zi*(DimX*DimY)+(yi+1)*DimX+xi+1],
//				   phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi], phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi+1],
//				   u, v, w);
	return trilerp(f[0], f[2], f[4], f[6],
				   f[1], f[3], f[5], f[7],
				   r, s, t);
#else
	return trilerp(phi[INDEX(xi, yi, zi)], phi[INDEX(xj, yi, zi)],
				   phi[INDEX(xi, yi, zj)], phi[INDEX(xj, yi, zj)],
				   phi[INDEX(xi, yj, zi)], phi[INDEX(xj, yj, zi)],
				   phi[INDEX(xi, yj, zj)], phi[INDEX(xj, yj, zj)],
				   r, s, t);
#endif
}

float ScalarField3D::TriInterpDebug(const Voxel &voxel, const Point &p) const{

	float delta = voxel.VoxelDelta();

	float x = p.x/delta - 0.5f;
	float y = p.y/delta - 0.5f;
	float z = p.z/delta - 0.5f;
	if(x < 0) x = 0;
	if(y < 0) y = 0;
	if(z < 0) z = 0;
	int xi = floor(x);
	int yi = floor(y);
	int zi = floor(z);
	//Point point = voxel.VoxelCenterPosition(xi,yi,zi);
	Point point = voxel.VoxelCornerPosition(xi,yi,zi,1);
	float r = x - point.x/delta;
	float s = y - point.y/delta;
	float t = z - point.z/delta;
	int xj, yj, zj;
	if( xi >= DimX-1 ){
		xi = DimX-1;
		xj = xi;
		r = 1.f;
	}
	else
		xj = xi + 1;
	if( yi >= DimY-1 ){
		yi = DimY-1;
		yj = yi;
		s = 1.f;
	}
	else
		yj = yi + 1;
	if( zi >= DimZ-1 ){
		zi = DimZ-1;
		zj = zi;
		t = 1.f;
	}
	else
		zj = zi + 1;

    printf("\npoint at (%f,%f,%f) \n", p.x, p.y, p.z);
    printf(" r = %f, s = %f, t = %f\n", r, s, t);
	printf("xi = %d, yi = %d, zi = %d \n", xi, yi, zi);
	printf("xj = %d, yj = %d, zj = %d \n", xj, yj, zj);
	printf("phi[%d,%d,%d] = %f, phi[%d,%d,%d] = %f, phi[%d,%d,%d] = %f, phi[%d,%d,%d] = %f,"
		   "phi[%d,%d,%d] = %f, phi[%d,%d,%d] = %f, phi[%d,%d,%d] = %f, phi[%d,%d,%d] = %f\n",
		    xi,yi,zi,phi[INDEX(xi, yi, zi)],
		    xj,yi,zi,phi[INDEX(xj, yi, zi)],
		    xj,yi,zj,phi[INDEX(xi, yi, zj)],
		    xj,yi,zj,phi[INDEX(xj, yi, zj)],
		    xi,yj,zi,phi[INDEX(xi, yj, zi)],
		    xj,yj,zi,phi[INDEX(xj, yj, zi)],
		    xi,yj,zj,phi[INDEX(xi, yj, zj)],
		    xj,yj,zj,phi[INDEX(xj, yj, zj)]);
//	return trilerp(phi[zi*(DimX*DimY)+yi*DimX+xi], phi[zi*(DimX*DimY)+yi*DimX+xi+1],
//				   phi[(zi+1)*(DimX*DimY)+yi*DimX+xi], phi[(zi+1)*(DimX*DimY)+yi*DimX+xi+1],
//				   phi[zi*(DimX*DimY)+(yi+1)*DimX+xi], phi[zi*(DimX*DimY)+(yi+1)*DimX+xi+1],
//				   phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi], phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi+1],
//				   u, v, w);
	float phii = trilerp(phi[INDEX(xi, yi, zi)], phi[INDEX(xj, yi, zi)],
				   phi[INDEX(xi, yi, zj)], phi[INDEX(xj, yi, zj)],
				   phi[INDEX(xi, yj, zi)], phi[INDEX(xj, yj, zi)],
				   phi[INDEX(xi, yj, zj)], phi[INDEX(xj, yj, zj)],
				   r, s, t);
	printf(" phii = %f \n\n", phii);
	return phii;
//    if(floor(p.x) == I && floor(p.y) == J && floor(p.z) == K){
//    printf("\npoint at (%f,%f,%f) \n", p.x, p.y, p.z);
//    printf(" u = %f, v = %f, w = %f\n", u, v, w);
//	printf("xi = %d, yi = %d, zi = %d \n\n", xi, yi, zi);
//    }

}


//void ProcessCoefficents(float phi_minus, float phi_plus, float minvalue,
//						float &phi){
//
//	   if( phi_minus < minvalue && phi_plus < minvalue){
////		   if(phi_minus < 0.f && phi_plus < 0.f){
//		   if(phi_minus < phi_plus)
//			   phi = minvalue - phi_minus;
//		   else
//			   phi = phi_plus - minvalue;
////		   }
////	    	if(fabs(phi_minus) < fabs(phi_plus)){
////	    		phi = minvalue - phi_minus;
////	    	}
////	    	else{
////	    		phi = phi_plus - minvalue;
////	    	}
//	    }
//	    else if(phi_minus < minvalue){
//	    	phi = minvalue - phi_minus;
//	    }
//	    else if(phi_plus < minvalue){
//	    	phi = phi_plus - minvalue;
//	    }
//	    else{
//	    	phi = 0.f;
//	    }
//
//	return;
//}

void ScalarField3D::GradPhi(const Voxel &voxel, const Point &p,  Vector &g) const{

	float delta = voxel.VoxelDelta();
//	float phix, phiy, phiz;
//	float phi_c = TriInterp(voxel, Point(p.x, p.y, p.z));
//	float phixp = TriInterp(voxel, Point(p.x+0.5f*delta, p.y, p.z));
//	float phixm = TriInterp(voxel, Point(p.x-0.5f*delta, p.y, p.z));
//	ProcessCoefficents(phixm, phixp, phi_c, phix);
//	float phiyp = TriInterp(voxel, Point(p.x, p.y+0.5f*delta, p.z));
//	float phiym = TriInterp(voxel, Point(p.x, p.y-0.5f*delta, p.z));
//	ProcessCoefficents(phiym, phiyp, phi_c, phiy);
//	float phizp = TriInterp(voxel, Point(p.x, p.y, p.z+0.5f*delta));
//	float phizm = TriInterp(voxel, Point(p.x, p.y, p.z-0.5f*delta));
//	ProcessCoefficents(phizm, phizp, phi_c, phiz);
//
//	g.x = phix;
//	g.y = phiy;
//	g.z = phiz;

	g.x = TriInterp(voxel, Point(p.x+0.5f*delta, p.y, p.z)) -
		  TriInterp(voxel, Point(p.x-0.5f*delta, p.y, p.z));
	g.y = TriInterp(voxel, Point(p.x, p.y+0.5f*delta, p.z)) -
		  TriInterp(voxel, Point(p.x, p.y-0.5f*delta, p.z));
	g.z = TriInterp(voxel, Point(p.x, p.y, p.z+0.5f*delta)) -
		  TriInterp(voxel, Point(p.x, p.y, p.z-0.5f*delta));
	g /= delta;
}

void ScalarField3D::NormalAt(const Voxel &voxel, const Point &p, Vector &g) const{
	GradPhi(voxel, p, g);
	if(g.Length() < 1.e-6){
//		printf("vanishing phi at (%f, %f, %f) with length = %16.14f and "
//				"g is (%16.12f, %16.12f, %16.12f)\n",
//				p.x, p.y, p.z, g.Length(), g.x, g.y, g.z);
		MTRand mt;
		double rn1 = mt();
		double rn2 = mt();
		double rn3 = mt();
		g = Vector((float)rn1, (float)rn2, (float)rn3);
	}
	g = Normalize(g);
}



void ScalarField3D::EvaluatePhi(const vector<Particle> &particles,
								 map<u_int, KnownPoint> &band,
								 const Voxel &v){

	// 1. Evaluating phi only on grid points near interface
//	map<u_int, KnownPoint>::iterator pos;
//	for(pos=band.begin(); pos!=band.end(); ++pos){
//		KnownPoint &p = pos->second;
//		u_int index = pos->first;
//		phi[index] = -0.5f;
//	}


//	2. An alternative way to evaluate phi, from Foster and Fedkiw (2001)

	FloatPair fp;
	KdTree<Particle, FloatPair> dist(particles);

	float vd = v.VoxelDelta();

	float startRadius = vd;
	float startRadiusSquared = vd*vd;

	for(int k=0;k<DimZ; k++){
		cout << "k= " << k << endl;
		for(int j=0;j<DimY; j++)
			for(int i=0;i<DimX; i++){
				Point p = v.VoxelCenterPosition(i,j,k);
//				int nn = 0;
				do{
					dist.Lookup(p, fp, startRadiusSquared);
					if( fp.Empty() ){
						startRadius += 5*vd;
						startRadiusSquared = startRadius * startRadius;
					}
//					nn++;
//					cout << "k = " << k << " j = " << j <<
//					 " i = " << i << " nn = " << nn << endl;
				}while(fp.nFoundParticles == 0);

				phi[INDEX(i,j,k)] = fp.value;
				if(phi[INDEX(i,j,k)] < 0.f){
					printf("found negative phi = %f \n", fp.value);
					phi[INDEX(i,j,k)] *= -1;
				}
				startRadius = vd;
				startRadiusSquared = vd*vd;
				fp.Clear();
			}
	}
//	int I = 10, J = 20, K = 10;
//	printf("A: phi = %f \n", phi[INDEX(I,J,K)]);
	return;
}

void ScalarField3D::EvaluatePhi(const Voxel &voxel, const Point &pos, float radius){

	float h = voxel.VoxelDelta();

	FOR_EACH_CELL

		Point p = voxel.VoxelCenterPosition(i,j,k);

        float val1 = ((pos - p).Length() - radius) * h;

//Double val3 = (pos - Vector(i,j)).Length() - (radius - (radius*0.2));
//val1 = max(-val3, val1);

        float leftWall = (pos[0] - radius/6);// * h;
        float rightWall = (pos[0] + radius/6);// * h;
        float topWall = (pos[1] + radius/1.5);// * h;
        float bottomWall = (pos[1] - radius);// * h;

        float val2 = INFINITY;

        if(bottomWall <= p[1] && p[1] <= topWall && leftWall <= p[0] && p[0] <= rightWall)
        {
                float top = topWall - p[1];
                float bottom = p[1] - bottomWall;
                float left = p[0] - leftWall;
                float right = rightWall - p[0];
                if(fabs(val2) > fabs(top))        val2 = top;
                if(fabs(val2) > fabs(bottom)) val2 = bottom;
                if(fabs(val2) > fabs(left))       val2 = left;
                if(fabs(val2) > fabs(right))      val2 = right;
        }
        else if(leftWall <= p[0] && p[0] <= rightWall)
        {
                if(fabs(val2) > fabs(p[1]-bottomWall))       val2 = p[1]-bottomWall;
                if(fabs(val2) > fabs(p[1]-topWall))          val2 = topWall-p[1];
        }
        else if(bottomWall <= p[1] && p[1] <= topWall)
        {
                if(fabs(val2) > fabs(p[0]-leftWall))         val2 = p[0]-leftWall;
                if(fabs(val2) > fabs(p[0]-rightWall))        val2 = rightWall-p[0];
        }
        else
        {
                float ul = sqrt((p[0]-leftWall)*(p[0]-leftWall)+(p[1]-topWall)*(p[1]-topWall));
                float ur = sqrt((p[0]-rightWall)*(p[0]-rightWall)+(p[1]-topWall)*(p[1]-topWall));
                float bl = sqrt((p[0]-leftWall)*(p[0]-leftWall)+(p[1]-bottomWall)*(p[1]-bottomWall));
                float br = sqrt((p[0]-rightWall)*(p[0]-rightWall)+(p[1]-bottomWall)*(p[1]-bottomWall));
                val2 = -min(min(min(ul,ur),bl),br);
        }
        val2 *= h;

        //val = 10*cos(i/Double(size[0])*10) + 10*sin(j/Double(size[1])*10) ;
        //gridPhi(i,j,k) = val1;
        phi[INDEX(i,j,k)] = max(val1,val2);
        if(i==I && j==J && k==K){
			printf("at cell (%d,%d,%d) phi = %f val1 = %f, val2 = %f\n",
					i,j,k, phi[INDEX(i,j,k)], val1, val2);
			printf(" p at (%f,%f,%f)\n", p.x, p.y, p.z);
			printf("leftwall = %f, rightwall = %f \n", leftWall, rightWall);
			printf("bottomwall = %f, topwall = %f \n", bottomWall, topWall);
		}

    END_FOR

}
void ScalarField3D::EvaluatePhi(const Voxel &voxel,
								int start_x0, int start_y0, int start_z0,
								int end_x0,  int end_y0,   int end_z0){

	float h = voxel.VoxelDelta();

	float leftWall = voxel.FacePostion(1, false, start_x0, start_y0, start_z0);
    float rightWall = voxel.FacePostion(1, true, end_x0, start_y0, start_z0);
    float backWall = voxel.FacePostion(2, false, start_x0, start_y0, start_z0);
    float frontWall = voxel.FacePostion(2, true, start_x0, end_y0, start_z0);
    float bottomWall =  voxel.FacePostion(3, false, start_x0, start_y0, start_z0);
    float topWall =  voxel.FacePostion(3, true, start_x0, start_y0, end_z0);

    printf("left wall at %f, right wall at %f \n",leftWall, rightWall);
    printf("back wall at %f, front wall at %f \n",backWall, frontWall);
    printf("bottom wall at %f, top wall at %f \n",bottomWall, topWall);

	FOR_EACH_CELL
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
	       phi[INDEX(i,j,k)] = d;
	    }
		else if(p[0] < leftWall)
			phi[INDEX(i,j,k)] = leftWall - p[0];
		else if(p[0] > rightWall)
			phi[INDEX(i,j,k)] = p[0] - rightWall;
		else if(p[1] < backWall)
			phi[INDEX(i,j,k)] = backWall - p[1];
		else if(p[1] > frontWall)
			phi[INDEX(i,j,k)] = p[1] - frontWall;
		else if(p[2] < bottomWall)
			phi[INDEX(i,j,k)] = bottomWall - p[2];
		else if(p[2] > topWall)
			phi[INDEX(i,j,k)] = p[2] - topWall;
		if( i == I && j == J && k == K)
			printf("p at (%f,%f,%f) with phi = %f \n", p[0], p[1], p[2], phi[INDEX(i,j,k)]);
	END_FOR
}

void ScalarField3D::EvaluatePhiAirWater(const Voxel &voxel,
								int start_x0, int start_y0, int start_z0,
								int end_x0,  int end_y0,   int end_z0){

	float h = voxel.VoxelDelta();

	float leftWall = voxel.FacePostion(1, false, start_x0, start_y0, start_z0);
    float rightWall = voxel.FacePostion(1, true, end_x0, start_y0, start_z0);
    float backWall = voxel.FacePostion(2, false, start_x0, start_y0, start_z0);
    float frontWall = voxel.FacePostion(2, true, start_x0, end_y0, start_z0);
    float bottomWall =  voxel.FacePostion(3, false, start_x0, start_y0, start_z0);
    float topWall =  voxel.FacePostion(3, true, start_x0, start_y0, end_z0);

    printf("left wall at %f, right wall at %f \n",leftWall, rightWall);
    printf("back wall at %f, front wall at %f \n",backWall, frontWall);
    printf("bottom wall at %f, top wall at %f \n",bottomWall, topWall);

	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k)){

			Point p = voxel.VoxelCenterPosition(i,j,k);

			if( bottomWall <= p[2] && p[2] <= topWall
				&& leftWall <= p[0] && p[0] <= rightWall
				&& backWall <= p[1] && p[1] <= frontWall) {
// 				falling water
//		       float dist[5];
//		       dist[0] = p[0] - leftWall;
//		       dist[1] = rightWall - p[0];
//		       dist[2] = p[1] - backWall;
//		       dist[3] = frontWall - p[1];
//		       dist[4] = topWall - p[2];
//		       float d = FindMinimum(dist, 5);
//		       phi[INDEX(i,j,k)] = d;

				// pouring water
		       float dist[6];
   		       dist[0] = p[0] - leftWall;
   		       dist[1] = rightWall - p[0];
   		       dist[2] = p[1] - backWall;
   		       dist[3] = frontWall - p[1];
   		       dist[4] = p[2] - bottomWall;
   		       dist[5] = topWall - p[2];
   		       float d = FindMinimum(dist, 6);
   		       phi[INDEX(i,j,k)] = -d;
		    }
			else if(p[0] < leftWall)
				phi[INDEX(i,j,k)] = leftWall - p[0];
			else if(p[0] > rightWall)
				phi[INDEX(i,j,k)] = p[0] - rightWall;
			else if(p[1] < backWall)
				phi[INDEX(i,j,k)] = backWall - p[1];
			else if(p[1] > frontWall)
				phi[INDEX(i,j,k)] = p[1] - frontWall;
			else if(p[2] < bottomWall)
				phi[INDEX(i,j,k)] = bottomWall - p[2];
			else if(p[2] > topWall)
				phi[INDEX(i,j,k)] = p[2] - topWall;
//		}
//		printf("p at (%f,%f,%f) with phi = %f \n", p[0], p[1], p[2], phi[INDEX(i,j,k)]);
	END_FOR
}
void ScalarField3D::EvaluatePhiAirWater(const Voxel &voxel,
								int start_x0, int start_y0, int start_z0,
								int end_x0,  int end_y0,   int end_z0,
								int start_x1, int start_y1, int start_z1,
								int end_x1,  int end_y1,   int end_z1){

	float h = voxel.VoxelDelta();

	float leftWall = voxel.FacePostion(1, false, start_x0, start_y0, start_z0);
    float rightWall = voxel.FacePostion(1, true, end_x0, start_y0, start_z0);
    float backWall = voxel.FacePostion(2, false, start_x0, start_y0, start_z0);
    float frontWall = voxel.FacePostion(2, true, start_x0, end_y0, start_z0);
    float bottomWall =  voxel.FacePostion(3, false, start_x0, start_y0, start_z0);
    float topWall =  voxel.FacePostion(3, true, start_x0, start_y0, end_z0);

    printf("left wall at %f, right wall at %f \n",leftWall, rightWall);
    printf("back wall at %f, front wall at %f \n",backWall, frontWall);
    printf("bottom wall at %f, top wall at %f \n",bottomWall, topWall);

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){

			Point p = voxel.VoxelCenterPosition(i,j,k);

			if( bottomWall <= p[2] && p[2] <= topWall
				&& leftWall <= p[0] && p[0] <= rightWall
				&& backWall <= p[1] && p[1] <= frontWall) {
// 				falling water
		       float dist[5];
		       dist[0] = p[0] - leftWall;
		       dist[1] = rightWall - p[0];
		       dist[2] = p[1] - backWall;
		       dist[3] = frontWall - p[1];
		       dist[4] = topWall - p[2];
		       float d = FindMinimum(dist, 5);
		       phi[INDEX(i,j,k)] = d;

				// pouring water
//		       float dist[6];
//   		       dist[0] = p[0] - leftWall;
//   		       dist[1] = rightWall - p[0];
//   		       dist[2] = p[1] - backWall;
//   		       dist[3] = frontWall - p[1];
//   		       dist[4] = p[2] - bottomWall;
//   		       dist[5] = topWall - p[2];
//   		       float d = FindMinimum(dist, 6);
//   		       phi[INDEX(i,j,k)] = d;
		    }
			else if(p[0] < leftWall)
				phi[INDEX(i,j,k)] = leftWall - p[0];
			else if(p[0] > rightWall)
				phi[INDEX(i,j,k)] = p[0] - rightWall;
			else if(p[1] < backWall)
				phi[INDEX(i,j,k)] = backWall - p[1];
			else if(p[1] > frontWall)
				phi[INDEX(i,j,k)] = p[1] - frontWall;
			else if(p[2] < bottomWall)
				phi[INDEX(i,j,k)] = bottomWall - p[2];
			else if(p[2] > topWall)
				phi[INDEX(i,j,k)] = p[2] - topWall;
		}
//		printf("p at (%f,%f,%f) with phi = %f \n", p[0], p[1], p[2], phi[INDEX(i,j,k)]);
	END_FOR

	leftWall = voxel.FacePostion(1, false, start_x1, start_y1, start_z1);
	rightWall = voxel.FacePostion(1, true, end_x1, start_y1, start_z1);
	backWall = voxel.FacePostion(2, false, start_x1, start_y1, start_z1);
	frontWall = voxel.FacePostion(2, true, start_x1, end_y1, start_z1);
	bottomWall =  voxel.FacePostion(3, false, start_x1, start_y1, start_z1);
	topWall =  voxel.FacePostion(3, true, start_x1, start_y1, end_z1);

	printf("left wall at %f, right wall at %f \n",leftWall, rightWall);
	printf("back wall at %f, front wall at %f \n",backWall, frontWall);
	printf("bottom wall at %f, top wall at %f \n",bottomWall, topWall);

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){

			Point p = voxel.VoxelCenterPosition(i,j,k);

			if( bottomWall <= p[2] && p[2] <= topWall
				&& leftWall <= p[0] && p[0] <= rightWall
				&& backWall <= p[1] && p[1] <= frontWall) {
// 				falling water
//			   float dist[5];
//			   dist[0] = p[0] - leftWall;
//			   dist[1] = rightWall - p[0];
//			   dist[2] = p[1] - backWall;
//			   dist[3] = frontWall - p[1];
//			   dist[4] = topWall - p[2];
//			   float d = FindMinimum(dist, 5);
//			   phi0[INDEX(i,j,k)] = d;

				// pouring water
		       float dist[6];
   		       dist[0] = p[0] - leftWall;
   		       dist[1] = rightWall - p[0];
   		       dist[2] = p[1] - backWall;
   		       dist[3] = frontWall - p[1];
   		       dist[4] = p[2] - bottomWall;
   		       dist[5] = topWall - p[2];
   		       float d = FindMinimum(dist, 6);
   		       phi0[INDEX(i,j,k)] = d;
			}
			else if(p[0] < leftWall)
				phi0[INDEX(i,j,k)] = leftWall - p[0];
			else if(p[0] > rightWall)
				phi0[INDEX(i,j,k)] = p[0] - rightWall;
			else if(p[1] < backWall)
				phi0[INDEX(i,j,k)] = backWall - p[1];
			else if(p[1] > frontWall)
				phi0[INDEX(i,j,k)] = p[1] - frontWall;
			else if(p[2] < bottomWall)
				phi0[INDEX(i,j,k)] = bottomWall - p[2];
			else if(p[2] > topWall)
				phi0[INDEX(i,j,k)] = p[2] - topWall;
		}
	END_FOR

	FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k)){
				phi[INDEX(i,j,k)] = phi[INDEX(i,j,k)] < phi0[INDEX(i,j,k)] ?
						phi[INDEX(i,j,k)] : phi0[INDEX(i,j,k)];
			}
	END_FOR

}

void ScalarField3D::EvaluatePhiAirWater(const Voxel &voxel,float d, float H){
	float h = voxel.VoxelDelta();
	int i,j,k;
	float *eta = new float[DimX];
	for(i=0; i<=DimX-1; ++i){
		Point p = voxel.VoxelCenterPosition(i,0,0);
		float tmp = 3*H/(4*d*d*d);
		eta[i] = d + H / (cosh(sqrt(tmp)*p.x)*cosh(sqrt(tmp)*p.x));
 	}
	// first along z-direction
	j = 0;
	for(i=0; i<=DimX-1; ++i){
		for(k=0; k<=DimZ-1; ++k){
			Point p = voxel.VoxelCenterPosition(i,j,k);
			if(p.z <= eta[i])
				phi[INDEX(i,j,k)] = p.z - eta[i];
			else
				phi[INDEX(i,j,k)] = eta[i] - p.z;
		}
	}
	// second do correction along x-direction
	for(k=0; k<=DimZ-1; ++k){
		float tmp = phi[INDEX(0,j,k)];
		for(i=1; i<=DimX-1; ++i){
			if(tmp * phi[INDEX(i,j,k)] < 0.f)
				break;
			tmp = phi[INDEX(i,j,k)];
		}
		float dtmp1 = h * fabs(tmp) / (fabs(tmp) + fabs(phi[INDEX(i,j,k)]));
		float dtmp2 = h * fabs(phi[INDEX(i,j,k)]) / (fabs(tmp) + fabs(phi[INDEX(i,j,k)]));
		if(fabs(phi[INDEX(i,j,k)]) < dtmp2){
			if(phi[INDEX(i,j,k)] > 0.f)
				phi[INDEX(i,j,k)] = dtmp2;
			else
				phi[INDEX(i,j,k)] = -dtmp2;
		}
		if(fabs(phi[INDEX(i-1,j,k)]) < dtmp1){
			if(phi[INDEX(i-1,j,k)] > 0.f)
				phi[INDEX(i-1,j,k)] = dtmp1;
			else
				phi[INDEX(i-1,j,k)] = -dtmp1;
		}
	}
	for(j=1; j<=DimY-1; ++j){
		for(i=0; i<=DimX-1; ++i){
			for(k=0; k<=DimZ-1; ++k){
					phi[INDEX(i,j,k)] = phi[INDEX(i,0,k)];
			}
		}
	}
	delete [] eta;
}

void ScalarField3D::OutputGridData(const Voxel &voxel, char *filename) const{
	FILE *fp=fopen(filename, "w");
	if(!fp){
      printf("Couldn't open file to write \n");
      exit(-1);
   	}
   	fprintf(fp,"%d %d %d \n",DimX, DimY, DimZ);
   	Point p0 = voxel.GetP0();
   	fprintf(fp,"%f %f %f \n",p0.x, p0.y, p0.z);
   	Point p1 = voxel.GetP1();
    fprintf(fp,"%f %f %f \n",p1.x, p1.y, p1.z);
    fprintf(fp,"%f \n",voxel.VoxelDelta());
    for(int k=0;k<DimZ;k++)
		for(int j=0;j<DimY;j++)
			for(int i=0;i<DimX;i++)
				fprintf(fp,"%f\n", phi[k*(DimX*DimY)+j*DimX+i]);
	fclose(fp);

}

void ScalarField3D::OutputBinaryGridData(const Voxel &voxel, char *filename) const{
	FILE *fp = fopen(filename, "wb");
	if(!fp){
      printf("Couldn't open file [ %s ] to write \n", filename);
      exit(-1);
   	}
//   	fprintf(fp,"%d %d %d \n", DimX, DimY, DimZ);
   	fwrite(&DimX, sizeof(int), 1, fp);
   	fwrite(&DimY, sizeof(int), 1, fp);
   	fwrite(&DimZ, sizeof(int), 1, fp);
   	Point p0 = voxel.GetP0();
//   	fprintf(fp,"%f %f %f \n",p0.x, p0.y, p0.z);
   	fwrite(&(p0.x), sizeof(float), 1, fp);
   	fwrite(&(p0.y), sizeof(float), 1, fp);
   	fwrite(&(p0.z), sizeof(float), 1, fp);
   	Point p1 = voxel.GetP1();
//    fprintf(fp,"%f %f %f \n",p1.x, p1.y, p1.z);
    fwrite(&(p1.x), sizeof(float), 1, fp);
    fwrite(&(p1.y), sizeof(float), 1, fp);
    fwrite(&(p1.z), sizeof(float), 1, fp);
    float delta = voxel.VoxelDelta();
    fwrite(&delta, sizeof(float), 1, fp);
//    fprintf(fp,"%f \n",voxel.VoxelDelta());

    for(int k=0;k<DimZ;k++)
		for(int j=0;j<DimY;j++)
			for(int i=0;i<DimX;i++)
				fwrite(&phi[INDEX(i,j,k)], sizeof(float), 1, fp);
//				fprintf(fp,"%f\n", phi[k*(DimX*DimY)+j*DimX+i]);
	fclose(fp);

}
void ScalarField3D::ReadBinaryGridData(char *filename){
	FILE *fp=fopen(filename, "rb");
   if(!fp){
	  printf("Couldn't open gridimplicit file \"%s\" for reading\n", filename);
	  return;
   }

   // read in dimensions
   /*if(3!=fscanf(fp, "%d %d %d", dim, dim+1, dim+2)){
	  Error("Problem reading gridimplicit file \"%s\" (dimensions)\n", filename);
	  dim[0]=dim[1]=dim[2]=0;
	  return;
   }*/
   fread(&DimX, sizeof(int), 1, fp);
   fread(&DimY, sizeof(int), 1, fp);
   fread(&DimZ, sizeof(int), 1, fp);


   Point p0;
   /*if(3!=fscanf(fp, "%f %f %f", &p0.x, &p0.y, &p0.z)){
		 Error("Problem reading gridimplicit file \"%s\" (dimensions)\n", filename);
		 p0.x=p0.y=p0.z=0;
		 return;
   }*/
   fread(&p0.x, sizeof(float), 1, fp);
   fread(&p0.y, sizeof(float), 1, fp);
   fread(&p0.z, sizeof(float), 1, fp);
   Point p1;
   /*if(3!=fscanf(fp, "%f %f %f", &p1.x, &p1.y, &p1.z)){
	   Error("Problem reading gridimplicit file \"%s\" (dimensions)\n", filename);
	   p1.x=p1.y=p1.z=0;
	   return;
   }*/
   fread(&p1.x, sizeof(float), 1, fp);
   fread(&p1.y, sizeof(float), 1, fp);
   fread(&p1.z, sizeof(float), 1, fp);
   /*if(1!=fscanf(fp, "%f", &delta)){
	   Error("Problem reading gridimplicit file \"%s\" (dimensions)\n", filename);
	   delta = 1.f;
	   return;
	}*/
   float delta;
   fread(&delta, sizeof(float), 1, fp);
   printf("p0 at (%f,%f,%f) p1 at (%f,%f,%f) delta = %f\n",
		   p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, delta);

   // read in samples
   //for(int i=0; i<dim[0]*dim[1]*dim[2]; ++i){
   FOR_EACH_CELL
	  /*if(1!=fscanf(fp, "%f", phi+i)){
		 Error("Problem reading gridimplicit file \"%s\" (at value %d)\n", filename, i);
		 dim[0]=dim[1]=dim[2]=0;
		 delete[] phi;
		 phi=0;
		 return;
	  }*/
	  //if(1!=fread(phi+i, sizeof(float), 1, fp)){
		if(1!=fread(&phi[INDEX(i,j,k)], sizeof(float), 1, fp)){
		  printf("Problem reading gridimplicit file \"%s\" (at value %d, %d, %d)\n",
				  filename, i, j, k);
		   delete[] phi;
		  phi=0;
		  return;
	  }
//		      assert(phi[i]==phi[i]);
	END_FOR
	fclose(fp);
}

void ScalarField3D::ReadinGridData(char *filename) {
	FILE *fp=fopen(filename, "r");
	if(!fp){
      printf("Couldn't open file to read \n");
      exit(-1);
   	}
   	printf("file opened \n");
   	int dim1, dim2, dim3;
   	fscanf(fp, "%d %d %d", &dim1, &dim2, &dim3);
    printf("dimension readin \n");
    float x,y,z;
    fscanf(fp, "%f %f %f", &x, &y, &z);
    fscanf(fp, "%f %f %f", &x, &y, &z);
    fscanf(fp, "%f", &x);
    for(int i=0; i<dim1*dim2*dim3; ++i){
      if(1!=fscanf(fp, "%f", phi+i)){
         printf("Problem reading gridimplicit file \"%s\" (at value %d)\n", filename, i);
         DimX=DimY=DimZ=0;
         return;
      }
      assert(phi[i]==phi[i]);
    }
    fclose(fp);
}


void ScalarField3D::WriteRestart(FILE* fp) const{

	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++)
				fprintf(fp,"%f ", phi[INDEX(i,j,k)]);
			fprintf(fp,"\n");
		}
	}

}

bool ScalarField3D::ReadRestart(FILE* fp) {
	char c;

	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++){
				if(1!=fscanf(fp, "%f", phi+INDEX(i,j,k))){
			         printf("Problem reading gridimplicit file \n");
			         DimX=DimY=DimZ=0;
			         return false;
			      }
				fscanf(fp, "%c", &c);
			}
		}
	}
	return true;
}



void ScalarField3D::SeedParticlesAtSource(const Voxel &voxel,
						list<Particle> &negParticles,
						int start_x0, int start_y0,  int start_z0,
						int end_x0,   int end_y0,   int  end_z0,
						int nPerCell) const{

	float delta = voxel.VoxelDelta();
	MTRand mt1, mt2;

	int i, j, k;

	if(negParticles.size() > 500000)
		return;

	for(k = start_z0; k <= end_z0; ++k)
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i){
				for(int m=0; m<nPerCell; ++m){
					float phi_i = 0.f;
					Particle ps;
					Point pos;
					float radius = 0.f;
					do{
						voxel.PlaceParticles(i, j, k, ps, 0);
						pos = ps.Position();
						phi_i = TriInterp(voxel, pos);
					}while(phi_i >= 0.f);

					if( -phi_i > 0.5f*delta ) radius = 0.5f*delta;
					else if( -phi_i < 0.1f*delta)  radius = 0.1f*delta;
					else radius = -phi_i;

					//if(i==I && j==J && k == K)
		//						printf("m=%d, (%f, %f,%f) with %f  and phi_i = %f phi = %f\n",
		//							m, pos.x, pos.y, pos.z, radius, phi_i, phi[INDEX(i,j,k)]);
					ps.SetParticle(pos, radius);
					negParticles.push_back(ps);
				}
		}

}

/*void ScalarField3D::SeedParticles(const Voxel &voxel,
						list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
						list<Particle> &posParticles,
#endif
						int nPerCell) const{

	float delta = voxel.VoxelDelta();
	MTRand mt1, mt2;

//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k)){
//			if(phi[INDEX(i,j,k)] < 0.f && phi[INDEX(i,j,k)] > -3*delta){
//				for(int m=0; m<nPerCell; ++m){
//					float phi_i = 0.f;
//					Particle ps;
//					Point pos;
//					float radius = 0.f;
//					do{
//						voxel.PlaceParticles(i, j, k, ps, 0);
//						pos = ps.Position();
//						phi_i = TriInterp(voxel, pos);
//					}while(phi_i >= 0.f);
//					if( -phi_i > 0.5f*delta ) radius = 0.5f*delta;
//					else if( -phi_i < 0.1f*delta)  radius = 0.1f*delta;
//					else radius = -phi_i;
//
//					//if(i==I && j==J && k == K)
////						printf("m=%d, (%f, %f,%f) with %f  and phi_i = %f phi = %f\n",
////							m, pos.x, pos.y, pos.z, radius, phi_i, phi[INDEX(i,j,k)]);
//					ps.SetParticle(pos, radius);
//					negParticles.push_back(ps);
//				}
//			}
//		}
//
//	END_FOR

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){

			bool seed = false;

			for(int m=0; m<8; ++m){
				Point corner = voxel.VoxelCornerPosition(i,j,k,m+1);
				float phi_i = TriInterp(voxel, corner);
				if(fabs(phi_i) < 3*delta){
					seed = true;
				    break;
				}
			}
			if(seed){
#ifdef POSITIVE_PARTICLES
				for(int m=0; m<nPerCell/2; ++m){
					Particle ps;
					voxel.PlaceParticles(i, j, k, ps, 0);
					negParticles.push_back(ps);
					voxel.PlaceParticles(i, j, k, ps, 0);
					posParticles.push_back(ps);
				}
#else
				for(int m=0; m<nPerCell; ++m){
					Particle ps;
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
//					if(k < 10)
//						printf("oldpos = (%f, %f, %f) at (%d, %d, %d) with phi = %f\n",
//							pos.x, pos.y, pos.z, i,j,k, phi[INDEX(i,j,k)]);
					negParticles.push_back(ps);
//					voxel.PlaceParticles(i, j, k, ps, 0);
//					posParticles.push_back(ps);
				}
#endif
			}

		}
	END_FOR

	printf("now there are %ld negative particles \n", negParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("now there are %ld positive particles \n", posParticles.size());
#endif

	list<Particle>::iterator iter_particle;

	for(iter_particle = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		//printf("A:\n");
		Vector N;
		NormalAt(voxel, pos, N);
		//printf("B:\n");
		float phi_i = TriInterp(voxel, pos);
		//printf("C:\n");
		double rn = mt1();
		float target_phi = -(RADIUS_MIN + (float)rn * 2.9f * delta);
//		printf("oldpos = (%f, %f, %f), newphi= %f, old_phi = %f \n",
//					pos.x, pos.y, pos.z,  new_phi, phi_i);
		Point new_pos = pos + (target_phi - phi_i) * N;
		float new_phi = TriInterp(voxel, new_pos);
	//	printf("D:\n");
//		printf("newpos= (%f, %f, %f),              newphi= %f \n",
//			 new_pos.y, new_pos.z, new_phi);
		int cycles = 0;
		if(new_phi < 0.f){
			TentativeGridPoint tpp = voxel.ContainsPoint(new_pos);
			if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk)){
				float lamda = 1.0f;
				do{
				repeat1: lamda /= 2;
				    cycles++;
					new_pos = pos + lamda * (target_phi - phi_i) * N;
					new_phi = TriInterp(voxel, new_pos);
					tpp = voxel.ContainsPoint(new_pos);
				}while(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk));


				if(new_phi < 0.f){
					if( -new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( -new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = -new_phi;
					ps.SetParticle(new_pos, radius);
					++iter_particle;
				}
				else{
				repeat2:	lamda /= 2 ;
					cycles++;
					new_pos = pos + lamda * (target_phi - phi_i) * N;
					pos = new_pos;
					lamda = 1.f;
					phi_i = TriInterp(voxel, pos);
					new_pos = pos + lamda * (target_phi - phi_i) * N;
					new_phi = TriInterp(voxel, new_pos);
					tpp = voxel.ContainsPoint(new_pos);
					if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk)){
						if(cycles > 15){
							iter_particle = negParticles.erase(iter_particle);
							continue;
						}
						goto repeat1;
					}
					else{
						if(new_phi > 0.f){
							if(cycles > 15){
								iter_particle = negParticles.erase(iter_particle);
								continue;
							}
							goto repeat2;
						}
					}
				}
//				printf("oldpos = (%f, %f, %f), target phi = %f, old_phi = %f \n",
//									pos.x, pos.y, pos.z,  target_phi, phi_i);
//				printf("newpos= (%f, %f, %f),              newphi= %f \n",
//						 new_pos.x, new_pos.y, new_pos.z, new_phi);
//				printf("targetphi-phi = %f, N.length = %f, N = (%f, %f, %f) \n",
//						target_phi - phi_i, N.Length(), N.x, N.y, N.z);

			}
			else{
				if( -new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( -new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = -new_phi;
					ps.SetParticle(new_pos, radius);
					++iter_particle;
			}
		}
		else
			iter_particle = negParticles.erase(iter_particle);
	}

	printf("now there are %ld negative particles left \n", negParticles.size());

#ifdef POSITIVE_PARTICLES
	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		//printf("A:\n");
		Vector N;
		NormalAt(voxel, pos, N);
		//printf("B:\n");
		float phi_i = TriInterp(voxel, pos);
		//printf("C:\n");
		double rn = mt2();
		float target_phi = RADIUS_MIN + (float)rn * 2.9f * delta;
//		printf("oldpos = (%f, %f, %f), newphi= %f, old_phi = %f \n",
//					pos.x, pos.y, pos.z,  new_phi, phi_i);
		Point new_pos = pos + (target_phi - phi_i) * N;
		float new_phi = TriInterp(voxel, new_pos);
	//	printf("D:\n");
//		printf("newpos= (%f, %f, %f),              newphi= %f \n",
//			 new_pos.y, new_pos.z, new_phi);
		int cycles = 0;
		if(new_phi > 0.f){
			TentativeGridPoint tpp = voxel.ContainsPoint(new_pos);
			if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk)){
//				printf("oldpos = (%f, %f, %f), target phi = %f, old_phi = %f \n",
//									pos.x, pos.y, pos.z,  target_phi, phi_i);
//				printf("newpos= (%f, %f, %f),              newphi= %f \n",
//						 new_pos.x, new_pos.y, new_pos.z, new_phi);
//				printf("targetphi-phi = %f, N.length = %f, N = (%f, %f, %f) \n",
//						target_phi - phi_i, N.Length(), N.x, N.y, N.z);

				float lamda = 1.0f;
				do{
				repeat3: lamda /= 2;
				    cycles++;
					new_pos = pos + lamda * (target_phi - phi_i) * N;
					new_phi = TriInterp(voxel, new_pos);
					tpp = voxel.ContainsPoint(new_pos);
				}while(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk));


				if(new_phi > 0.f){
					if( new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = new_phi;
					ps.SetParticle(new_pos, radius);
					++iter_particle;
				}
				else{
				repeat4:	lamda /= 2 ;
					cycles++;
					new_pos = pos + lamda * (target_phi - phi_i) * N;
					pos = new_pos;
					lamda = 1.f;
					phi_i = TriInterp(voxel, pos);
					new_pos = pos + lamda * (target_phi - phi_i) * N;
					new_phi = TriInterp(voxel, new_pos);
					tpp = voxel.ContainsPoint(new_pos);
					if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk)){
						if(cycles > 15){
							iter_particle = posParticles.erase(iter_particle);
							continue;
						}
						goto repeat3;
					}
					else{
						if(new_phi < 0.f){
							if(cycles > 15){
								iter_particle = posParticles.erase(iter_particle);
								continue;
							}
							goto repeat4;
						}
					}
				}
//				printf("oldpos = (%f, %f, %f), target phi = %f, old_phi = %f \n",
//									pos.x, pos.y, pos.z,  target_phi, phi_i);
//				printf("newpos= (%f, %f, %f),              newphi= %f \n",
//						 new_pos.x, new_pos.y, new_pos.z, new_phi);
//				printf("targetphi-phi = %f, N.length = %f, N = (%f, %f, %f) \n",
//						target_phi - phi_i, N.Length(), N.x, N.y, N.z);

			}
			else{
				if( new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = new_phi;
					ps.SetParticle(new_pos, radius);
					++iter_particle;
			}
		}
		else
			iter_particle = posParticles.erase(iter_particle);

	}
	printf("now there are %ld positive particles left \n", posParticles.size());


	printf("Checking positive particles... \n");

	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		float phi_i= TriInterp(voxel, pos);
		if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk) || phi_i < 0.f)
			printf("pos at (%f, %f,%f) in cell (%d, %d, %d) with phi = %f\n",
				pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
	}
#endif

	printf("Checking negative particles... \n");

	for(iter_particle = negParticles.begin();
		iter_particle != negParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		float phi_i= TriInterp(voxel, pos);
		if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk) || phi_i > 0.f ||
				radius > RADIUS_MAX || radius < RADIUS_MIN  )
			printf("pos at (%f, %f,%f) in cell (%d, %d, %d) with phi = %f\n",
				pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
	}



	printf("(ScalarField3D::SeedParticles) particles are seeded \n");

}*/

void ScalarField3D::NumOfNonEscapedParticles(const Voxel &voxel, int *cellParticles, int s,
			list<Particle> &particles) const{
	float delta = voxel.VoxelDelta();
	list<Particle>::iterator iter_particle;
	memset(cellParticles, 0, DimX*DimY*DimZ * sizeof(int));
	if(s < 0){
		for(iter_particle  = particles.begin();
			iter_particle != particles.end();
			){
			Particle &ps = *iter_particle;
			Point pos = ps.Position();
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			float phi_i = TriInterp(voxel, pos);
			if(phi_i <= 0.f){
				if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
					++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				}
			}
			++iter_particle;
		}
	}
	else{
		for(iter_particle  = particles.begin();
			iter_particle != particles.end();
			){
			Particle &ps = *iter_particle;
			Point pos = ps.Position();
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			float phi_i = TriInterp(voxel, pos);
			if(phi_i > 0.f){
				if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
					++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				}
			}
			++iter_particle;
		}
	}
}

/*void ScalarField3D::ReSeedParticles(const Voxel &voxel,
				const ScalarField3D &object,
				list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
				list<Particle> &posParticles,
#endif
				int nPerCell) const{


	float delta = voxel.VoxelDelta();

	typedef priority_queue<float> PBUFFER;

	map<u_int, PBUFFER> forDelete;
	map<u_int, PBUFFER>::iterator iter_cell;

	int *cellNegParticles = new int[DimX*DimY*DimZ];
	memset(cellNegParticles, 0, DimX*DimY*DimZ * sizeof(int));
#ifdef POSITIVE_PARTICLES
	int *cellPosParticles = new int[DimX*DimY*DimZ];
	memset(cellPosParticles, 0, DimX*DimY*DimZ * sizeof(int));
#endif
	list<Particle> tempNegParticles;
#ifdef POSITIVE_PARTICLES
	list<Particle> tempPosParticles;
#endif
	list<Particle>::iterator iter_particle;

	float *phiobj = object.getScalar();

	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellNegParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else{
				iter_particle = negParticles.erase(iter_particle);
			}
		}
		else
			++iter_particle;
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellPosParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = posParticles.erase(iter_particle);
		}
		else
			++iter_particle;

	}
#endif

	printf(" (1) at cell [%d, %d, %d] there are %d negative non-escaped particles \n",
			  I, J, K, cellNegParticles[INDEX(I,J,K)]);
#ifdef POSITIVE_PARTICLES
	printf(" (1) at cell [%d, %d, %d] there are %d positive non-escaped particles \n",
			I, J, K, cellPosParticles[INDEX(I,J,K)]);
#endif

	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
		int n;
		if( (fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT && cellNegParticles[INDEX(i,j,k)] < nPerCell) &&
				(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp,n)) ){
			int N = nPerCell - cellNegParticles[INDEX(i,j,k)];
			for(int m=0; m<N; ++m){
				MTRand mt;
				double rn = mt();
				Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
				voxel.PlaceParticles(i, j, k, ps, 0);
				tempNegParticles.push_back(ps);
			}
		}
#ifdef POSITIVE_PARTICLES
		if( (fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT && cellPosParticles[INDEX(i,j,k)] < nPerCell) &&
			(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp,n)) ){
				int N = nPerCell - cellPosParticles[INDEX(i,j,k)];
				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, 1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					tempPosParticles.push_back(ps);
				}
			}
#endif
	END_FOR
//	NumOfNonEscapedParticles(voxel,cellNegParticles,-1,tempNegParticles);
//	printf(" (3) at cell [%d, %d, %d] there are %d negative non-escaped particles \n",
//				  I, J, K, cellNegParticles[INDEX(I,J,K)]);
//#ifdef POSITIVE_PARTICLES
//	NumOfNonEscapedParticles(voxel,cellPosParticles,1,tempPosParticles);
//	printf(" (3) at cell [%d, %d, %d] there are %d positive non-escaped particles \n",
//			I, J, K, cellPosParticles[INDEX(I,J,K)]);
//#endif

	AttractParticles(voxel, object, tempNegParticles,
#ifdef POSITIVE_PARTICLES
				tempPosParticles,
#endif
				nPerCell);
//	NumOfNonEscapedParticles(voxel,cellNegParticles,-1,tempNegParticles);
//	printf(" (4) at cell [%d, %d, %d] there are %d negative non-escaped particles \n",
//				  I, J, K, cellNegParticles[INDEX(I,J,K)]);
//#ifdef POSITIVE_PARTICLES
//	NumOfNonEscapedParticles(voxel,cellPosParticles,1,tempPosParticles);
//	printf(" (4) at cell [%d, %d, %d] there are %d positive non-escaped particles \n",
//			I, J, K, cellPosParticles[INDEX(I,J,K)]);
//#endif

	// find cells with more non-escaped particles than necessary
	FOR_EACH_CELL
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellNegParticles[INDEX(i,j,k)] > nPerCell){
			PBUFFER thisheap;
			forDelete.insert(make_pair(INDEX(i,j,k), thisheap));
		}
	END_FOR


	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float radius = ps.Radius();
				float value = -phi_i - radius;
				thisheap.push(value);
			}
		}
	}

// compute threshold for these cells
	for(iter_cell = forDelete.begin();
		iter_cell != forDelete.end();
		++iter_cell){
		PBUFFER &thisheap = iter_cell->second;
		u_int index = iter_cell->first;
		int N = nPerCell;
		do{
			thisheap.pop();
		}while(thisheap.size() > N);
		// after this do-while loop, thisheap.top() will
		// return the nPerCell-th big distance of the particles in
		// this cell. any particle with distance larger will be
		// deleted
	}

	// get rid of particles beyond threshold
	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float radius = ps.Radius();
				float value = -phi_i - radius;
				float th = thisheap.top();
				if(value > th)
					iter_particle = negParticles.erase(iter_particle);
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			++iter_particle;
	}

	NumOfNonEscapedParticles(voxel,cellNegParticles,-1,negParticles);
	printf(" (2) at cell [%d, %d, %d] there are %d negative non-escaped particles \n",
				  I, J, K, cellNegParticles[INDEX(I,J,K)]);

#ifdef POSITIVE_PARTICLES
	forDelete.clear();
// find cells with more non-escaped particles than necessary
	FOR_EACH_CELL
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellPosParticles[INDEX(i,j,k)] > nPerCell){
			PBUFFER thisheap;
			forDelete.insert(make_pair(INDEX(i,j,k), thisheap));
		}
	END_FOR

	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float radius = ps.Radius();
				float value = phi_i - radius;
				thisheap.push(value);
			}
		}
	}

	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float radius = ps.Radius();
				float value = phi_i - radius;
				float th = thisheap.top();
				if(value > th)
					iter_particle = posParticles.erase(iter_particle);
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			++iter_particle;
	}
	NumOfNonEscapedParticles(voxel,cellPosParticles,1,posParticles);
	printf(" (2) at cell [%d, %d, %d] there are %d positive non-escaped particles \n",
				I, J, K, cellPosParticles[INDEX(I,J,K)]);
#endif

#ifdef POSITIVE_PARTICLES
	printf("Checking positive particles... \n");

	for(iter_particle = tempPosParticles.begin();
		iter_particle != tempPosParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = object.TriInterp(voxel, pos);
		if(phi_i > 0.f){
//			printf("positive particle pos at (%f, %f, %f) with phi_obj = %f\n",
//				pos.x, pos.y, pos.z, phi_i);
			iter_particle = tempPosParticles.erase(iter_particle);
		}
		else{
			phi_i = TriInterp(voxel, pos);
			TentativeGridPoint tpp = voxel.ContainsPoint(pos);
			if(phi_i < RADIUS_MIN || phi_i > PARTICLE_SEED_LIMIT)
				printf("positive particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
					pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
			++iter_particle;
		}
	}
#endif



	printf("Checking negative particles... \n");
	for(iter_particle = tempNegParticles.begin();
		iter_particle != tempNegParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = object.TriInterp(voxel, pos);
		if(phi_i > 0.f){
//			printf("negative particle pos at (%f, %f, %f) with phi_obj = %f\n",
//				pos.x, pos.y, pos.z, phi_i);
			iter_particle = tempNegParticles.erase(iter_particle);
		}
		else{
			phi_i = TriInterp(voxel, pos);
			TentativeGridPoint tpp = voxel.ContainsPoint(pos);
			if(phi_i > -RADIUS_MIN || phi_i < -PARTICLE_SEED_LIMIT)
				printf("negative particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
					pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
			++iter_particle;
		}
	}

	printf("There are %ld newly added negative particles \n ", tempNegParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("There are %ld newly added positive particles \n ", tempPosParticles.size());
#endif

	for(iter_particle  = tempNegParticles.begin();
		iter_particle != tempNegParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		negParticles.push_back(ps);
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = tempPosParticles.begin();
		iter_particle != tempPosParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		posParticles.push_back(ps);
	}
#endif



//	NumOfNonEscapedParticles(voxel,cellNegParticles,-1,negParticles);
//#ifdef POSITIVE_PARTICLES
//	NumOfNonEscapedParticles(voxel,cellPosParticles,1,posParticles);
//#endif
//	FOR_EACH_CELL
//		TentativeGridPoint tp(0.f, i, j, k);
//		int n;
//		if( (phi[INDEX(i,j,k)] <= 0.f && fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
//			cellNegParticles[INDEX(i,j,k)] < nPerCell) &&
//			(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp,n)) ){
//			printf(" at cell [%d, %d, %d] with phi = %f and phiobj = %f, there are %d, not enough negative particles \n",
//					i, j, k, phi[INDEX(i,j,k)], phiobj[INDEX(i,j,k)], cellNegParticles[INDEX(i,j,k)]);
//		}
//
//#ifdef POSITIVE_PARTICLES
//		if( (phi[INDEX(i,j,k)] > 0.f && fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
//			cellPosParticles[INDEX(i,j,k)] < nPerCell) &&
//			(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp,n)) ){
//			printf(" at cell [%d, %d, %d] with phi = %f and phiobj = %f, there are %d, not enough positive particles \n",
//						i, j, k, phi[INDEX(i,j,k)], phiobj[INDEX(i,j,k)], cellPosParticles[INDEX(i,j,k)]);
//			}
//#endif
//	END_FOR

	delete [] cellNegParticles;
	delete [] cellPosParticles;

	printf("after reseeding now there are %ld negative particles \n", negParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("after reseeding now there are %ld positive particles \n", posParticles.size());
#endif
	printf("(ScalarField3D::ReSeedParticles) particles are reseeded \n");
}*/

void ScalarField3D::ReSeedParticlesNew(const Voxel &voxel, const ScalarField3D &object,
				list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
				list<Particle> &posParticles,
#endif
				int nPerCell) const{

	float delta = voxel.VoxelDelta();

	typedef priority_queue<float> PBUFFER;

	map<u_int, PBUFFER> forDelete;
	map<u_int, PBUFFER>::iterator iter_cell;

	int *cellParticles = new int[DimX * DimY *DimZ];
	memset(cellParticles, 0, DimX * DimY *DimZ * sizeof(int));
	list<Particle>::iterator iter_particle;

	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = negParticles.erase(iter_particle);
		}
		else
			++iter_particle;
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = posParticles.erase(iter_particle);
		}
		else
			++iter_particle;

	}
#endif


	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && ((fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellParticles[INDEX(i,j,k)] < nPerCell))){
//				(fabs(phi[INDEX(i,j,k)]) < 0.5 * delta &&
//				cellParticles[INDEX(i,j,k)] < 2 * nPerCell)) ){
			int N;
//			if(fabs(phi[INDEX(i,j,k)]) <= 0.5f * delta)
//				N = 2 * nPerCell - cellParticles[INDEX(i,j,k)];
//			else
				N = nPerCell - cellParticles[INDEX(i,j,k)];

			for(int m=0; m<N; ++m){
				MTRand mt;
				double rn = mt();
				Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
				voxel.PlaceParticles(i, j, k, ps, 0);
				Point pos = ps.Position();
				float phi_i = TriInterp(voxel, pos);
				if(phi_i <= 0.f){
					float radius;
					if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = -phi_i;
					ps.SetParticle(pos, radius);
					negParticles.push_back(ps);
				}
				else{
					float radius;
					if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = phi_i;
					ps.SetParticle(pos, radius);
#ifdef POSITIVE_PARTICLES
					posParticles.push_back(ps);
#endif
				}
			}
		}
	END_FOR


	// find cells with more non-escaped particles than necessary
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
//			if(fabs(phi[INDEX(i,j,k)]) < 0.5 * delta &&
//			   cellParticles[INDEX(i,j,k)] > 2 * nPerCell){
//				PBUFFER thisheap;
//				forDelete.insert(make_pair(INDEX(i,j,k), thisheap));
//			}
//			else
			if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
					cellParticles[INDEX(i,j,k)] > nPerCell){
				PBUFFER thisheap;
				forDelete.insert(make_pair(INDEX(i,j,k), thisheap));
			}
		}
	END_FOR


	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float radius = ps.Radius();
				float value = -phi_i - radius;
				thisheap.push(value);
			}
		}
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float radius = ps.Radius();
				float value = phi_i - radius;
				thisheap.push(value);
			}
		}
	}
#endif

	// compute threshold for these cells
	for(iter_cell = forDelete.begin();
	    iter_cell != forDelete.end();
	    ++iter_cell){
		PBUFFER &thisheap = iter_cell->second;
		u_int index = iter_cell->first;
		int N;
//		if(fabs(phi[index]) <= 0.5 * delta)
//			N = 2 * nPerCell;
//		else
			N = nPerCell;
		do{
			thisheap.pop();
		}while(thisheap.size() > N);
		// after this do-while loop, thisheap.top() will
		// return the nPerCell-th big distance of the particles in
		// this cell. any particle with distance larger will be
		// deleted
	}

	// get rid of particles beyond threshold
	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		float phi_obj = object.TriInterp(voxel, pos);
		if(phi_obj <= 0.f){ // seeded particles are outside of solid object
			if(phi_i <= 0.f){
				TentativeGridPoint tp = voxel.ContainsPoint(pos);
				iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
				if(iter_cell != forDelete.end()){
					PBUFFER &thisheap = iter_cell->second;
					float radius = ps.Radius();
					float value = -phi_i - radius;
					float th = thisheap.top();
					if(value > th)
						iter_particle = negParticles.erase(iter_particle);
					else
						++iter_particle;
				}
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			iter_particle = negParticles.erase(iter_particle);
	}
#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = TriInterp(voxel, pos);
		float phi_obj = object.TriInterp(voxel, pos);
		if(phi_obj <= 0.f){ // seeded particles are outside of solid object
			if(phi_i > 0.f){
				TentativeGridPoint tp = voxel.ContainsPoint(pos);
				iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
				if(iter_cell != forDelete.end()){
					PBUFFER &thisheap = iter_cell->second;
					float radius = ps.Radius();
					float value = phi_i - radius;
					float th = thisheap.top();
					if(value > th)
						iter_particle = posParticles.erase(iter_particle);
					else
						++iter_particle;
				}
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			iter_particle = posParticles.erase(iter_particle);
	}
#endif


	// we also seed particles in those solid cells which are close to non-solid cells
	FOR_EACH_CELL
	    TentativeGridPoint tp(0.f, i, j, k);
		int n;
		if(voxel.InSolid(i,j,k) && fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT
				&& voxel.CloseToNonSolid(tp, n)){

			int N;
//			if(fabs(phi[INDEX(i,j,k)]) <= 0.5f * delta)
//				N = 2 * nPerCell;
//			else
				N = nPerCell;

			for(int m=0; m<N; ++m){
				MTRand mt;
				double rn = mt();
				Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
				voxel.PlaceParticles(i, j, k, ps, 0);
				Point pos = ps.Position();
				float phi_obj = object.TriInterp(voxel, pos);
				if(phi_obj <= 0.f){
					float phi_i = TriInterp(voxel, pos);
					if(phi_i <= 0.f){
						float radius;
						if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = -phi_i;
						ps.SetParticle(pos, radius);
						negParticles.push_back(ps);
					}
					else{
						float radius;
						if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = phi_i;
						ps.SetParticle(pos, radius);
#ifdef POSITIVE_PARTICLES
						posParticles.push_back(ps);
#endif
					}
				}
			}
		}
	END_FOR

	delete [] cellParticles;
	printf("(ScalarField3D::ReSeedParticlesNew) particles are reseeded \n");
}

void ScalarField3D::ReSeedSourceParticles(const Voxel &voxel, const ScalarField3D &object,
				list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
				list<Particle> &posParticles,
#endif
				int nPerCell) const{

	float delta = voxel.VoxelDelta();

	float *phi_obj = object.getScalar();

	int *cellParticles = new int[DimX * DimY * DimZ];
	memset(cellParticles, 0, DimX * DimY * DimZ * sizeof(int));

	list<Particle> tempNegParticles;
	list<Particle> tempPosParticles;
	list<Particle>::iterator iter_particle;

	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= radius){ // non-escaped
			if( voxel.InSource(tp.ii, tp.jj, tp.kk) && fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT ){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
			}
			++iter_particle;
		}
		else
			++iter_particle;
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i >= -radius){
			if(voxel.InSource(tp.ii, tp.jj, tp.kk) && fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
			}
			++iter_particle;

		}
		else
			++iter_particle;

	}
#endif

//	if(i == I && j == J && k == K)
	printf(" A. at (%d, %d, %d) there are %d nonescaped particles \n",
			I, J, K, cellParticles[INDEX(I, J, K)]);

	/*FOR_EACH_CELL
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellParticles[INDEX(i,j,k)] < nPerCell){
			TentativeGridPoint tp(0.f, i, j, k);
			int n;
			if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){

				int N = nPerCell - cellParticles[INDEX(i,j,k)];

				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
					float phi_i = TriInterp(voxel, pos);
					if(phi_i <= 0.f){
						float radius;
						if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = -phi_i;
						ps.SetParticle(pos, radius);
						tempNegParticles.push_back(ps);
					}
#ifdef POSITIVE_PARTICLES
					else{
						float radius;
						if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = phi_i;
						ps.SetParticle(pos, radius);
						tempPosParticles.push_back(ps);
					}
#endif
				}
			}
		}
	END_FOR*/

//	u_int N1 = 0, N2 = 0;
	FOR_EACH_CELL
	  if(voxel.InSource(i, j, k) && fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellParticles[INDEX(i,j,k)] < nPerCell){


				int NP = cellParticles[INDEX(i,j,k)];

				int N = nPerCell - NP;
				int iters = 0;
				if(i == I && j == J && k == K)
					printf(" at (%d, %d, %d) there are %d nonescaped particles %d will be added \n",
							i, j, k, NP, N);
				do{
				int MP = 0;
				++iters;
				int SP = 0;
				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
					float phi_i = TriInterp(voxel, pos);
					if(phi_i <= 0.f){
						float radius;
						if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
//						else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else if( -phi_i < RADIUS_MIN ){
//								++N1;
							++MP;
							continue;
						}
						else radius = -phi_i;
						ps.SetParticle(pos, radius);
						tempNegParticles.push_back(ps);
						++SP;
					}
#ifdef POSITIVE_PARTICLES
					else{
						float radius;
						if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
//						else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else if( phi_i < RADIUS_MIN ){
//								++N2;
							++MP;
							continue;
						}
						else radius = phi_i;
						ps.SetParticle(pos, radius);
						tempPosParticles.push_back(ps);
						++SP;
					}
#endif
				}
				if(i == I && j == J && k == K)
					printf(" at (%d, %d, %d) iter = %d added %d particles, there are %d particles too close to interface \n",
							i, j, k, iters, SP, MP);
				N = MP;
			}while (N > 0 && iters < 15);

		}
		END_FOR

#ifdef POSITIVE_PARTICLES
	printf("Checking positive particles... \n");

	for(iter_particle = tempPosParticles.begin();
		iter_particle != tempPosParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){
			iter_particle = tempPosParticles.erase(iter_particle);
		}
		else{
			float phi_i = object.TriInterp(voxel, pos);
			if(phi_i > 0.f){
	//			printf("positive particle pos at (%f, %f, %f) with phi_obj = %f\n",
	//				pos.x, pos.y, pos.z, phi_i);
				iter_particle = tempPosParticles.erase(iter_particle);
			}
			else{
//				phi_i = TriInterp(voxel, pos);
//				if(phi_i < RADIUS_MIN || phi_i > PARTICLE_SEED_LIMIT)
//					printf("positive particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
//						pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
				++iter_particle;
			}
		}
	}
#endif



	printf("Checking negative particles... \n");
	for(iter_particle = tempNegParticles.begin();
		iter_particle != tempNegParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){
			iter_particle = tempNegParticles.erase(iter_particle);
		}
		else{
			float phi_i = object.TriInterp(voxel, pos);
			if(phi_i > 0.f){
	//			printf("negative particle pos at (%f, %f, %f) with phi_obj = %f\n",
	//				pos.x, pos.y, pos.z, phi_i);
				iter_particle = tempNegParticles.erase(iter_particle);
			}
			else{
//				phi_i = TriInterp(voxel, pos);
//				TentativeGridPoint tpp = voxel.ContainsPoint(pos);
//				if(phi_i > -RADIUS_MIN || phi_i < -PARTICLE_SEED_LIMIT)
//					printf("negative particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
//						pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
				++iter_particle;
			}
		}
	}

	printf("There are %u newly added negative particles \n ", tempNegParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("There are %u newly added positive particles \n ", tempPosParticles.size());
#endif

	for(iter_particle  = tempNegParticles.begin();
		iter_particle != tempNegParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		if(tp.ii == I && tp.jj ==J && tp.kk == K)
			printf("CM. at (%d, %d, %d) negative particle (%f, %f, %f) added \n",
					tp.ii, tp.jj, tp.kk, pos.x, pos.y, pos.z);
		negParticles.push_back(ps);
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = tempPosParticles.begin();
		iter_particle != tempPosParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		if(tp.ii == I && tp.jj ==J && tp.kk == K)
			printf("CP. at (%d, %d, %d) positive particle (%f, %f, %f) added \n",
					tp.ii, tp.jj, tp.kk, pos.x, pos.y, pos.z);
		posParticles.push_back(ps);
	}
#endif

		delete [] cellParticles;
}

void ScalarField3D::ReSeedParticles(const Voxel &voxel, const ScalarField3D &object,
				ParticleList &negParticles,
#ifdef POSITIVE_PARTICLES
				ParticleList &posParticles,
#endif
				int nPerCell) const{

	float delta = voxel.VoxelDelta();

	float *phi_obj = object.getScalar();

	typedef priority_queue<float> PBUFFER;
#ifdef WIN32
	hash_map<u_int, PBUFFER> forDelete;
	hash_map<u_int, PBUFFER>::iterator iter_cell;
#else
	unordered_map<u_int, PBUFFER> forDelete;
	unordered_map<u_int, PBUFFER>::iterator iter_cell;
#endif

	int *cellParticles = new int[DimX * DimY * DimZ];
	memset(cellParticles, 0, DimX * DimY * DimZ * sizeof(int));

//	ParticleVector tempNegParticles;
//	ParticleVector tempPosParticles;
//	ParticleVector::iterator iter_tempParticle;
//	ParticleList tempNegParticles;
//	ParticleList tempPosParticles;

	ParticleList::iterator iter_particle;

	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= radius){ // non-escaped
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = negParticles.erase(iter_particle);
		}
		else
			++iter_particle;
	}


#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i >= -radius){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = posParticles.erase(iter_particle);
		}
		else
			++iter_particle;

	}
#endif

	/*FOR_EACH_CELL
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellParticles[INDEX(i,j,k)] < nPerCell){
			TentativeGridPoint tp(0.f, i, j, k);
			int n;
			if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){

				int N = nPerCell - cellParticles[INDEX(i,j,k)];

				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
					float phi_i = TriInterp(voxel, pos);
					if(phi_i <= 0.f){
						float radius;
						if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = -phi_i;
						ps.SetParticle(pos, radius);
						tempNegParticles.push_back(ps);
					}
#ifdef POSITIVE_PARTICLES
					else{
						float radius;
						if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = phi_i;
						ps.SetParticle(pos, radius);
						tempPosParticles.push_back(ps);
					}
#endif
				}
			}
		}
	END_FOR*/

//	u_int N1 = 0, N2 = 0;


//	printf(" there are %u negative particles not added due to their approximity to interface\n ", N1);
//#ifdef POSITIVE_PARTICLES
//	printf(" there are %u positive particles not added due to their approximity to interface\n ", N2);
//#endif

	/*AttractParticles(voxel, object, tempNegParticles,
#ifdef POSITIVE_PARTICLES
					tempPosParticles,
#endif
					nPerCell);*/


	// find cells with more non-escaped particles than necessary
	FOR_EACH_CELL
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellParticles[INDEX(i,j,k)] > nPerCell){
			PBUFFER thisheap;
			forDelete.insert(make_pair(INDEX(i,j,k), thisheap));
		}
	END_FOR


	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = -phi_i - radius;
				thisheap.push(value);
			}
		}
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i >= -radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = phi_i - radius;
				thisheap.push(value);
			}
		}
	}
#endif

	// compute threshold for these cells
	for(iter_cell = forDelete.begin();
	    iter_cell != forDelete.end();
	    ++iter_cell){
		PBUFFER &thisheap = iter_cell->second;
//		u_int index = iter_cell->first;
		int N = nPerCell;
		do{
			thisheap.pop();
		}while(thisheap.size() > N);
		// after this do-while loop, thisheap.top() will
		// return the nPerCell-th big distance of the particles in
		// this cell. any particle with distance larger will be
		// deleted
	}

	// get rid of particles beyond threshold
	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = -phi_i - radius;
				float th = thisheap.top();
				if(value > th)
					iter_particle = negParticles.erase(iter_particle);
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			++iter_particle;
	}
#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i >= -radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = phi_i - radius;
				float th = thisheap.top();
				if(value > th)
					iter_particle = posParticles.erase(iter_particle);
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			++iter_particle;
	}
#endif

	printf("After removing unused particles, there are %u  negative particles \n ", negParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("After removing unused particles, there are %u  positive particles \n ", posParticles.size());
#endif

	u_int NL = 0;
	u_int NTP = 0;
	u_int NTM = 0;
	FOR_EACH_CELL
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellParticles[INDEX(i,j,k)] < nPerCell){
			TentativeGridPoint tp(0.f, i, j, k);
			int n;
			if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){
				++NL;
				int NP = cellParticles[INDEX(i,j,k)];
				int N = nPerCell - NP;
				int iters = 0;
				do{
				int MP = 0;
				++iters;
				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
					float phi_i = TriInterp(voxel, pos);
					float phi_i_obj = object.TriInterp(voxel, pos);
					if(phi_i <= 0.f){
						float radius;
						if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
//						else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else if( -phi_i < RADIUS_MIN ){
//							++N1;
							++MP;
							continue;
						}
						else radius = -phi_i;
						ps.SetParticle(pos, radius);
						if( phi_i_obj <= 0.f && phi_i >= (-PARTICLE_SEED_LIMIT) ){
							negParticles.push_back(ps);
							++NTM;
						}
					}
#ifdef POSITIVE_PARTICLES
					else{
						float radius;
						if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
//						else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else if( phi_i < RADIUS_MIN ){
//							++N2;
							++MP;
							continue;
						}
						else radius = phi_i;
						ps.SetParticle(pos, radius);
						if( phi_i_obj <= 0.f && phi_i <= PARTICLE_SEED_LIMIT ){
							posParticles.push_back(ps);
							++NTP;
						}
					}
#endif
				}
				N = MP;
			}while (N > 0 && iters < 15);
			}
		}
	END_FOR

	printf("\nthere are %u grid cells lacking particles \n\n ", NL);

//#ifdef POSITIVE_PARTICLES
//	printf("Checking positive particles... \n");
//
//	for(iter_particle = tempPosParticles.begin();
//		iter_particle != tempPosParticles.end();
//		){
//		Particle &ps = *iter_particle;
//		Point pos = ps.Position();
//		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
//		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){
//			iter_particle = tempPosParticles.erase(iter_particle);
//		}
//		else{
//			float phi_i = object.TriInterp(voxel, pos);
//			if(phi_i > 0.f){
//	//			printf("positive particle pos at (%f, %f, %f) with phi_obj = %f\n",
//	//				pos.x, pos.y, pos.z, phi_i);
//				iter_particle = tempPosParticles.erase(iter_particle);
//			}
//			else{
////				phi_i = TriInterp(voxel, pos);
////				if(phi_i < RADIUS_MIN || phi_i > PARTICLE_SEED_LIMIT)
////					printf("positive particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
////						pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
//				++iter_particle;
//			}
//		}
//	}
//#endif



//	printf("Checking negative particles... \n");
//	for(iter_particle = tempNegParticles.begin();
//		iter_particle != tempNegParticles.end();
//		){
//		Particle &ps = *iter_particle;
//		Point pos = ps.Position();
//		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
//		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){
//			iter_particle = tempNegParticles.erase(iter_particle);
//		}
//		else{
//			float phi_i = object.TriInterp(voxel, pos);
//			if(phi_i > 0.f){
//	//			printf("negative particle pos at (%f, %f, %f) with phi_obj = %f\n",
//	//				pos.x, pos.y, pos.z, phi_i);
//				iter_particle = tempNegParticles.erase(iter_particle);
//			}
//			else{
////				phi_i = TriInterp(voxel, pos);
////				TentativeGridPoint tpp = voxel.ContainsPoint(pos);
////				if(phi_i > -RADIUS_MIN || phi_i < -PARTICLE_SEED_LIMIT)
////					printf("negative particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
////						pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
//				++iter_particle;
//			}
//		}
//	}

	printf("There are %u newly added negative particles \n ", NTM);
#ifdef POSITIVE_PARTICLES
	printf("There are %u newly added positive particles \n ", NTP);
#endif
//
//	for(iter_tempParticle  = tempNegParticles.begin();
//		iter_TempParticle != tempNegParticles.end();
//		++iter_particle){
//		Particle &ps = *iter_tempParticle;
//		negParticles.push_back(ps);
//	}
//	tempNegParticles.clear();
//
//#ifdef POSITIVE_PARTICLES
//	for(iter_tempParticle  = tempPosParticles.begin();
//		iter_tempParticle != tempPosParticles.end();
//		++iter_particle){
//		Particle &ps = *iter_tempParticle;
//		posParticles.push_back(ps);
//	}
//	tempPosParticles.clear();
//#endif


	// check the density of particles near interface
	memset(cellParticles, 0, DimX * DimY * DimZ * sizeof(int));

	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= radius){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = negParticles.erase(iter_particle);
		}
		else
			++iter_particle;
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i >= -radius){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = posParticles.erase(iter_particle);
		}
		else
			++iter_particle;

	}
#endif
	FOR_EACH_CELL
		if( phi_obj[INDEX(i,j,k)] < -1.5f * delta && !voxel.InSource(i,j,k) &&
				fabs(phi[INDEX(i,j,k)]) <= 1.5f * delta &&
				cellParticles[INDEX(i,j,k)] > 0 ){
			if(cellParticles[INDEX(i,j,k)] != nPerCell)
//		if(i == I && j == J && k == K)
			printf(" at (%d, %d, %d) phi = %f there are %d particles \n",
						i,j,k,phi[INDEX(i,j,k)],cellParticles[INDEX(i,j,k)]);
		}
	END_FOR


	delete [] cellParticles;
	printf("(ScalarField3D::ReSeedParticles) particles are reseeded \n");
}

void ScalarField3D::ReSeedParticles2(const Voxel &voxel, const ScalarField3D &object,
				list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
				list<Particle> &posParticles,
#endif
				int nPerCell) const{

	float delta = voxel.VoxelDelta();

	float *phi_obj = object.getScalar();

	typedef priority_queue<float> PBUFFER;

	map<u_int, PBUFFER> forDelete;
	map<u_int, PBUFFER>::iterator iter_cell;

	int *cellParticles = new int[DimX * DimY * DimZ];
	memset(cellParticles, 0, DimX * DimY * DimZ * sizeof(int));
	int *cellPosParticles = new int[DimX * DimY * DimZ];
	memset(cellPosParticles, 0, DimX * DimY * DimZ * sizeof(int));
	int *cellNegParticles = new int[DimX * DimY * DimZ];
	memset(cellNegParticles, 0, DimX * DimY * DimZ * sizeof(int));

	list<Particle> tempNegParticles;
	list<Particle> tempPosParticles;
	list<Particle>::iterator iter_particle;

	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f || phi_i <= radius){ // non-escaped
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellNegParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = negParticles.erase(iter_particle);
		}
		else
			++iter_particle;
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f || phi_i >= -radius){  //non-escaped
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellPosParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = posParticles.erase(iter_particle);
		}
		else
			++iter_particle;

	}
#endif

	FOR_EACH_CELL
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellNegParticles[INDEX(i,j,k)] < nPerCell){
			TentativeGridPoint tp(0.f, i, j, k);
			int n;
			if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){

				int N = nPerCell - cellNegParticles[INDEX(i,j,k)];

				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					tempNegParticles.push_back(ps);
				}
			}
		}
#ifdef POSITIVE_PARTICLES
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
				cellPosParticles[INDEX(i,j,k)] < nPerCell){
			TentativeGridPoint tp(0.f, i, j, k);
			int n;
			if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){

				int N = nPerCell - cellPosParticles[INDEX(i,j,k)];

				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					tempPosParticles.push_back(ps);
				}
			}
		}
#endif
	END_FOR

	AttractParticles(voxel, object, tempNegParticles,
#ifdef POSITIVE_PARTICLES
					tempPosParticles,
#endif
					nPerCell);

#ifdef POSITIVE_PARTICLES
	printf("Checking positive particles... \n");

	for(iter_particle = tempPosParticles.begin();
		iter_particle != tempPosParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		if(radius < 0.f){
			printf("Error! radius is negative \n");
			exit(1);
		}
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){
			iter_particle = tempPosParticles.erase(iter_particle);
		}
		else{
			float phi_i = object.TriInterp(voxel, pos);
			if(phi_i > 0.f){
	//			printf("positive particle pos at (%f, %f, %f) with phi_obj = %f\n",
	//				pos.x, pos.y, pos.z, phi_i);
				iter_particle = tempPosParticles.erase(iter_particle);
			}
			else{
				phi_i = TriInterp(voxel, pos);
				if(phi_i < RADIUS_MIN || phi_i > PARTICLE_SEED_LIMIT)
					printf("positive particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
						pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
				++iter_particle;
			}
		}
	}
#endif



	printf("Checking negative particles... \n");
	for(iter_particle = tempNegParticles.begin();
		iter_particle != tempNegParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		if(radius < 0.f){
			printf("Error! radius is negative \n");
			exit(1);
		}
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){
			iter_particle = tempNegParticles.erase(iter_particle);
		}
		else{
			float phi_i = object.TriInterp(voxel, pos);
			if(phi_i > 0.f){
	//			printf("negative particle pos at (%f, %f, %f) with phi_obj = %f\n",
	//				pos.x, pos.y, pos.z, phi_i);
				iter_particle = tempNegParticles.erase(iter_particle);
			}
			else{
				phi_i = TriInterp(voxel, pos);
				TentativeGridPoint tpp = voxel.ContainsPoint(pos);
				if(phi_i > -RADIUS_MIN || phi_i < -PARTICLE_SEED_LIMIT)
					printf("negative particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
						pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
				++iter_particle;
			}
		}
	}

	printf("There are %u newly added negative particles \n ", tempNegParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("There are %u newly added positive particles \n ", tempPosParticles.size());
#endif


	for(iter_particle  = tempNegParticles.begin();
		iter_particle != tempNegParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		negParticles.push_back(ps);
	}

#ifdef POSITIVE_PARTICLES
		for(iter_particle  = tempPosParticles.begin();
			iter_particle != tempPosParticles.end();
			++iter_particle){
			Particle &ps = *iter_particle;
			posParticles.push_back(ps);
		}
#endif

	tempNegParticles.clear();
	tempPosParticles.clear();

	for(iter_particle  = negParticles.begin();
	    iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f || phi_i <= radius){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = negParticles.erase(iter_particle);
		}
		else
			++iter_particle;
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f || phi_i >= -radius){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = posParticles.erase(iter_particle);
		}
		else
			++iter_particle;

	}
#endif

	// find cells with more non-escaped particles than necessary
	FOR_EACH_CELL
			if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
					cellParticles[INDEX(i,j,k)] > nPerCell){
				PBUFFER thisheap;
				forDelete.insert(make_pair(INDEX(i,j,k), thisheap));
			}
	END_FOR


	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f || phi_i <= radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = -phi_i - radius;
				thisheap.push(value);
			}
		}
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f || phi_i >= -radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = phi_i - radius;
				thisheap.push(value);
			}
		}
	}
#endif

	// compute threshold for these cells
	for(iter_cell = forDelete.begin();
	    iter_cell != forDelete.end();
	    ++iter_cell){
		PBUFFER &thisheap = iter_cell->second;
		u_int index = iter_cell->first;
		int N = nPerCell;
		do{
			thisheap.pop();
		}while(thisheap.size() > N);
		// after this do-while loop, thisheap.top() will
		// return the nPerCell-th big distance of the particles in
		// this cell. any particle with distance larger will be
		// deleted
	}

	// get rid of particles beyond threshold
	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f || phi_i <= radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = -phi_i - radius;
				float th = thisheap.top();
				if(value > th)
					iter_particle = negParticles.erase(iter_particle);
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			++iter_particle;
	}
#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f || phi_i >= -radius){
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			iter_cell = forDelete.find(INDEX(tp.ii, tp.jj, tp.kk));
			if(iter_cell != forDelete.end()){
				PBUFFER &thisheap = iter_cell->second;
				float value = phi_i - radius;
				float th = thisheap.top();
				if(value > th)
					iter_particle = posParticles.erase(iter_particle);
				else
					++iter_particle;
			}
			else
				++iter_particle;
		}
		else
			++iter_particle;
	}
#endif

	printf("After removing unused particles, there are %u  negative particles \n ", negParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("After removing unused particles, there are %u  positive particles \n ", posParticles.size());
#endif



	// check the density of particles near interface
	memset(cellParticles, 0, DimX * DimY * DimZ * sizeof(int));

	for(iter_particle  = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f || phi_i <= radius){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = negParticles.erase(iter_particle);
		}
		else
			++iter_particle;
	}

#ifdef POSITIVE_PARTICLES
	for(iter_particle  = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(phi_i > 0.f || phi_i >= -radius){
			if(fabs(phi[INDEX(tp.ii, tp.jj, tp.kk)]) <= PARTICLE_SEED_LIMIT){
				++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
				++iter_particle;
			}
			else
				iter_particle = posParticles.erase(iter_particle);
		}
		else
			++iter_particle;

	}
#endif
	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && phi_obj[INDEX(i,j,k)] < -delta &&  !voxel.InSource(i,j,k) &&
//				fabs(phi[INDEX(i,j,k)]) <= 1.5f * delta &&
//				cellParticles[INDEX(i,j,k)] > 0){
//			if(cellParticles[INDEX(i,j,k)] < nPerCell)
//				printf(" at (%d, %d, %d) phi = %f there are not enought particles %d \n",
//						i,j,k,phi[INDEX(i,j,k)],cellParticles[INDEX(i,j,k)]);
//		}
		if(fabs(phi[INDEX(i,j,k)]) <= PARTICLE_SEED_LIMIT &&
			cellParticles[INDEX(i,j,k)] < nPerCell){
			TentativeGridPoint tp(0.f, i, j, k);
			int n;
			if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){
				int N = nPerCell - cellParticles[INDEX(i,j,k)];
				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
					float phi_i_obj = object.TriInterp(voxel, pos);
					if( phi_i_obj <= 0.f){
						float phi_i = TriInterp(voxel, pos);
						if(phi_i <= 0.f){
							float radius;
							if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
	//						else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
							else if( -phi_i < RADIUS_MIN ) continue;
							else radius = -phi_i;
							ps.SetParticle(pos, radius);
							negParticles.push_back(ps);
						}
#ifdef POSITIVE_PARTICLES
						else{
							float radius;
							if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
	//						else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
							else if( phi_i < RADIUS_MIN ) continue;
							else radius = phi_i;
							ps.SetParticle(pos, radius);
							posParticles.push_back(ps);
						}
#endif
					}
				}
			}
		}
	END_FOR


	delete [] cellParticles;
	delete [] cellNegParticles;
	delete [] cellPosParticles;
	forDelete.clear();
	printf("(ScalarField3D::ReSeedParticles2) particles are reseeded \n");
}

void ScalarField3D::AdjustWaterParticleRadii(const Voxel &voxel,
						list<WaterParticle> &currentWaterParticles,
						BlendParticles *water) const{

	float delta = voxel.VoxelDelta();
	float volume = delta * delta * delta;
	list<WaterParticle>::iterator iter_waterparticle;
	int *cellParticles = new int[DimX * DimY *DimZ];
	memset(cellParticles, 0, DimX * DimY *DimZ * sizeof(int));

	typedef list<WaterParticle *> PBUFFER;
	map<u_int, PBUFFER> vCells;
	map<u_int, PBUFFER>::iterator iter_cell;
	list<WaterParticle *>::iterator iter_escape;

	for(iter_waterparticle = currentWaterParticles.begin();
		iter_waterparticle != currentWaterParticles.end();
		){
		WaterParticle &wps = *iter_waterparticle;
		Point pos = wps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
		++iter_waterparticle;
	}

	FOR_EACH_CELL
		if(cellParticles[INDEX(i,j,k)] > 0){
			PBUFFER thisheap;
			vCells.insert(make_pair(INDEX(i,j,k), thisheap));
		}
	END_FOR

	for(iter_waterparticle = currentWaterParticles.begin();
		iter_waterparticle != currentWaterParticles.end();
		){
		WaterParticle &wps = *iter_waterparticle;
		Point pos = wps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		iter_cell = vCells.find(INDEX(tp.ii, tp.jj, tp.kk));
		if(iter_cell != vCells.end()){
			PBUFFER &thisheap = iter_cell->second;
			thisheap.push_back(&wps);
		}
		++iter_waterparticle;
	}

#define MONTE_CARLO_ITERATIONS 2048
	FOR_EACH_CELL
		iter_cell = vCells.find(INDEX(i, j, k));
		if(iter_cell != vCells.end()){
			PBUFFER &thisheap = iter_cell->second;
			// estimating volume loss within the given cell
			// using the Monte Carlo method
			int hit = 0;
			for(int n=0; n < MONTE_CARLO_ITERATIONS; ++n){
				Point p;
				voxel.PlacePoints(i,j,k,p,0);
				for(iter_escape = thisheap.begin();
					iter_escape != thisheap.end();
					++iter_escape){
					WaterParticle *wps = *iter_escape;
					if(wps->PointInside(p)){
						++hit;
						break;
					}
				}
			}
			float V = volume * hit / MONTE_CARLO_ITERATIONS;
			// distribute the lost volume onto particles
			float new_rad;
			if(hit > 0){
				float tmp = 3 * V / (4 * M_PI * cellParticles[INDEX(i,j,k)]);
				new_rad = pow(tmp, 1.f/3);
			}
			else
				new_rad = RADIUS_MIN;
			for(iter_escape = thisheap.begin();
				iter_escape != thisheap.end();
				++iter_escape){
				WaterParticle *wps = *iter_escape;
				Point pos = wps->Position();
				wps->SetParticle(pos, new_rad);
			}
		}
	END_FOR

	delete [] cellParticles;

	for(iter_waterparticle = currentWaterParticles.begin();
		iter_waterparticle != currentWaterParticles.end();
		){
		WaterParticle &wps = *iter_waterparticle;
		water->AddWaterParticles(wps);
		++iter_waterparticle;
	}
}


void ScalarField3D::AdjustParticleRadii(const Voxel &voxel, const VectorField3D &vel,
				int niter, int iters,
				ParticleList &negParticles
#ifdef POSITIVE_PARTICLES
				,ParticleList &posParticles
#endif
#ifdef COUPLED_SPH
				, SPH *water
#else
#ifdef COUPLED_FLIP
				, WaterParticleList &wpl
#else
				, BlendParticles *water
#endif
#endif
				, WaterParticleList &absorbed
				) const{

	float delta = voxel.VoxelDelta();
//	printf(" at here(1) \n");

	// timeToDelete controls when to delete escaped particles
	// true: every subcycling time step
	// niter == iters-1: every velocity time step
	// may have other options later
//	bool timeToDelete = (niter == iters - 1);
	bool timeToDelete = true;

	ParticleList::iterator iter_particle;

	u_int escapedParticles = 0, waterParticles = 0;
	if(absorbed.size() > 0)
		absorbed.clear();
//	printf(" at here(2) \n");
	list<WaterParticle> currentWaterParticles;
	WaterParticleList::iterator iter_waterparticle;

#ifdef COUPLED_FLIP
	for(iter_waterparticle  = wpl.begin();
		iter_waterparticle != wpl.end();
#else
	for(iter_waterparticle  = water->waterParticles.begin();
		iter_waterparticle != water->waterParticles.end();
#endif
		){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
//		float r = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		if(phi_i <= 0.f){ // get rid of particles back into the liquid
			absorbed.push_back(ps);
#ifdef COUPLED_FLIP
			iter_waterparticle = wpl.erase(iter_waterparticle);
#else
			iter_waterparticle = water->waterParticles.erase(iter_waterparticle);
#endif
		}
		else
			++iter_waterparticle;
	}
	printf("\n%u water particles are absorbed back into water in this time step\n\n", absorbed.size());

	for(iter_particle = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float rad = ps.Radius();
		float cd  = ps.CollisionDistance();
		float phi_i = TriInterp(voxel, pos);
//		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
//		if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk))
//			iter_particle = negParticles.erase(iter_particle);
//		else{
			if(phi_i <= 0.f){
//				if(fabs(phi_i) > 5*delta){
//					iter_particle = negParticles.erase(iter_particle);
//				}
//				else{
					if( -phi_i > RADIUS_MAX ) rad = RADIUS_MAX;
					else if( -phi_i < RADIUS_MIN ) rad = RADIUS_MIN;
					else rad = -phi_i;
					ps.SetParticle(pos, rad);
					++iter_particle;
//				}
			}
			else{
				// we delete particles which are escaped more than radius
				// in the original PLS paper, Enright suggests not deleting these
				// particles, but that approach is not suited here for fluid simulation
				//
				// later, Enright et al. (2004) mentioned deleting them if
				// they are escaped more than 1.5 times their radii.
				//
				// Bridson(2008) just propose to delete them all, since these particles
				// are no more useful in representing information on the grid. A more
				// serious problem is they may produce grid-aligned artifacts when they
				// come across close enough to a grid in one step, and then disappear
				// again in the next.
				if( phi_i > 1.f*rad ){
//				if( phi_i > PARTICLE_SEED_LIMIT ){
					if(timeToDelete){
						iter_particle = negParticles.erase(iter_particle);
						++waterParticles;
	//					radius = phi_i;
	//					radius = RADIUS_MIN;
						float u, v, w;
						vel.GetVelocityAtPos(voxel, pos, delta, u, v, w);
//						MTRand mt;
//						double rn = mt();
						WaterParticle wps(pos, rad, -1, cd, u, v, w);
	//					currentWaterParticles.push_back(wps);
#ifdef COUPLED_FLIP
						wpl.push_back(wps);
#else
						water->AddWaterParticles(wps);
#endif
					}
					else{
						rad = RADIUS_MIN;
						ps.SetParticle(pos, rad);
						++iter_particle;
					}
				}
				else{ // escaped particles have their radius set to minimum
					rad = RADIUS_MIN;
					ps.SetParticle(pos, rad);
					++iter_particle;
				}
				++escapedParticles;
			}
//		}

	}

//	AdjustWaterParticleRadii(voxel, currentWaterParticles, water);

#ifdef POSITIVE_PARTICLES
	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
//		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		if(phi_i > 0.f){
//			if(phi_i > 5*delta){
//				iter_particle = posParticles.erase(iter_particle);
//			}
//			else{
				if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
				else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
				else radius = phi_i;
				ps.SetParticle(pos, radius);
				++iter_particle;
//			}
		}
		else{
			// the same reason as for negative particles, see above
			if( phi_i < -1.f*radius ){
//			if( phi_i < -PARTICLE_SEED_LIMIT){
				if(timeToDelete)
					iter_particle = posParticles.erase(iter_particle);
				else{
					radius = RADIUS_MIN;
					ps.SetParticle(pos, radius);
					++iter_particle;
				}
			}
			else{// escaped particles have their radius set to minimum
				radius = RADIUS_MIN;
				ps.SetParticle(pos, radius);
				++iter_particle;
			}
			++escapedParticles;
		}
	}
#endif




	printf("(ScalarField3D::AdjustParticleRadii) there are %u escaped particles, %u water particles in this time step\n",
			  escapedParticles, waterParticles);
	printf("(ScalarField3D::AdjustParticleRadii) particle radius adjusted \n");

//	if(niter == iters - 1){
//		u_int N = 0;
//		FOR_EACH_CELL
//			if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f)
//				++N;
//		END_FOR
//		printf("\nat this step, there are %u liquid points while last step with %u difference is %d \n\n",
//				N, LiquidLastStep, N-LiquidLastStep);
//		LiquidLastStep = N;
//	}

}

void ScalarField3D::SeedParticles(const Voxel &voxel,
						ParticleList &negParticles,
#ifdef POSITIVE_PARTICLES
						ParticleList &posParticles,
#endif
						int nPerCell) const{

	float delta = voxel.VoxelDelta();

   ParticleList::iterator iter_particle;

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){

			bool seed = false;

			bool seedtwice = false;

			for(int m=0; m<8; ++m){
				Point corner = voxel.VoxelCornerPosition(i,j,k,m+1);
				float phi_i = TriInterp(voxel, corner);
				if(fabs(phi_i) < 3*delta){
					seed = true;
					break;
				}
			}
			for(int m=0; m<8; ++m){
				Point corner = voxel.VoxelCornerPosition(i,j,k,m+1);
				float phi_i = TriInterp(voxel, corner);
				if(fabs(phi_i) <= 0.5f*delta){
					seedtwice = true;
					break;
				}
			}
			if(seed){
				int N;
#ifdef POSITIVE_PARTICLES
//				if( seedtwice)
//					N = 2*nPerCell;
//				else
				N = nPerCell;
				int NP = 0, NM = 0;
				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
					float phi_i = TriInterp(voxel, pos);
					if(phi_i <= 0.f){
						float radius;
						if( -phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( -phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = -phi_i;
						ps.SetParticle(pos, radius);
						negParticles.push_back(ps);
						++NM;
					}
					else{
						float radius;
						if( phi_i > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( phi_i < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = phi_i;
						ps.SetParticle(pos, radius);
						posParticles.push_back(ps);
						++NP;
					}
//					if(i==I && j==J && k==K){
//						printf("m = %d, pos = (%f,%f,%f) phi_i = %f r = %f\n",
//								m, pos.x, pos.y, pos.z, phi_i,ps.Radius());
//					}
				}
				if(i==I && j==J && k==K){
					printf("at cell (%d,%d,%d) phi = %f partciels seeded \n",
							i,j,k, phi[INDEX(i,j,k)]);
					printf("there are total %d particles with %d negative and %d positive \n",
							N, NM, NP);
				}
#else
				for(int m=0; m<nPerCell; ++m){
					Particle ps;
					voxel.PlaceParticles(i, j, k, ps, 0);
					Point pos = ps.Position();
					negParticles.push_back(ps);
				}
#endif
			}

		}
	END_FOR

	printf("now there are %ld negative particles \n", negParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("now there are %ld positive particles \n", posParticles.size());


	printf("Checking positive particles... \n");

	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk) || phi_i < 0.f)
			printf("pos at (%f, %f,%f) in cell (%d, %d, %d) with phi = %f\n",
				pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
	}
#endif

	printf("Checking negative particles... \n");

	for(iter_particle = negParticles.begin();
		iter_particle != negParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos);
		if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk) || phi_i > 0.f ||
				radius > RADIUS_MAX || radius < RADIUS_MIN  )
			printf("pos at (%f, %f,%f) in cell (%d, %d, %d) with phi = %f\n",
				pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
	}



	printf("(ScalarField3D::SeedParticles) particles are seeded \n");

}

void ScalarField3D::AttractParticles(const Voxel &voxel, const ScalarField3D &object,
								list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
								list<Particle> &posParticles,
#endif
								int nPerCell) const{

	list<Particle>::iterator iter_particle;
	float delta = voxel.VoxelDelta();

	u_int  N1 = 0, N2 = 0, N3 = 0, N4 = 0, N5 = 0;

	for(iter_particle = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		//printf("A:\n");
		Vector N;
		NormalAt(voxel, pos, N);
		//printf("B:\n");
		float phi_i = TriInterp(voxel, pos);
		//printf("C:\n");
		MTRand mt;
		double rn = mt();
		float target_phi = -RADIUS_MIN - (float)rn * (PARTICLE_SEED_LIMIT - RADIUS_MIN);
//		printf("oldpos = (%f, %f, %f), newphi= %f, old_phi = %f \n",
//					pos.x, pos.y, pos.z,  new_phi, phi_i);
		Point new_pos = pos + (target_phi - phi_i) * N;
		Point ori_pos = pos;
		Point ori_new_pos = new_pos;
		float ori_phi = phi_i;
		Vector ori_N = N;
	//	printf("D:\n");
//		printf("newpos= (%f, %f, %f),              newphi= %f \n",
//			 new_pos.y, new_pos.z, new_phi);

		TentativeGridPoint tpp = voxel.ContainsPoint(new_pos);
		TentativeGridPoint ori_tpp = voxel.ContainsPoint(pos);;
		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){  //  out of the computational domain
			int cycles = 0;
			float lamda = 1.0f;
			do{
				lamda /= 2;
				cycles++;
				new_pos = pos + lamda * (target_phi - phi_i) * N;
				tpp = voxel.ContainsPoint(new_pos);
			}while(cycles < SEED_PARTICLE_ITERATIONS && !CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk));
			if(cycles >= SEED_PARTICLE_ITERATIONS){
				iter_particle = negParticles.erase(iter_particle);
				++N1;
			}
			else{
				float phi_obj = object.TriInterp(voxel, new_pos);
				if(phi_obj > 0.f){
					iter_particle = negParticles.erase(iter_particle);
					++N2;
				}
				else{
					float new_phi = TriInterp(voxel, new_pos);
					if(new_phi < -RADIUS_MIN && new_phi > -PARTICLE_SEED_LIMIT){
						if( -new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( -new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = -new_phi;
						ps.SetParticle(new_pos, radius);
						++iter_particle;
					}
					else{
						lamda /= 2;
						new_pos = pos + lamda * (target_phi - phi_i) * N;
						lamda = 1.f;
						cycles = 0;
						do{
							cycles++;
							pos = new_pos;
							phi_i = TriInterp(voxel, pos);
							NormalAt(voxel, pos, N);
							new_pos = pos + lamda * (target_phi - phi_i) * N;
							lamda /= 2;
						}while(cycles < SEED_PARTICLE_ITERATIONS &&
							(new_phi <= -PARTICLE_SEED_LIMIT || new_phi >= -RADIUS_MIN));
						if(cycles >= SEED_PARTICLE_ITERATIONS){
							iter_particle = negParticles.erase(iter_particle);
							++N3;
						}
						else{
							if( -new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
							else if( -new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
							else radius = -new_phi;
							ps.SetParticle(new_pos, radius);
							++iter_particle;
						}
					}
				}
			}
		}
		else{
			float phi_obj = object.TriInterp(voxel, new_pos);
			if(phi_obj > 0.f){
				iter_particle = negParticles.erase(iter_particle);
				++N4;
//				printf("oldpos = (%f, %f, %f), target_phi = %f, old_phi = %f \n",
//						pos.x, pos.y, pos.z,  target_phi, phi_i);
//				float new_phi = TriInterp(voxel, new_pos);
//				printf("newpos= (%f, %f, %f),  newphi = %f phi_obj = %f \n",
//							new_pos.x, new_pos.y, new_pos.z, new_phi, phi_obj);
			}
			else{
				float new_phi = TriInterp(voxel, new_pos);
				if(new_phi < -RADIUS_MIN && new_phi > -PARTICLE_SEED_LIMIT){
					if( -new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( -new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = -new_phi;
					ps.SetParticle(new_pos, radius);
					++iter_particle;
//					if(fabsf(new_phi-target_phi) >  1.e-3f){
//						printf("\nnew_phi = %f, target_phi = %f \n", new_phi, target_phi);
//						TentativeGridPoint spp = voxel.ContainsPoint(new_pos);
//						printf("target_phi = %f, old_phi = %f in cell (%d, %d, %d) N = (%f, %f, %f) \n",
//								target_phi, phi_i, spp.ii, spp.jj, spp.kk, N.x, N.y, N.z);
//						phi_obj = object.TriInterp(voxel, new_pos);
//						printf("newpos = (%f, %f, %f),  newphi = %f phi_obj = %f \n",
//									new_pos.x, new_pos.y, new_pos.z, new_phi, phi_obj);
//						printf("ori_pos = (%f, %f, %f), ori_new_pos = (%f, %f, %f)\n",
//							ori_pos.x, ori_pos.y, ori_pos.z, ori_new_pos.x, ori_new_pos.y, ori_new_pos.z);
//						printf("target_phi = %f, ori_phi = %f in ori_cell (%d, %d, %d) ori_N = (%f, %f, %f) \n",
//							target_phi, ori_phi, ori_tpp.ii, ori_tpp.jj, ori_tpp.kk, ori_N.x, ori_N.y, ori_N.z);
//						printf("phi at ori_cell = %f \n\n", phi[INDEX(ori_tpp.ii, ori_tpp.jj, ori_tpp.kk)]);
//					}
				}
				else{
					float lamda = 1.0f;
					int cycles = 0;
					do{
						cycles++;
						pos = new_pos;
						phi_i = TriInterp(voxel, pos);
						NormalAt(voxel, pos, N);
						new_pos = pos + lamda * (target_phi - phi_i) * N;
						new_phi = TriInterp(voxel, new_pos);
						if(new_phi < -RADIUS_MIN && new_phi > -PARTICLE_SEED_LIMIT)
							break;
						//lamda /= 2;
					}while(cycles < SEED_PARTICLE_ITERATIONS);
					if(cycles >= SEED_PARTICLE_ITERATIONS){
						if(ori_phi < -RADIUS_MIN && ori_phi > -PARTICLE_SEED_LIMIT){
							if( -ori_phi > RADIUS_MAX ) radius = RADIUS_MAX;
							else if( -ori_phi < RADIUS_MIN ) radius = RADIUS_MIN;
							else radius = -ori_phi;
							ps.SetParticle(ori_pos, radius);
							++iter_particle;
						}
						else{
							iter_particle = negParticles.erase(iter_particle);
							++N5;
						}
						TentativeGridPoint spp = voxel.ContainsPoint(new_pos);
//						printf("\ntarget_phi = %f, old_phi = %f in cell (%d, %d, %d) N = (%f, %f, %f) \n",
//								target_phi, phi_i, spp.ii, spp.jj, spp.kk, N.x, N.y, N.z);
//						phi_obj = object.TriInterp(voxel, new_pos);
//						printf("newpos = (%f, %f, %f),  newphi = %f phi_obj = %f \n",
//									new_pos.x, new_pos.y, new_pos.z, new_phi, phi_obj);
//						printf("ori_pos = (%f, %f, %f), ori_new_pos = (%f, %f, %f)\n",
//							ori_pos.x, ori_pos.y, ori_pos.z, ori_new_pos.x, ori_new_pos.y, ori_new_pos.z);
//						printf("target_phi = %f, ori_phi = %f in ori_cell (%d, %d, %d) ori_N = (%f, %f, %f) \n",
//							target_phi, ori_phi, ori_tpp.ii, ori_tpp.jj, ori_tpp.kk, ori_N.x, ori_N.y, ori_N.z);
//						printf("phi at ori_cell = %f \n\n", phi[INDEX(ori_tpp.ii, ori_tpp.jj, ori_tpp.kk)]);
					}
					else{
						if( -new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( -new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = -new_phi;
						ps.SetParticle(new_pos, radius);
						++iter_particle;
					}
				}
			}
		}
	}

	printf(" after attraction of negative N1 = %u, N2 = %u, N3 = %u, N4 = %u, N5 = %u \n",
			N1, N2, N3, N4, N5);
#ifdef POSITIVE_PARTICLES
	N1 = 0; N2 = 0; N3 = 0; N4 = 0; N5 = 0;
	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		//printf("A:\n");
		Vector N;
		NormalAt(voxel, pos, N);
		//printf("B:\n");
		float phi_i = TriInterp(voxel, pos);
		//printf("C:\n");
		MTRand mt;
		double rn = mt();
		float target_phi = RADIUS_MIN + (float)rn * (PARTICLE_SEED_LIMIT - RADIUS_MIN);
//		printf("oldpos = (%f, %f, %f), newphi= %f, old_phi = %f \n",
//					pos.x, pos.y, pos.z,  new_phi, phi_i);
		Point new_pos = pos + (target_phi - phi_i) * N;
		Point ori_pos = pos;
		Point ori_new_pos = new_pos;
		float ori_phi = phi_i;
		Vector ori_N = N;
	//	printf("D:\n");
//		printf("newpos= (%f, %f, %f),              newphi= %f \n",
//			 new_pos.y, new_pos.z, new_phi);

		TentativeGridPoint tpp = voxel.ContainsPoint(new_pos);
		if(!CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk)){  //  out of the computational domain
			int cycles = 0;
			float lamda = 1.0f;
			do{
				lamda /= 2;
				cycles++;
				new_pos = pos + lamda * (target_phi - phi_i) * N;
				tpp = voxel.ContainsPoint(new_pos);
			}while(cycles < SEED_PARTICLE_ITERATIONS && !CheckIndexNoExit(tpp.ii, tpp.jj, tpp.kk));
			if(cycles >= SEED_PARTICLE_ITERATIONS){
				iter_particle = posParticles.erase(iter_particle);
				++N1;
			}
			else{
				float phi_obj = object.TriInterp(voxel, new_pos);
				if(phi_obj > 0.f){
					iter_particle = posParticles.erase(iter_particle);
					++N2;
				}
				else{
					float new_phi = TriInterp(voxel, new_pos);
					if(new_phi > RADIUS_MIN && new_phi < PARTICLE_SEED_LIMIT){
						if( new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = new_phi;
						ps.SetParticle(new_pos, radius);
						++iter_particle;
					}
					else{
						lamda /= 2;
						new_pos = pos + lamda * (target_phi - phi_i) * N;
						lamda = 1.f;
						cycles = 0;
						do{
							cycles++;
							pos = new_pos;
							phi_i = TriInterp(voxel, pos);
							NormalAt(voxel, pos, N);
							new_pos = pos + lamda * (target_phi - phi_i) * N;
//							tpp = voxel.ContainsPoint(new_pos);
//							if(!CheckIndexNoExit(i,j,k)){
//								iter_particle = posParticles.erase(iter_particle);
//								continue;
//							}
							lamda /= 2;
						}while(cycles < SEED_PARTICLE_ITERATIONS &&
							(new_phi >= PARTICLE_SEED_LIMIT || new_phi <= RADIUS_MIN));
						if(cycles >= SEED_PARTICLE_ITERATIONS){
							iter_particle = posParticles.erase(iter_particle);
							++N3;
						}
						else{
							if( new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
							else if( new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
							else radius = new_phi;
							ps.SetParticle(new_pos, radius);
							++iter_particle;
						}
					}
				}
			}
		}
		else{
			float phi_obj = object.TriInterp(voxel, new_pos);
			if(phi_obj > 0.f){
				iter_particle = posParticles.erase(iter_particle);
				++N4;
			}
			else{
				float new_phi = TriInterp(voxel, new_pos);
				if(new_phi > RADIUS_MIN && new_phi < PARTICLE_SEED_LIMIT){
					if( new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
					else if( new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
					else radius = new_phi;
					ps.SetParticle(new_pos, radius);
					++iter_particle;
				}
				else{
					float lamda = 1.0f;
					int cycles = 0;
					do{
						cycles++;
						pos = new_pos;
						phi_i = TriInterp(voxel, pos);
						NormalAt(voxel, pos, N);
						new_pos = pos + lamda * (target_phi - phi_i) * N;
						new_phi = TriInterp(voxel, new_pos);
						if(new_phi > RADIUS_MIN && new_phi < PARTICLE_SEED_LIMIT)
							break;
						//lamda /= 2;
					}while(cycles < SEED_PARTICLE_ITERATIONS);
					if(cycles >= SEED_PARTICLE_ITERATIONS){
						if(ori_phi > RADIUS_MIN && ori_phi < PARTICLE_SEED_LIMIT){
							if( ori_phi > RADIUS_MAX ) radius = RADIUS_MAX;
							else if( ori_phi < RADIUS_MIN ) radius = RADIUS_MIN;
							else radius = ori_phi;
							ps.SetParticle(ori_pos, radius);
							++iter_particle;
						}
						else{
							iter_particle = posParticles.erase(iter_particle);
							++N5;
						}
					}
					else{
						if( new_phi > RADIUS_MAX ) radius = RADIUS_MAX;
						else if( new_phi < RADIUS_MIN ) radius = RADIUS_MIN;
						else radius = new_phi;
						ps.SetParticle(new_pos, radius);
						++iter_particle;
					}
				}
			}
		}
	}
	printf(" after attraction of positive N1 = %u, N2 = %u, N3 = %u, N4 = %u, N5 = %u \n",
				N1, N2, N3, N4, N5);
#endif
}



void ScalarField3D::SeedParticles(const Voxel &voxel, const ScalarField3D &object,
						list<Particle> &negParticles,
#ifdef POSITIVE_PARTICLES
						list<Particle> &posParticles,
#endif
						int nPerCell) const{

	float delta = voxel.VoxelDelta();

	list<Particle>::iterator iter_particle;
	FOR_EACH_CELL
	if(i == I && j == J && k == K){
		printf(" SeedParticles: phi = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
					 phi[INDEX(i,j,k)],
					phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)] );
		printf("SeedParticles: phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f \n",
					phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
					phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)] );
	}
   END_FOR

	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
	    int n;
		if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){
			bool seed = false;

			for(int m=0; m<8; ++m){
				Point corner = voxel.VoxelCornerPosition(i,j,k,m+1);
				float phi_i = TriInterp(voxel, corner);
				if(fabs(phi_i) < PARTICLE_SEED_LIMIT){
					seed = true;
					break;
				}
			}

			if(seed){
				int N = nPerCell;
				int NP = 0, NM = 0;
				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, -1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					negParticles.push_back(ps);
					++NM;
				}
#ifdef POSITIVE_PARTICLES
				for(int m=0; m<N; ++m){
					MTRand mt;
					double rn = mt();
					Particle ps(Point(0,0,0), 0, 1, (0.1f+(float)rn*0.9f)*delta);
					voxel.PlaceParticles(i, j, k, ps, 0);
					posParticles.push_back(ps);
					++NP;
				}
#endif
				if(i==I && j==J && k==K){
					printf("at cell (%d,%d,%d) phi = %f partciels seeded \n",
							i,j,k, phi[INDEX(i,j,k)]);
					printf("there are total %d particles with %d negative and %d positive \n",
							N, NM, NP);
				}
			}

		}
	END_FOR

	printf("now there are %u negative particles \n", negParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("now there are %u positive particles \n", posParticles.size());
#endif

	AttractParticles(voxel, object, negParticles,
#ifdef POSITIVE_PARTICLES
			posParticles,
#endif
			nPerCell);


#ifdef POSITIVE_PARTICLES

	printf("Checking positive particles... \n");

	u_int NP = 0;
	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = object.TriInterp(voxel, pos);
		if(phi_i > 0.f){
//			printf("positive particle pos at (%f, %f, %f) with phi_obj = %f\n",
//				pos.x, pos.y, pos.z, phi_i);
			iter_particle = posParticles.erase(iter_particle);
			++NP;
		}
		else{
			phi_i = TriInterp(voxel, pos);
			TentativeGridPoint tpp = voxel.ContainsPoint(pos);
			if(phi_i < RADIUS_MIN || phi_i > PARTICLE_SEED_LIMIT)
				printf("positive particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
					pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
			++iter_particle;
		}
	}
	printf("there are %u positive particles removed due to its location in solid \n", NP);
#endif



	printf("Checking negative particles... \n");
	u_int NM = 0;
	for(iter_particle = negParticles.begin();
		iter_particle != negParticles.end();
		){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float phi_i = object.TriInterp(voxel, pos);
		if(phi_i > 0.f){
//			printf("negative particle pos at (%f, %f, %f) with phi_obj = %f\n",
//				pos.x, pos.y, pos.z, phi_i);
			iter_particle = negParticles.erase(iter_particle);
			++NM;
		}
		else{
			phi_i = TriInterp(voxel, pos);
			TentativeGridPoint tpp = voxel.ContainsPoint(pos);
			if(phi_i > -RADIUS_MIN || phi_i < -PARTICLE_SEED_LIMIT)
				printf("negative particle pos at (%f, %f, %f) in cell (%d, %d, %d) with phi_i = %f\n",
					pos.x, pos.y, pos.z, tpp.ii, tpp.jj, tpp.kk, phi_i);
			++iter_particle;
		}
	}
	printf("there are %u negative particles removed due to its location in solid \n", NM);

	printf("after seeding now there are %u negative particles \n", negParticles.size());
#ifdef POSITIVE_PARTICLES
	printf("after seeding now there are %u positive particles \n", posParticles.size());
#endif

	printf("(ScalarField3D::SeedParticles) particles are seeded \n");

}

void ScalarField3D::AdvectEscapedParticles(const Voxel &voxel,
		          list<Particle> &posParticles,
		          list<Particle> &negParticles,
		          float dt) const {

	list<Particle>::iterator iter_particle;
	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		){
		Particle &ps = *iter_particle;
		ps.MoveUnderGravity(dt);
		Point pos = ps.Position();
		float radius = ps.Radius();
		TentativeGridPoint tpp = voxel.ContainsPoint(pos);
		if(voxel.InSolid(tpp.ii, tpp.jj, tpp.kk))
			iter_particle = posParticles.erase(iter_particle);
		else{
			float phi_i = TriInterp(voxel, pos);
			if(phi_i < 1.5*radius ){
				negParticles.push_back(ps);
				iter_particle = posParticles.erase(iter_particle);
			}
			else
				++iter_particle;
		}
	}
}


static void ProcessCoefficents(float phi_minus, float phi_plus, float minvalue,
						float vel_minus, float vel_plus,
						float &phi, float &X, float &x){
//	if( phi_minus < minvalue && phi_plus < minvalue){
//    	if(fabs(phi_minus) < fabs(phi_plus)){
//    		phi = minvalue - phi_minus;
//    		X = phi;
//    		x = -phi*vel_minus;
//    	}
//    	else{
//    		phi = phi_plus - minvalue;
//    		X = -phi;
//    		x = phi*vel_plus;
//    	}
//    }
//    else if(phi_minus < minvalue){
//    	phi = minvalue - phi_minus;
//    	X = phi;
//		x = -phi*vel_minus;
//    }
//    else if(phi_plus < minvalue){
//    	phi = phi_plus - minvalue;
//    	X = -phi;
//		x = phi*vel_plus;
//    }
//    else{
//    	phi = 0.f;
//    	X = 0.f;
//    	x = 0.f;
//    }
	   if( phi_minus < minvalue && phi_plus < minvalue){
		   if(phi_minus < 0.f && phi_plus < 0.f){
			   // information should come from closest point at interface
		    	if(fabs(phi_minus) > fabs(phi_plus))
		    		phi = minvalue - phi_minus;
		    	else
		    		phi = phi_plus - minvalue;
		   }
		   else{
		    	if(phi_minus < phi_plus)
		    		phi = minvalue - phi_minus;
		    	else
		    		phi = phi_plus - minvalue;
		   }
	    }
	    else if(phi_minus == minvalue &&  phi_plus > minvalue)
	  	    phi = minvalue - phi_minus;
	  	else if(phi_plus == minvalue &&  phi_minus > minvalue)
	  	   	phi = phi_plus - minvalue;
	    else if(phi_minus < minvalue){
	    	phi = minvalue - phi_minus;
	    }
	    else if(phi_plus < minvalue){
	    	phi = phi_plus - minvalue;
	    }
	    else{
	    	phi = 0.f;
	    }

	   if(phi > 0.f){
		   X = phi;
		   x = phi*vel_minus;
	   }
	   else if (phi < 0.f){
		   X = -phi;
  		   x = -phi*vel_plus;
	   }
	   else{
		   X = 0.f;
       	   x = 0.f;
	   }


	return;
}

#define EXTRAPOLATE_PHI_INTO_OBJECT_LIMIT 2*delta

#ifdef SPMD
static void MarchingOut(const TentativeGridPoint &minTentative, int DimX, int DimY, int DimZ,
				float dtau, float inv_delta2, float *phi_tmp, float *vel, float *vel0){

	float phi_minus_x = INFINITY, phi_plus_x = INFINITY;
	float phi_minus_y = INFINITY, phi_plus_y = INFINITY;
	float phi_minus_z = INFINITY, phi_plus_z = INFINITY;
	float vel_minus_x, vel_minus_y, vel_minus_z;
	float vel_plus_x, vel_plus_y, vel_plus_z;
	float phi_x, phi_y, phi_z;
	float vel_x, vel_y, vel_z;
	float A=0.f, B=0.f, C=0.f;
	float a=0.f, b=0.f, c=0.f;
	float minValue = minTentative.value;
	int i = minTentative.ii;
	int j = minTentative.jj;
	int k = minTentative.kk;
	u_int pos = INDEX(i,j,k);
	vector<TentativeGridPoint> neighbors;
	minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
	for(int m=0; m<neighbors.size(); m++){
		if(minTentative.LeftNeighbor(neighbors[m])){
			phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			vel_minus_x = vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		}
		if(minTentative.RightNeighbor(neighbors[m])){
			phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			vel_plus_x = vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		}
		if(minTentative.BackNeighbor(neighbors[m])){
			phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			vel_minus_y = vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		}
		if(minTentative.FrontNeighbor(neighbors[m])){
			phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			vel_plus_y = vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		}
		if(minTentative.BottomNeighbor(neighbors[m])){
			phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			vel_minus_z = vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		}
		if(minTentative.TopNeighbor(neighbors[m])){
			phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			vel_plus_z = vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		}
	}
	ProcessCoefficents(phi_minus_x, phi_plus_x, minValue,
					vel_minus_x, vel_plus_x, phi_x, A, a);
	ProcessCoefficents(phi_minus_y, phi_plus_y, minValue,
					vel_minus_y, vel_plus_y, phi_y, B, b);
	ProcessCoefficents(phi_minus_z, phi_plus_z, minValue,
					vel_minus_z, vel_plus_z, phi_z, C, c);

	if(phi_x != 0.f || phi_y != 0.f || phi_z != 0.f){
		if(phi_x > 0.f)
			vel_x = vel0[pos] - vel_minus_x;
		else if(phi_x < 0.f)
			vel_x = vel_plus_x - vel0[pos];
		else
			vel_x = 0.f;
		if(phi_y > 0.f)
			vel_y = vel0[pos] - vel_minus_y;
		else if(phi_y < 0.f)
			vel_y = vel_plus_y - vel0[pos];
		else
			vel_y = 0.f;
		if(phi_z > 0.f)
			vel_z = vel0[pos] - vel_minus_z;
		else if(phi_z < 0.f)
			vel_z = vel_plus_z - vel0[pos];
		else
			vel_z = 0.f;
		float mag = sqrtf(inv_delta2 * (phi_x * phi_x + phi_y * phi_y + phi_z * phi_z));
		vel[pos] = vel0[pos] - dtau * inv_delta2 / mag * (phi_x * vel_x +  phi_y * vel_y + phi_z * vel_z);
	}

}

#ifdef TBB
struct ExtrapolateBody {
	float * const t1;
	float * const t2;
	float * const t3;
	const Voxel * const voxel;
	float limit;
	int DimX, DimY, DimZ;
	ExtrapolateBody(const Voxel * const vox, float *s1, float *s2, float *s3,
			     float l, int nx, int ny, int nz)
		: t1(s1), t2(s2), t3(s3), voxel(vox),
		   limit(l), DimX(nx), DimY(ny), DimZ(nz){
//		printf("DimX = %d, DimY = %d, DimZ = %d\n ", DimX, DimY, DimZ);
	}
    void operator()( const tbb::blocked_range<int>& range ) const {
    	float *phi = t1;
    	float *phi_tmp = t2;
    	float *phi_tmp0 = t3;


		float delta = voxel->VoxelDelta();
		float delta2 = delta * delta;
		float inv_delta2 = 1.f / delta2;
		float dtau = delta / 2 ;

//		printf("got here \n");
        int k_end = range.end();
        for( int k=range.begin(); k!=k_end; ++k ){
        	for(int j=0;j<DimY;j++) {
        		for(int i=0;i<DimX;i++){
        			u_int pos = INDEX(i,j,k);
        			if(voxel->InSolid(i,j,k) && phi[pos] < limit){
						KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi[pos], tmp);
						MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi, phi_tmp, phi_tmp0);
//						if(i == I && j == J && k == K){
//							printf(" iter = %d, i = %d, j = %d, k = %d, phi = %f, phi[j+1] = %f, phi[i-1] = %f\n",
//									n,i,j,k, phi_tmp[INDEX(i,j,k)], phi_tmp[INDEX(i,j+1,k)], phi_tmp[INDEX(i-1,j,k)]);
//						}
					}
        		}
        	}
        }
  }
};
#endif

void ScalarField3D::Extrapolate(const Voxel &voxel, float *phi_tmp, float *phi_tmp0, bool limit) const{

	float delta = voxel.VoxelDelta();

	printf("Start extrapolating phi into solid object\n ");
#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
	if(limit){
		int iterations = int(round(LS_ADV_LIMIT / delta * 2));
		for (int n = 0; n < iterations; ++n){
			SWAP(phi_tmp, phi_tmp0);
			ExtrapolateBody body(&voxel, phi, phi_tmp, phi_tmp0,
					LS_ADV_LIMIT, DimX, DimY, DimZ);
			tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
							  body);
		}
	}
	else{
		int iterations = int(round(EXTRAPOLATE_PHI_INTO_OBJECT_LIMIT / delta * 2));
		for (int n = 0; n < iterations; ++n){
			SWAP(phi_tmp, phi_tmp0);
			ExtrapolateBody body(&voxel, phi, phi_tmp, phi_tmp0,
					EXTRAPOLATE_PHI_INTO_OBJECT_LIMIT, DimX, DimY, DimZ);
			tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
							  body);
		}
	}
#else
	float inv_delta2 = 1.f / (delta * delta);
	float dtau = delta / 6;

#ifdef CUDA
	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k))
			obj[INDEX(i,j,k)] = 1;
		else
			obj[INDEX(i,j,k)] = 0;
	 END_FOR

	int dims[3];
	dims[0] = DimX; dims[1] = DimY; dims[2] = DimZ;
	float slimit;
	if(limit)
		slimit = LS_ADV_LIMIT;
	else
		slimit = EXTRAPOLATE_PHI_INTO_OBJECT_LIMIT;

	ExtrapolatePhiCUDA( phi_tmp, phi, obj, I, J, K,
	 				    delta, dtau, INFINITY, slimit, dims );

#else

	if(limit){
		int iterations = int(round(LS_ADV_LIMIT / delta * 2));
		for (int n = 0; n < iterations; ++n){
			SWAP(phi_tmp, phi_tmp0);
			FOR_EACH_CELL
				if(voxel.InSolid(i,j,k) && phi[INDEX(i,j,k)] < LS_ADV_LIMIT){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi, phi_tmp, phi_tmp0);
					if(i == I && j == J && k == K){
						printf(" iter = %d, i = %d, j = %d, k = %d, phi = %f, phi[j+1] = %f, phi[i-1] = %f\n",
								n,i,j,k, phi_tmp[INDEX(i,j,k)], phi_tmp[INDEX(i,j+1,k)], phi_tmp[INDEX(i-1,j,k)]);
					}
				}
			END_FOR
		}
	}
	else{
		int iterations = int(round(EXTRAPOLATE_PHI_INTO_OBJECT_LIMIT / delta * 2));
		for (int n = 0; n < iterations; ++n){
			SWAP(phi_tmp, phi_tmp0);
			FOR_EACH_CELL
				if(voxel.InSolid(i,j,k) && phi[INDEX(i,j,k)] < EXTRAPOLATE_PHI_INTO_OBJECT_LIMIT ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi, phi_tmp, phi_tmp0);
				}
			END_FOR
		}
	}
#endif
#endif
	printf("End extrapolating phi into solid object\n ");
}
#else
void ScalarField3D::Extrapolate(const Voxel &voxel, float *phi_tmp, bool limit) const{

	float delta = voxel.VoxelDelta();

	float minValue;
	MinHeap<TentativeGridPoint, float> heap;

	printf("Start extrapolating phi into solid object\n ");

	if(limit){
		FOR_EACH_CELL
			if(voxel.InSolid(i,j,k)){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi[INDEX(i,j,k)]);
			}
		END_FOR
	}
	else{
		FOR_EACH_CELL
			if(voxel.InSolid(i,j,k) && phi[INDEX(i,j,k)] < EXTRAPOLATE_PHI_INTO_OBJECT_LIMIT){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi[INDEX(i,j,k)]);
			}
		END_FOR
	}

//	printf("Extrapolating velocity: b = %d with %ld in heap \n ",
//			    vel_index, heap.size());

	do{
		float phi_minus_x = INFINITY, phi_plus_x = INFINITY;
		float phi_minus_y = INFINITY, phi_plus_y = INFINITY;
		float phi_minus_z = INFINITY, phi_plus_z = INFINITY;
		float vel_minus_x, vel_minus_y, vel_minus_z;
		float vel_plus_x, vel_plus_y, vel_plus_z;
		float phi_x, phi_y, phi_z;
		float A=0.f, B=0.f, C=0.f;
		float a=0.f, b=0.f, c=0.f;
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		vector<TentativeGridPoint> neighbors;
	    minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
	    for(int m=0; m<neighbors.size(); m++){
	    	if(minTentative.LeftNeighbor(neighbors[m])){
	    		phi_minus_x = phi[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.RightNeighbor(neighbors[m])){
	    		phi_plus_x = phi[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.BackNeighbor(neighbors[m])){
	    		phi_minus_y = phi[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.FrontNeighbor(neighbors[m])){
	    		phi_plus_y = phi[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.BottomNeighbor(neighbors[m])){
	    		phi_minus_z = phi[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.TopNeighbor(neighbors[m])){
	    		phi_plus_z = phi[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    }

	    ProcessCoefficents(phi_minus_x, phi_plus_x, minValue,
	    				vel_minus_x, vel_plus_x, phi_x, A, a);
	    ProcessCoefficents(phi_minus_y, phi_plus_y, minValue,
	    	    		vel_minus_y, vel_plus_y, phi_y, B, b);
	    ProcessCoefficents(phi_minus_z, phi_plus_z, minValue,
	    	    		vel_minus_z, vel_plus_z, phi_z, C, c);

//	    if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K && vel_index == 3){
////	    if(vel_index == 3){
//	    	printf("\n before interp: at i = %d, j = %d, k = %d vel = %f  b = %d\n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//	    			vel_index);
//	    	printf("before interp: at i = %d, j = %d, k = %d, phi = %f, phi-x = %f, phi+x = %f ,phi-y = %f, phi+y = %f, phi-z = %f, phi+z = %f \n ",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			minValue, phi_minus_x, phi_plus_x, phi_minus_y, phi_plus_y,
//	    			phi_minus_z, phi_plus_z);
//	    	printf("before interp: at i = %d, j = %d, k = %d, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
//	    			vel_minus_z, vel_plus_z);
//	    	printf("before interp: at i = %d, j = %d, k = %d, A = %f, B = %f, C = %f, a = %f, b = %f, c = %f \n ",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			A, B, C, a, b, c);
//
//	    }

	    if( (A+B+C) > 1.e-10 ){
//	    	printf("ii = %d, jj = %d, kk = %d, phi = %f, phi_x = %f, phi_y = %f, phi_z = %f \n ",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			minValue, phi_x, phi_y, phi_z);
	    	phi_tmp[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
	    }
	    else{
//	    	exit(1);
	    }
//	    if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
////	    if(vel_index == 3){
////	    if(fabs(vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]) > 40.f){
////	    if(phi_minus_x < 0.f && phi_plus_x < 0.f ||
////    	   phi_minus_y < 0.f && phi_plus_y < 0.f ||
////    	   phi_minus_z < 0.f && phi_plus_z < 0.f ) {
////	    if( phi_x == 0.f && A != 0.f || phi_y == 0.f && B != 0.f || phi_z == 0.f && C != 0.f ){
//	    	printf("(ExtrapolateOneVelocity) vel_index = %d\n", vel_index);
//	    	printf("after interp: at i = %d, j = %d, k = %d vel = %f \n\n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]);
//	    	printf("before interp: A = %f, B = %f, C = %f, a = %f, b = %f, c = %f \n ",
//	    			A, B, C, a, b, c);
//	    	printf(" phix = %f, phiy = %f, phiz = %f \n",
//	    			phi_x, phi_y, phi_z);
//	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
//	    			minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
//	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
//		    			minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
//	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
//	    		    	minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
////	    	exit(1);
//	    }

	} while(!heap.empty());

	printf("End extrapolating phi into solid object\n ");
}
#endif

/*void ScalarField3D::Initialize( Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band){
	float r;
	float *phi_tmp = new float[DimX*DimY*DimZ];
	memset(phi_tmp, 0, DimX*DimY*DimZ);

	FOR_EACH_CELL
		if(voxel.IsDone(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, DONE);
		if(voxel.IsClose(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, CLOSE);
	END_FOR

// comment the follwoing lines for Zalesak's disk
// uncomment for water fall and pourting water test case
	u_int air=0, liquid=0;
	FOR_EACH_CELL
		if(voxel.InLiquid(i,j,k)){
			phi[INDEX(i,j,k)] *= -1;
			liquid++;
		}
		else if(voxel.InAir(i,j,k))
			air++;

	END_FOR
	printf("before marching air points = %d, liquid points = %d \n", air, liquid);



	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			if(phi[INDEX(i,j,k)] == 0.f)
				phi_tmp[INDEX(i,j,k)] = 0.f;
			else{
				TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
				vector<TentativeGridPoint> neighbors;
				tp.Neigbor(neighbors, DimX, DimY, DimZ);
				Point p = voxel.VoxelCenterPosition(i,j,k);
				float *candidates = new float[neighbors.size()];
				for(int m=0; m<neighbors.size(); m++){
					candidates[m] = INFINITY;
					Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
														 neighbors[m].jj,
														 neighbors[m].kk);
					neighbors[m].value = phi[POS(neighbors[m])];
					if(voxel.InSolid(neighbors[m]))
						continue;
					else if(tp.value*neighbors[m].value == 0.f){
						candidates[m] = fabs(tp.value);
					}
					else if(tp.value*neighbors[m].value < 0.f){

						if(tp.LeftNeighbor(neighbors[m])){
							r = (tp.value * pn.x - neighbors[m].value * p.x) /
							       (tp.value - neighbors[m].value);
							candidates[m] = p.x - r;
						}
						else if(tp.RightNeighbor(neighbors[m])){
							r = (tp.value * pn.x - neighbors[m].value * p.x) /
							       (tp.value - neighbors[m].value);
							candidates[m] = r - p.x;
						}
						else if(tp.BackNeighbor(neighbors[m])){
							r = (tp.value * pn.y - neighbors[m].value * p.y) /
							       (tp.value - neighbors[m].value);
							candidates[m] = p.y - r;
						}
						else if(tp.FrontNeighbor(neighbors[m])){
							r = (tp.value * pn.y - neighbors[m].value * p.y) /
							       (tp.value - neighbors[m].value);
							candidates[m] = r - p.y;
						}
						else if(tp.BottomNeighbor(neighbors[m])){
							r = (tp.value * pn.z - neighbors[m].value * p.z) /
							       (tp.value - neighbors[m].value);
							candidates[m] = p.z - r;
						}
						else if(tp.TopNeighbor(neighbors[m])){
							r = (tp.value * pn.z - neighbors[m].value * p.z) /
							       (tp.value - neighbors[m].value);
							candidates[m] = r - p.z;
						}
						if(candidates[m] < 0.f){
							printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
									i,j,k, tp.value);
							printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
									candidates[m], r, pn.x, pn.y, pn.z);

						}
//						if(i == I && j == J && k == K){
//							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
//									i, j, k, tp.value, neighbors[m].value, candidates[m]);
//							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
//												neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
//												neighbors[m].value, phi[POS(neighbors[m])]);
//						}
					}
					if(i == I && j == J && k == K){
						printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
							i, j, k, tp.value, neighbors[m].value, candidates[m]);
						printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
								neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
								neighbors[m].value, phi[POS(neighbors[m])]);

					}
				}
				float phi_min = FindMinimum(candidates, neighbors.size());
		//		if(phi_min < 0.f){
		//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
		//					i,j,k, tp.value);
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", candidates[n]);
		//			printf("\n");
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", neighbors[n].value);
		//			printf("\n");
		//		}
				//if(phi_min == INFINITY)
				//	continue;
		//		printf("before i = %d, j= %d, k = %d, phi = %f \n",
		//						i, j, k, phi[INDEX(i,j,k)]);
				//KnownPoint tmp(i, j, k);
				if(tp.value > 0.f){
					phi_tmp[INDEX(i,j,k)] = phi_min;
					//if(phi_min != INFINITY)
						//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				else{
					phi_tmp[INDEX(i,j,k)] = -phi_min;
		//			if(phi_min != INFINITY)
		//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				if(candidates) delete [] candidates;
		//		printf("after i = %d, j= %d, k = %d, phi = %f \n",
		//				i, j, k, phi[INDEX(i,j,k)]);
			}
		}
		else
			phi_tmp[INDEX(i,j,k)] = INFINITY;
	END_FOR

	FOR_EACH_CELL
		KnownPoint tmp(i, j, k);
		if(phi_tmp[INDEX(i,j,k)] != INFINITY && phi_tmp[INDEX(i,j,k)] != -INFINITY){
			if(i == I && j == J && k == K)
				printf("grid point (%d,%d,%d) in band \n", i, j, k);
			if(phi_tmp[INDEX(i,j,k)] == 0.f){
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				if(phi[INDEX(i,j,k)] > 0.f)
					pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				else
					neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			voxel.UpdateCellFastMarchType(i, j, k, SET, DONE);
		}
	END_FOR

//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);
	if(phi_tmp) delete [] phi_tmp;


}*/
#define MAX_LOSSASO_ITERATIONS 10
#define EPSILON 1.e-5
float ScalarField3D::LossasoDistance(const Voxel &voxel, const Point &p, float phimin) const{
	char iter = 0;
	float dist = phimin;
	Point pp = p;
	Vector N;
	float phi_i;
//	TentativeGridPoint tp = voxel.ContainsPoint(p);
	do{
		phi_i = TriInterp(voxel, pp);
		NormalAt(voxel, pp, N);
		pp += (-phi_i) * N;
		dist = (pp - p).Length();
//		printf(" at (%d, %d, %d) and iter = %d,  dist = %f, phimin = %f \n", tp.ii, tp.jj, tp.kk, iter, dist, phimin);
		if(fabs(phi_i) < EPSILON)
			return dist;
		++iter;
	}while(iter < MAX_LOSSASO_ITERATIONS);

	return phimin;
}

void ScalarField3D::Initialize( Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band,
								const ScalarField3D &solid){
	float r;
	float delta = voxel.VoxelDelta();
//	float *phi_tmp = new float[DimX*DimY*DimZ];
	float *phi_tmp = phi0;
//	memset(phi_tmp, 0, DimX*DimY*DimZ);
	SetZero(phi_tmp);
    pos_band.clear();
    neg_band.clear();
	float *phi_obj = solid.getScalar();

	FOR_EACH_CELL
		if(voxel.IsDone(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, DONE);
		if(voxel.IsClose(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, CLOSE);
	END_FOR

// comment the following lines for Zalesak's disk
// uncomment for water fall and pourting water test case
	u_int air=0, liquid=0;
	FOR_EACH_CELL
		if(voxel.InLiquid(i,j,k)){
			if(phi[INDEX(i,j,k)] > 0.f){
				phi[INDEX(i,j,k)] *= -1;
				liquid++;
				if(i == I && j == J && k == K)
							printf(" i = %d, j = %d, k = %d, phi = %f \n",
								     i,j,k, phi[INDEX(i,j,k)]);
			}
		}
		else if(voxel.InAir(i,j,k))
			air++;
//	if(!voxel.InSolid(i,j,k)){
//		if(phi[INDEX(i,j,k)] <= 0.f)
//			liquid++;
//		else
//			air++;
//	}

	phi_tmp[INDEX(i,j,k)] = INFINITY;
	END_FOR
	printf("before marching air points = %d, liquid points = %d \n", air, liquid);




	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			if(phi[INDEX(i,j,k)] == 0.f)
				phi_tmp[INDEX(i,j,k)] = 0.f;
			else{
				TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
				vector<TentativeGridPoint> neighbors;
				tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
				Point p = voxel.VoxelCenterPosition(i,j,k);
				float *candidates = new float[6];
				for(int m=0; m<6; m++)
					candidates[m] = INFINITY;
//				Point *candidatesP = new Point[6];
				for(int m=0; m<6; m++){
					if(tp.NeighborExists(m, DimX, DimY, DimZ)){
						Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
															 neighbors[m].jj,
															 neighbors[m].kk);
						neighbors[m].value = phi[POS(neighbors[m])];
						if(voxel.InSolid(neighbors[m])){
							if(tp.value <= 0.f && voxel.InMovingSolid(neighbors[m])){
								if(tp.LeftNeighbor(neighbors[m])){
									r = delta * phi_obj[POS(neighbors[m])] /
									       (phi_obj[POS(neighbors[m])] - phi_obj[INDEX(i,j,k)]);
									candidates[m] = min(fabsf(tp.value), (delta - r));
								}
								else if(tp.RightNeighbor(neighbors[m])){
									r = delta * phi_obj[INDEX(i,j,k)] /
										(phi_obj[INDEX(i,j,k)] - phi_obj[POS(neighbors[m])]);
									candidates[m] = min(fabsf(tp.value), r);
								}
								else if(tp.BackNeighbor(neighbors[m])){
									r = delta * phi_obj[POS(neighbors[m])] /
									       (phi_obj[POS(neighbors[m])] - phi_obj[INDEX(i,j,k)]);
									candidates[m] = min(fabsf(tp.value), (delta - r));
								}
								else if(tp.FrontNeighbor(neighbors[m])){
									r = delta * phi_obj[INDEX(i,j,k)] /
										(phi_obj[INDEX(i,j,k)] - phi_obj[POS(neighbors[m])]);
									candidates[m] = min(fabsf(tp.value), r);
								}
								else if(tp.BottomNeighbor(neighbors[m])){
									r = delta * phi_obj[POS(neighbors[m])] /
									       (phi_obj[POS(neighbors[m])] - phi_obj[INDEX(i,j,k)]);
									candidates[m] = min(fabsf(tp.value), (delta - r));
								}
								else if(tp.TopNeighbor(neighbors[m])){
									r = delta * phi_obj[INDEX(i,j,k)] /
										(phi_obj[INDEX(i,j,k)] - phi_obj[POS(neighbors[m])]);
									candidates[m] = min(fabsf(tp.value), r);
								}

							}
							else
								continue;
						}
//						if(voxel.InSolid(neighbors[m]))
//							continue;
//						else if(tp.value*neighbors[m].value == 0.f){
//							candidates[m] = delta;
//							candidatesP[m] = pn;
//						}
						else if(tp.value*neighbors[m].value <= 0.f){

							if(tp.LeftNeighbor(neighbors[m])){
								r = (tp.value * pn.x - neighbors[m].value * p.x) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.x - r;
//								candidatesP[m] = Point(r, p.y, p.z);
							}
							else if(tp.RightNeighbor(neighbors[m])){
								r = (tp.value * pn.x - neighbors[m].value * p.x) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.x;
//								candidatesP[m] = Point(r, p.y, p.z);
							}
							else if(tp.BackNeighbor(neighbors[m])){
								r = (tp.value * pn.y - neighbors[m].value * p.y) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.y - r;
//								candidatesP[m] = Point(p.x, r, p.z);
							}
							else if(tp.FrontNeighbor(neighbors[m])){
								r = (tp.value * pn.y - neighbors[m].value * p.y) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.y;
//								candidatesP[m] = Point(p.x, r, p.z);
							}
							else if(tp.BottomNeighbor(neighbors[m])){
								r = (tp.value * pn.z - neighbors[m].value * p.z) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.z - r;
//								candidatesP[m] = Point(p.x, p.y, r);
							}
							else if(tp.TopNeighbor(neighbors[m])){
								r = (tp.value * pn.z - neighbors[m].value * p.z) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.z;
//								candidatesP[m] = Point(p.x, p.y, r);
							}
							if(candidates[m] < 0.f){
								printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
										i,j,k, tp.value);
								printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
										candidates[m], r, pn.x, pn.y, pn.z);
								candidates[m] = fabs(candidates[m]);
							}
	//						if(i == I && j == J && k == K){
	//							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//									i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//												neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//												neighbors[m].value, phi[POS(neighbors[m])]);
	//						}
						}
						if(i == I && j == J && k == K){
							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
								i, j, k, tp.value, neighbors[m].value, candidates[m]);
							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
									neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
									neighbors[m].value, phi[POS(neighbors[m])]);

						}
					}
				}
				float phi_min[3], phimin = 0.f;
				Point phi_min_p[3];
				phi_min[0]= min(candidates[0], candidates[1]);
//				phi_min_p[0] = (candidates[0] < candidates[1]) ? candidatesP[0]: candidatesP[1];
				phi_min[1]= min(candidates[2], candidates[3]);
//				phi_min_p[1] = (candidates[2] < candidates[3]) ? candidatesP[2]: candidatesP[3];
				phi_min[2]= min(candidates[4], candidates[5]);
//				phi_min_p[2] = (candidates[4] < candidates[5]) ? candidatesP[4]: candidatesP[5];
				if( phi_min[0] != INFINITY &&
					phi_min[1] != INFINITY &&
					phi_min[2] != INFINITY ){  // 3 points forming a plane
//					Vector v1 = phi_min_p[0] - phi_min_p[1];
//					Vector v2 = phi_min_p[0] - phi_min_p[2];
//					Vector v3 =  Cross(v1, v2);
//					v3 = Normalize(v3);
//					Vector v4 = p - phi_min_p[0];
//					phimin = AbsDot(v3, v4);
					phimin = phi_min[0]*phi_min[1]*phi_min[2] /
					       sqrt( phi_min[0]*phi_min[0]*phi_min[1]*phi_min[1] +
					    		 phi_min[0]*phi_min[0]*phi_min[2]*phi_min[2] +
					    		 phi_min[2]*phi_min[2]*phi_min[1]*phi_min[1]);
				}
				else if( phi_min[0] != INFINITY && phi_min[1] != INFINITY ){
					// 2 points forming a line segment
//					phimin = DistanceToLine(phi_min_p[0], phi_min_p[1], p, tp.value);
					phimin = TwoPointDistance(phi_min[0], phi_min[1]);
				}
				else if( phi_min[1] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
//					phimin = DistanceToLine(phi_min_p[1], phi_min_p[2], p, tp.value);
					phimin = TwoPointDistance(phi_min[1], phi_min[2]);
				}
				else if( phi_min[0] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
//					phimin = DistanceToLine(phi_min_p[0], phi_min_p[2], p, tp.value);
					phimin = TwoPointDistance(phi_min[0], phi_min[2]);
				}
				else if(phi_min[0] != INFINITY)
					phimin = phi_min[0];
				else if(phi_min[1] != INFINITY)
					phimin = phi_min[1];
				else if(phi_min[2] != INFINITY)
					phimin = phi_min[2];
				else{
//					printf("wrong! at least 1 point should be valid at [%d, %d, %d] with phi = %f \n",
//						i,j,k, tp.value);
//					for(int n=0; n < neighbors.size(); ++n)
//						printf("%f ", candidates[n]);
//					printf("\n");
//					for(int n=0; n < neighbors.size(); ++n)
//						printf("%f ", neighbors[n].value);
//					printf("\n");
//					exit(1);
					phimin = INFINITY;
				}
		//		if(phi_min < 0.f){
		//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
		//					i,j,k, tp.value);
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", candidates[n]);
		//			printf("\n");
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", neighbors[n].value);
		//			printf("\n");
		//		}
				//if(phi_min == INFINITY)
				//	continue;
		//		printf("before i = %d, j= %d, k = %d, phi = %f \n",
		//						i, j, k, phi[INDEX(i,j,k)]);
				//KnownPoint tmp(i, j, k);
//				if(phimin != INFINITY && !voxel.CloseToSolid(tp)){
//					float lossdist = LossasoDistance(voxel, p, phimin);
//					phimin = min(phimin, lossdist);
//				}
				if(tp.value > 0.f){
					phi_tmp[INDEX(i,j,k)] = phimin;
					//if(phi_min != INFINITY)
						//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				else{
//					if(voxel.CloseToSolid(tp))
//						phimin = min(phimin, -phi_obj[INDEX(i,j,k)]);
					phi_tmp[INDEX(i,j,k)] = -phimin;
		//			if(phi_min != INFINITY)
		//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				if(candidates) delete [] candidates;
//				if(candidatesP) delete [] candidatesP;
		//		printf("after i = %d, j= %d, k = %d, phi = %f \n",
		//				i, j, k, phi[INDEX(i,j,k)]);
			}
		}
		else
			phi_tmp[INDEX(i,j,k)] = INFINITY;
	END_FOR

	FOR_EACH_CELL
		KnownPoint tmp(i, j, k);
		if(phi_tmp[INDEX(i,j,k)] != INFINITY && phi_tmp[INDEX(i,j,k)] != -INFINITY){
			if(i == I && j == J && k == K)
				printf("grid point (%d,%d,%d) in band \n", i, j, k);
			if(phi_tmp[INDEX(i,j,k)] == 0.f){
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				if(phi[INDEX(i,j,k)] > 0.f)
					pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				else
					neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			voxel.UpdateCellFastMarchType(i, j, k, SET, DONE);
		}
	END_FOR


	printf("Positive band size = %u, negtive band size = %u \n", pos_band.size(),neg_band.size());
//	if(phi_tmp) delete [] phi_tmp;


}

void ScalarField3D::InitializeInterface( Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band,
								const ScalarField3D &solid, float *value){
	float r;
	float delta = voxel.VoxelDelta();
//	float *phi_tmp = new float[DimX*DimY*DimZ];
	float *phi_tmp = phiNeg;
//	memset(phi_tmp, 0, DimX*DimY*DimZ);
//	SetZero(phi_tmp);
	SetInfinity(phi_tmp);

	// obtain the solid levelset
	float *phi_obj = solid.getScalar();

//	FOR_EACH_CELL
////		if(voxel.InSolid(i,j,k))
//			phi_tmp[INDEX(i,j,k)] = INFINITY;
//		if(voxel.IsDone(i,j,k))
//			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, DONE);
//		if(voxel.IsClose(i,j,k))
//			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, CLOSE);
//	END_FOR


	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			if(fabs(value[INDEX(i,j,k)]) < 1.e-6f)
				phi_tmp[INDEX(i,j,k)] = value[INDEX(i,j,k)];
			else{
				TentativeGridPoint tp(value[INDEX(i,j,k)], i, j, k);
				vector<TentativeGridPoint> neighbors;
				tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
				Point p = voxel.VoxelCenterPosition(i,j,k);
				float *candidates = new float[6];
				for(int m=0; m<6; m++)
					candidates[m] = INFINITY;
//				Point *candidatesP = new Point[6];
				for(int m=0; m<6; m++){
					if(tp.NeighborExists(m, DimX, DimY, DimZ)){
						Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
															 neighbors[m].jj,
															 neighbors[m].kk);
						neighbors[m].value = phi[POS(neighbors[m])];
						if(voxel.InSolid(neighbors[m])){
							if(tp.value <= 0.f && voxel.InMovingSolid(neighbors[m])){
								if(tp.LeftNeighbor(neighbors[m])){
									r = delta * phi_obj[POS(neighbors[m])] /
									       (phi_obj[POS(neighbors[m])] - phi_obj[INDEX(i,j,k)]);
									candidates[m] = min(fabsf(tp.value), (delta - r));
								}
								else if(tp.RightNeighbor(neighbors[m])){
									r = delta * phi_obj[INDEX(i,j,k)] /
										(phi_obj[INDEX(i,j,k)] - phi_obj[POS(neighbors[m])]);
									candidates[m] = min(fabsf(tp.value), r);
								}
								else if(tp.BackNeighbor(neighbors[m])){
									r = delta * phi_obj[POS(neighbors[m])] /
									       (phi_obj[POS(neighbors[m])] - phi_obj[INDEX(i,j,k)]);
									candidates[m] = min(fabsf(tp.value), (delta - r));
								}
								else if(tp.FrontNeighbor(neighbors[m])){
									r = delta * phi_obj[INDEX(i,j,k)] /
										(phi_obj[INDEX(i,j,k)] - phi_obj[POS(neighbors[m])]);
									candidates[m] = min(fabsf(tp.value), r);
								}
								else if(tp.BottomNeighbor(neighbors[m])){
									r = delta * phi_obj[POS(neighbors[m])] /
									       (phi_obj[POS(neighbors[m])] - phi_obj[INDEX(i,j,k)]);
									candidates[m] = min(fabsf(tp.value), (delta - r));
								}
								else if(tp.TopNeighbor(neighbors[m])){
									r = delta * phi_obj[INDEX(i,j,k)] /
										(phi_obj[INDEX(i,j,k)] - phi_obj[POS(neighbors[m])]);
									candidates[m] = min(fabsf(tp.value), r);
								}

							}
							else
								continue;
						}
//						else if(tp.value*neighbors[m].value == 0.f){
//							candidates[m] = delta;
//							candidatesP[m] = pn;
//						}
//						if(voxel.InSolid(neighbors[m]))
//							continue;
						else if(tp.value*neighbors[m].value <= 0.f){

							if(tp.LeftNeighbor(neighbors[m])){
								r = (tp.value * pn.x - neighbors[m].value * p.x) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.x - r;
//								candidatesP[m] = Point(r, p.y, p.z);
							}
							else if(tp.RightNeighbor(neighbors[m])){
								r = (tp.value * pn.x - neighbors[m].value * p.x) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.x;
//								candidatesP[m] = Point(r, p.y, p.z);
							}
							else if(tp.BackNeighbor(neighbors[m])){
								r = (tp.value * pn.y - neighbors[m].value * p.y) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.y - r;
//								candidatesP[m] = Point(p.x, r, p.z);
							}
							else if(tp.FrontNeighbor(neighbors[m])){
								r = (tp.value * pn.y - neighbors[m].value * p.y) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.y;
//								candidatesP[m] = Point(p.x, r, p.z);
							}
							else if(tp.BottomNeighbor(neighbors[m])){
								r = (tp.value * pn.z - neighbors[m].value * p.z) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.z - r;
//								candidatesP[m] = Point(p.x, p.y, r);
							}
							else if(tp.TopNeighbor(neighbors[m])){
								r = (tp.value * pn.z - neighbors[m].value * p.z) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.z;
//								candidatesP[m] = Point(p.x, p.y, r);
							}
							if(candidates[m] < 0.f){
								printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
										i,j,k, tp.value);
								printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
										candidates[m], r, pn.x, pn.y, pn.z);
								candidates[m] = fabs(candidates[m]);
							}
	//						if(i == I && j == J && k == K){
	//							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//									i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//												neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//												neighbors[m].value, phi[POS(neighbors[m])]);
	//						}
						}
	//					if(i == I && j == J && k == K){
	//						printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//							i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//						printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//								neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//								neighbors[m].value, phi[POS(neighbors[m])]);
	//
	//					}
					}
				}
				float phi_min[3], phimin = 0.f;
//				Point phi_min_p[3];
				phi_min[0]= min(candidates[0], candidates[1]);
//				phi_min_p[0] = (candidates[0] < candidates[1]) ? candidatesP[0]: candidatesP[1];
				phi_min[1]= min(candidates[2], candidates[3]);
//				phi_min_p[1] = (candidates[2] < candidates[3]) ? candidatesP[2]: candidatesP[3];
				phi_min[2]= min(candidates[4], candidates[5]);
//				phi_min_p[2] = (candidates[4] < candidates[5]) ? candidatesP[4]: candidatesP[5];
				if( phi_min[0] != INFINITY &&
					phi_min[1] != INFINITY &&
					phi_min[2] != INFINITY ){  // 3 points forming a plane
//					Vector v1 = phi_min_p[0] - phi_min_p[1];
//					Vector v2 = phi_min_p[0] - phi_min_p[2];
//					Vector v3 =  Cross(v1, v2);
//					if(v3.Length() == 0.f){
//						printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
//								i,j,k, tp.value);
//						for(int n=0; n < neighbors.size(); ++n)
//							printf("%f ", candidates[n]);
//						printf("\n");
//						for(int n=0; n < neighbors.size(); ++n)
//							printf("%f ", neighbors[n].value);
//						printf("\n");
//						printf(" p0 = (%f, %f, %f), p1 = (%f, %f, %f), p2 = (%f, %f, %f)\n",
//								phi_min_p[0].x, phi_min_p[0].y, phi_min_p[0].z,
//								phi_min_p[1].x, phi_min_p[1].y, phi_min_p[1].z,
//								phi_min_p[2].x, phi_min_p[2].y, phi_min_p[2].z);
//						phimin = fabs(tp.value);
////						exit(1);
//					}
//					else{
//						v3 = Normalize(v3);
//						Vector v4 = p - phi_min_p[0];
//						phimin = AbsDot(v3, v4);
//					}
					phimin = phi_min[0]*phi_min[1]*phi_min[2] /
					       sqrt( phi_min[0]*phi_min[0]*phi_min[1]*phi_min[1] +
					    		 phi_min[0]*phi_min[0]*phi_min[2]*phi_min[2] +
					    		 phi_min[2]*phi_min[2]*phi_min[1]*phi_min[1]);
				}
				else if( phi_min[0] != INFINITY && phi_min[1] != INFINITY ){
					// 2 points forming a line segment
//					phimin = DistanceToLine(phi_min_p[0], phi_min_p[1], p, tp.value);
					phimin = TwoPointDistance(phi_min[0], phi_min[1]);
				}
				else if( phi_min[1] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
//					phimin = DistanceToLine(phi_min_p[1], phi_min_p[2], p, tp.value);
					phimin = TwoPointDistance(phi_min[1], phi_min[2]);
				}
				else if( phi_min[0] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
//					phimin = DistanceToLine(phi_min_p[0], phi_min_p[2], p, tp.value);
					phimin = TwoPointDistance(phi_min[0], phi_min[2]);
				}
				else if(phi_min[0] != INFINITY)
					phimin = phi_min[0];
				else if(phi_min[1] != INFINITY)
					phimin = phi_min[1];
				else if(phi_min[2] != INFINITY)
					phimin = phi_min[2];
				else{
//					printf("wrong! at least 1 point should be valid at [%d, %d, %d] with phi = %f \n",
//						i,j,k, tp.value);
//					for(int n=0; n < neighbors.size(); ++n)
//						printf("%f ", candidates[n]);
//					printf("\n");
//					for(int n=0; n < neighbors.size(); ++n)
//						printf("%f ", neighbors[n].value);
//					printf("\n");
//					exit(1);
					phimin = INFINITY;
				}
		//		if(phi_min < 0.f){
		//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
		//					i,j,k, tp.value);
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", candidates[n]);
		//			printf("\n");
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", neighbors[n].value);
		//			printf("\n");
		//		}
				//if(phi_min == INFINITY)
				//	continue;
		//		printf("before i = %d, j= %d, k = %d, phi = %f \n",
		//						i, j, k, phi[INDEX(i,j,k)]);
				//KnownPoint tmp(i, j, k);
//				if(phimin != INFINITY && !voxel.CloseToSolid(tp)){
//					float lossdist = LossasoDistance(voxel, p, phimin);
//					phimin = min(phimin, lossdist);
//				}
				if(tp.value > 0.f){
					phi_tmp[INDEX(i,j,k)] = phimin;
					//if(phi_min != INFINITY)
						//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				else{
//					if(voxel.CloseToSolid(tp))
//						phimin = min(phimin, -phi_obj[INDEX(i,j,k)]);
					phi_tmp[INDEX(i,j,k)] = -phimin;
		//			if(phi_min != INFINITY)
		//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				if(candidates) delete [] candidates;
//				if(candidatesP) delete [] candidatesP;
		//		printf("after i = %d, j= %d, k = %d, phi = %f \n",
		//				i, j, k, phi[INDEX(i,j,k)]);
			}
		}
//		else{
////			phi_tmp[INDEX(i,j,k)] = INFINITY;
//			if(phi[INDEX(i,j,k)] == 0.f)
//				phi_tmp[INDEX(i,j,k)] = 0.f;
//			else{
//				TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
//				vector<TentativeGridPoint> neighbors;
//				tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
//				Point p = voxel.VoxelCenterPosition(i,j,k);
//				float *candidates = new float[6];
//				for(int m=0; m<6; m++)
//					candidates[m] = INFINITY;
//				Point *candidatesP = new Point[6];
//				for(int m=0; m<6; m++){
//					if(tp.NeighborExists(m, DimX, DimY, DimZ)){
//						Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
//															 neighbors[m].jj,
//															 neighbors[m].kk);
//						neighbors[m].value = phi[POS(neighbors[m])];
//						if(neighbors[m].value <= -0.5f*delta)
//							candidates[m] = 0.5f*delta;
//						else if(tp.value*neighbors[m].value < 0.f){
//
//							if(tp.LeftNeighbor(neighbors[m])){
//								r = (tp.value * pn.x - neighbors[m].value * p.x) /
//								       (tp.value - neighbors[m].value);
//								candidates[m] = p.x - r;
//								candidatesP[m] = Point(r, p.y, p.z);
//							}
//							else if(tp.RightNeighbor(neighbors[m])){
//								r = (tp.value * pn.x - neighbors[m].value * p.x) /
//								       (tp.value - neighbors[m].value);
//								candidates[m] = r - p.x;
//								candidatesP[m] = Point(r, p.y, p.z);
//							}
//							else if(tp.BackNeighbor(neighbors[m])){
//								r = (tp.value * pn.y - neighbors[m].value * p.y) /
//								       (tp.value - neighbors[m].value);
//								candidates[m] = p.y - r;
//								candidatesP[m] = Point(p.x, r, p.z);
//							}
//							else if(tp.FrontNeighbor(neighbors[m])){
//								r = (tp.value * pn.y - neighbors[m].value * p.y) /
//								       (tp.value - neighbors[m].value);
//								candidates[m] = r - p.y;
//								candidatesP[m] = Point(p.x, r, p.z);
//							}
//							else if(tp.BottomNeighbor(neighbors[m])){
//								r = (tp.value * pn.z - neighbors[m].value * p.z) /
//								       (tp.value - neighbors[m].value);
//								candidates[m] = p.z - r;
//								candidatesP[m] = Point(p.x, p.y, r);
//							}
//							else if(tp.TopNeighbor(neighbors[m])){
//								r = (tp.value * pn.z - neighbors[m].value * p.z) /
//								       (tp.value - neighbors[m].value);
//								candidates[m] = r - p.z;
//								candidatesP[m] = Point(p.x, p.y, r);
//							}
//							if(candidates[m] < 0.f){
//								printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
//										i,j,k, tp.value);
//								printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
//										candidates[m], r, pn.x, pn.y, pn.z);
//								candidates[m] = fabs(candidates[m]);
//							}
//	//						if(i == I && j == J && k == K){
//	//							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
//	//									i, j, k, tp.value, neighbors[m].value, candidates[m]);
//	//							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
//	//												neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
//	//												neighbors[m].value, phi[POS(neighbors[m])]);
//	//						}
//						}
//	//					if(i == I && j == J && k == K){
//	//						printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
//	//							i, j, k, tp.value, neighbors[m].value, candidates[m]);
//	//						printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
//	//								neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
//	//								neighbors[m].value, phi[POS(neighbors[m])]);
//	//
//	//					}
//					}
//				}
//				float phi_min[3], phimin = 0.f;
//				Point phi_min_p[3];
//				phi_min[0]= min(candidates[0], candidates[1]);
//				phi_min_p[0] = (candidates[0] < candidates[1]) ? candidatesP[0]: candidatesP[1];
//				phi_min[1]= min(candidates[2], candidates[3]);
//				phi_min_p[1] = (candidates[2] < candidates[3]) ? candidatesP[2]: candidatesP[3];
//				phi_min[2]= min(candidates[4], candidates[5]);
//				phi_min_p[2] = (candidates[4] < candidates[5]) ? candidatesP[4]: candidatesP[5];
//				if( phi_min[0] != INFINITY &&
//					phi_min[1] != INFINITY &&
//					phi_min[2] != INFINITY ){  // 3 points forming a plane
////					Vector v1 = phi_min_p[0] - phi_min_p[1];
////					Vector v2 = phi_min_p[0] - phi_min_p[2];
////					Vector v3 =  Cross(v1, v2);
////					if(v3.Length() == 0.f){
////						printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
////								i,j,k, tp.value);
////						for(int n=0; n < neighbors.size(); ++n)
////							printf("%f ", candidates[n]);
////						printf("\n");
////						for(int n=0; n < neighbors.size(); ++n)
////							printf("%f ", neighbors[n].value);
////						printf("\n");
////						printf(" p0 = (%f, %f, %f), p1 = (%f, %f, %f), p2 = (%f, %f, %f)\n",
////								phi_min_p[0].x, phi_min_p[0].y, phi_min_p[0].z,
////								phi_min_p[1].x, phi_min_p[1].y, phi_min_p[1].z,
////								phi_min_p[2].x, phi_min_p[2].y, phi_min_p[2].z);
////						phimin = fabs(tp.value);
//////						exit(1);
////					}
////					else{
////						v3 = Normalize(v3);
////						Vector v4 = p - phi_min_p[0];
////						phimin = AbsDot(v3, v4);
////					}
//					phimin = phi_min[0]*phi_min[1]*phi_min[2] /
//					       sqrt( phi_min[0]*phi_min[0]*phi_min[1]*phi_min[1] +
//					    		 phi_min[0]*phi_min[0]*phi_min[2]*phi_min[2] +
//					    		 phi_min[2]*phi_min[2]*phi_min[1]*phi_min[1]);
//				}
//				else if( phi_min[0] != INFINITY && phi_min[1] != INFINITY ){
//					// 2 points forming a line segment
////					phimin = DistanceToLine(phi_min_p[0], phi_min_p[1], p, tp.value);
//					phimin = TwoPointDistance(phi_min[0], phi_min[1]);
//				}
//				else if( phi_min[1] != INFINITY && phi_min[2] != INFINITY ){
//					// 2 points forming a line segment
////					phimin = DistanceToLine(phi_min_p[1], phi_min_p[2], p, tp.value);
//					phimin = TwoPointDistance(phi_min[1], phi_min[2]);
//				}
//				else if( phi_min[0] != INFINITY && phi_min[2] != INFINITY ){
//					// 2 points forming a line segment
////					phimin = DistanceToLine(phi_min_p[0], phi_min_p[2], p, tp.value);
//					phimin = TwoPointDistance(phi_min[0], phi_min[2]);
//				}
//				else if(phi_min[0] != INFINITY)
//					phimin = phi_min[0];
//				else if(phi_min[1] != INFINITY)
//					phimin = phi_min[1];
//				else if(phi_min[2] != INFINITY)
//					phimin = phi_min[2];
//				else{
////					printf("wrong! at least 1 point should be valid at [%d, %d, %d] with phi = %f \n",
////						i,j,k, tp.value);
////					for(int n=0; n < neighbors.size(); ++n)
////						printf("%f ", candidates[n]);
////					printf("\n");
////					for(int n=0; n < neighbors.size(); ++n)
////						printf("%f ", neighbors[n].value);
////					printf("\n");
////					exit(1);
//					phimin = INFINITY;
//				}
//		//		if(phi_min < 0.f){
//		//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
//		//					i,j,k, tp.value);
//		//			for(int n=0; n < neighbors.size(); ++n)
//		//				printf("%f ", candidates[n]);
//		//			printf("\n");
//		//			for(int n=0; n < neighbors.size(); ++n)
//		//				printf("%f ", neighbors[n].value);
//		//			printf("\n");
//		//		}
//				//if(phi_min == INFINITY)
//				//	continue;
//		//		printf("before i = %d, j= %d, k = %d, phi = %f \n",
//		//						i, j, k, phi[INDEX(i,j,k)]);
//				//KnownPoint tmp(i, j, k);
//				if(tp.value > 0.f){
//					phi_tmp[INDEX(i,j,k)] = phimin;
//					//if(phi_min != INFINITY)
//						//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
//				}
//				else{
//					phi_tmp[INDEX(i,j,k)] = -phimin;
//		//			if(phi_min != INFINITY)
//		//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
//				}
//				if(candidates) delete [] candidates;
//				if(candidatesP) delete [] candidatesP;
//		//		printf("after i = %d, j= %d, k = %d, phi = %f \n",
//		//				i, j, k, phi[INDEX(i,j,k)]);
//			}
//		}
	END_FOR

	int zero_points = 0;
	FOR_EACH_CELL
		KnownPoint tmp(i, j, k);
		if(phi_tmp[INDEX(i,j,k)] != INFINITY && phi_tmp[INDEX(i,j,k)] != -INFINITY){
			if(i == I && j == J && k == K)
				printf("grid point (%d,%d,%d) in band \n", i, j, k);
			if(phi_tmp[INDEX(i,j,k)] == 0.f){
				value[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
				++zero_points;
			}
			else{
				value[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				if(value[INDEX(i,j,k)] > 0.f)
					pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				else
					neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			voxel.UpdateCellFastMarchType(i, j, k, SET, DONE);
		}
	END_FOR

	printf("there are %d duplicated points in both pos and neg band \n",zero_points);
//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);
//	if(phi_tmp) delete [] phi_tmp;


}

void ScalarField3D::InitializeObject( Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band){
	float r;
//	float *phi_tmp = new float[DimX*DimY*DimZ];
	float *phi_tmp = phi0;
//	memset(phi_tmp, 0, DimX*DimY*DimZ);
	SetZero(phi_tmp);

	pos_band.clear();
	neg_band.clear();

//	FOR_EACH_CELL
//		if(voxel.IsDone(i,j,k))
//			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, DONE);
//		if(voxel.IsClose(i,j,k))
//			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, CLOSE);
//	END_FOR

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k))
			phi[INDEX(i,j,k)] *= -1;
	END_FOR


	FOR_EACH_CELL
			if(phi[INDEX(i,j,k)] == 0.f)
				phi_tmp[INDEX(i,j,k)] = 0.f;
			else{
				TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
				vector<TentativeGridPoint> neighbors;
				tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
				Point p = voxel.VoxelCenterPosition(i,j,k);
				float *candidates = new float[6];
				for(int m=0; m<6; m++)
					candidates[m] = INFINITY;
				Point *candidatesP = new Point[6];
				for(int m=0; m<6; m++){
					if(tp.NeighborExists(m, DimX, DimY, DimZ)){
						Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
															 neighbors[m].jj,
															 neighbors[m].kk);
						neighbors[m].value = phi[POS(neighbors[m])];
						if(tp.value*neighbors[m].value == 0.f){
							candidates[m] = fabs(tp.value);
							candidatesP[m] = pn;
						}
						else if(tp.value*neighbors[m].value < 0.f){

							if(tp.LeftNeighbor(neighbors[m])){
								r = (tp.value * pn.x - neighbors[m].value * p.x) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.x - r;
								candidatesP[m] = Point(r, p.y, p.z);
							}
							else if(tp.RightNeighbor(neighbors[m])){
								r = (tp.value * pn.x - neighbors[m].value * p.x) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.x;
								candidatesP[m] = Point(r, p.y, p.z);
							}
							else if(tp.BackNeighbor(neighbors[m])){
								r = (tp.value * pn.y - neighbors[m].value * p.y) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.y - r;
								candidatesP[m] = Point(p.x, r, p.z);
							}
							else if(tp.FrontNeighbor(neighbors[m])){
								r = (tp.value * pn.y - neighbors[m].value * p.y) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.y;
								candidatesP[m] = Point(p.x, r, p.z);
							}
							else if(tp.BottomNeighbor(neighbors[m])){
								r = (tp.value * pn.z - neighbors[m].value * p.z) /
								       (tp.value - neighbors[m].value);
								candidates[m] = p.z - r;
								candidatesP[m] = Point(p.x, p.y, r);
							}
							else if(tp.TopNeighbor(neighbors[m])){
								r = (tp.value * pn.z - neighbors[m].value * p.z) /
								       (tp.value - neighbors[m].value);
								candidates[m] = r - p.z;
								candidatesP[m] = Point(p.x, p.y, r);
							}
							if(candidates[m] < 0.f){
								printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
										i,j,k, tp.value);
								printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
										candidates[m], r, pn.x, pn.y, pn.z);
								candidates[m] = fabs(candidates[m]);
							}
	//						if(i == I && j == J && k == K){
	//							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//									i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//												neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//												neighbors[m].value, phi[POS(neighbors[m])]);
	//						}
						}
	//					if(i == I && j == J && k == K){
	//						printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//							i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//						printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//								neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//								neighbors[m].value, phi[POS(neighbors[m])]);
	//
	//					}
					}
				}
				float phi_min[3], phimin = 0.f;
				Point phi_min_p[3];
				phi_min[0]= min(candidates[0], candidates[1]);
				phi_min_p[0] = (candidates[0] < candidates[1]) ? candidatesP[0]: candidatesP[1];
				phi_min[1]= min(candidates[2], candidates[3]);
				phi_min_p[1] = (candidates[2] < candidates[3]) ? candidatesP[2]: candidatesP[3];
				phi_min[2]= min(candidates[4], candidates[5]);
				phi_min_p[2] = (candidates[4] < candidates[5]) ? candidatesP[4]: candidatesP[5];
				if( phi_min[0] != INFINITY &&
					phi_min[1] != INFINITY &&
					phi_min[2] != INFINITY ){  // 3 points forming a plane
					Vector v1 = phi_min_p[0] - phi_min_p[1];
					Vector v2 = phi_min_p[0] - phi_min_p[2];
					Vector v3 =  Cross(v1, v2);
					v3 = Normalize(v3);
					Vector v4 = p - phi_min_p[0];
					phimin = AbsDot(v3, v4);
				}
				else if( phi_min[0] != INFINITY && phi_min[1] != INFINITY ){
					// 2 points forming a line segment
					phimin = DistanceToLine(phi_min_p[0], phi_min_p[1], p, tp.value);
				}
				else if( phi_min[1] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
					phimin = DistanceToLine(phi_min_p[1], phi_min_p[2], p, tp.value);
				}
				else if( phi_min[0] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
					phimin = DistanceToLine(phi_min_p[0], phi_min_p[2], p, tp.value);
				}
				else if(phi_min[0] != INFINITY)
					phimin = phi_min[0];
				else if(phi_min[1] != INFINITY)
					phimin = phi_min[1];
				else if(phi_min[2] != INFINITY)
					phimin = phi_min[2];
				else{
//					printf("wrong! at least 1 point should be valid at [%d, %d, %d] with phi = %f \n",
//						i,j,k, tp.value);
//					for(int n=0; n < neighbors.size(); ++n)
//						printf("%f ", candidates[n]);
//					printf("\n");
//					for(int n=0; n < neighbors.size(); ++n)
//						printf("%f ", neighbors[n].value);
//					printf("\n");
//					exit(1);
					phimin = INFINITY;
				}
		//		if(phi_min < 0.f){
		//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
		//					i,j,k, tp.value);
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", candidates[n]);
		//			printf("\n");
		//			for(int n=0; n < neighbors.size(); ++n)
		//				printf("%f ", neighbors[n].value);
		//			printf("\n");
		//		}
				//if(phi_min == INFINITY)
				//	continue;
		//		printf("before i = %d, j= %d, k = %d, phi = %f \n",
		//						i, j, k, phi[INDEX(i,j,k)]);
				//KnownPoint tmp(i, j, k);
				if(tp.value > 0.f){
					phi_tmp[INDEX(i,j,k)] = phimin;
					//if(phi_min != INFINITY)
						//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				else{
					phi_tmp[INDEX(i,j,k)] = -phimin;
		//			if(phi_min != INFINITY)
		//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
				}
				if(candidates) delete [] candidates;
				if(candidatesP) delete [] candidatesP;
			}
		if(i == I && j == J && k == K)
					printf("(InitializeObject) point (%d,%d,%d) phi_tmp = %f \n",
							i, j, k, phi_tmp[INDEX(i,j,k)]);
	END_FOR

	FOR_EACH_CELL
		KnownPoint tmp(i, j, k);
		if(phi_tmp[INDEX(i,j,k)] != INFINITY && phi_tmp[INDEX(i,j,k)] != -INFINITY){
//			if(i == I && j == J && k == K)
//				printf("grid point (%d,%d,%d) in band \n", i, j, k);
			if(phi_tmp[INDEX(i,j,k)] == 0.f){
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				if(phi[INDEX(i,j,k)] > 0.f)
					pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				else
					neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
//			voxel.UpdateCellFastMarchType(i, j, k, SET, DONE);
		}
	END_FOR

	printf("InitializeObject Finished! \n");
//	if(phi_tmp) delete [] phi_tmp;


}

/*void ScalarField3D::InitializeObject( Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band){


	float r;
	float *phi_tmp = new float[DimX*DimY*DimZ];
	memset(phi_tmp, 0, DimX*DimY*DimZ);

	pos_band.clear();
	neg_band.clear();

	FOR_EACH_CELL
		if(voxel.IsDone(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, DONE);
		if(voxel.IsClose(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, CLOSE);
	END_FOR

// comment the follwoing lines for Zalesak's disk
// uncomment for water fall and pourting water test case

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k))
			phi[INDEX(i,j,k)] *= -1;
	END_FOR


	FOR_EACH_CELL
		if(phi[INDEX(i,j,k)] == 0.f)
			phi_tmp[INDEX(i,j,k)] = 0.f;
		else{
			TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
			vector<TentativeGridPoint> neighbors;
			tp.Neigbor(neighbors, DimX, DimY, DimZ);
			Point p = voxel.VoxelCenterPosition(i,j,k);
			float *candidates = new float[neighbors.size()];
			for(int m=0; m<neighbors.size(); m++){
				candidates[m] = INFINITY;
				Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
													 neighbors[m].jj,
													 neighbors[m].kk);
				neighbors[m].value = phi[POS(neighbors[m])];
				if(tp.value*neighbors[m].value == 0.f){
					candidates[m] = fabs(tp.value);
				}
				else if(tp.value*neighbors[m].value < 0.f){

					if(tp.LeftNeighbor(neighbors[m])){
						r = (tp.value * pn.x - neighbors[m].value * p.x) /
						       (tp.value - neighbors[m].value);
						candidates[m] = p.x - r;
					}
					else if(tp.RightNeighbor(neighbors[m])){
						r = (tp.value * pn.x - neighbors[m].value * p.x) /
						       (tp.value - neighbors[m].value);
						candidates[m] = r - p.x;
					}
					else if(tp.BackNeighbor(neighbors[m])){
						r = (tp.value * pn.y - neighbors[m].value * p.y) /
						       (tp.value - neighbors[m].value);
						candidates[m] = p.y - r;
					}
					else if(tp.FrontNeighbor(neighbors[m])){
						r = (tp.value * pn.y - neighbors[m].value * p.y) /
						       (tp.value - neighbors[m].value);
						candidates[m] = r - p.y;
					}
					else if(tp.BottomNeighbor(neighbors[m])){
						r = (tp.value * pn.z - neighbors[m].value * p.z) /
						       (tp.value - neighbors[m].value);
						candidates[m] = p.z - r;
					}
					else if(tp.TopNeighbor(neighbors[m])){
						r = (tp.value * pn.z - neighbors[m].value * p.z) /
						       (tp.value - neighbors[m].value);
						candidates[m] = r - p.z;
					}
					if(candidates[m] < 0.f){
						printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
								i,j,k, tp.value);
						printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
								candidates[m], r, pn.x, pn.y, pn.z);

					}
	//				printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//							i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//				printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//										neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//										neighbors[m].value, phi[POS(neighbors[m])]);
				}
			}
			float phi_min = FindMinimum(candidates, neighbors.size());
	//		if(phi_min < 0.f){
	//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
	//					i,j,k, tp.value);
	//			for(int n=0; n < neighbors.size(); ++n)
	//				printf("%f ", candidates[n]);
	//			printf("\n");
	//			for(int n=0; n < neighbors.size(); ++n)
	//				printf("%f ", neighbors[n].value);
	//			printf("\n");
	//		}
			//if(phi_min == INFINITY)
			//	continue;
	//		printf("before i = %d, j= %d, k = %d, phi = %f \n",
	//						i, j, k, phi[INDEX(i,j,k)]);
			//KnownPoint tmp(i, j, k);
			if(tp.value > 0.f){
				phi_tmp[INDEX(i,j,k)] = phi_min;
				//if(phi_min != INFINITY)
					//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				phi_tmp[INDEX(i,j,k)] = -phi_min;
	//			if(phi_min != INFINITY)
	//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			if(candidates) delete [] candidates;
	//		printf("after i = %d, j= %d, k = %d, phi = %f \n",
	//				i, j, k, phi[INDEX(i,j,k)]);
		}
	END_FOR

	FOR_EACH_CELL
		KnownPoint tmp(i, j, k);
		if(phi_tmp[INDEX(i,j,k)] != INFINITY && phi_tmp[INDEX(i,j,k)] != -INFINITY){
			if(phi_tmp[INDEX(i,j,k)] == 0.f){
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				if(phi[INDEX(i,j,k)] > 0.f)
					pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				else
					neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			voxel.UpdateCellFastMarchType(i, j, k, SET, DONE);
		}
	END_FOR

//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);
	if(phi_tmp) delete [] phi_tmp;


}*/

/*void ScalarField3D::InitializeInterface(Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band){
	float r;
	float *phi_tmp = new float[DimX*DimY*DimZ];
	memset(phi_tmp, 0, DimX*DimY*DimZ);

	FOR_EACH_CELL
		if(voxel.IsDone(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, DONE);
		if(voxel.IsClose(i,j,k))
			voxel.UpdateCellFastMarchType(i, j, k, CLEAR, CLOSE);
	END_FOR

	FOR_EACH_CELL
	if(!voxel.InSolid(i,j,k)){
		if(phi[INDEX(i,j,k)] == 0.f)
			phi_tmp[INDEX(i,j,k)] = 0.f;
		else{
			TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
			vector<TentativeGridPoint> neighbors;
			tp.Neigbor(neighbors, DimX, DimY, DimZ);
			Point p = voxel.VoxelCenterPosition(i,j,k);
			float *candidates = new float[neighbors.size()];
			for(int m=0; m<neighbors.size(); m++){
				candidates[m] = INFINITY;
				Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
													 neighbors[m].jj,
													 neighbors[m].kk);
				neighbors[m].value = phi[POS(neighbors[m])];
				if(voxel.InSolid(neighbors[m]))
					continue;
				else if(tp.value*neighbors[m].value == 0.f){
					candidates[m] = fabs(tp.value);
				}
				else if(!voxel.InSolid(neighbors[m]) &&
						tp.value*neighbors[m].value < 0.f){

					if(tp.LeftNeighbor(neighbors[m])){
						r = (tp.value * pn.x - neighbors[m].value * p.x) /
						       (tp.value - neighbors[m].value);
						candidates[m] = p.x - r;
					}
					else if(tp.RightNeighbor(neighbors[m])){
						r = (tp.value * pn.x - neighbors[m].value * p.x) /
						       (tp.value - neighbors[m].value);
						candidates[m] = r - p.x;
					}
					else if(tp.BackNeighbor(neighbors[m])){
						r = (tp.value * pn.y - neighbors[m].value * p.y) /
						       (tp.value - neighbors[m].value);
						candidates[m] = p.y - r;
					}
					else if(tp.FrontNeighbor(neighbors[m])){
						r = (tp.value * pn.y - neighbors[m].value * p.y) /
						       (tp.value - neighbors[m].value);
						candidates[m] = r - p.y;
					}
					else if(tp.BottomNeighbor(neighbors[m])){
						r = (tp.value * pn.z - neighbors[m].value * p.z) /
						       (tp.value - neighbors[m].value);
						candidates[m] = p.z - r;
					}
					else if(tp.TopNeighbor(neighbors[m])){
						r = (tp.value * pn.z - neighbors[m].value * p.z) /
						       (tp.value - neighbors[m].value);
						candidates[m] = r - p.z;
					}
					if(candidates[m] < 0.f){
						printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
								i,j,k, tp.value);
						printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
								candidates[m], r, pn.x, pn.y, pn.z);

					}
	//				printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//							i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//				printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//										neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//										neighbors[m].value, phi[POS(neighbors[m])]);
				}
			}
			float phi_min = FindMinimum(candidates, neighbors.size());
	//		if(phi_min < 0.f){
	//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
	//					i,j,k, tp.value);
	//			for(int n=0; n < neighbors.size(); ++n)
	//				printf("%f ", candidates[n]);
	//			printf("\n");
	//			for(int n=0; n < neighbors.size(); ++n)
	//				printf("%f ", neighbors[n].value);
	//			printf("\n");
	//		}
			//if(phi_min == INFINITY)
			//	continue;
	//		printf("before i = %d, j= %d, k = %d, phi = %f \n",
	//						i, j, k, phi[INDEX(i,j,k)]);
			//KnownPoint tmp(i, j, k);
			if(tp.value > 0.f){
				phi_tmp[INDEX(i,j,k)] = phi_min;
				//if(phi_min != INFINITY)
					//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				phi_tmp[INDEX(i,j,k)] = -phi_min;
	//			if(phi_min != INFINITY)
	//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			if(candidates) delete [] candidates;
	//		printf("after i = %d, j= %d, k = %d, phi = %f \n",
	//				i, j, k, phi[INDEX(i,j,k)]);
		}
	}
	else
		phi_tmp[INDEX(i,j,k)] = INFINITY;
	END_FOR

	FOR_EACH_CELL
		KnownPoint tmp(i, j, k);
		if(phi_tmp[INDEX(i,j,k)] != INFINITY && phi_tmp[INDEX(i,j,k)] != -INFINITY){
			if(phi_tmp[INDEX(i,j,k)] == 0.f){
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				phi[INDEX(i,j,k)] = phi_tmp[INDEX(i,j,k)];
				if(phi[INDEX(i,j,k)] > 0.f)
					pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				else
					neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			voxel.UpdateCellFastMarchType(i, j, k, SET, DONE);
		}
	END_FOR

//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);
	if(phi_tmp) delete [] phi_tmp;


}*/

void ScalarField3D::FastMarching(map<u_int, KnownPoint> &posband,
 								 map<u_int, KnownPoint> &negband,
								 Voxel &v, bool initial){
	// initialize grid points within the band of the interface
	// already initialized?
	float delta = v.VoxelDelta();

	int num_points=0;
	MinHeap<TentativeGridPoint, float> heap;
	map<u_int, KnownPoint> posadj;
	map<u_int, KnownPoint> negadj;
	map<u_int, KnownPoint>::iterator pos;


	//Process grid points at phi > 0
	for(pos=posband.begin(); pos!=posband.end(); ++pos){
		KnownPoint &p = pos->second;
//			u_int index = pos->first;
		TentativeGridPoint tp(0.f, p);
		vector<TentativeGridPoint> neigbor;
	    tp.Neigbor(neigbor, DimX, DimY, DimZ);
	    	map<u_int, KnownPoint>::iterator foundpos;
	    	map<u_int, KnownPoint>::iterator foundpos1;
		for(int i=0;i<neigbor.size();i++){
			foundpos1 = posband.find(Index(neigbor[i]));
			if(foundpos1==posband.end()){ // if this neighbor is not in accepted band
//			if(!v.IsDone(neigbor[i])){
				foundpos = posadj.find(Index(neigbor[i]));
				if(foundpos==posadj.end()){ // if this neighbor is not in close band
//				if(!v.IsClose(neigbor[i])){
					if(phi[POS(neigbor[i])] > 0.f){
						KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
						posadj.insert(make_pair(Index(neigbor[i]),tmp)); // add it to close band
						v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
					}
				}
			}
		}
//			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
//			i++;
	}
	printf("pos band size = %ld, adj_points = %ld \n", posband.size(), posadj.size());

	for(pos=posadj.begin(); pos!=posadj.end(); ++pos){
			KnownPoint &p = pos->second;
			u_int index = pos->first;
			float tentativeValue = UpdateTentativeValue(p, index, posband, v, phi);
//			float tentativeValue = UpdateTentativeValue(p, index, v);
//			printf("pos band phi = %f \n", tentativeValue);
			TentativeGridPoint tp(tentativeValue, p);
			heap.insert(tp, tentativeValue);
	}
	float minValue;
	do{
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		++num_points;
//		printf("i = %d \n", i);
		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
		phi[POS(t)] = minValue;
//		if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K)
//			printf("\n point(%d, %d, %d) phi = %f \n\n",
//					minTentative.ii,minTentative.jj,minTentative.kk, phi[POS(t)]);
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
		posband.insert( make_pair(Index(minTentative), t) );
		posadj.erase(Index(minTentative));
		if(!v.IsDone(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, SET, DONE);
		if(v.IsClose(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, CLEAR, CLOSE);
		if(!initial)
			if(minValue > FASTMARCH_LIMIT)
				break;
//		printf("heap size = %d, pos_adj size = %d\n", heap.size(), posadj.size());
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
	    vector<TentativeGridPoint> neigbor;
	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
	    for(int i=0;i<neigbor.size();i++){
	    	map<u_int, KnownPoint>::iterator foundinband = posband.find(Index(neigbor[i]));
	    	if(foundinband == posband.end()){ // not currently in the band
//	    	if(!v.IsDone(neigbor[i])){
	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), posband, v, phi);
//	    		if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
//	    			printf("\nneigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
//	    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    					minTentative.ii, minTentative.jj, minTentative.kk );
//	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), v);
//	    		if(neigbor[i].value == INFINITY){
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
//	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    		}
		    	map<u_int, KnownPoint>::iterator found = posadj.find(Index(neigbor[i]));
	    		if(found == posadj.end()){ // not currently in the heap
//	    		if(!v.IsClose(neigbor[i])){ // not currently in the heap
	    			if(phi[Index(neigbor[i])] > 0.f){
			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
			    		posadj.insert(make_pair(Index(neigbor[i]), tmp));
		//		    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
			    		heap.insert(neigbor[i], neigbor[i].value);
			    		v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
//			    		if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
//			    			printf("\n add close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
//			    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
	    			}
	    			else{
//	    				printf("\nfound negative phi in postive band \n");
//	    				printf("neigbor i = %d, value = %f, phi = %16.13f, id = (%d, %d, %d) \n",
//	    						i, neigbor[i].value,  phi[Index(neigbor[i])],
//			    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//			    		printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n\n", i, minTentative.value,
//			    			minTentative.ii, minTentative.jj, minTentative.kk );
	    			}
	    		}
	    		else{											// already in heap,
	    			heap.update(neigbor[i], neigbor[i].value);	// just update values
													    		// according to newly added point
//	    			if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
//		    			printf("\n update close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
//		    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
//					printf("update neigbor i = %d, value = %f \n",
//	    							i, neigbor[i].value );
	    		}
	    	}
	    }
	    neigbor.clear();
//	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
	} while(!heap.empty());

	printf("pos band size = %ld, num_points=%d \n", posband.size(), num_points);

	// in case when marching stops at a certain distance, the heap is not empty
	// make it empty for negative band
	if(!heap.empty())
		heap.clear();

	//Process grid points at phi < 0
	for(pos=negband.begin(); pos!=negband.end(); ++pos){
			KnownPoint &p = pos->second;
//			u_int index = pos->first;
			TentativeGridPoint tp(0.f, p);
			vector<TentativeGridPoint> neigbor;
		    tp.Neigbor(neigbor, DimX, DimY, DimZ);
	    	map<u_int, KnownPoint>::iterator foundpos;
	    	map<u_int, KnownPoint>::iterator foundpos1;
			for(int i=0;i<neigbor.size();i++){
				foundpos1 = negband.find(Index(neigbor[i]));
				if(foundpos1==negband.end()){
//				if(!v.IsDone(neigbor[i])){
					foundpos = negadj.find(Index(neigbor[i]));
					if(foundpos==negadj.end()){
//					if(!v.IsClose(neigbor[i])){
//						if((neigbor[i].ii < x1 && neigbor[i].ii > x0 &&
//						    neigbor[i].jj < y1 && neigbor[i].jj > y0 &&
//					   	    neigbor[i].kk < z1 && neigbor[i].kk > z0)){
						if(phi[POS(neigbor[i])] < 0.f){
							KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
							negadj.insert(make_pair(Index(neigbor[i]),tmp));
							v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
						}
					}
				}
			}
//			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
//			i++;
	}
	num_points=0;
	printf("neg band size = %ld, num_points=%d, adjacent points = %d \n",
			negband.size(), num_points, negadj.size());
	if(negadj.size() == 0){
		return;
	}
	for(pos=negband.begin(); pos!=negband.end(); ++pos){
		KnownPoint &p = pos->second;
		//printf("i = %d, j = %d, k = %d, phi = %f\n", p.ii, p.jj, p.kk, phi[POS(p)]);
		phi[POS(p)] *= -1;
	}
	for(pos=negadj.begin(); pos!=negadj.end(); ++pos){
			KnownPoint &p = pos->second;
			//printf("i = %d, j = %d, k = %d, phi = %f\n", p.ii, p.jj, p.kk, phi[POS(p)]);
			u_int index = pos->first;
			float tatitiveValue = UpdateTentativeValue(p, index, negband, v, phi);
//			float tatitiveValue = UpdateTentativeValue(p, index, v);
			TentativeGridPoint tp(tatitiveValue, p);
			heap.insert(tp, tatitiveValue);
//			++i;
//			printf("i = %d, index = %d \n", i, index);
	}
	do{
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		++num_points;
//		printf("i = %d \n", i);
		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
		phi[POS(t)] = minValue;
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
	    negband.insert( make_pair(Index(minTentative), t) );
		negadj.erase(Index(minTentative));
		if(!v.IsDone(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, SET, DONE);
		if(v.IsClose(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, CLEAR, CLOSE);
		if(!initial)
			if(minValue > FASTMARCH_LIMIT)
				break;
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
	    vector<TentativeGridPoint> neigbor;
	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
	    for(int i=0;i<neigbor.size();i++){
	    	map<u_int, KnownPoint>::iterator foundinband = negband.find(Index(neigbor[i]));
	    	if(foundinband == negband.end()){ // not currently in the band
//	    	if(!v.IsDone(neigbor[i])){
	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), negband, v, phi);
//	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), v);
//	    		if(neigbor[i].value == INFINITY){
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
//	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    		}
		    	map<u_int, KnownPoint>::iterator found = negadj.find(Index(neigbor[i]));
	    		if(found == negadj.end()){ // not currently in the heap
//	    		if(!v.IsClose(neigbor[i])){ // not currently in the heap
	    			if(phi[Index(neigbor[i])] < 0.f){
			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
			    		negadj.insert(make_pair(Index(neigbor[i]), tmp));
		//	    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
			    		heap.insert(neigbor[i], neigbor[i].value);
			    		v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
	    			}
	    			else{
//	    				printf("\n found positive phi in negative band \n");
//	    				printf("neigbor i = %d, value = %f, phi = %16.13f, id = (%d, %d, %d) \n",
//	    						i, neigbor[i].value, phi[Index(neigbor[i])],
//			    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//			    		printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n\n", i, minTentative.value,
//			    			minTentative.ii, minTentative.jj, minTentative.kk );
	    			}
	    		}
	    		else{											// already in heap,
//	    		printf("update neigbor i = %d, value = %f \n", i, neigbor[i].value );
	    			heap.update(neigbor[i], neigbor[i].value);	// just update values
													    		// according to newly added point
	    		}
	    	}
	    }
	    neigbor.clear();
//	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
	} while(!heap.empty());
	printf("neg band size = %ld, num_points=%d \n", negband.size(), num_points);

	// multiply by -1 to get negative distance in phi < 0 area
	for(pos=negband.begin(); pos!=negband.end(); ++pos){
		KnownPoint &p = pos->second;
		phi[POS(p)] *= -1;
	}

//	u_int air=0, liquid=0;
//	FOR_EACH_CELL
//		if(phi[INDEX(i,j,k)] < 0.f)
//			liquid++;
//		else
//			air++;
//	END_FOR

//	printf("D: phi = %f \n", phi[INDEX(I,J,K)]);
//	printf("after marching air points = %d, liquid points = %d \n", air, liquid);

}

void ScalarField3D::FastMarchingAirWater(map<u_int, KnownPoint> &posband,
										 map<u_int, KnownPoint> &negband,
										 Voxel &v, bool initial, float *phi_tmp){
	// initialize grid points within the band of the interface
	// already initialized?
	float delta = v.VoxelDelta();

	int num_points=0;
	MinHeap<TentativeGridPoint, float> heap;
	map<u_int, KnownPoint> posadj;
	map<u_int, KnownPoint> negadj;
	map<u_int, KnownPoint>::iterator pos;

	if(posband.size() == 0)
		return;

	//Process grid points at phi > 0
	for(pos=posband.begin(); pos!=posband.end(); ++pos){
		KnownPoint &p = pos->second;
//			u_int index = pos->first;
		TentativeGridPoint tp(0.f, p);
		vector<TentativeGridPoint> neigbor;
	    tp.Neigbor(neigbor, DimX, DimY, DimZ);
	    map<u_int, KnownPoint>::iterator foundpos;
	    map<u_int, KnownPoint>::iterator foundpos1;
		for(int i=0;i<neigbor.size();i++){
			foundpos1 = posband.find(Index(neigbor[i]));
			if(foundpos1==posband.end()){ // if this neighbor is not in accepted band
//			if(!v.IsDone(neigbor[i])){
				foundpos = posadj.find(Index(neigbor[i]));
				if(foundpos==posadj.end()){ // if this neighbor is not in close band
//				if(!v.IsClose(neigbor[i])){
					if(phi_tmp[POS(neigbor[i])] > 0.f && !v.InSolid(neigbor[i])){
						KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
						posadj.insert(make_pair(Index(neigbor[i]),tmp)); // add it to close band
						v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
					}
				}
			}
		}
//			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
//			i++;
	}
	printf("pos band size = %ld, adj_points = %ld \n", posband.size(), posadj.size());

	for(pos=posadj.begin(); pos!=posadj.end(); ++pos){
			KnownPoint &p = pos->second;
			u_int index = pos->first;
			float tentativeValue = UpdateTentativeValue(p, index, posband, v, phi_tmp);
//			float tentativeValue = UpdateTentativeValue(p, index, v);
//			printf("pos band phi = %f \n", tentativeValue);
			TentativeGridPoint tp(tentativeValue, p);
			heap.insert(tp, tentativeValue);
	}
	float minValue;
	do{
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		++num_points;
//		printf("i = %d \n", i);
		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
		phi_tmp[POS(t)] = minValue;
//		if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K)
//			printf("\n point(%d, %d, %d) phi = %f \n\n",
//					minTentative.ii,minTentative.jj,minTentative.kk, phi[POS(t)]);
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
		posband.insert( make_pair(Index(minTentative), t) );
		posadj.erase(Index(minTentative));
		if(!v.IsDone(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, SET, DONE);
		if(v.IsClose(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, CLEAR, CLOSE);
		if(!initial)
			if(minValue > FASTMARCH_LIMIT)
				break;
//		printf("heap size = %d, pos_adj size = %d\n", heap.size(), posadj.size());
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
	    vector<TentativeGridPoint> neigbor;
	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
	    for(int i=0;i<neigbor.size();i++){
	    	map<u_int, KnownPoint>::iterator foundinband = posband.find(Index(neigbor[i]));
	    	if(foundinband == posband.end()){ // not currently in the band
//	    	if(!v.IsDone(neigbor[i])){
	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), posband, v, phi_tmp);
//	    		if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
//	    			printf("\nneigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
//	    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    					minTentative.ii, minTentative.jj, minTentative.kk );
//	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), v);
//	    		if(neigbor[i].value == INFINITY){
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
//	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    		}
		    	map<u_int, KnownPoint>::iterator found = posadj.find(Index(neigbor[i]));
	    		if(found == posadj.end()){ // not currently in the heap
//	    		if(!v.IsClose(neigbor[i])){ // not currently in the heap
	    			if(phi_tmp[Index(neigbor[i])] > 0.f && !v.InSolid(neigbor[i])){
			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
			    		posadj.insert(make_pair(Index(neigbor[i]), tmp));
		//		    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
			    		heap.insert(neigbor[i], neigbor[i].value);
			    		v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
//			    		if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
//			    			printf("\n add close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
//			    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
	    			}
	    			else{
//	    				printf("\nfound negative phi in postive band \n");
//	    				printf("neigbor i = %d, value = %f, phi = %16.13f, id = (%d, %d, %d) \n",
//	    						i, neigbor[i].value,  phi[Index(neigbor[i])],
//			    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//			    		printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n\n", i, minTentative.value,
//			    			minTentative.ii, minTentative.jj, minTentative.kk );
	    			}
	    		}
	    		else{											// already in heap,
	    			heap.update(neigbor[i], neigbor[i].value);	// just update values
													    		// according to newly added point
//	    			if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
//		    			printf("\n update close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
//		    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
//					printf("update neigbor i = %d, value = %f \n",
//	    							i, neigbor[i].value );
	    		}
	    	}
	    }
	    neigbor.clear();
//	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
	} while(!heap.empty());

	printf("pos band size = %ld, num_points=%d \n", posband.size(), num_points);

	// in case when marching stops at a certain distance, the heap is not empty
	// make it empty for negative band
	if(!heap.empty())
		heap.clear();

	//Process grid points at phi < 0
	for(pos=negband.begin(); pos!=negband.end(); ++pos){
			KnownPoint &p = pos->second;
//			u_int index = pos->first;
			TentativeGridPoint tp(0.f, p);
			vector<TentativeGridPoint> neigbor;
		    tp.Neigbor(neigbor, DimX, DimY, DimZ);
	    	map<u_int, KnownPoint>::iterator foundpos;
	    	map<u_int, KnownPoint>::iterator foundpos1;
			for(int i=0;i<neigbor.size();i++){
				foundpos1 = negband.find(Index(neigbor[i]));
				if(foundpos1==negband.end()){
//				if(!v.IsDone(neigbor[i])){
					foundpos = negadj.find(Index(neigbor[i]));
					if(foundpos==negadj.end()){
//					if(!v.IsClose(neigbor[i])){
//						if((neigbor[i].ii < x1 && neigbor[i].ii > x0 &&
//						    neigbor[i].jj < y1 && neigbor[i].jj > y0 &&
//					   	    neigbor[i].kk < z1 && neigbor[i].kk > z0)){
						if(phi_tmp[POS(neigbor[i])] < 0.f && !v.InSolid(neigbor[i])){
							KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
							negadj.insert(make_pair(Index(neigbor[i]),tmp));
							v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
						}
					}
				}
			}
//			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
//			i++;
	}
	num_points=0;
	printf("neg band size = %ld, num_points=%d, adjacent points = %d \n",
			negband.size(), num_points, negadj.size());
	if(negadj.size() == 0){
		return;
	}
	for(pos=negband.begin(); pos!=negband.end(); ++pos){
		KnownPoint &p = pos->second;
		//printf("i = %d, j = %d, k = %d, phi = %f\n", p.ii, p.jj, p.kk, phi[POS(p)]);
		phi_tmp[POS(p)] *= -1;
	}
	for(pos=negadj.begin(); pos!=negadj.end(); ++pos){
			KnownPoint &p = pos->second;
			//printf("i = %d, j = %d, k = %d, phi = %f\n", p.ii, p.jj, p.kk, phi[POS(p)]);
			u_int index = pos->first;
			float tatitiveValue = UpdateTentativeValue(p, index, negband, v, phi_tmp);
//			float tatitiveValue = UpdateTentativeValue(p, index, v);
			TentativeGridPoint tp(tatitiveValue, p);
			heap.insert(tp, tatitiveValue);
//			++i;
//			printf("i = %d, index = %d \n", i, index);
	}
	do{
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		++num_points;
//		printf("i = %d \n", i);
		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
		phi_tmp[POS(t)] = minValue;
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
	    negband.insert( make_pair(Index(minTentative), t) );
		negadj.erase(Index(minTentative));
		if(!v.IsDone(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, SET, DONE);
		if(v.IsClose(minTentative))
			v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, CLEAR, CLOSE);
		if(!initial)
			if(minValue > FASTMARCH_LIMIT)
				break;
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
	    vector<TentativeGridPoint> neigbor;
	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
	    for(int i=0;i<neigbor.size();i++){
	    	map<u_int, KnownPoint>::iterator foundinband = negband.find(Index(neigbor[i]));
	    	if(foundinband == negband.end()){ // not currently in the band
//	    	if(!v.IsDone(neigbor[i])){
	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), negband, v, phi_tmp);
//	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), v);
//	    		if(neigbor[i].value == INFINITY){
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
//	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    		}
		    	map<u_int, KnownPoint>::iterator found = negadj.find(Index(neigbor[i]));
	    		if(found == negadj.end()){ // not currently in the heap
//	    		if(!v.IsClose(neigbor[i])){ // not currently in the heap
	    			if(phi_tmp[Index(neigbor[i])] < 0.f && !v.InSolid(neigbor[i])){
			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
			    		negadj.insert(make_pair(Index(neigbor[i]), tmp));
		//	    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
			    		heap.insert(neigbor[i], neigbor[i].value);
			    		v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
	    			}
	    			else{
//	    				printf("\n found positive phi in negative band \n");
//	    				printf("neigbor i = %d, value = %f, phi = %16.13f, id = (%d, %d, %d) \n",
//	    						i, neigbor[i].value, phi[Index(neigbor[i])],
//			    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//			    		printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n\n", i, minTentative.value,
//			    			minTentative.ii, minTentative.jj, minTentative.kk );
	    			}
	    		}
	    		else{											// already in heap,
//	    		printf("update neigbor i = %d, value = %f \n", i, neigbor[i].value );
	    			heap.update(neigbor[i], neigbor[i].value);	// just update values
													    		// according to newly added point
	    		}
	    	}
	    }
	    neigbor.clear();
//	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
	} while(!heap.empty());
	printf("neg band size = %ld, num_points=%d \n", negband.size(), num_points);

	// multiply by -1 to get negative distance in phi < 0 area
	for(pos=negband.begin(); pos!=negband.end(); ++pos){
		KnownPoint &p = pos->second;
		phi_tmp[POS(p)] *= -1;
	}

//	u_int air=0, liquid=0;
//	FOR_EACH_CELL
//		if(phi[INDEX(i,j,k)] < 0.f)
//			liquid++;
//		else
//			air++;
//	END_FOR

//	printf("D: phi = %f \n", phi[INDEX(I,J,K)]);
//	printf("after marching air points = %d, liquid points = %d \n", air, liquid);

}

void ScalarField3D::print_grad_phi(int I, int J, int K) const {
	float grad_x = 0.5f*(phi[INDEX(I+1,J,K)] - phi[INDEX(I-1,J,K)]);
	float grad_y = 0.5f*(phi[INDEX(I,J+1,K)] - phi[INDEX(I,J-1,K)]);
	float grad_z = 0.5f*(phi[INDEX(I,J,K+1)] - phi[INDEX(I,J,K-1)]);
	float grad = sqrtf(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z);
	printf("after marching, grad(phi) at i = %d, j = %d, k = %d is %f \n",
			I,J,K, grad);
}

void ScalarField3D::UpdateVoxel(Voxel &voxel) {

	float cell_phi[4];

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			cell_phi[0] = phi[INDEX(i,j,k)];
			cell_phi[1] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i+1,j,k)]);
			cell_phi[2] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j+1,k)]);
			cell_phi[3] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j,k+1)]);
//			cell_phi[1] = 0.5f*(phi[INDEX(i-1,j,k)]+phi[INDEX(i,j,k)]);
//			cell_phi[2] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i+1,j,k)]);
//			cell_phi[3] = 0.5f*(phi[INDEX(i,j-1,k)]+phi[INDEX(i,j,k)]);
//			cell_phi[4] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j+1,k)]);
//			cell_phi[5] = 0.5f*(phi[INDEX(i,j,k-1)]+phi[INDEX(i,j,k)]);
//			cell_phi[6] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j,k+1)]);
			int negative = 0;
			for(int m=0; m<4; ++m){
				if(cell_phi[m] < 0.f)
					negative++;
			}
			if(negative == 4){ // liquid cell
				if(voxel.InAir(i,j,k))
					voxel.UpdateCellType(i, j, k, CLEAR, EMPTY);
				if(voxel.InSurface(i,j,k))
					voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
				if(!voxel.InLiquid(i,j,k))
					voxel.UpdateCellType(i, j, k, SET, LIQUID);
			}
			else if(negative == 0){  // air cell
				if(voxel.InLiquid(i,j,k))
					voxel.UpdateCellType(i, j, k, CLEAR, LIQUID);
				if(voxel.InSurface(i,j,k))
					voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
				if(!voxel.InAir(i,j,k))
					voxel.UpdateCellType(i, j, k, SET, EMPTY);
			}
			else{
				if(cell_phi[0] < 0.f){   // surface cell
					if(voxel.InAir(i,j,k))
						voxel.UpdateCellType(i, j, k, CLEAR, EMPTY);
					if(!voxel.InSurface(i,j,k))
						voxel.UpdateCellType(i, j, k, SET, SURFACE);
					if(!voxel.InLiquid(i,j,k))
						voxel.UpdateCellType(i, j, k, SET, LIQUID);
				}
				else{   // air cell
					if(!voxel.InAir(i,j,k))
						voxel.UpdateCellType(i, j, k, SET, EMPTY);
					if(voxel.InSurface(i,j,k))
						voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
					if(voxel.InLiquid(i,j,k))
						voxel.UpdateCellType(i, j, k, CLEAR, LIQUID);
				}
			}
			 if(i==I && j==J && k==K){
				if(voxel.InLiquid(I,J,K))
				  printf(" \n (%d, %d, %d) in liquid \n\n", I, J, K);
			    if(voxel.InAir(I,J,K))
				  printf(" \n (%d, %d, %d) in air \n\n", I, J, K);
			    if(voxel.InSurface(I,J,K))
				  printf("\n  (%d, %d, %d) in surface \n\n", I, J, K);

		    	printf(" \n negative = %d, cell phi = (%f, %f, %f, %f) \n\n",
		    		negative, cell_phi[0], cell_phi[1], cell_phi[2], cell_phi[3]);
			 }
		}

	END_FOR

	FOR_EACH_CELL
		if(voxel.InSurface(i,j,k)){
			bool realSurfaceCell = false;
			TentativeGridPoint p(0.f,i,j,k);
			vector<TentativeGridPoint> neighbors;
			p.Neigbor(neighbors, DimX, DimY, DimZ);
			for(int m=0; m<neighbors.size(); m++){
				if(voxel.InAir(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk))
					realSurfaceCell = true;
			}
			if(!realSurfaceCell)
				voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
		}
//		if(voxel.InAir(i,j,k)){
//			float face_phi[7];
//			face_phi[1] = 0.5f*(phi[INDEX(i-1,j,k)]+phi[INDEX(i,j,k)]);
//			face_phi[2] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i+1,j,k)]);
//			face_phi[3] = 0.5f*(phi[INDEX(i,j-1,k)]+phi[INDEX(i,j,k)]);
//			face_phi[4] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j+1,k)]);
//			face_phi[5] = 0.5f*(phi[INDEX(i,j,k-1)]+phi[INDEX(i,j,k)]);
//			face_phi[6] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j,k+1)]);
//			int negative = 0;
//			float tmp = 0;
//			for(int m=1; m<7; ++m){
//				if(face_phi[m] < 0.f){
//					negative++;
//					tmp += face_phi[m];
//				}
//			}
//			if(negative == 6){
//				phi[INDEX(i,j,k)] = tmp / 6.f;
//				if(voxel.InAir(i,j,k))
//					voxel.UpdateCellType(i, j, k, CLEAR, EMPTY);
//				if(voxel.InSurface(i,j,k))
//					voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
//				if(!voxel.InLiquid(i,j,k))
//					voxel.UpdateCellType(i, j, k, SET, LIQUID);
//			}
//		}
//		if(voxel.InSurface(i,j,k)){
//			float face_phi[7];
//			face_phi[1] = 0.5f*(phi[INDEX(i-1,j,k)]+phi[INDEX(i,j,k)]);
//			face_phi[2] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i+1,j,k)]);
//			face_phi[3] = 0.5f*(phi[INDEX(i,j-1,k)]+phi[INDEX(i,j,k)]);
//			face_phi[4] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j+1,k)]);
//			face_phi[5] = 0.5f*(phi[INDEX(i,j,k-1)]+phi[INDEX(i,j,k)]);
//			face_phi[6] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j,k+1)]);
//			int positive = 0;
//			float tmp = 0;
//			for(int m=1; m<7; ++m){
//				if(face_phi[m] > 0.f){
//					positive++;
//					tmp += face_phi[m];
//				}
//			}
//			if(positive == 6){
//				phi[INDEX(i,j,k)] = tmp / 6.f;
//				if(!voxel.InAir(i,j,k))
//					voxel.UpdateCellType(i, j, k, SET, EMPTY);
//				if(voxel.InSurface(i,j,k))
//					voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
//				if(voxel.InLiquid(i,j,k))
//					voxel.UpdateCellType(i, j, k, CLEAR, LIQUID);
//			}
//		}
	END_FOR



}

void ScalarField3D::CorrectNeighborLevelSet(const Voxel &voxel, const Point &pos,
		                                    const Particle &ps, float delta, char s){
	float x = pos.x/delta - 0.5f;
    float y = pos.y/delta - 0.5f;
    float z = pos.z/delta - 0.5f;
    if(x < 0) x = 0;
    if(y < 0) y = 0;
    if(z < 0) z = 0;
	int xi = floor(x);
	int yi = floor(y);
	int zi = floor(z);

	double r,t;

//	printf("\nxi = %d, yi = %d, zi = %d \n", xi, yi, zi);
//	printf(" particle at (%f, %f, %f) \n", pos.x, pos.y, pos.z);
	for(int i=xi; i<=xi+1; ++i){
		for(int j=yi; j<=yi+1; ++j){
			for(int k=zi; k<=zi+1; ++k){

				if(!IndexOutofBounds(i, j, k)){

				if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//				if(tp.ii == I && tp.jj == J && tp.kk == K){
//					for(int n=0; n < 8; ++n){
//						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
//						save_corner[n] = ps.EvaluatePhiP(corner);
//					}
//					printf(" negative particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
//							pos.x, pos.y, pos.z, phi_i, radius);
//					for(int n=0; n < 8; ++n)
//						printf(" %f", save_corner[n]);
//					printf("\n");
//				}
					Point center = voxel.VoxelCenterPosition(i, j, k);
#ifdef USE_MESH
					// we only use visible particles to correct phi
					if(!mPhysWorld->SegmentIntersectsMesh(center, pos, &r, &t)){
#endif
//					if(s < 0)
//						printf("at (%d, %d, %d), negPhi = %f\n",
//								i, j, k,phiNeg[INDEX(i,j,k)]);
//					else
//						printf("at (%d, %d, %d), posPhi = %f \n",
//								i, j, k,phiPos[INDEX(i,j,k)]);
						float phi_center;
						if(s < 0){
							phi_center = ps.EvaluatePhiP(center);
							phiNeg[INDEX(i,j,k)] = min(phi_center, phiNeg[INDEX(i,j,k)]);
						}
						else{
							phi_center = -ps.EvaluatePhiP(center);
							phiPos[INDEX(i,j,k)] = max(phi_center, phiPos[INDEX(i,j,k)]);
						}
#ifdef USE_MESH
					}
#endif
//					if(s < 0)
//						printf("at (%d, %d, %d), negPhi = %f, phip = %f \n",
//								i, j, k,phiNeg[INDEX(i,j,k)], phi_center);
//					else
//						printf("at (%d, %d, %d), posPhi = %f, phip = %f \n",
//								i, j, k,phiPos[INDEX(i,j,k)], phi_center);
				}
				}
			}
		}
	}
//	printf("\n\n");

//		float x = pos.x/delta;
//	    float y = pos.y/delta;
//	    float z = pos.z/delta;
//	    if(x < 0) x = 0;
//	    if(y < 0) y = 0;
//	    if(z < 0) z = 0;
//		int xi = floor(x);
//		int yi = floor(y);
//		int zi = floor(z);
//
////		printf("xi = %d, yi = %d, zi = %d \n", xi, yi, zi);
////		printf(" particle at (%f, %f, %f) \n", pos.x, pos.y, pos.z);
//		int i = xi;
//		int j = yi;
//		int k = zi;
//
//	//float phi_p = TriInterp(voxel, pos);
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
////				if(tp.ii == I && tp.jj == J && tp.kk == K){
////					for(int n=0; n < 8; ++n){
////						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
////						save_corner[n] = ps.EvaluatePhiP(corner);
////					}
////					printf(" negative particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
////							pos.x, pos.y, pos.z, phi_i, radius);
////					for(int n=0; n < 8; ++n)
////						printf(" %f", save_corner[n]);
////					printf("\n");
////				}
//			Point center = voxel.VoxelCenterPosition(i, j, k);
////			if(s < 0)
////				printf("at (%d, %d, %d), negPhi = %f\n",
////						i, j, k,phiNeg[INDEX(i,j,k)]);
////			else
////				printf("at (%d, %d, %d), posPhi = %f \n",
////						i, j, k,phiPos[INDEX(i,j,k)]);
//			float phi_center;
//			if(s < 0){
//				phi_center = ps.EvaluatePhiP(center);
//				phiNeg[INDEX(i,j,k)] = min(phi_center, phiNeg[INDEX(i,j,k)]);
//			}
//			else{
//				phi_center = -ps.EvaluatePhiP(center);
//				phiPos[INDEX(i,j,k)] = max(phi_center, phiPos[INDEX(i,j,k)]);
//			}
////			if(s < 0)
////				printf("at (%d, %d, %d), negPhi = %f, phip = %f \n",
////						i, j, k,phiNeg[INDEX(i,j,k)], phi_center);
////			else
////				printf("at (%d, %d, %d), posPhi = %f, phip = %f \n",
////						i, j, k,phiPos[INDEX(i,j,k)], phi_center);
//		}

//		printf("\n\n");
}


void ScalarField3D::ErrorCorrection(const Voxel &voxel,
							ParticleList &particles
//							list<Particle> &particles
#ifdef POSITIVE_PARTICLES
//							,list<Particle> &posParticles
							,ParticleList &posParticles
#endif
							){

	float delta = voxel.VoxelDelta();

//	list<Particle>::iterator iter_particle;
	ParticleList::iterator iter_particle;

//	float *x = new float[DimX*DimY*DimZ];
//	memset(x, 0, DimX*DimY*DimZ);
	SetEqual(phi, phi0);

	SetEqual(phi, phiPos);
	SetEqual(phi, phiNeg);

	u_int N = 0, NM = 0, NP = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f)
			++N;
	END_FOR

	printf("before error correction, there are %u liquid points \n", N);

	printf("start error correction ... \n");

	//for(int m=0; m<particles.size(); ++m){
	for(iter_particle  = particles.begin();
		iter_particle != particles.end();
		++iter_particle){
//		Point pos = particles[m].Position();
//		float radius = particles[m].Radius();
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
//		Point center = voxel.VoxelCenterPosition(tp.ii, tp.jj, tp.kk);
//		if(tp.ii == I && tp.jj == J && tp.kk == K){
//		if(tp.kk == K && pos.z < center.z && phi_i > 0.f){
//			printf(" negative particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
//					pos.x, pos.y, pos.z, phi_i, radius);
//			printf("at (%d,%d,%d), phi = %f, phi[i+1] = %f, phi[i-1] = %f \n",
//					tp.ii, tp.jj, tp.kk, phi[INDEX(tp.ii,tp.jj,tp.kk)],phi[INDEX(tp.ii+1,tp.jj,tp.kk)], phi[INDEX(tp.ii-1,tp.jj,tp.kk)]);
//			printf("at (%d,%d,%d), phi = %f, phi[j+1] = %f, phi[j-1] = %f \n",
//					tp.ii, tp.jj, tp.kk,phi[INDEX(tp.ii,tp.jj,tp.kk)],phi[INDEX(tp.ii,tp.jj+1,tp.kk)], phi[INDEX(tp.ii,tp.jj-1,tp.kk)]);
//		}
		if( !IndexOutofBounds(tp.ii, tp.jj, tp.kk) && phi_i > radius ){
//			printf("\nphi_i = %f, radius = %f \n", phi_i, radius);
			++NM;
			CorrectNeighborLevelSet(voxel, pos, ps, delta, -1);
//			Point center = voxel.VoxelCenterPosition(tp.ii, tp.jj, tp.kk);
//			float phi_center = ps.EvaluatePhiP(center);
//			phiNeg[INDEX(tp.ii,tp.jj,tp.kk)] = min(phi_center, phiNeg[INDEX(tp.ii,tp.jj,tp.kk)]);
		}

	}
	printf("%u negative particles participating in error correction  \n", NM);

#ifdef POSITIVE_PARTICLES
		for(iter_particle = posParticles.begin();
			iter_particle != posParticles.end();
			++iter_particle){
	//		Point pos = particles[m].Position();
	//		float radius = particles[m].Radius();
			Particle &ps = *iter_particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			float phi_i = TriInterp(voxel, pos);
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			if( !IndexOutofBounds(tp.ii, tp.jj, tp.kk) && -phi_i > radius ){
//				if(tp.ii == I && tp.jj == J && tp.kk == K){
//					printf(" positive particle in (%d, %d, %d) at (%f, %f, %f) phi_i = %f, radius = %f \n",
//								tp.ii, tp.jj, tp.kk, pos.x, pos.y, pos.z, phi_i, radius);
//				}
//				printf("\nphi_i = %f, radius = %f \n", phi_i, radius);
//				printf(" positive particle in (%d, %d, %d) at (%f, %f, %f) phi_i = %f, radius = %f \n",
//							tp.ii, tp.jj, tp.kk, pos.x, pos.y, pos.z, phi_i, radius);
				++NP;
				CorrectNeighborLevelSet(voxel, pos, ps, delta, 1);
//				Point center = voxel.VoxelCenterPosition(tp.ii, tp.jj, tp.kk);
//				float phi_center = -ps.EvaluatePhiP(center);
//				phiPos[INDEX(tp.ii,tp.jj,tp.kk)] = max(phi_center, phiPos[INDEX(tp.ii,tp.jj,tp.kk)]);
			}

		}
		printf("%u positive particles participating in error correction  \n", NP);
#endif

		FOR_EACH_CELL
			if(i == I && j == J && k == K){
			    printf("before correction: at i = %d, j = %d, k = %d,  phi = %f, "
			   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
				 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
				 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
				 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
			}
		END_FOR


#ifdef POSITIVE_PARTICLES
		FOR_EACH_CELL
		    u_int pos = INDEX(i,j,k);
			phi[pos] = fabs(phiPos[pos]) <= fabs(phiNeg[pos]) ?	 phiPos[pos] : phiNeg[pos];
		END_FOR
#endif

	N = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f)
			++N;
	END_FOR

	printf("after error correction, there are %u liquid points \n", N);

	NM = 0;
	NP = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
			if(phi0[INDEX(i,j,k)] <= 0.f && phi[INDEX(i,j,k)] > 0.f){
	//			printf("phi changes from liquid to air at (%d, %d, %d) old phi = %f, new phi = %f\n",
	//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
				++NM;
			}
			if(phi0[INDEX(i,j,k)] > 0.f && phi[INDEX(i,j,k)] <= 0.f){
	//			printf("phi changes from air to liquid at (%d, %d, %d) old phi = %f, new phi = %f\n",
	//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
				++NP;
			}
		}
	END_FOR
	printf("\nafter error correction there are %d points liquid->air, %d points air->liquid \n\n", NM, NP);
//	delete [] x;
	FOR_EACH_CELL
		if(i == I && j == J && k == K){
		    printf("after correction: at i = %d, j = %d, k = %d,  phi = %20.16f, "
		   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
			 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
			 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
			 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
		}
	END_FOR
}

/*void ScalarField3D::ErrorCorrection(const Voxel &voxel,
							list<Particle> &particles
#ifdef POSITIVE_PARTICLES
							,list<Particle> &posParticles
#endif
							){

	map<u_int, float> corners;
#ifdef POSITIVE_PARTICLES
	map<u_int, float> corners_pos;
#endif
	map<u_int, float>::iterator iter_corner, found;
	list<Particle>::iterator iter_particle;

	float *x = new float[DimX*DimY*DimZ];
	memset(x, 0, DimX*DimY*DimZ);
	FOR_EACH_CELL
		x[INDEX(i,j,k)] = phi[INDEX(i,j,k)];
	END_FOR

	float save_corner[8];
	Point save_point[8];

	//for(int m=0; m<particles.size(); ++m){
	for(iter_particle = particles.begin();
		iter_particle != particles.end();
		++iter_particle){
//		Point pos = particles[m].Position();
//		float radius = particles[m].Radius();
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
//		if(tp.ii == I && tp.jj == J && tp.kk == K){
//			if(source->IsSourceCell(tp.ii, tp.jj, tp.kk))
//				printf(" (%d,%d,%d) is source cell \n", tp.ii, tp.jj, tp.kk);
//			else
//				printf(" (%d,%d,%d) is not source cell \n", tp.ii, tp.jj, tp.kk);
//		}
		if(tp.ii == I && tp.jj == J && tp.kk == K){
			printf(" negative particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
					pos.x, pos.y, pos.z, phi_i, radius);
		}
		if(!voxel.InSolid(tp.ii, tp.jj, tp.kk) && !source->IsSourceCell(tp.ii, tp.jj, tp.kk)){
			//float phi_p = TriInterp(voxel, pos);
			if( phi_i > radius ){
//				Point center  = voxel.VoxelCenterPosition(tp.ii, tp.jj, tp.kk);
//				float newphi  = ps.EvaluatePhiP(center);
//				float curv = Curvature(voxel, tp.ii, tp.jj, tp.kk);
//				if(tp.ii == I && tp.jj == J && tp.kk == K){
//					for(int n=0; n < 8; ++n){
//						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
//						save_corner[n] = ps.EvaluatePhiP(corner);
//					}
//					printf(" negative particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
//							pos.x, pos.y, pos.z, phi_i, radius);
//					for(int n=0; n < 8; ++n)
//						printf(" %f", save_corner[n]);
//					printf("\n");
//				}
				float local_phi = phi[INDEX(tp.ii, tp.jj, tp.kk)];
				found = corners.find(INDEX(tp.ii, tp.jj, tp.kk));

				if(found == corners.end()){
					float phi_corner;
					//float phi_tmp = -INFINITY;

						Point corner = voxel.VoxelCenterPosition(tp.ii,tp.jj,tp.kk);
						phi_corner = ps.EvaluatePhiP(corner);
						if(local_phi > phi_corner )
							local_phi = phi_corner;
	//					if(fabs(phi_corner[n]) < fabs(phi[INDEX(tp.ii,tp.jj, tp.kk)] ))
	//						phi[INDEX(tp.ii,tp.jj, tp.kk)] = phi_corner[n];

					corners.insert(make_pair(INDEX(tp.ii, tp.jj, tp.kk), local_phi));
				}
				else{
						float &lphi = found->second;
						Point corner = voxel.VoxelCenterPosition(tp.ii,tp.jj,tp.kk);
//						save_point[n] = corner;
						float phi_corner = ps.EvaluatePhiP(corner);
//						save_corner[n] = phi_corner;
//						if(phi_corner < pp[n])
						if(lphi > phi_corner)
							lphi = phi_corner;

				}
//				if(tp.ii == I && tp.jj == J && tp.kk == K){
//					printf("(Corretion): phi_i = %f, newphi = %f,  oldphi = %f, curvature = %f \n",
//							phi_i, newphi, phi[INDEX(tp.ii,tp.jj, tp.kk)], curv);
//					found = corners.find(INDEX(tp.ii, tp.jj, tp.kk));
//					if(found != corners.end()){
//						float *pp = found->second;
////						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
////								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
////						    			pos.x, pos.y, pos.z, radius,
////						    			pp[0], pp[1], pp[2], pp[3],
////						    			pp[4], pp[5], pp[6], pp[7]);
//						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
//								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
//						    			pos.x, pos.y, pos.z, radius,
//						    			save_corner[0], save_corner[1], save_corner[2], save_corner[3],
//						    			save_corner[4], save_corner[5], save_corner[6], save_corner[7]);
////						printf("(Corretion): cell corners at (%f, %f, %f) and (%f,%f,%f) \n",
////						    		save_point[0].x, save_point[0].y, save_point[0].z,
////						    		save_point[6].x, save_point[6].y, save_point[6].z);
//					}
//				}


			}

		}

	}
#ifdef POSITIVE_PARTICLES
		for(iter_particle = posParticles.begin();
			iter_particle != posParticles.end();
			++iter_particle){
	//		Point pos = particles[m].Position();
	//		float radius = particles[m].Radius();
			Particle &ps = *iter_particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			float phi_i = TriInterp(voxel, pos);
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			if(tp.ii == I && tp.jj == J && tp.kk == K){
				printf(" positive particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
						pos.x, pos.y, pos.z, phi_i, radius);
			}
			if(!voxel.InSolid(tp.ii, tp.jj, tp.kk) && !source->IsSourceCell(tp.ii, tp.jj, tp.kk)){
				//float phi_p = TriInterp(voxel, pos);
				if( -phi_i > radius ){
//					if(tp.ii == I && tp.jj == J && tp.kk == K){
//						printf(" positive particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
//								pos.x, pos.y, pos.z, phi_i, radius);
//					}
//					Point center  = voxel.VoxelCenterPosition(tp.ii, tp.jj, tp.kk);
//					float newphi  = ps.EvaluatePhiP(center);
//					float curv = Curvature(voxel, tp.ii, tp.jj, tp.kk);

					float local_phi = phi[INDEX(tp.ii, tp.jj, tp.kk)];
					found = corners_pos.find(INDEX(tp.ii, tp.jj, tp.kk));

					if(found == corners_pos.end()){
						float phi_corner;
						//float phi_tmp = -INFINITY;

							Point corner = voxel.VoxelCenterPosition(tp.ii,tp.jj,tp.kk);
							phi_corner = -ps.EvaluatePhiP(corner);
							if(local_phi < phi_corner )
								local_phi = phi_corner;
		//					if(fabs(phi_corner[n]) < fabs(phi[INDEX(tp.ii,tp.jj, tp.kk)] ))
		//						phi[INDEX(tp.ii,tp.jj, tp.kk)] = phi_corner[n];

						corners_pos.insert(make_pair(INDEX(tp.ii, tp.jj, tp.kk), local_phi));
					}
					else{
							float &lphi = found->second;

							Point corner = voxel.VoxelCenterPosition(tp.ii,tp.jj,tp.kk);
	//						save_point[n] = corner;
							float phi_corner = -ps.EvaluatePhiP(corner);
	//						save_corner[n] = phi_corner;
	//						if(phi_corner < pp[n])
							if(lphi < phi_corner)
								lphi = phi_corner;

					}
//				if(tp.ii == I && tp.jj == J && tp.kk == K){
//					printf("(Corretion): phi_i = %f, newphi = %f,  oldphi = %f, curvature = %f \n",
//							phi_i, newphi, phi[INDEX(tp.ii,tp.jj, tp.kk)], curv);
//					found = corners.find(INDEX(tp.ii, tp.jj, tp.kk));
//					if(found != corners.end()){
//						float *pp = found->second;
////						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
////								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
////						    			pos.x, pos.y, pos.z, radius,
////						    			pp[0], pp[1], pp[2], pp[3],
////						    			pp[4], pp[5], pp[6], pp[7]);
//						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
//								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
//						    			pos.x, pos.y, pos.z, radius,
//						    			save_corner[0], save_corner[1], save_corner[2], save_corner[3],
//						    			save_corner[4], save_corner[5], save_corner[6], save_corner[7]);
////						printf("(Corretion): cell corners at (%f, %f, %f) and (%f,%f,%f) \n",
////						    		save_point[0].x, save_point[0].y, save_point[0].z,
////						    		save_point[6].x, save_point[6].y, save_point[6].z);
//					}
//				}


				}

			}

		}
#endif

		FOR_EACH_CELL
			if(i == I & j == J && k == K){
			    printf("before correction: at i = %d, j = %d, k = %d,  phi = %f, "
			   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
				 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
				 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
				 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
			}
		END_FOR


#ifdef POSITIVE_PARTICLES
	for(iter_corner=corners.begin(); iter_corner!=corners.end(); ++iter_corner){
			float &pp = iter_corner->second;
			u_int ind = iter_corner->first;
			found = corners_pos.find(ind);
			float oldphi = phi[ind];
			if(found != corners_pos.end()){
				float &ppp = found->second;
				if(fabs(pp) < fabs(ppp))
					phi[ind] = pp;
				else
					phi[ind] = ppp;
			}
			else{
//				printf("no postive particel correction old_phi = %f, new_phi = %f\n ",
//						phi[ind], pp);
				phi[ind] = pp;
			}
			if(ind == INDEX(I,J,K))
				printf("oldphi = %f, newphi = %f \n",
					oldphi, phi[ind]);
		}
	for(iter_corner=corners_pos.begin(); iter_corner!=corners_pos.end(); ++iter_corner){
		float &pp = iter_corner->second;
		u_int ind = iter_corner->first;
		found = corners.find(ind);
		float oldphi = phi[ind];
		if(found != corners.end()){
			float &ppp = found->second;
			if(fabs(pp) < fabs(ppp))
				phi[ind] = pp;
			else{
//				printf("no negative particel correction old_phi = %f, new_phi = %f\n ",
//						phi[ind], pp);
				phi[ind] = ppp;
			}
		}
		else
			phi[ind] = pp;
//		if(oldphi > phi[ind])
		if(ind == INDEX(I,J,K))
			printf("oldphi = %f, newphi = %f \n",
					oldphi, phi[ind]);
	}
#else
	for(iter_corner=corners.begin(); iter_corner!=corners.end(); ++iter_corner){
		float &pp = iter_corner->second;
		u_int ind = iter_corner->first;

		float oldphi = phi[ind];
		phi[ind] = pp;
//		if(oldphi > phi[ind])
//			printf("oldphi = %f, newphi = %f \n",
//				oldphi, phi[ind]);
	}
#endif
	u_int N = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && phi[INDEX(i,j,k)] < 0.f)
			++N;
	END_FOR

	printf("after error correction, there are %ld liquid points \n", N);

	u_int NM = 0, NP = 0;
	FOR_EACH_CELL
		if(x[INDEX(i,j,k)] < 0.f && phi[INDEX(i,j,k)] > 0.f){
//			printf("phi changes from liquid to air at (%d, %d, %d) old phi = %f, new phi = %f\n",
//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
			++NM;
		}
		if(x[INDEX(i,j,k)] > 0.f && phi[INDEX(i,j,k)] < 0.f){
//			printf("phi changes from air to liquid at (%d, %d, %d) old phi = %f, new phi = %f\n",
//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
			++NP;
		}
	END_FOR
	printf("\nafter error correction there are %d points liquid->air, %d points air->liquid \n\n", NM, NP);
	delete [] x;
	FOR_EACH_CELL
		if(i == I & j == J && k == K){
		    printf("after correction: at i = %d, j = %d, k = %d,  phi = %f, "
		   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
			 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
			 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
			 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
		}
	END_FOR
}*/

struct DualFloat{
	DualFloat(){
		pos = INFINITY;
		neg = INFINITY;
	}

	float pos;
	float neg;
};

/*void ScalarField3D::ErrorCorrection(const Voxel &voxel,
							list<Particle> &negParticles,
							list<Particle> &posParticles){

	map<u_int, DualFloat*> corners;
	map<u_int, DualFloat*>::iterator iter_corner, found;
	list<Particle>::iterator iter_particle;

	float save_corner[8];
	Point save_point[8];

	for(iter_particle = negParticles.begin();
		iter_particle != negParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		if(!voxel.InSolid(tp.ii, tp.jj, tp.kk)){
			//float phi_p = TriInterp(voxel, pos);
			if( phi_i > radius ){
//				Point center  = voxel.VoxelCenterPosition(tp.ii, tp.jj, tp.kk);
//				float newphi  = ps.EvaluatePhiP(center);
//				float curv = Curvature(voxel, tp.ii, tp.jj, tp.kk);

				found = corners.find(INDEX(tp.ii, tp.jj, tp.kk));

				if(found == corners.end()){
					DualFloat *phi_corner = new DualFloat[8];
					//float phi_tmp = -INFINITY;
					for(int n=0; n < 8; ++n){
						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
						phi_corner[n].neg = -ps.EvaluatePhiP(corner);
						save_corner[n] = phi_corner[n].neg;

	//					if(fabs(phi_corner[n]) < fabs(phi[INDEX(tp.ii,tp.jj, tp.kk)] ))
	//						phi[INDEX(tp.ii,tp.jj, tp.kk)] = phi_corner[n];
					}
					for(int n=0; n < 8; ++n){
						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
						save_point[n] = corner;
						float corner_phi = TriInterp(voxel, corner);
//						save_corner[n] = corner_phi;
						if(corner_phi < phi_corner[n].neg)
							phi_corner[n].neg = corner_phi;
					}
					corners.insert(make_pair(INDEX(tp.ii, tp.jj, tp.kk), phi_corner));
				}
				else{
					DualFloat *pp = found->second;
					for(int n=0; n < 8; ++n){
						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
						save_point[n] = corner;
						float phi_corner = -ps.EvaluatePhiP(corner);
						save_corner[n] = phi_corner;
						if(phi_corner < pp[n].neg)
							pp[n].neg = phi_corner;
					}
				}
//				if(tp.ii == I && tp.jj == J && tp.kk == K){
//					printf("(Corretion): phi_i = %f, oldphi = %f,  \n",
//							phi_i, phi[INDEX(tp.ii,tp.jj, tp.kk)]);
//					found = corners.find(INDEX(tp.ii, tp.jj, tp.kk));
//					if(found != corners.end()){
//						DualFloat *pp = found->second;
//						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
//								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
//						    			pos.x, pos.y, pos.z, radius,
//						    			pp[0].neg, pp[1].neg, pp[2].neg, pp[3].neg,
//						    			pp[4].neg, pp[5].neg, pp[6].neg, pp[7].neg);
//						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
//								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
//						    			pos.x, pos.y, pos.z, radius,
//						    			save_corner[0], save_corner[1], save_corner[2], save_corner[3],
//						    			save_corner[4], save_corner[5], save_corner[6], save_corner[7]);
////						printf("(Corretion): cell corners at (%f, %f, %f) and (%f,%f,%f) \n",
////						    		save_point[0].x, save_point[0].y, save_point[0].z,
////						    		save_point[6].x, save_point[6].y, save_point[6].z);
//					}
//				}


			}

		}
	}

	for(iter_particle = posParticles.begin();
		iter_particle != posParticles.end();
		++iter_particle){
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float phi_i = TriInterp(voxel, pos);
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		if(!voxel.InSolid(tp.ii, tp.jj, tp.kk)){
			//float phi_p = TriInterp(voxel, pos);
			if( -phi_i > radius ){

				found = corners.find(INDEX(tp.ii, tp.jj, tp.kk));

				if(found == corners.end()){
					DualFloat *phi_corner = new DualFloat[8];
					//float phi_tmp = -INFINITY;
					for(int n=0; n < 8; ++n){
						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
						phi_corner[n].pos = -ps.EvaluatePhiP(corner);
						save_corner[n] = phi_corner[n].pos;

	//					if(fabs(phi_corner[n]) < fabs(phi[INDEX(tp.ii,tp.jj, tp.kk)] ))
	//						phi[INDEX(tp.ii,tp.jj, tp.kk)] = phi_corner[n];
					}
					for(int n=0; n < 8; ++n){
						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
						save_point[n] = corner;
						float corner_phi = TriInterp(voxel, corner);
//						save_corner[n] = corner_phi;
						if(corner_phi > phi_corner[n].pos)
							phi_corner[n].pos = corner_phi;
					}
					corners.insert(make_pair(INDEX(tp.ii, tp.jj, tp.kk), phi_corner));
				}
				else{
					DualFloat *pp = found->second;
					for(int n=0; n < 8; ++n){
						Point corner = voxel.VoxelCornerPosition(tp.ii,tp.jj,tp.kk,n+1);
						save_point[n] = corner;
						float phi_corner = -ps.EvaluatePhiP(corner);
						save_corner[n] = phi_corner;
						if(phi_corner > pp[n].pos)
							pp[n].pos = phi_corner;
					}
				}
//				if(tp.ii == I && tp.jj == J && tp.kk == K){
//					printf("(Corretion): phi_i = %f, newphi = %f,  oldphi = %f, curvature = %f \n",
//							phi_i, newphi, phi[INDEX(tp.ii,tp.jj, tp.kk)], curv);
//					found = corners.find(INDEX(tp.ii, tp.jj, tp.kk));
//					if(found != corners.end()){
//						float *pp = found->second;
////						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
////								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
////						    			pos.x, pos.y, pos.z, radius,
////						    			pp[0], pp[1], pp[2], pp[3],
////						    			pp[4], pp[5], pp[6], pp[7]);
//						printf("(Corretion): cell contains particle at (%f, %f, %f) wtih radius %f,"
//								" corner value = (%f, %f, %f, %f,%f, %f, %f, %f)\n",
//						    			pos.x, pos.y, pos.z, radius,
//						    			save_corner[0], save_corner[1], save_corner[2], save_corner[3],
//						    			save_corner[4], save_corner[5], save_corner[6], save_corner[7]);
////						printf("(Corretion): cell corners at (%f, %f, %f) and (%f,%f,%f) \n",
////						    		save_point[0].x, save_point[0].y, save_point[0].z,
////						    		save_point[6].x, save_point[6].y, save_point[6].z);
//					}
//				}


			}

		}
	}


	for(iter_corner=corners.begin(); iter_corner!=corners.end(); ++iter_corner){
		DualFloat *pp = iter_corner->second;
		u_int ind = iter_corner->first;
		float *fvalue = new float[8];
		for(int n=0; n < 8; ++n){
			if(pp[n].pos != INFINITY && pp[n].neg != INFINITY){
				fvalue[n] = fabs(pp[n].pos) < fabs(pp[n].neg) ? pp[n].pos : pp[n].neg;
			}
			else if (pp[n].pos != INFINITY)
				fvalue[n] = pp[n].pos;
			else if (pp[n].neg != INFINITY)
				fvalue[n] = pp[n].neg;
			else{
				printf("Error! Either one should not be INFINITY \n ");
				exit(1);
			}
		}

		float phi_t = fabs(fvalue[0]);
		int min_index = 0;
		for(int n=1; n < 8; ++n){
			if(phi_t > fabs(fvalue[n])){
				phi_t = fabs(fvalue[n]);
				min_index = n;
			}
		}
		//delete [] fvalue;
		//phi_t *= 0.125f;
		phi[ind] = fvalue[min_index];
		delete [] fvalue;
		if(pp != NULL ) delete [] pp;
	}

//	FOR_EACH_CELL
//	if(i == I & j == J && k == K){
//	    printf("after correction: at i = %d, j = %d, k = %d,  phi = %f, "
//			   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
//				 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
//				 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
//				 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
//	}
//	END_FOR
}*/

float ScalarField3D::Curvature(const Voxel &voxel, int i, int j, int k) const{

	float h = voxel.VoxelDelta();

	float laplac_phi = ((phi[INDEX(i+1,j,k)] - 2.f * phi[INDEX(i,j,k)] + phi[INDEX(i-1,j,k)]) +
					    (phi[INDEX(i,j+1,k)] - 2.f * phi[INDEX(i,j,k)] + phi[INDEX(i,j-1,k)]) +
					    (phi[INDEX(i,j,k+1)] - 2.f * phi[INDEX(i,j,k)] + phi[INDEX(i,j,k-1)]))/(h*h);

	return laplac_phi;
}

#ifdef SPMD
void ScalarField3D::ReInitialize(Voxel &voxel, const ScalarField3D &solid,
								bool initial){
	float delta = voxel.VoxelDelta();
	FOR_EACH_CELL
	if(i == I && j == J && k == K){
//		   if(voxel.InSurface(i,j,k)) printf("SURFACE CELL\n");
//		   else printf("NOT SURFACE CELL\n");
//		   if(voxel.InAir(i,j,k)) printf(" AIR CELL\n");
//		   else printf("NOT AIR CELL\n");
		   printf("before reinit: at i = %d, j = %d, k = %d,  phi = %f, "
						"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
							I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
							phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
							phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
//		   printf("before reinit: at i = %d, j = %d, k = %d,  "
//			   		"phi[i+1,j,k+1] = %f, phi[i+1,j+1,k+1] = %f, phi[i,j+1,k+1]= %f, "
//			   		"phi[i+1,j+1,k] = %f \n",
//				 I,J,K, phi[INDEX(i+1,j,k+1)], phi[INDEX(i+1,j+1,k+1)],
//				        phi[INDEX(i,j+1,k+1)], phi[INDEX(i+1,j+1,k)]);
		}
	END_FOR
	u_int N = 0, NM = 0, NP = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f)
			++N;
	END_FOR
	printf("\nbefore reinitialization, there are %u liquid points \n\n", N);

	set_time_base();
	printf(" Reinitializing Phi Begins ... \n");

	SetEqual(phi, phi0);
	float *phitmp = new float[DimX*DimY*DimZ];
	SetEqual(phi, phitmp);
	ReInitialize(voxel, solid, initial, phi, phi0);

	printf(" Reinitializing Phi Finished \n");
	printf("Time spent in Reinitializing phi = %f \n\n", get_time_in_seconds());

	FOR_EACH_CELL
	if(i == I && j == J && k == K){
//	   if(voxel.InSurface(i,j,k)) printf("SURFACE CELL\n");
//	   else printf("NOT SURFACE CELL\n");
//	   if(voxel.InAir(i,j,k)) printf(" AIR CELL\n");
//	   else printf("NOT AIR CELL\n");
	   printf("after reinit: at i = %d, j = %d, k = %d,  phi = %f, "
					"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
						I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
						phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
						phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
//	   printf("after reinit: at i = %d, j = %d, k = %d, "
//		   		"phi[i+1,j,k+1] = %f, phi[i+1,j+1,k+1] = %f, phi[i,j+1,k+1]= %f, "
//		   		"phi[i+1,j+1,k] = %f \n",
//			 I,J,K, phi[INDEX(i+1,j,k+1)], phi[INDEX(i+1,j+1,k+1)],
//			        phi[INDEX(i,j+1,k+1)], phi[INDEX(i+1,j+1,k)]);
//	   float dx, dy, dz;
//	   dx = (phi[INDEX(i+1,j,k)] - phi[INDEX(i-1,j,k)])/(2 * delta);
//	   dy = (phi[INDEX(i,j+1,k)] - phi[INDEX(i,j-1,k)])/(2 * delta);
//	   dz = (phi[INDEX(i,j,k+1)] - phi[INDEX(i,j,k-1)])/(2 * delta);
//	   Vector N(dx, dy, dz);
//	   printf("gradient (%f, %f, %f) with length = %f \n", dx, dy, dz, N.Length());
	}
//	if(k == 51 && i < DimX/2 && phi[INDEX(i,j,k)] < 0.f){
//		printf("at (%d, %d, %d) phinew = %f phiold = %f \n",
//				i, j, k, phi[INDEX(i,j,k)], phitmp[INDEX(i,j,k)] );
//		exit(1);
//	}
	END_FOR

	N = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi[INDEX(i,j,k)] <= 0.f)
			++N;
	END_FOR

	printf("\nafter reinitialization, there are %u liquid points \n", N);

	NM = 0;
	NP = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
			if(phitmp[INDEX(i,j,k)] <= 0.f && phi[INDEX(i,j,k)] > 0.f){
//				printf("phi changes from liquid to air at (%d, %d, %d) old phi = %f, new phi = %f\n",
//						i, j, k, phitmp[INDEX(i,j,k)], phi[INDEX(i,j,k)]);
				++NM;
			}
			if(phitmp[INDEX(i,j,k)] > 0.f && phi[INDEX(i,j,k)] <= 0.f){
//				printf("phi changes from air to liquid at (%d, %d, %d) old phi = %f, new phi = %f\n",
//						i, j, k, phitmp[INDEX(i,j,k)], phi[INDEX(i,j,k)]);
				++NP;
			}
		}
	END_FOR
	printf("after reinitialization there are %d points liquid->air, %d points air->liquid \n\n", NM, NP);

	delete [] phitmp;

}

#else
void ScalarField3D::ReInitialize(Voxel &voxel,
					map<u_int, KnownPoint> &pos_band,
					map<u_int, KnownPoint> &neg_band,
					const ScalarField3D &solid){

//	FOR_EACH_CELL
//	if(i == I & j == J && k == K){
//		   if(voxel.InSurface(i,j,k)) printf("SURFACE CELL\n");
//		   else printf("NOT SURFACE CELL\n");
//		   if(voxel.InAir(i,j,k)) printf(" AIR CELL\n");
//		   else printf("NOT AIR CELL\n");
//		   printf("before reinit: at i = %d, j = %d, k = %d,  phi = %f, "
//				   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
//					 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
//					 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
//					 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
//		   printf("before reinit: at i = %d, j = %d, k = %d,  "
//			   		"phi[i+1,j,k+1] = %f, phi[i+1,j+1,k+1] = %f, phi[i,j+1,k+1]= %f, "
//			   		"phi[i+1,j+1,k] = %f \n",
//				 I,J,K, phi[INDEX(i+1,j,k+1)], phi[INDEX(i+1,j+1,k+1)],
//				        phi[INDEX(i,j+1,k+1)], phi[INDEX(i+1,j+1,k)]);
//		}
//	END_FOR
	pos_band.clear();
	neg_band.clear();
	InitializeInterface(voxel, pos_band, neg_band, solid, phi);
//	FastMarching(pos_band, neg_band, voxel, false);
	FOR_EACH_CELL
	//	printf("E: phi = %f \n", phi[INDEX(I,J,K)]);
		if(i == I & j == J && k == K){
	//	   if(voxel.InSurface(i,j,k)) printf("SURFACE CELL\n");
	//	   else printf("NOT SURFACE CELL\n");
	//	   if(voxel.InAir(i,j,k)) printf(" AIR CELL\n");
	//	   else printf("NOT AIR CELL\n");
		   printf("after InitializeInterface: at i = %d, j = %d, k = %d,  phi = %f, "
				   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
					 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
					 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
					 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
	//	   printf("after reinit: at i = %d, j = %d, k = %d, "
	//		   		"phi[i+1,j,k+1] = %f, phi[i+1,j+1,k+1] = %f, phi[i,j+1,k+1]= %f, "
	//		   		"phi[i+1,j+1,k] = %f \n",
	//			 I,J,K, phi[INDEX(i+1,j,k+1)], phi[INDEX(i+1,j+1,k+1)],
	//			        phi[INDEX(i,j+1,k+1)], phi[INDEX(i+1,j+1,k)]);
		}
	END_FOR
	FastMarchingAirWater(pos_band, neg_band, voxel, false, phi);

	FOR_EACH_CELL
//	printf("E: phi = %f \n", phi[INDEX(I,J,K)]);
	if(i == I & j == J && k == K){
//	   if(voxel.InSurface(i,j,k)) printf("SURFACE CELL\n");
//	   else printf("NOT SURFACE CELL\n");
//	   if(voxel.InAir(i,j,k)) printf(" AIR CELL\n");
//	   else printf("NOT AIR CELL\n");
	   printf("after reinit: at i = %d, j = %d, k = %d,  phi = %f, "
			   		"phi[i-1] = %f, phi[i+1] = %f, phi[j-1]= %f, phi[j+1] = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
				 		I,J,K, phi[INDEX(i,j,k)], phi[INDEX(i-1,j,k)], phi[INDEX(i+1,j,k)],
				 		phi[INDEX(i,j-1,k)], phi[INDEX(i,j+1,k)],
				 		phi[INDEX(i,j,k-1)], phi[INDEX(i,j,k+1)]);
//	   printf("after reinit: at i = %d, j = %d, k = %d, "
//		   		"phi[i+1,j,k+1] = %f, phi[i+1,j+1,k+1] = %f, phi[i,j+1,k+1]= %f, "
//		   		"phi[i+1,j+1,k] = %f \n",
//			 I,J,K, phi[INDEX(i+1,j,k+1)], phi[INDEX(i+1,j+1,k+1)],
//			        phi[INDEX(i,j+1,k+1)], phi[INDEX(i+1,j+1,k)]);
	}
	END_FOR
//	 printf("reinit levelset finished \n");
}

#endif


#ifdef SPMD


#define EXTENT 5

bool ScalarField3D::NearInterface(int ii, int jj, int kk, float e, float limit, float *tmp) const{

	for(int i = ii-EXTENT; i <= ii+EXTENT; ++i)
		for(int j = jj-EXTENT; j <= jj+EXTENT; ++j)
			for(int k = kk-EXTENT; k <= kk+EXTENT; ++k){
				if(fabsf(tmp[INDEX(i,j,k)])+e < limit)
					return true;
			}
	return false;

}

#ifdef WIN32
#else
typedef float v4sf __attribute__ ((vector_size (16)));
union f4vector
{
  v4sf v;
  float f[4];
};
#endif

static inline float HJ_WENO_Coefficients(float q1, float q2, float q3, float q4, float q5){

//	float epsil = 1.e-6*max(q1,max(q2,max(q3,max(q4,q5))))+1.e-99;
//#ifdef WIN32
	float t1 = q1-2*q2+q3;
	float t2 = q1-4*q2+3*q3;
	float s1 = thirteen_twelve*t1*t1 + 0.25f*t2*t2 + E_EPSIL;
	t1 = q2-2*q3+q4;
	t2 = q2-q4;
	float s2 = thirteen_twelve*t1*t1 + 0.25f*t2*t2 + E_EPSIL;
	t1 = q3-2*q4+q5;
	t2 = 3*q3-4*q4+q5;
	float s3 = thirteen_twelve*t1*t1 + 0.25f*t2*t2 + E_EPSIL;
//#else
//	union f4vector values0, values1, values2;
//	union f4vector v0, v1, v2;
//
//	values0.f[0] = q1-2*q2+q3;
//	values1.f[0] = q1-4*q2+3*q3;
//	values0.f[1] = q2-2*q3+q4;
//	values1.f[1] = q2-q4;
//	values0.f[2] = q3-2*q4+q5;
//	values1.f[2] = 3*q3-4*q4+q5;
//	values0.f[3] = 0.f;
//	values1.f[3] = 0.f;
//	for(int n=0; n<4;++n){
//		v0.f[n] =  thirteen_twelve;
//		v1.f[n] =  0.25f;
//		v2.f[n] =  E_EPSIL;
//	}
//
//	values2.v = v0.v * values0.v * values0.v + v1.v * values1.v * values1.v  + v2.v;
//
//	float s1 = values2.f[0];
//	float s2 = values2.f[1];
//	float s3 = values2.f[2];
//
//#endif
	float a1 = 0.1f/(s1*s1);
	float a2 = 0.6f/(s2*s2);
	float a3 = 0.3f/(s3*s3);
	float at = 1.f / (a1 + a2 + a3);
	return at*(a1*( q1 * one_third - seven_sixth * q2 + eleven_sixth * q3) +
	           a2*(-q2 * one_sixth + five_sixth  * q3 + q4 * one_third) +
	           a3*( q3 * one_third + five_sixth  * q4 - q5 * one_sixth));
}


static inline float Godunov(float s_p, float s_m, float phi_m, float phi_p,
 					   float S){
 	if(s_p > 0.f && s_m > 0.f)
 		return phi_m;
 	else if(s_p < 0.f && s_m < 0.f)
 		return phi_p;
 	else if(s_p >= 0.f && s_m <= 0.f)
 		return 0.f;
 	else if(s_p <= 0.f && s_m >= 0.f)
// 		float s = S *(fabsf(phi_p) - fabsf(phi_m))/(phi_p - phi_m);
// 		return s > 0.f ? phi_m : phi_p;
 		return fabsf(s_p) > fabsf(s_m)? phi_p : phi_m;
 	else{
 		printf("Error! Shouldn't reach here \n");
 		exit(1);
 	}

}

#ifdef TBB
struct FindSmoothFactorBody {
	float * const t1;
	float * const t2;
	char * const t3;
	const Voxel * const voxel;
	int DimX, DimY, DimZ;
	FindSmoothFactorBody(const Voxel * const vox, float *s1, float *s2, char *n,
			       int nx, int ny, int nz)
		: t1(s1), t2(s2), t3(n), voxel(vox),
		   DimX(nx), DimY(ny), DimZ(nz){
//		printf("DimX = %d, DimY = %d, DimZ = %d\n ", DimX, DimY, DimZ);
	}
    void operator()( const tbb::blocked_range<int>& range ) const {
    	float *x0 = t1;
    	float *value0 = t1;
    	float *S = t2;
    	char *needed = t3;

    	float s_xp, s_yp, s_zp;
		float s_xm, s_ym, s_zm;
		float vel_minus_x, vel_minus_y, vel_minus_z;
		float vel_plus_x, vel_plus_y, vel_plus_z;
		float phi_x, phi_y, phi_z;
		float phi_xm, phi_ym, phi_zm;
		float phi_xp, phi_yp, phi_zp;

		float delta = voxel->VoxelDelta();
		float delta2 = delta * delta;
		float inv_delta = 1.f / delta;
		float q1, q2, q3, q4, q5;
//		printf("got here \n");
        int k_end = range.end();
//        printf("range at [%d  %d] \n", range.begin(), k_end);
        for( int k=range.begin(); k!=k_end; ++k ){
        	for(int j=0;j<DimY;j++) {
        		for(int i=0;i<DimX;i++){
        		u_int pos = INDEX(i,j,k);

        		if(needed[pos]){
        			bool solid_xm = false, solid_ym = false, solid_zm = false;
        			bool solid_xp = false, solid_yp = false, solid_zp = false;
        			TentativeGridPoint minTentative(value0[pos], i, j, k);
        			vector<TentativeGridPoint> neighbors;
        			minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
        			if(neighbors.size() == 6){
        				vel_minus_x = value0[INDEX(i-1, j, k)];
        				if(voxel->InSolid(i-1, j, k) && !voxel->InMovingSolid(i-1, j, k))
        					solid_xm = true;
        				vel_plus_x = value0[INDEX(i+1, j, k)];
						if(voxel->InSolid(i+1, j, k) && !voxel->InMovingSolid(i+1, j, k))
							solid_xp = true;
						vel_minus_y = value0[INDEX(i, j-1, k)];
						if(voxel->InSolid(i, j-1, k) && !voxel->InMovingSolid(i, j-1, k))
							solid_ym = true;
						vel_plus_y = value0[INDEX(i, j+1, k)];
						if(voxel->InSolid(i, j+1, k) && !voxel->InMovingSolid(i, j+1, k))
							solid_yp = true;
						vel_minus_z = value0[INDEX(i, j, k-1)];
						if(voxel->InSolid(i, j, k-1) && !voxel->InMovingSolid(i, j, k-1))
							solid_zm = true;
						vel_plus_z = value0[INDEX(i, j, k+1)];
						if(voxel->InSolid(i, j, k+1) && !voxel->InMovingSolid(i, j, k+1))
							solid_zp = true;
        			}
        			else{
						for(int m=0; m<neighbors.size(); m++){
							if(minTentative.LeftNeighbor(neighbors[m])){
								vel_minus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_xm = true;
							}
							if(minTentative.RightNeighbor(neighbors[m])){
								vel_plus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_xp = true;
							}
							if(minTentative.BackNeighbor(neighbors[m])){
								vel_minus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_ym = true;
							}
							if(minTentative.FrontNeighbor(neighbors[m])){
								vel_plus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_yp = true;
							}
							if(minTentative.BottomNeighbor(neighbors[m])){
								vel_minus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_zm = true;
							}
							if(minTentative.TopNeighbor(neighbors[m])){
								vel_plus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_zp = true;
							}
						}
        			}

        			float s = value0[pos];

        			if( (voxel->InSolid(i-3,j,k) && !voxel->InMovingSolid(i-3, j, k)) ||
						(voxel->InSolid(i-2,j,k) && !voxel->InMovingSolid(i-2, j, k)) ||
						(voxel->InSolid(i-1,j,k) && !voxel->InMovingSolid(i-1, j, k)) ||
					    (voxel->InSolid(i+1,j,k) && !voxel->InMovingSolid(i+1, j, k)) ||
					    (voxel->InSolid(i+2,j,k) && !voxel->InMovingSolid(i+2, j, k)) ){
						if(solid_xm)
							vel_minus_x = value0[pos] + sign(value0[pos]) * delta;
						phi_xm = (value0[pos] - vel_minus_x) * inv_delta ;
					}
					else{
						q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
						q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
						phi_xm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}

					if( (voxel->InSolid(i-2,j,k) && !voxel->InMovingSolid(i-2, j, k)) ||
						(voxel->InSolid(i-1,j,k) && !voxel->InMovingSolid(i-1, j, k)) ||
						(voxel->InSolid(i+3,j,k) && !voxel->InMovingSolid(i+3, j, k)) ||
						(voxel->InSolid(i+1,j,k) && !voxel->InMovingSolid(i+1, j, k)) ||
						(voxel->InSolid(i+2,j,k) && !voxel->InMovingSolid(i+2, j, k)) ){
						if(solid_xp)
							vel_plus_x = value0[pos] + sign(value0[pos]) * delta;
						phi_xp = (vel_plus_x - value0[pos]) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
						q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
						q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
						phi_xp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}

					if( (voxel->InSolid(i,j-3,k) && !voxel->InMovingSolid(i, j-3, k)) ||
						(voxel->InSolid(i,j-2,k) && !voxel->InMovingSolid(i, j-2, k)) ||
						(voxel->InSolid(i,j-1,k) && !voxel->InMovingSolid(i, j-1, k)) ||
					    (voxel->InSolid(i,j+1,k) && !voxel->InMovingSolid(i, j+1, k)) ||
					    (voxel->InSolid(i,j+2,k) && !voxel->InMovingSolid(i, j+2, k)) ){
						if(solid_ym)
							vel_minus_y = value0[pos] + sign(value0[pos]) * delta;
						phi_ym = (value0[pos] - vel_minus_y) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
						q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
						phi_ym = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}
					if( (voxel->InSolid(i,j-2,k) && !voxel->InMovingSolid(i, j-2, k)) ||
						(voxel->InSolid(i,j-1,k) && !voxel->InMovingSolid(i, j-1, k)) ||
						(voxel->InSolid(i,j+3,k) && !voxel->InMovingSolid(i, j+3, k)) ||
					    (voxel->InSolid(i,j+1,k) && !voxel->InMovingSolid(i, j+1, k)) ||
					    (voxel->InSolid(i,j+2,k) && !voxel->InMovingSolid(i, j+2, k)) ){
						if(solid_yp)
							vel_plus_y = value0[pos] + sign(value0[pos]) * delta;
						phi_yp = (vel_plus_y - value0[pos]) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
						q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
						q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
						q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
						phi_yp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}

					if( (voxel->InSolid(i,j,k-3) && !voxel->InMovingSolid(i, j, k-3)) ||
					    (voxel->InSolid(i,j,k-2) && !voxel->InMovingSolid(i, j, k-2)) ||
					    (voxel->InSolid(i,j,k-1) && !voxel->InMovingSolid(i, j, k-1)) ||
					    (voxel->InSolid(i,j,k+1) && !voxel->InMovingSolid(i, j, k+1)) ||
					    (voxel->InSolid(i,j,k+2) && !voxel->InMovingSolid(i, j, k+2)) ){
						if(solid_zm)
							vel_minus_z = value0[pos] + sign(value0[pos]) * delta;
						phi_zm = (value0[pos] - vel_minus_z) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
						q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
						phi_zm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}
					if( (voxel->InSolid(i,j,k-2) && !voxel->InMovingSolid(i, j, k-2)) ||
						(voxel->InSolid(i,j,k-1) && !voxel->InMovingSolid(i, j, k-1)) ||
						(voxel->InSolid(i,j,k+3) && !voxel->InMovingSolid(i, j, k+3)) ||
					    (voxel->InSolid(i,j,k+1) && !voxel->InMovingSolid(i, j, k+1)) ||
					    (voxel->InSolid(i,j,k+2) && !voxel->InMovingSolid(i, j, k+2)) ){
						if(solid_zp)
							vel_plus_z = value0[pos] + sign(value0[pos]) * delta;
						phi_zp = (vel_plus_z - value0[pos]) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
						q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
						q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
						phi_zp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}

        //				if(i == I && j == J && k == K)
        //					printf("S = %f, phi_ym = %f, phi_yp = %f value0 = %f, vel_minus_y = %f\n ",
        //							S, phi_ym, phi_yp, value0[pos], vel_minus_y);

        			s_xp = s * phi_xp;
        			s_xm = s * phi_xm;
        			s_yp = s * phi_yp;
        			s_ym = s * phi_ym;
        			s_zp = s * phi_zp;
        			s_zm = s * phi_zm;

        			phi_x = Godunov(s_xp, s_xm, phi_xm, phi_xp, s);
        			phi_y = Godunov(s_yp, s_ym, phi_ym, phi_yp, s);
        			phi_z = Godunov(s_zp, s_zm, phi_zm, phi_zp, s);

        			S[pos] = value0[pos] / sqrtf(value0[pos] * value0[pos] +
        					(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) * delta2);

        		}
				}
			}
		}
  }
};

struct EulerStepBody {
	float * const t1;
	float * const t2;
	float * const t3;
	char * const t4;
	const Voxel * const voxel;
	int DimX, DimY, DimZ;
	EulerStepBody(const Voxel * const vox, float *s1, float *s2, float *s3, char *n,
			       int nx, int ny, int nz)
		: t1(s1), t2(s2), t3(s3), t4(n), voxel(vox),
		   DimX(nx), DimY(ny), DimZ(nz){
//		printf("DimX = %d, DimY = %d, DimZ = %d\n ", DimX, DimY, DimZ);
	}
    void operator()( const tbb::blocked_range<int>& range ) const {
    	float *x0 = t1;
    	float *value0 = t1;
    	float *S = t2;
    	float *value = t3;
    	char *needed = t4;

    	float s_xp, s_yp, s_zp;
		float s_xm, s_ym, s_zm;
		float vel_minus_x, vel_minus_y, vel_minus_z;
		float vel_plus_x, vel_plus_y, vel_plus_z;
		float phi_x, phi_y, phi_z;
		float phi_xm, phi_ym, phi_zm;
		float phi_xp, phi_yp, phi_zp;

		float delta = voxel->VoxelDelta();
		float delta2 = delta * delta;
		float inv_delta = 1.f / delta;
		float dtau = delta / 2 ;

		float q1, q2, q3, q4, q5;
//		printf("got here \n");
        int k_end = range.end();
        for( int k=range.begin(); k!=k_end; ++k ){
        	for(int j=0;j<DimY;j++) {
        		for(int i=0;i<DimX;i++){
        		u_int pos = INDEX(i,j,k);

        		if(needed[pos]){
        			bool solid_xm = false, solid_ym = false, solid_zm = false;
        			bool solid_xp = false, solid_yp = false, solid_zp = false;
        			TentativeGridPoint minTentative(value0[pos], i, j, k);
        			vector<TentativeGridPoint> neighbors;
        			minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
        			if(neighbors.size() == 6){
						vel_minus_x = value0[INDEX(i-1, j, k)];
						if(voxel->InSolid(i-1, j, k) && !voxel->InMovingSolid(i-1, j, k))
							solid_xm = true;
						vel_plus_x = value0[INDEX(i+1, j, k)];
						if(voxel->InSolid(i+1, j, k) && !voxel->InMovingSolid(i+1, j, k))
							solid_xp = true;
						vel_minus_y = value0[INDEX(i, j-1, k)];
						if(voxel->InSolid(i, j-1, k) && !voxel->InMovingSolid(i, j-1, k))
							solid_ym = true;
						vel_plus_y = value0[INDEX(i, j+1, k)];
						if(voxel->InSolid(i, j+1, k) && !voxel->InMovingSolid(i, j+1, k))
							solid_yp = true;
						vel_minus_z = value0[INDEX(i, j, k-1)];
						if(voxel->InSolid(i, j, k-1) && !voxel->InMovingSolid(i, j, k-1))
							solid_zm = true;
						vel_plus_z = value0[INDEX(i, j, k+1)];
						if(voxel->InSolid(i, j, k+1) && !voxel->InMovingSolid(i, j, k+1))
							solid_zp = true;
					}
					else{
						for(int m=0; m<neighbors.size(); m++){
							if(minTentative.LeftNeighbor(neighbors[m])){
								vel_minus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_xm = true;
							}
							if(minTentative.RightNeighbor(neighbors[m])){
								vel_plus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_xp = true;
							}
							if(minTentative.BackNeighbor(neighbors[m])){
								vel_minus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_ym = true;
							}
							if(minTentative.FrontNeighbor(neighbors[m])){
								vel_plus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_yp = true;
							}
							if(minTentative.BottomNeighbor(neighbors[m])){
								vel_minus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_zm = true;
							}
							if(minTentative.TopNeighbor(neighbors[m])){
								vel_plus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
								if(voxel->InSolid(neighbors[m]) && !voxel->InMovingSolid(neighbors[m]))
									solid_zp = true;
							}
						}
					}

        			float s = S[pos];

        			if( (voxel->InSolid(i-3,j,k) && !voxel->InMovingSolid(i-3, j, k)) ||
						(voxel->InSolid(i-2,j,k) && !voxel->InMovingSolid(i-2, j, k)) ||
						(voxel->InSolid(i-1,j,k) && !voxel->InMovingSolid(i-1, j, k)) ||
					    (voxel->InSolid(i+1,j,k) && !voxel->InMovingSolid(i+1, j, k)) ||
					    (voxel->InSolid(i+2,j,k) && !voxel->InMovingSolid(i+2, j, k)) ){
						if(solid_xm)
							vel_minus_x = value0[pos] + sign(value0[pos]) * delta;
						phi_xm = (value0[pos] - vel_minus_x) * inv_delta ;
					}
					else{
						q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
						q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
						phi_xm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}

					if( (voxel->InSolid(i-2,j,k) && !voxel->InMovingSolid(i-2, j, k)) ||
						(voxel->InSolid(i-1,j,k) && !voxel->InMovingSolid(i-1, j, k)) ||
						(voxel->InSolid(i+3,j,k) && !voxel->InMovingSolid(i+3, j, k)) ||
						(voxel->InSolid(i+1,j,k) && !voxel->InMovingSolid(i+1, j, k)) ||
						(voxel->InSolid(i+2,j,k) && !voxel->InMovingSolid(i+2, j, k)) ){
						if(solid_xp)
							vel_plus_x = value0[pos] + sign(value0[pos]) * delta;
						phi_xp = (vel_plus_x - value0[pos]) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
						q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
						q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
						phi_xp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}

					if( (voxel->InSolid(i,j-3,k) && !voxel->InMovingSolid(i, j-3, k)) ||
						(voxel->InSolid(i,j-2,k) && !voxel->InMovingSolid(i, j-2, k)) ||
						(voxel->InSolid(i,j-1,k) && !voxel->InMovingSolid(i, j-1, k)) ||
					    (voxel->InSolid(i,j+1,k) && !voxel->InMovingSolid(i, j+1, k)) ||
					    (voxel->InSolid(i,j+2,k) && !voxel->InMovingSolid(i, j+2, k)) ){
						if(solid_ym)
							vel_minus_y = value0[pos] + sign(value0[pos]) * delta;
						phi_ym = (value0[pos] - vel_minus_y) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
						q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
						phi_ym = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}
					if( (voxel->InSolid(i,j-2,k) && !voxel->InMovingSolid(i, j-2, k)) ||
						(voxel->InSolid(i,j-1,k) && !voxel->InMovingSolid(i, j-1, k)) ||
						(voxel->InSolid(i,j+3,k) && !voxel->InMovingSolid(i, j+3, k)) ||
					    (voxel->InSolid(i,j+1,k) && !voxel->InMovingSolid(i, j+1, k)) ||
					    (voxel->InSolid(i,j+2,k) && !voxel->InMovingSolid(i, j+2, k)) ){
						if(solid_yp)
							vel_plus_y = value0[pos] + sign(value0[pos]) * delta;
						phi_yp = (vel_plus_y - value0[pos]) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
						q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
						q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
						q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
						phi_yp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}

					if( (voxel->InSolid(i,j,k-3) && !voxel->InMovingSolid(i, j, k-3)) ||
					    (voxel->InSolid(i,j,k-2) && !voxel->InMovingSolid(i, j, k-2)) ||
					    (voxel->InSolid(i,j,k-1) && !voxel->InMovingSolid(i, j, k-1)) ||
					    (voxel->InSolid(i,j,k+1) && !voxel->InMovingSolid(i, j, k+1)) ||
					    (voxel->InSolid(i,j,k+2) && !voxel->InMovingSolid(i, j, k+2)) ){
						if(solid_zm)
							vel_minus_z = value0[pos] + sign(value0[pos]) * delta;
						phi_zm = (value0[pos] - vel_minus_z) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
						q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
						q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
						phi_zm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}
					if( (voxel->InSolid(i,j,k-2) && !voxel->InMovingSolid(i, j, k-2)) ||
						(voxel->InSolid(i,j,k-1) && !voxel->InMovingSolid(i, j, k-1)) ||
						(voxel->InSolid(i,j,k+3) && !voxel->InMovingSolid(i, j, k+3)) ||
					    (voxel->InSolid(i,j,k+1) && !voxel->InMovingSolid(i, j, k+1)) ||
					    (voxel->InSolid(i,j,k+2) && !voxel->InMovingSolid(i, j, k+2)) ){
						if(solid_zp)
							vel_plus_z = value0[pos] + sign(value0[pos]) * delta;
						phi_zp = (vel_plus_z - value0[pos]) * inv_delta;
					}
					else{
						q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
						q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
						q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
						q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
						q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
						phi_zp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
					}



        //				if(i == I && j == J && k == K)
        //					printf("S = %f, phi_ym = %f, phi_yp = %f value0 = %f, vel_minus_y = %f\n ",
        //							S, phi_ym, phi_yp, value0[pos], vel_minus_y);

        			s_xp = s * phi_xp;
        			s_xm = s * phi_xm;
        			s_yp = s * phi_yp;
        			s_ym = s * phi_ym;
        			s_zp = s * phi_zp;
        			s_zm = s * phi_zm;

        			phi_x = Godunov(s_xp, s_xm, phi_xm, phi_xp, s);
        			phi_y = Godunov(s_yp, s_ym, phi_ym, phi_yp, s);
        			phi_z = Godunov(s_zp, s_zm, phi_zm, phi_zp, s);


        			value[pos] = value0[pos] - dtau * s * ( sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) - 1.f );
        		}
				}
			}
		}
  }
};

struct AddTwoBody {
	float * const t1;
	float * const t2;
	float * const t3;
	char * const mask;
	float r1, r2;
	int DimX, DimY, DimZ;
	AddTwoBody(char *n, float *s1, float *s2, float *s3,
				float R1, float R2,
			       int nx, int ny, int nz)
		: t1(s1), t2(s2), t3(s3), mask(n), r1(R1), r2(R2),
		   DimX(nx), DimY(ny), DimZ(nz){
//		printf("DimX = %d, DimY = %d, DimZ = %d\n ", DimX, DimY, DimZ);
	}
    void operator()( const tbb::blocked_range<int>& range ) const {

        int k_end = range.end();
        for( int k=range.begin(); k!=k_end; ++k ){
        	for(int j=0;j<DimY;j++) {
        		for(int i=0;i<DimX;i++){
        		u_int pos = INDEX(i,j,k);
        		if(mask[pos])
        			t3[pos] = r1 * t1[pos] + r2 * t2[pos];
				}
			}
		}
  }
};

#endif

void ScalarField3D::FindSmoothFactor(const Voxel &voxel, char *needed, float *value0, float *S){


#ifdef TBB
	FindSmoothFactorBody body(&voxel, value0, S, needed, DimX, DimY, DimZ);
//	printf("FindSmoothFactor body created \n");
//	tbb::parallel_for(tbb::blocked_range<int>(0, DimZ),
//			body, tbb::auto_partitioner( ));
	tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
				body);
#else
	float s_xp, s_yp, s_zp;
	float s_xm, s_ym, s_zm;
	float vel_minus_x, vel_minus_y, vel_minus_z;
	float vel_plus_x, vel_plus_y, vel_plus_z;
	float phi_x, phi_y, phi_z;
	float phi_xm, phi_ym, phi_zm;
	float phi_xp, phi_yp, phi_zp;

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float inv_delta = 1.f / delta;
	float dtau = delta / 2 ;

 	float q1, q2, q3, q4, q5;
	float *x0 = value0;

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);

		if(needed[pos]){
			bool solid_xm = false, solid_ym = false, solid_zm = false;
			bool solid_xp = false, solid_yp = false, solid_zp = false;
			TentativeGridPoint minTentative(value0[pos], i, j, k);
			vector<TentativeGridPoint> neighbors;
			minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
			if(neighbors.size() == 6){
				vel_minus_x = value0[INDEX(i-1, j, k)];
				if(voxel.InSolid(i-1, j, k) && !voxel.InMovingSolid(i-1, j, k))
					solid_xm = true;
				vel_plus_x = value0[INDEX(i+1, j, k)];
				if(voxel.InSolid(i+1, j, k) && !voxel.InMovingSolid(i+1, j, k))
					solid_xp = true;
				vel_minus_y = value0[INDEX(i, j-1, k)];
				if(voxel.InSolid(i, j-1, k) && !voxel.InMovingSolid(i, j-1, k))
					solid_ym = true;
				vel_plus_y = value0[INDEX(i, j+1, k)];
				if(voxel.InSolid(i, j+1, k) && !voxel.InMovingSolid(i, j+1, k))
					solid_yp = true;
				vel_minus_z = value0[INDEX(i, j, k-1)];
				if(voxel.InSolid(i, j, k-1) && !voxel.InMovingSolid(i, j, k-1))
					solid_zm = true;
				vel_plus_z = value0[INDEX(i, j, k+1)];
				if(voxel.InSolid(i, j, k+1) && !voxel.InMovingSolid(i, j, k+1))
					solid_zp = true;
			}
			else{
				for(int m=0; m<neighbors.size(); m++){
					if(minTentative.LeftNeighbor(neighbors[m])){
						vel_minus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_xm = true;
					}
					if(minTentative.RightNeighbor(neighbors[m])){
						vel_plus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_xp = true;
					}
					if(minTentative.BackNeighbor(neighbors[m])){
						vel_minus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_ym = true;
					}
					if(minTentative.FrontNeighbor(neighbors[m])){
						vel_plus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_yp = true;
					}
					if(minTentative.BottomNeighbor(neighbors[m])){
						vel_minus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_zm = true;
					}
					if(minTentative.TopNeighbor(neighbors[m])){
						vel_plus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_zp = true;
					}
				}
			}

//				float S = value0[pos] / sqrtf(value0[pos] * value0[pos] + delta2);
			float s = value0[pos];

			if( (voxel.InSolid(i-3,j,k) && !voxel.InMovingSolid(i-3, j, k)) ||
				(voxel.InSolid(i-2,j,k) && !voxel.InMovingSolid(i-2, j, k)) ||
				(voxel.InSolid(i-1,j,k) && !voxel.InMovingSolid(i-1, j, k)) ||
			    (voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1, j, k)) ||
			    (voxel.InSolid(i+2,j,k) && !voxel.InMovingSolid(i+2, j, k)) ){
				if(solid_xm)
					vel_minus_x = value0[pos] + sign(value0[pos]) * delta;
				phi_xm = (value0[pos] - vel_minus_x) * inv_delta ;
			}
			else{
				q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
				q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
				phi_xm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}

			if( (voxel.InSolid(i-2,j,k) && !voxel.InMovingSolid(i-2, j, k))	||
				(voxel.InSolid(i-1,j,k) && !voxel.InMovingSolid(i-1, j, k)) ||
				(voxel.InSolid(i+3,j,k) && !voxel.InMovingSolid(i+3, j, k)) ||
				(voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1, j, k)) ||
				(voxel.InSolid(i+2,j,k) && !voxel.InMovingSolid(i+2, j, k)) ){
				if(solid_xp)
					vel_plus_x = value0[pos] + sign(value0[pos]) * delta;
				phi_xp = (vel_plus_x - value0[pos]) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
				q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
				q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
				phi_xp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}

			if( (voxel.InSolid(i,j-3,k) && !voxel.InMovingSolid(i, j-3, k)) ||
				(voxel.InSolid(i,j-2,k) && !voxel.InMovingSolid(i, j-2, k)) ||
				(voxel.InSolid(i,j-1,k) && !voxel.InMovingSolid(i, j-1, k)) ||
			    (voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i, j+1, k)) ||
			    (voxel.InSolid(i,j+2,k) && !voxel.InMovingSolid(i, j+2, k)) ){
				if(solid_ym)
					vel_minus_y = value0[pos] + sign(value0[pos]) * delta;
				phi_ym = (value0[pos] - vel_minus_y) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
				q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
				phi_ym = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}
			if( (voxel.InSolid(i,j-2,k) && !voxel.InMovingSolid(i, j-2, k)) ||
				(voxel.InSolid(i,j-1,k) && !voxel.InMovingSolid(i, j-1, k)) ||
				(voxel.InSolid(i,j+3,k) && !voxel.InMovingSolid(i, j+3, k)) ||
			    (voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i, j+1, k)) ||
			    (voxel.InSolid(i,j+2,k) && !voxel.InMovingSolid(i, j+2, k)) ){
				if(solid_yp)
					vel_plus_y = value0[pos] + sign(value0[pos]) * delta;
				phi_yp = (vel_plus_y - value0[pos]) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
				q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
				q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
				q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
				phi_yp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}

			if( (voxel.InSolid(i,j,k-3) && !voxel.InMovingSolid(i, j, k-3)) ||
			    (voxel.InSolid(i,j,k-2) && !voxel.InMovingSolid(i, j, k-2)) ||
			    (voxel.InSolid(i,j,k-1) && !voxel.InMovingSolid(i, j, k-1)) ||
			    (voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i, j, k+1)) ||
			    (voxel.InSolid(i,j,k+2) && !voxel.InMovingSolid(i, j, k+2)) ){
				if(solid_zm)
					vel_minus_z = value0[pos] + sign(value0[pos]) * delta;
				phi_zm = (value0[pos] - vel_minus_z) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
				q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
				phi_zm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}
			if( (voxel.InSolid(i,j,k-2) && !voxel.InMovingSolid(i, j, k-2)) ||
				(voxel.InSolid(i,j,k-1) && !voxel.InMovingSolid(i, j, k-1)) ||
				(voxel.InSolid(i,j,k+3) && !voxel.InMovingSolid(i, j, k+3)) ||
			    (voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i, j, k+1)) ||
			    (voxel.InSolid(i,j,k+2) && !voxel.InMovingSolid(i, j, k+2)) ){
				if(solid_zp)
					vel_plus_z = value0[pos] + sign(value0[pos]) * delta;
				phi_zp = (vel_plus_z - value0[pos]) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
				q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
				q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
				phi_zp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}



//				if(i == I && j == J && k == K)
//					printf("S = %f, phi_ym = %f, phi_yp = %f value0 = %f, vel_minus_y = %f\n ",
//							S, phi_ym, phi_yp, value0[pos], vel_minus_y);

			s_xp = s * phi_xp;
			s_xm = s * phi_xm;
			s_yp = s * phi_yp;
			s_ym = s * phi_ym;
			s_zp = s * phi_zp;
			s_zm = s * phi_zm;

			phi_x = Godunov(s_xp, s_xm, phi_xm, phi_xp, s);
			phi_y = Godunov(s_yp, s_ym, phi_ym, phi_yp, s);
			phi_z = Godunov(s_zp, s_zm, phi_zm, phi_zp, s);

			S[pos] = value0[pos] / sqrtf(value0[pos] * value0[pos] +
					(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) * delta2);

		}
	END_FOR
#endif
}


/*void ScalarField3D::EulerStep(const Voxel &voxel, char *needed, float *S, float *value, float *value0){


#ifdef TBB
	EulerStepBody body(&voxel, value0, S, value, needed, DimX, DimY, DimZ);
//	printf("Eulerstep body created \n");
//	tbb::parallel_for(tbb::blocked_range<int>(0, DimZ),
//			body, tbb::auto_partitioner( ));
	tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
					body);
#else
	float s_xp, s_yp, s_zp;
	float s_xm, s_ym, s_zm;
	float vel_minus_x, vel_minus_y, vel_minus_z;
	float vel_plus_x, vel_plus_y, vel_plus_z;
	float phi_x, phi_y, phi_z;
	float phi_xm, phi_ym, phi_zm;
	float phi_xp, phi_yp, phi_zp;

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float inv_delta = 1.f / delta;
	float dtau = delta / 2 ;

 	float q1, q2, q3, q4, q5;
	float *x0 = value0;

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);

		if(needed[pos]){
			bool solid_xm = false, solid_ym = false, solid_zm = false;
			bool solid_xp = false, solid_yp = false, solid_zp = false;
			TentativeGridPoint minTentative(value0[pos], i, j, k);
			vector<TentativeGridPoint> neighbors;
			minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
			if(neighbors.size() == 6){
				vel_minus_x = value0[INDEX(i-1, j, k)];
				if(voxel.InSolid(i-1, j, k) && !voxel.InMovingSolid(i-1, j, k))
					solid_xm = true;
				vel_plus_x = value0[INDEX(i+1, j, k)];
				if(voxel.InSolid(i+1, j, k) && !voxel.InMovingSolid(i+1, j, k))
					solid_xp = true;
				vel_minus_y = value0[INDEX(i, j-1, k)];
				if(voxel.InSolid(i, j-1, k) && !voxel.InMovingSolid(i, j-1, k))
					solid_ym = true;
				vel_plus_y = value0[INDEX(i, j+1, k)];
				if(voxel.InSolid(i, j+1, k) && !voxel.InMovingSolid(i, j+1, k))
					solid_yp = true;
				vel_minus_z = value0[INDEX(i, j, k-1)];
				if(voxel.InSolid(i, j, k-1) && !voxel.InMovingSolid(i, j, k-1))
					solid_zm = true;
				vel_plus_z = value0[INDEX(i, j, k+1)];
				if(voxel.InSolid(i, j, k+1) && !voxel.InMovingSolid(i, j, k+1))
					solid_zp = true;
			}
			else{
				for(int m=0; m<neighbors.size(); m++){
					if(minTentative.LeftNeighbor(neighbors[m])){
						vel_minus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_xm = true;
					}
					if(minTentative.RightNeighbor(neighbors[m])){
						vel_plus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_xp = true;
					}
					if(minTentative.BackNeighbor(neighbors[m])){
						vel_minus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_ym = true;
					}
					if(minTentative.FrontNeighbor(neighbors[m])){
						vel_plus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_yp = true;
					}
					if(minTentative.BottomNeighbor(neighbors[m])){
						vel_minus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_zm = true;
					}
					if(minTentative.TopNeighbor(neighbors[m])){
						vel_plus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
							solid_zp = true;
					}
				}
			}

//				float S = value0[pos] / sqrtf(value0[pos] * value0[pos] + delta2);
			float s = S[pos];

			if( (voxel.InSolid(i-3,j,k) && !voxel.InMovingSolid(i-3, j, k)) ||
				(voxel.InSolid(i-2,j,k) && !voxel.InMovingSolid(i-2, j, k)) ||
				(voxel.InSolid(i-1,j,k) && !voxel.InMovingSolid(i-1, j, k)) ||
			    (voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1, j, k)) ||
			    (voxel.InSolid(i+2,j,k) && !voxel.InMovingSolid(i+2, j, k)) ){
				if(solid_xm)
					vel_minus_x = value0[pos] + sign(value0[pos]) * delta;
				phi_xm = (value0[pos] - vel_minus_x) * inv_delta ;
			}
			else{
				q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
				q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
				phi_xm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}

			if( (voxel.InSolid(i-2,j,k) && !voxel.InMovingSolid(i-2, j, k))	||
				(voxel.InSolid(i-1,j,k) && !voxel.InMovingSolid(i-1, j, k)) ||
				(voxel.InSolid(i+3,j,k) && !voxel.InMovingSolid(i+3, j, k)) ||
				(voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1, j, k)) ||
				(voxel.InSolid(i+2,j,k) && !voxel.InMovingSolid(i+2, j, k)) ){
				if(solid_xp)
					vel_plus_x = value0[pos] + sign(value0[pos]) * delta;
				phi_xp = (vel_plus_x - value0[pos]) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
				q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
				q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
				phi_xp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}

			if( (voxel.InSolid(i,j-3,k) && !voxel.InMovingSolid(i, j-3, k)) ||
				(voxel.InSolid(i,j-2,k) && !voxel.InMovingSolid(i, j-2, k)) ||
				(voxel.InSolid(i,j-1,k) && !voxel.InMovingSolid(i, j-1, k)) ||
			    (voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i, j+1, k)) ||
			    (voxel.InSolid(i,j+2,k) && !voxel.InMovingSolid(i, j+2, k)) ){
				if(solid_ym)
					vel_minus_y = value0[pos] + sign(value0[pos]) * delta;
				phi_ym = (value0[pos] - vel_minus_y) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
				q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
				phi_ym = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}
			if( (voxel.InSolid(i,j-2,k) && !voxel.InMovingSolid(i, j-2, k)) ||
				(voxel.InSolid(i,j-1,k) && !voxel.InMovingSolid(i, j-1, k)) ||
				(voxel.InSolid(i,j+3,k) && !voxel.InMovingSolid(i, j+3, k)) ||
			    (voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i, j+1, k)) ||
			    (voxel.InSolid(i,j+2,k) && !voxel.InMovingSolid(i, j+2, k)) ){
				if(solid_yp)
					vel_plus_y = value0[pos] + sign(value0[pos]) * delta;
				phi_yp = (vel_plus_y - value0[pos]) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
				q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
				q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
				q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
				phi_yp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}

			if( (voxel.InSolid(i,j,k-3) && !voxel.InMovingSolid(i, j, k-3)) ||
			    (voxel.InSolid(i,j,k-2) && !voxel.InMovingSolid(i, j, k-2)) ||
			    (voxel.InSolid(i,j,k-1) && !voxel.InMovingSolid(i, j, k-1)) ||
			    (voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i, j, k+1)) ||
			    (voxel.InSolid(i,j,k+2) && !voxel.InMovingSolid(i, j, k+2)) ){
				if(solid_zm)
					vel_minus_z = value0[pos] + sign(value0[pos]) * delta;
				phi_zm = (value0[pos] - vel_minus_z) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
				q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
				q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
				phi_zm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}
			if( (voxel.InSolid(i,j,k-2) && !voxel.InMovingSolid(i, j, k-2)) ||
				(voxel.InSolid(i,j,k-1) && !voxel.InMovingSolid(i, j, k-1)) ||
				(voxel.InSolid(i,j,k+3) && !voxel.InMovingSolid(i, j, k+3)) ||
			    (voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i, j, k+1)) ||
			    (voxel.InSolid(i,j,k+2) && !voxel.InMovingSolid(i, j, k+2)) ){
				if(solid_zp)
					vel_plus_z = value0[pos] + sign(value0[pos]) * delta;
				phi_zp = (vel_plus_z - value0[pos]) * inv_delta;
			}
			else{
				q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
				q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
				q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
				q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
				q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
				phi_zp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
			}



//				if(i == I && j == J && k == K)
//					printf("S = %f, phi_ym = %f, phi_yp = %f value0 = %f, vel_minus_y = %f\n ",
//							S, phi_ym, phi_yp, value0[pos], vel_minus_y);

			s_xp = s * phi_xp;
			s_xm = s * phi_xm;
			s_yp = s * phi_yp;
			s_ym = s * phi_ym;
			s_zp = s * phi_zp;
			s_zm = s * phi_zm;

			phi_x = Godunov(s_xp, s_xm, phi_xm, phi_xp, s);
			phi_y = Godunov(s_yp, s_ym, phi_ym, phi_yp, s);
			phi_z = Godunov(s_zp, s_zm, phi_zm, phi_zp, s);


			value[pos] = value0[pos] - dtau * s * ( sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) - 1.f );

//			if(old_phi*value[pos]< 0.f)
//				if(i == I && j == J && k == K){
//					printf(" (reinit) phi changes sign at iter = %d i=%d, j=%d, k=%d, old = %f, new = %f \n",
//							n, i, j, k, old_phi, value[pos]);
//					printf("S = %f, phi_x = %f, phi_y = %f, phi_z = %f \n ", S, phi_x, phi_y, phi_z);
//				}
		}
	END_FOR
#endif

}*/


static inline float WENOCoeffsX(int i, int j, int k, int DimX, int DimY, int DimZ,
						const float *x0, float inv_delta, char s){
	float q1, q2, q3, q4, q5;
	if(s < 0){
		q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
		q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
		q3 = (x0[INDEX(i,j,k)]   - x0[INDEX(i-1,j,k)]) * inv_delta;
		q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
		q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
	}
	else{
		q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
		q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
		q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
		q4 = (x0[INDEX(i,j,k)]   - x0[INDEX(i-1,j,k)]) * inv_delta;
		q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
	}
	return HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
}
static inline float WENOCoeffsY(int i, int j, int k, int DimX, int DimY, int DimZ,
						const float *x0, float inv_delta, char s){
	float q1, q2, q3, q4, q5;
	if(s < 0){
		q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
		q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
		q3 = (x0[INDEX(i,j,k)]   - x0[INDEX(i,j-1,k)]) * inv_delta;
		q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
		q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
	}
	else{
		q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
		q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
		q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
		q4 = (x0[INDEX(i,j,k)]   - x0[INDEX(i,j-1,k)]) * inv_delta;
		q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
	}
	return HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
}
static inline float WENOCoeffsZ(int i, int j, int k, int DimX, int DimY, int DimZ,
						const float *x0, float inv_delta, char s){
	float q1, q2, q3, q4, q5;
	if(s < 0){
		q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
		q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
		q3 = (x0[INDEX(i,j,k)]   - x0[INDEX(i,j,k-1)]) * inv_delta;
		q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
		q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
	}
	else{
		q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
		q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
		q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
		q4 = (x0[INDEX(i,j,k)]   - x0[INDEX(i,j,k-1)]) * inv_delta;
		q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
	}
	return HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
}


#define SOLID_MINUS   0x01
#define SOLID_PLUS    0x02
#define NOWENO_MINUS  0x04
#define NOWENO_PLUS   0x08

void ScalarField3D::FindSolidNeighbors(const Voxel &voxel,
		const char *needed, const float *x0,
		char *axis_x, char *axis_y, char *axis_z) const{

	FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);

			if(needed[pos]){

				if(x0[pos] <= 0.f){   // liquid point
					TentativeGridPoint minTentative(x0[pos], i, j, k);
					vector<TentativeGridPoint> neighbors;
					minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
					if(neighbors.size() == 6){
						if(voxel.InSolid(i-1, j, k) && !voxel.InMovingSolid(i-1, j, k))
							axis_x[pos] |= SOLID_MINUS;
						if(voxel.InSolid(i+1, j, k) && !voxel.InMovingSolid(i+1, j, k))
							axis_x[pos] |= SOLID_PLUS;
						if(voxel.InSolid(i, j-1, k) && !voxel.InMovingSolid(i, j-1, k))
							axis_y[pos] |= SOLID_MINUS;
						if(voxel.InSolid(i, j+1, k) && !voxel.InMovingSolid(i, j+1, k))
							axis_y[pos] |= SOLID_PLUS;
						if(voxel.InSolid(i, j, k-1) && !voxel.InMovingSolid(i, j, k-1))
							axis_z[pos] |= SOLID_MINUS;
						if(voxel.InSolid(i, j, k+1) && !voxel.InMovingSolid(i, j, k+1))
							axis_z[pos] |= SOLID_PLUS;
					}
					else{
						for(int m=0; m<neighbors.size(); m++){
							if(minTentative.LeftNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
									axis_x[pos] |= SOLID_MINUS;
							}
							if(minTentative.RightNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
									axis_x[pos] |= SOLID_PLUS;
							}
							if(minTentative.BackNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
									axis_y[pos] |= SOLID_MINUS;
							}
							if(minTentative.FrontNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
									axis_y[pos] |= SOLID_PLUS;
							}
							if(minTentative.BottomNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
									axis_z[pos] |= SOLID_MINUS;
							}
							if(minTentative.TopNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]) && !voxel.InMovingSolid(neighbors[m]))
									axis_z[pos] |= SOLID_PLUS;
							}
						}
					}

					if( (voxel.InSolid(i-3,j,k) && !voxel.InMovingSolid(i-3, j, k)) ||
						(voxel.InSolid(i-2,j,k) && !voxel.InMovingSolid(i-2, j, k)) ||
						(voxel.InSolid(i-1,j,k) && !voxel.InMovingSolid(i-1, j, k)) ||
						(voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1, j, k)) ||
						(voxel.InSolid(i+2,j,k) && !voxel.InMovingSolid(i+2, j, k)) )
							axis_x[pos] |= NOWENO_MINUS;

					if( (voxel.InSolid(i-2,j,k) && !voxel.InMovingSolid(i-2, j, k))	||
						(voxel.InSolid(i-1,j,k) && !voxel.InMovingSolid(i-1, j, k)) ||
						(voxel.InSolid(i+3,j,k) && !voxel.InMovingSolid(i+3, j, k)) ||
						(voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1, j, k)) ||
						(voxel.InSolid(i+2,j,k) && !voxel.InMovingSolid(i+2, j, k)) )
							axis_x[pos] |= NOWENO_PLUS;

					if( (voxel.InSolid(i,j-3,k) && !voxel.InMovingSolid(i, j-3, k)) ||
						(voxel.InSolid(i,j-2,k) && !voxel.InMovingSolid(i, j-2, k)) ||
						(voxel.InSolid(i,j-1,k) && !voxel.InMovingSolid(i, j-1, k)) ||
						(voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i, j+1, k)) ||
						(voxel.InSolid(i,j+2,k) && !voxel.InMovingSolid(i, j+2, k)) )
							axis_y[pos] |= NOWENO_MINUS;
					if( (voxel.InSolid(i,j-2,k) && !voxel.InMovingSolid(i, j-2, k)) ||
						(voxel.InSolid(i,j-1,k) && !voxel.InMovingSolid(i, j-1, k)) ||
						(voxel.InSolid(i,j+3,k) && !voxel.InMovingSolid(i, j+3, k)) ||
						(voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i, j+1, k)) ||
						(voxel.InSolid(i,j+2,k) && !voxel.InMovingSolid(i, j+2, k)) )
							axis_y[pos] |= NOWENO_PLUS;
					if( (voxel.InSolid(i,j,k-3) && !voxel.InMovingSolid(i, j, k-3)) ||
						(voxel.InSolid(i,j,k-2) && !voxel.InMovingSolid(i, j, k-2)) ||
						(voxel.InSolid(i,j,k-1) && !voxel.InMovingSolid(i, j, k-1)) ||
						(voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i, j, k+1)) ||
						(voxel.InSolid(i,j,k+2) && !voxel.InMovingSolid(i, j, k+2)) )
							axis_z[pos] |= NOWENO_MINUS;
					if( (voxel.InSolid(i,j,k-2) && !voxel.InMovingSolid(i, j, k-2)) ||
						(voxel.InSolid(i,j,k-1) && !voxel.InMovingSolid(i, j, k-1)) ||
						(voxel.InSolid(i,j,k+3) && !voxel.InMovingSolid(i, j, k+3)) ||
						(voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i, j, k+1)) ||
						(voxel.InSolid(i,j,k+2) && !voxel.InMovingSolid(i, j, k+2)) )
							axis_z[pos] |= NOWENO_PLUS;
				}
				else{  // air point
					TentativeGridPoint minTentative(x0[pos], i, j, k);
					vector<TentativeGridPoint> neighbors;
					minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
					if(neighbors.size() == 6){
						if(voxel.InSolid(i-1, j, k))
							axis_x[pos] |= SOLID_MINUS;
						if(voxel.InSolid(i+1, j, k))
							axis_x[pos] |= SOLID_PLUS;
						if(voxel.InSolid(i, j-1, k))
							axis_y[pos] |= SOLID_MINUS;
						if(voxel.InSolid(i, j+1, k))
							axis_y[pos] |= SOLID_PLUS;
						if(voxel.InSolid(i, j, k-1))
							axis_z[pos] |= SOLID_MINUS;
						if(voxel.InSolid(i, j, k+1))
							axis_z[pos] |= SOLID_PLUS;
					}
					else{
						for(int m=0; m<neighbors.size(); m++){
							if(minTentative.LeftNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]))
									axis_x[pos] |= SOLID_MINUS;
							}
							if(minTentative.RightNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]))
									axis_x[pos] |= SOLID_PLUS;
							}
							if(minTentative.BackNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]))
									axis_y[pos] |= SOLID_MINUS;
							}
							if(minTentative.FrontNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]))
									axis_y[pos] |= SOLID_PLUS;
							}
							if(minTentative.BottomNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]))
									axis_z[pos] |= SOLID_MINUS;
							}
							if(minTentative.TopNeighbor(neighbors[m])){
								if(voxel.InSolid(neighbors[m]))
									axis_z[pos] |= SOLID_PLUS;
							}
						}
					}

					if( voxel.InSolid(i-3,j,k) ||
						voxel.InSolid(i-2,j,k) ||
						voxel.InSolid(i-1,j,k) ||
						voxel.InSolid(i+1,j,k) ||
						voxel.InSolid(i+2,j,k)  )
							axis_x[pos] |= NOWENO_MINUS;
					if( voxel.InSolid(i-2,j,k)	||
						voxel.InSolid(i-1,j,k)  ||
						voxel.InSolid(i+3,j,k)  ||
						voxel.InSolid(i+1,j,k)  ||
						voxel.InSolid(i+2,j,k)  )
							axis_x[pos] |= NOWENO_PLUS;
					if( voxel.InSolid(i,j-3,k)  ||
						voxel.InSolid(i,j-2,k)  ||
						voxel.InSolid(i,j-1,k)  ||
						voxel.InSolid(i,j+1,k)  ||
						voxel.InSolid(i,j+2,k)  )
							axis_y[pos] |= NOWENO_MINUS;
					if( voxel.InSolid(i,j-2,k)  ||
						voxel.InSolid(i,j-1,k)  ||
						voxel.InSolid(i,j+3,k)  ||
						voxel.InSolid(i,j+1,k)  ||
						voxel.InSolid(i,j+2,k)  )
							axis_y[pos] |= NOWENO_PLUS;
					if( voxel.InSolid(i,j,k-3)  ||
						voxel.InSolid(i,j,k-2)  ||
						voxel.InSolid(i,j,k-1)  ||
						voxel.InSolid(i,j,k+1)  ||
						voxel.InSolid(i,j,k+2)  )
							axis_z[pos] |= NOWENO_MINUS;
					if( voxel.InSolid(i,j,k-2)  ||
						voxel.InSolid(i,j,k-1)  ||
						voxel.InSolid(i,j,k+3)  ||
						voxel.InSolid(i,j,k+1)  ||
						voxel.InSolid(i,j,k+2)  )
							axis_z[pos] |= NOWENO_PLUS;
				}
			}
	END_FOR

}


void ScalarField3D::EulerStep(const Voxel &voxel, float dtau, const char *needed,
		const float *S, const char *axis_x, const char *axis_y, const char *axis_z,
		float *value, const float *value0){


#ifdef TBB
	EulerStepBody body(&voxel, value0, S, value, needed, DimX, DimY, DimZ);
//	printf("Eulerstep body created \n");
//	tbb::parallel_for(tbb::blocked_range<int>(0, DimZ),
//			body, tbb::auto_partitioner( ));
	tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
					body);
#else
	float s_xp, s_yp, s_zp;
	float s_xm, s_ym, s_zm;
	float vel_minus_x, vel_minus_y, vel_minus_z;
	float vel_plus_x, vel_plus_y, vel_plus_z;
	float phi_x, phi_y, phi_z;
	float phi_xm, phi_ym, phi_zm;
	float phi_xp, phi_yp, phi_zp;

	float delta = voxel.VoxelDelta();

	float inv_delta = 1.f / delta;

	const float *x0 = value0;

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);

		if(needed[pos]){

			float s = S[pos];

			if( axis_x[pos] & NOWENO_MINUS ){
				if( axis_x[pos] & SOLID_MINUS )
					vel_minus_x = value0[pos] + sign(value0[pos]) * delta;
				else
					vel_minus_x = value0[INDEX(i-1,j,k)];
				phi_xm = (value0[pos] - vel_minus_x) * inv_delta ;
			}
			else
				phi_xm = WENOCoeffsX(i, j, k, DimX, DimY, DimZ, x0, inv_delta, -1);

			if( axis_x[pos] & NOWENO_PLUS ){
				if( axis_x[pos] & SOLID_PLUS )
					vel_plus_x = value0[pos] + sign(value0[pos]) * delta;
				else
					vel_plus_x = value0[INDEX(i+1,j,k)];
				phi_xp = (vel_plus_x - value0[pos]) * inv_delta;
			}
			else
				phi_xp = WENOCoeffsX(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1);

			if( axis_y[pos] & NOWENO_MINUS ){
				if( axis_y[pos] & SOLID_MINUS )
					vel_minus_y = value0[pos] + sign(value0[pos]) * delta;
				else
					vel_minus_y = value0[INDEX(i,j-1,k)];
				phi_ym = (value0[pos] - vel_minus_y) * inv_delta;
			}
			else
				phi_ym = WENOCoeffsY(i, j, k, DimX, DimY, DimZ, x0, inv_delta, -1);

			if( axis_y[pos] & NOWENO_PLUS ){
				if( axis_y[pos] & SOLID_PLUS )
					vel_plus_y = value0[pos] + sign(value0[pos]) * delta;
				else
					vel_plus_y = value0[INDEX(i,j+1,k)];
				phi_yp = (vel_plus_y - value0[pos]) * inv_delta;
			}
			else
				phi_yp = WENOCoeffsY(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1);

			if( axis_z[pos] & NOWENO_MINUS ){
				if( axis_z[pos] & SOLID_MINUS )
					vel_minus_z = value0[pos] + sign(value0[pos]) * delta;
				else
					vel_minus_z = value0[INDEX(i,j,k-1)];
				phi_zm = (value0[pos] - vel_minus_z) * inv_delta;
			}
			else
				phi_zm = WENOCoeffsZ(i, j, k, DimX, DimY, DimZ, x0, inv_delta, -1);

			if( axis_z[pos] & NOWENO_PLUS ){
				if( axis_z[pos] & SOLID_PLUS )
					vel_plus_z = value0[pos] + sign(value0[pos]) * delta;
				else
					vel_plus_z = value0[INDEX(i,j,k+1)];
				phi_zp = (vel_plus_z - value0[pos]) * inv_delta;
			}
			else
				phi_zp = WENOCoeffsZ(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1);


//				if(i == I && j == J && k == K){
//					printf("S = %f, phi_xm = %f, phi_xp = %f value0 = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
//						s, phi_xm, phi_xp, value0[pos], vel_minus_x, vel_plus_x);
//					printf("S = %f, phi_ym = %f, phi_yp = %f value0 = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
//							s, phi_ym, phi_yp, value0[pos], vel_minus_y, vel_plus_y);
//					printf("S = %f, phi_zm = %f, phi_zp = %f value0 = %f, vel_minus_z = %f, vel_plus_z = %f\n ",
//						s, phi_zm, phi_zp, value0[pos], vel_minus_z, vel_plus_z);
//				}

			s_xp = s * phi_xp;
			s_xm = s * phi_xm;
			s_yp = s * phi_yp;
			s_ym = s * phi_ym;
			s_zp = s * phi_zp;
			s_zm = s * phi_zm;

			phi_x = Godunov(s_xp, s_xm, phi_xm, phi_xp, s);
			phi_y = Godunov(s_yp, s_ym, phi_ym, phi_yp, s);
			phi_z = Godunov(s_zp, s_zm, phi_zm, phi_zp, s);

//			s = value0[pos] / sqrtf(value0[pos] * value0[pos] +
//								(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) * delta2);
//			if(fabsf(value0[pos]) < E_EPSIL)
//				s = value0[pos] / sqrtf(value0[pos] * value0[pos] + delta2);
//			float cfl = dtau / delta * (fabsf(s*phi_x)+fabsf(s*phi_y)+fabsf(s*phi_z))/sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z);
//			if(cfl > 1.f){
////			if(old_phi*value[pos]< 0.f)
//			if(i == I && j == J && k == K){
////				printf(" (reinit) phi changes sign at i=%d, j=%d, k=%d, old = %f, new = %f cfl = %f\n",
////						i, j, k, value0[pos], value[pos], cfl);
//				printf("S = %f, phi_x = %f, phi_y = %f, phi_z = %f, sqrt = %f \n ",
//						s, phi_x, phi_y, phi_z, sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z));
////				printf(" phi_xm = %f, phi_xp = %f axis_x = %d \n", phi_xm, phi_xp, axis_x[pos]);
////				printf(" phi_ym = %f, phi_yp = %f axis_y = %d \n", phi_ym, phi_yp, axis_y[pos]);
////				printf(" phi_zm = %f, phi_zp = %f axis_z = %d \n", phi_zm, phi_zp, axis_z[pos]);
////				if(axis_x[pos] & NOWENO_MINUS)
////					printf("X: NOWENO_MINUS is true \n");
////				else
////					printf("X: NOWENO_MINUS is false \n");
////				if(axis_x[pos] & NOWENO_PLUS)
////					printf("X: NOWENO_PLUS is true \n");
////				else
////					printf("X: NOWENO_PLUS is false \n");
////				if(axis_x[pos] & SOLID_MINUS)
////					printf("X: SOLID_MINUS is true \n");
////				else
////					printf("X: SOLID_MINUS is false \n");
////				if(axis_x[pos] & SOLID_PLUS)
////					printf("X: SOLID_PLUS is true \n");
////				else
////					printf("X: SOLID_PLUS is false \n");
////				if(axis_y[pos] & NOWENO_MINUS)
////					printf("Y: NOWENO_MINUS is true \n");
////				else
////					printf("Y: NOWENO_MINUS is false \n");
////				if(axis_y[pos] & NOWENO_PLUS)
////					printf("Y: NOWENO_PLUS is true \n");
////				else
////					printf("Y: NOWENO_PLUS is false \n");
////				if(axis_y[pos] & SOLID_MINUS)
////					printf("Y: SOLID_MINUS is true \n");
////				else
////					printf("Y: SOLID_MINUS is false \n");
////				if(axis_y[pos] & SOLID_PLUS)
////					printf("Y: SOLID_PLUS is true \n");
////				else
////					printf("Y: SOLID_PLUS is false \n");
////				if(axis_z[pos] & NOWENO_MINUS)
////					printf("Z: NOWENO_MINUS is true \n");
////				else
////					printf("Z: NOWENO_MINUS is false \n");
////				if(axis_z[pos] & NOWENO_PLUS)
////					printf("Z: NOWENO_PLUS is true \n");
////				else
////					printf("Z: NOWENO_PLUS is false \n");
////				if(axis_z[pos] & SOLID_MINUS)
////					printf("Z: SOLID_MINUS is true \n");
////				else
////					printf("Z: SOLID_MINUS is false \n");
////				if(axis_z[pos] & SOLID_PLUS)
////					printf("Z: SOLID_PLUS is true \n");
////				else
////					printf("Z: SOLID_PLUS is false \n");
//
////				printf(" reinit cfl > 1, exiting \n ");
////				exit(1);
//			}
			value[pos] = value0[pos] - dtau * s * ( sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) - 1.f );

//			if(old_phi*value[pos]< 0.f)
//				if(i == I && j == J && k == K){
//					printf(" (reinit) phi changes sign at i=%d, j=%d, k=%d, old = %f, new = %f \n",
//							i, j, k, value0[pos], value[pos]);
//					printf("S = %f, phi_x = %f, phi_y = %f, phi_z = %f \n ", s, phi_x, phi_y, phi_z);
//				}
//			if(isnan(value[pos])){
//				printf(" (reinit) phi changes sign at i=%d, j=%d, k=%d, old = %f, new = %f s = %f\n",
//								i, j, k, value0[pos], value[pos], s);
//				exit(1);
//			}
		}
	END_FOR
#endif

}

void ScalarField3D::ReInitialize(Voxel &voxel, const ScalarField3D &solid,
							bool initial, float *value, float *value0){

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float delta3 = delta2 * delta;
	float inv_delta = 1.f / delta;
	float dtau = delta / 10;
	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;

	float eps = 1.e-4f * delta;

	int iterations = int(round(FASTMARCH_LIMIT / dtau));

	// obtain the solid levelset
	//float *phi_obj = solid.getScalar();

	char *needed = new char[DimX*DimY*DimZ];
	//	memset(valid, 0, DimX*DimY*DimZ);
//	float *S = new float[DimX*DimY*DimZ];

	SetZero(needed);
//	SetZero(S);

	u_int N = 0;

	float old_phi = value[INDEX(I, J, K)];
#ifdef CUDA
//	printf("A. iteraations = %d \n", iterations);


	char *source = needed;
	FOR_EACH_CELL
			if(voxel.InSolid(i,j,k))
				obj[INDEX(i,j,k)] = 1;
			else
				obj[INDEX(i,j,k)] = 0;
			if(voxel.InMovingSolid(i,j,k))
				mov_obj[INDEX(i,j,k)] = 1;
			else
				mov_obj[INDEX(i,j,k)] = 0;
			if(voxel.InSource(i,j,k))
				source[INDEX(i,j,k)] = 1;
			else
				source[INDEX(i,j,k)] = 0;
    END_FOR
//    printf("B. iteraations = %d \n", iterations);
//    int dims[3];
//    dims[0] = DimX; dims[1] = DimY; dims[2] = DimZ;

    ReinitializeCUDA(value, obj, mov_obj, source, initial, iterations, I, J, K,
    					delta, dtau, eps, FASTMARCH_LIMIT, dims);


//	float *cudaphi = new float[DimX*DimY*DimZ];
//	SetEqual(value, cudaphi);
//	SetEqual(value0, value);
//	u_int N9 = 0;
//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && cudaphi[INDEX(i,j,k)] <= 0.f)
//			++N9;
//	END_FOR
//	printf("\nafter cuda reinit, there are %u liquid points \n\n", N9);

#else
	if(initial){
		// if called for the first time, we do reinitializtion over
		// the whole domain except for the solid points

		FOR_EACH_CELL
			if( !voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) ){
				needed[INDEX(i,j,k)] = 1;
				++N;
			}
		END_FOR
	}
	else{
		FOR_EACH_CELL
		if( !voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)
			&& NearInterface(i, j, k, eps, FASTMARCH_LIMIT, value0) ){
			needed[INDEX(i,j,k)] = 1;
			++N;
		}
		END_FOR
	}

	printf("Grid cells needing reinitialization determined ! \n");
	printf("dtau = %f, iterations = %d \n", dtau, iterations);

	SetEqual(value, phiPos);
	SetEqual(value, phiNeg);
	float *S = new float[DimX*DimY*DimZ];
	SetEqual(value, S);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(needed[pos])
			S[pos] = value[pos]/sqrtf(value[pos] * value[pos] + delta2);
	END_FOR
	char *axis_x = new char[DimX*DimY*DimZ];
	memset(axis_x, 0, DimX*DimY*DimZ*sizeof(char));
	char *axis_y = new char[DimX*DimY*DimZ];
	memset(axis_y, 0, DimX*DimY*DimZ*sizeof(char));
	char *axis_z = new char[DimX*DimY*DimZ];
	memset(axis_z, 0, DimX*DimY*DimZ*sizeof(char));
	FindSolidNeighbors(voxel, needed, value, axis_x, axis_y, axis_z);

#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif

	for(int n = 0; n < iterations; ++n){

		SWAP(value, value0);
//		FindSmoothFactor(voxel, needed, value0, S);
		EulerStep(voxel, dtau, needed, S, axis_x, axis_y, axis_z, phiPos, value0);
		EulerStep(voxel, dtau, needed, S, axis_x, axis_y, axis_z, phiNeg, phiPos);
#ifdef TBB
		AddTwoBody body1(needed, value0, phiNeg, phiPos, 0.75f, 0.25f, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
			body1);
#else
//		FOR_EACH_CELL
//			u_int pos = INDEX(i,j,k);
//			if(needed[pos])
//				phiPos[pos] = 0.75f*value0[pos] + 0.25f * phiNeg[pos];
//		END_FOR
#endif
//		EulerStep(voxel, dtau, needed, S, phiNeg, phiPos);
#ifdef TBB
		AddTwoBody body2(needed, value0, phiNeg, value, one_third, two_third, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
			body2);
#else
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
//				value[pos] = one_third * value0[pos] + two_third * phiNeg[pos];
				value[pos] = 0.5f  * (value0[pos] + phiNeg[pos]);
		END_FOR
#endif

//		printf(" (reinit) phi changes sign at iter = %d i=%d, j=%d, k=%d, "
//				"needed = %d, old = %f, new = %f \n",
//			n, I, J, K, needed[INDEX(I,J,K)], value0[INDEX(I,J,K)], value[INDEX(I,J,K)]);
	}

	delete [] S;
	delete [] axis_x;
	delete [] axis_y;
	delete [] axis_z;

//	FOR_EACH_CELL
////		if(i == I && j == J && k == K){
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) &&
//			value[INDEX(i,j,k)] <= 0.f && cudaphi[INDEX(i,j,k)] > 0.f){
//			printf("at (%d, %d, %d) cpu_phi = %f, cudaphi = %f \n",
//				 i, j, k, value[INDEX(i,j,k)], cudaphi[INDEX(i,j,k)]);
//		}
//	END_FOR
//
//	delete [] cudaphi;

#endif
	delete [] needed;


	// check reinitialization results
//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//			if(fabs(value[INDEX(i,j,k)]) < delta){
//				Point p = voxel.VoxelCenterPosition(i,j,k);
//				Vector N;
//				GradPhi(voxel, p, N);
//				printf("at (%d, %d, %d) phi = %f, grad(phi) = %f \n", i,j,k, value[INDEX(i,j,k)], N.Length());
//			}
//		}
//	END_FOR


	if(!initial){
		u_int N1 = 0;
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && fabs(value[INDEX(i,j,k)]) > FASTMARCH_LIMIT){
				value[INDEX(i,j,k)] = value[INDEX(i,j,k)] > 0.f ? FASTMARCH_LIMIT : -FASTMARCH_LIMIT;
				++N1;
			}
			else if(!voxel.InSolid(i,j,k) && FASTMARCH_LIMIT - fabs(value[INDEX(i,j,k)]) < delta){
				value[INDEX(i,j,k)] = value[INDEX(i,j,k)] > 0.f ? FASTMARCH_LIMIT : -FASTMARCH_LIMIT;
				++N1;
			}
		END_FOR
		printf("there are %u points set back to  FASTMARCH_LIMIT\n", N1);
	}

	printf("\nthere are %u points reinitialized \n\n", N);

}

#if 0
void ScalarField3D::ReInitializeForRendering(const Voxel &voxel, float *value, float *value0){

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float delta3 = delta2 * delta;
	float inv_delta = 1.f / delta;
	float dtau = delta / 2 ;
	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;

	int iterations = 10;

	float s_xp, s_yp, s_zp;
	float s_xm, s_ym, s_zm;
	float vel_minus_x, vel_minus_y, vel_minus_z;
	float vel_plus_x, vel_plus_y, vel_plus_z;
	float phi_x, phi_y, phi_z;
	float phi_xm, phi_ym, phi_zm;
	float phi_xp, phi_yp, phi_zp;


	char *needed = new char[DimX*DimY*DimZ];
	//	memset(valid, 0, DimX*DimY*DimZ);
//	float *S = new float[DimX*DimY*DimZ];

	SetZero(needed);
//	SetZero(S);

	u_int N = 0;

	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
		int n;
		if(voxel.InSolid(i,j,k)){
			if(voxel.CloseToNonSolid(tp, n) && fabs(value[INDEX(i,j,k)]) < 2 * delta){
				needed[INDEX(i,j,k)] = 1;
				++N;
			}
		}
		else if(voxel.CloseToSolid(tp) && fabs(value[INDEX(i,j,k)]) < 2 * delta){
			needed[INDEX(i,j,k)] = 1;
			++N;
		}
	END_FOR

	SetEqual(value, phiPos);
	SetEqual(value, phiNeg);

	float old_phi = value[INDEX(I, J, K)];
#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif
	for(int n = 0; n < iterations; ++n){

		SWAP(value, value0);
//		FindSmoothFactor(voxel, needed, value0, S);
		EulerStep(voxel, dtau, needed, phiPos, value0);
		EulerStep(voxel, dtau, needed, phiNeg, phiPos);

#ifdef TBB
		AddTwoBody body1(needed, value0, phiNeg, phiPos, 0.75f, 0.25f, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ),
			body1, tbb::auto_partitioner( ));
#else
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
				phiPos[pos] = 0.75f*value0[pos] + 0.25f * phiNeg[pos];
		END_FOR
#endif
		EulerStep(voxel, dtau, needed, phiNeg, phiPos);
#ifdef TBB
		AddTwoBody body2(needed, value0, phiNeg, value, one_third, two_third, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ),
			body2, tbb::auto_partitioner( ));
#else
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
				value[pos] = one_third*value0[pos] + two_third * phiNeg[pos];
		END_FOR
#endif

	}

	delete [] needed;
//	delete [] S;

	printf("\nthere are %u points reinitialized \n\n", N);

}
#endif

#endif


float ScalarField3D::
	UpdateTentativeValue(const KnownPoint &p, u_int index,
						 map<u_int, KnownPoint> &band,
						 const Voxel &v, float *phi_tmp){
	//float tatitiveValue;
	float t[6];
	bool axis[3]={true, true, true};
	for(int i=0;i<6;i++){
		u_int indadj;
		if(i==0) {	 // 0,1 ==> x-axis
			if(p.ii==0){
				t[0] = INFINITY;
				continue;
			}
		 	indadj = p.kk*(DimX*DimY) +
				     p.jj*DimX +
				    (p.ii - 1);
		}
		else if(i==1) {	 // 0,1 ==> x-axis
			if(p.ii==DimX-1){
				t[1] = INFINITY;
				continue;
			}
		 	indadj = p.kk*(DimX*DimY) +
				     p.jj*DimX +
				    (p.ii + 1);
		}
		else if(i==4){ // 4,5 ==> z-axis
			if(p.kk==0){
				t[4] = INFINITY;
				continue;
			}
		    indadj = (p.kk-1)*(DimX*DimY) +
		   			 p.jj*DimX +
		   			 p.ii;
//		   	printf("i=%d, p.kk = %d, p.kk+-1 = %d \n", i, p.kk, p.kk+((i%2==0) ? (-1) : 1) );
		}
		else if(i==5){ // 4,5 ==> z-axis
			if(p.kk==DimZ-1){
				t[5] = INFINITY;
				continue;
			}
		    indadj = (p.kk+1)*(DimX*DimY) +
		   			 p.jj*DimX +
		   			 p.ii;
//		   	printf("i=%d, p.kk = %d, p.kk+-1 = %d \n", i, p.kk, p.kk+((i%2==0) ? (-1) : 1) );
		}
		else if(i==2){         // 2,3 ==> y-axis
			if(p.jj==0){
				t[2] = INFINITY;
				continue;
			}
			indadj = p.kk*(DimX*DimY) +
		   			 (p.jj-1)*DimX +
		   			 p.ii;
		}
		else{         // 2,3 ==> y-axis
			if(p.jj==DimY-1){
				t[3] = INFINITY;
				continue;
			}
		    indadj = p.kk*(DimX*DimY) +
		   			 (p.jj+1)*DimX +
		   			 p.ii;
		}
//		printf("trying to find in band i=%d ii=%d, jj=%d, kk=%d, %d \n",
//					i,	p.ii, p.jj, p.kk, indadj);
		map<u_int, KnownPoint>::iterator found = band.find(indadj);
		if(found != band.end()){
//		if(v.IsDone(indadj)){
//			printf("found in band %d \n", indadj);
			t[i] = phi_tmp[indadj];
			if(t[i] < 0.f)
				printf("i = %d, ii=%d, jj=%d, kk=%d, t[%d]=%f \n",
						i, p.ii, p.jj, p.kk, i, t[i]);
		}
		else
			t[i] = INFINITY;
	}
	float phi1 = min(t[0], t[1]);
	float phi2 = min(t[2], t[3]);
	float phi3 = min(t[4], t[5]);
	if(phi1 == INFINITY)
		axis[0] = false;
	if(phi2 == INFINITY)
		axis[1] = false;
	if(phi3 == INFINITY)
		axis[2] = false;

//	if(p.ii == I && p.jj == J && p.kk == K)
//		printf("p.ii = %d, p.jj = %d, p.kk = %d ,   "
//				" t0=%f, t1=%f, t2=%f, t3=%f, t4=%f, t5=%f \n",
//			p.ii, p.jj, p.kk, t[0], t[1], t[2], t[3], t[4], t[5]);

	if(axis[0] && axis[1] && axis[2]){ // 3 non-zero terms
		if(t[0] <= t[1]){     // quardrants: 0,1,2,3
			if(t[2] <= t[3]) { // quardrants: 0,2
				if(t[4] <= t[5]){ // quardrants: 0
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[2], t[4], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//			 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 2
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[2], t[5], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//								 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
			}
			else{               // quardrants: 1,3
				if(t[4] <= t[5]){ // quardrants: 1
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[3], t[4], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 3
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[3], t[5], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
			}
		}
		else{			// quardrants: 4,5,6,7
			if(t[2] <= t[3]) { // quardrants: 4,6
				if(t[4] <= t[5]){ // quardrants: 4
					float solution = ProcessThreeTerms(0, 1, 2, t[1], t[2], t[4], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 6
					float solution =  ProcessThreeTerms(0, 1, 2, t[1], t[2], t[5], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
			}
			else{               // quardrants: 5,7
				if(t[4] <= t[5]){ // quardrants: 5
					float solution = ProcessThreeTerms(0, 1, 2, t[1], t[3], t[4], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 7
					float solution =  ProcessThreeTerms(0, 1, 2, t[1], t[3], t[5], v);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//								t[0], t[1], t[2], t[3], t[4], t[5]);
				return solution;
				}
			}
		}
	}
	else if(axis[0] && axis[1] && !axis[2]){ // 2 non-zero terms
		if( t[0] <= t[1]){
			if(t[2] <= t[3])
				return ProcessTwoTerms(0, 1, t[0], t[2], v);
			else
				return ProcessTwoTerms(0, 1, t[0], t[3], v);
		}
		else{
			if(t[2] <= t[3])
				return ProcessTwoTerms(0, 1, t[1], t[2], v);
			else
				return ProcessTwoTerms(0, 1, t[1], t[3], v);
		}
	}
	else if(axis[0] && !axis[1] && axis[2]){ // 2 non-zero terms
		if( t[0] <= t[1]){
			if(t[4] <= t[5])
				return ProcessTwoTerms(0, 2, t[0], t[4], v);
			else
				return ProcessTwoTerms(0, 2, t[0], t[5], v);
		}
		else{
			if(t[4] <= t[5])
				return ProcessTwoTerms(0, 2, t[1], t[4], v);
			else
				return ProcessTwoTerms(0, 2, t[1], t[5], v);
		}
	}
	else if(!axis[0] && axis[1] && axis[2]){ // 2 non-zero terms
		if( t[2] <= t[3]){
			if(t[4] <= t[5])
				return ProcessTwoTerms(1, 2, t[2], t[4], v);
			else
				return ProcessTwoTerms(1, 2, t[2], t[5], v);
		}
		else{
			if(t[4] <= t[5])
				return ProcessTwoTerms(1, 2, t[3], t[4], v);
			else
				return ProcessTwoTerms(1, 2, t[3], t[5], v);
		}
	}
	else if(axis[0] && !axis[1] && !axis[2]){ // 1 non-zero term
		if(t[0] <= t[1])
			return ProcessOneTerm(0, t[0], v);
		else
			return ProcessOneTerm(0, t[1], v);
	}
	else if(!axis[0] && axis[1] && !axis[2]){ // 1 non-zero term
		if(t[2] <= t[3])
			return ProcessOneTerm(1, t[2], v);
		else
			return ProcessOneTerm(1, t[3], v);
	}
	else if(!axis[0] && !axis[1] && axis[2]){ // 1 non-zero term
		if(t[4] <= t[5])
			return ProcessOneTerm(2, t[4], v);
		else
			return ProcessOneTerm(2, t[5], v);
	}
	else
		return INFINITY;

}

float  ScalarField3D::ProcessOneTerm(int axis, float f, const Voxel &v){
	float delta = v.VoxelDelta();
	if(f + delta < 0.f)
		printf("Wrong! phi < 0 occurs in One Term, phi = %f \n", f + delta);
	return f + delta;
}

float  ScalarField3D::ProcessTwoTerms(int axis1, int axis2,
									  float phi1, float phi2,
									  const Voxel &v){
	float delta1 = v.VoxelDelta();
	float delta2 = v.VoxelDelta();
	float phimax = max(phi1, phi2);
	float phi = phi1 <= phi2 ? phi1 : phi2;
	int axis = phi1 <= phi2 ? axis1 : axis2;
	float P = ((phimax - phi1)/delta1) * ((phimax - phi1)/delta1) +
	 		  ((phimax - phi2)/delta2) * ((phimax - phi2)/delta2);
	if(P > 1)
	 	return ProcessOneTerm(axis, phi, v);
	else{
		float A, B, C, t0, t1;
	 	A = 1.f/(delta1*delta1) + 1.f/(delta2*delta2);
	 	B = -2*(phi1/(delta1*delta1) + phi2/(delta2*delta2));
	 	C = phi1*phi1/(delta1*delta1) +
	 		phi2*phi2/(delta2*delta2) - 1.f;

	 	if(Quadratic(A, B, C, &t0, &t1)){
//	 		if(t1 < 0.f)
//	 			printf("Wrong! phi < 0 occurs in Two Terms t0 = %f, t1 = %f, phi1 = %f, phi2 = %f \n",
//	 					t0, t1, phi1, phi2);
	 		return t1;
	 	}
	 	else
	 		return INFINITY;
	}
}

float  ScalarField3D::ProcessThreeTerms(int axis1, int axis2, int axis3,
									    float phi1, float phi2, float phi3,
									    const Voxel &v){
	float delta1 = v.VoxelDelta();
	float delta2 = v.VoxelDelta();
	float delta3 = v.VoxelDelta();
	float phimax = max(phi1, phi2);
	phimax = max(phimax, phi3);
	float P = ((phimax - phi1)/delta1) * ((phimax - phi1)/delta1) +
	 		  ((phimax - phi2)/delta2) * ((phimax - phi2)/delta2) +
	 		  ((phimax - phi3)/delta3) * ((phimax - phi3)/delta3);
	if(P > 1){
		if(fabsf(phi1-phimax) < E_EPSIL)
	 		return ProcessTwoTerms(axis2, axis3, phi2, phi3, v);
	 	else if(fabsf(phi2-phimax) < E_EPSIL)
	 		return ProcessTwoTerms(axis1, axis3, phi1, phi3, v);
	 	else
	 		return ProcessTwoTerms(axis1, axis2, phi1, phi2, v);
	}
	else{
		float A, B, C, t0, t1;
	 	A = delta2*delta2*delta3*delta3 +
	 	    delta1*delta1*delta3*delta3 +
	 	    delta1*delta1*delta2*delta2;
	 	B = -2*(phi1*delta2*delta2*delta3*delta3
	 	      + phi2*delta1*delta1*delta3*delta3
	 	      + phi3*delta1*delta1*delta2*delta2);
	 	C = phi1*phi1*delta2*delta2*delta3*delta3 +
	 		phi2*phi2*delta1*delta1*delta3*delta3 +
	 		phi3*phi3*delta1*delta1*delta2*delta2 -
	 		delta1*delta1*delta2*delta2*delta3*delta3 ;
	 	if(Quadratic(A, B, C, &t0, &t1)){
//	 		if(t1 < 0.f)
//	 			printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, phi1 = %f, phi2 = %f, phi3 = %f \n",
//	 				 					t0, t1, phi1, phi2, phi3);

	 		return t1;
	 	}
	 	else
	 		return INFINITY;
	}
}

