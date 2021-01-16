
#include <queue>
using std::queue;

#include "vectorfield3d.h"
#include "mtrand.h"
#include "pcg_solver.h"
#include "container.h"
#include "wallclocktime.h"
#include "physicsworld.h"

#ifdef CUDA
#include "reinit_cuda.h"
#endif

#define EPSIL 0.5*delta
#define POS(x) (x).kk*(DimX*DimY)+(x).jj*DimX+(x).ii

#define AIRUVELSET 0x01
#define AIRVVELSET 0x02
#define AIRWVELSET 0x04



//Vector VectorField3D::TriInterp(const Voxel &Vox, const Point &p,
//								float deltax, float deltay, float deltaz) const{
//	float x = p.x/deltax;
//	float y = p.y/deltay;
//	float z = p.z/deltaz;
//	int xi = floor(x);
//	int yi = floor(y);
//	int zi = floor(z);
//	float u = x - xi;
//	float v = y - yi;
//	float w = z - zi;
//
//	float vel_u = lerp( U[zi*(DimX*DimY)+yi*DimX+xi].x,
//			        U[zi*(DimX*DimY)+yi*DimX+(xi+1)].x,
//			        u);
//	float vel_v = lerp(U[zi*(DimX*DimY)+yi*DimX+xi].y, U[zi*(DimX*DimY)+(yi+1)*DimX+xi].y, v);
//	float vel_w = lerp(U[zi*(DimX*DimY)+yi*DimX+xi].z, U[(zi+1)*(DimX*DimY)+yi*DimX+xi].z, w);
//	return Vector(vel_u, vel_v, vel_w);
//}

void VectorField3D::OutputBindaryData(const Voxel &voxel, char *filename) const{
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
			for(int i=0;i<DimX;i++){
				fwrite(&u[INDEX(i,j,k)], sizeof(float), 1, fp);
				fwrite(&v[INDEX(i,j,k)], sizeof(float), 1, fp);
				fwrite(&w[INDEX(i,j,k)], sizeof(float), 1, fp);
			}
//				fprintf(fp,"%f\n", phi[k*(DimX*DimY)+j*DimX+i]);
	fclose(fp);

}

void VectorField3D::WriteRestart(FILE* fp) const{

	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++)
				fprintf(fp,"%f ", u[INDEX(i,j,k)]);
			fprintf(fp,"\n");
		}
	}
	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++)
				fprintf(fp,"%f ", v[INDEX(i,j,k)]);
			fprintf(fp,"\n");
		}
	}
	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++)
				fprintf(fp,"%f ", w[INDEX(i,j,k)]);
			fprintf(fp,"\n");
		}
	}
}


bool VectorField3D::ReadRestart(FILE* fp) {
	char c;

	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++){
				if(1!=fscanf(fp, "%f", u+INDEX(i,j,k))){
			         printf("Problem reading gridimplicit file \n");
			         DimX=DimY=DimZ=0;
			         return false;
			      }
				fscanf(fp, "%c", &c);
			}
		}
	}
	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++){
				if(1!=fscanf(fp, "%f", v+INDEX(i,j,k))){
			         printf("Problem reading gridimplicit file \n");
			         DimX=DimY=DimZ=0;
			         return false;
			      }
				fscanf(fp, "%c", &c);
			}
		}
	}
	for(int k=0;k<DimZ;k++){
		for(int j=0;j<DimY;j++){
			for(int i=0;i<DimX;i++){
				if(1!=fscanf(fp, "%f", w+INDEX(i,j,k))){
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


float VectorField3D::TriInterp(const Voxel &voxel, const Point &p, float *phi) const{

	float delta = voxel.VoxelDelta();

	float x = p.x/delta - 0.5f;
    float y = p.y/delta - 0.5f;
    float z = p.z/delta - 0.5f;
    if(x < 0) x = 0.f;
    if(y < 0) y = 0.f;
    if(z < 0) z = 0.f;
	int xi = floor(x);
	int yi = floor(y);
	int zi = floor(z);
//	printf(" xi = %d, yi = %d, zi = %d \n", xi, yi, zi);
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

//	printf(" r = %f,  s = %f,  t = %f \n",  r, s, t);
	return trilerp(phi[INDEX(xi, yi, zi)], phi[INDEX(xj, yi, zi)],
				   phi[INDEX(xi, yi, zj)], phi[INDEX(xj, yi, zj)],
				   phi[INDEX(xi, yj, zi)], phi[INDEX(xj, yj, zi)],
				   phi[INDEX(xi, yj, zj)], phi[INDEX(xj, yj, zj)],
				   r, s, t);
}

float VectorField3D::TriInterp(const Voxel &voxel, const Point &p, float *phi,
				unsigned char index) const{

	float delta = voxel.VoxelDelta();
	float x, y, z;
	if(index == 1){
		x = p.x/delta - 1.f;
		y = p.y/delta - 0.5f;
		z = p.z/delta - 0.5f;
	}
	else if(index == 2){
		x = p.x/delta - 0.5f;
		y = p.y/delta - 1.f;
		z = p.z/delta - 0.5f;
	}
	else if(index == 3){
		x = p.x/delta - 0.5f;
		y = p.y/delta - 0.5f;
		z = p.z/delta - 1.f;
	}
    if(x < 0) x = 0.f;
    if(y < 0) y = 0.f;
    if(z < 0) z = 0.f;
	int xi = floor(x);
	int yi = floor(y);
	int zi = floor(z);
//	printf(" xi = %d, yi = %d, zi = %d \n", xi, yi, zi);
	//Point point = voxel.VoxelCenterPosition(xi,yi,zi);
//	Point point = voxel.VoxelCornerPosition(xi,yi,zi,1);
//	float r = x - point.x/delta;
//    float s = y - point.y/delta;
//    float t = z - point.z/delta;
	float r = x - xi;
	float s = y - yi;
	float t = z - zi;
//	printf(" r = %f,  s = %f,  t = %f \n",  r, s, t);
	return trilerp(phi[zi*(DimX*DimY)+yi*DimX+xi], phi[zi*(DimX*DimY)+yi*DimX+xi+1],
				   phi[(zi+1)*(DimX*DimY)+yi*DimX+xi], phi[(zi+1)*(DimX*DimY)+yi*DimX+xi+1],
				   phi[zi*(DimX*DimY)+(yi+1)*DimX+xi], phi[zi*(DimX*DimY)+(yi+1)*DimX+xi+1],
				   phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi], phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi+1],
				   r, s, t);
}

void VectorField3D::AddSource(const Voxel &voxel){
	//int I = 20, J = 20, K = 10;
//	FOR_EACH_CELL
//		if(i == I & j == J && k == K){
//		   printf("(Addsource) Before addsource at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
//			 			I,J,K, u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)]);
//		   printf("(Addsource) Before addsource at i = %d, j = %d, k = %d  "
//				   "u0 = %f, v0 = %f, w0 = %f \n",
//				   I, J, K, u0[INDEX(i,j,k)], v0[INDEX(i,j,k)], w0[INDEX(i,j,k)]);
//
//		}
//	END_FOR
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//		u[INDEX(i,j,k)] += dt*0.5f*(u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
//
//		v[INDEX(i,j,k)] += dt*0.5f*(v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
//
//		w[INDEX(i,j,k)] += dt*0.5f*(w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
		u[pos] += dt*u0[pos];

		v[pos] += dt*v0[pos];

		w[pos] += dt*w0[pos];
	END_FOR
//	SetSolidBoundaryNoExtrapolate(voxel,1,u);
//	SetSolidBoundaryNoExtrapolate(voxel,2,v);
//	SetSolidBoundaryNoExtrapolate(voxel,3,w);

//	FOR_EACH_CELL
//	if(i == I & j == J && k == K){
//	   printf("(Addsource) after addsource: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
//		 			I,J,K, u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)]);
//	   printf("(Addsource) after addsource: at i = %d, j = %d, k = %d  "
//			   "u0 = %f, v0 = %f, w0 = %f \n",
//			   I, J, K, u0[INDEX(i,j,k)], v0[INDEX(i,j,k)], w0[INDEX(i,j,k)]);
//
//	}
//	END_FOR
//	SetBoundary(voxel, 1, u);
//	SetBoundary(voxel, 2, v);
//	SetBoundary(voxel, 3, w);
//	Extrapolate(voxel);
//	FOR_EACH_CELL
//	if(i == I & j == J && k == K){
//	   printf("(Addsource) After Extrapolation: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
//		 			I,J,K, u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)],
//		 			u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)]);
//	}
//	END_FOR
}

#ifdef SPMD
static inline float Godunov(float s_p, float s_m, float phi_m, float phi_p,
 					   float S){
 	if(s_p >= 0.f && s_m >= 0.f){
 		return phi_m;
 	}
 	else if(s_p <= 0.f && s_m <= 0.f){
 		return phi_p;
 	}
 	else if(s_p > 0.f && s_m < 0.f)
 		return 0.f;
 	else if(s_p < 0.f && s_m > 0.f){
 		float s = S *(fabsf(phi_p) - fabsf(phi_m))/(phi_p - phi_m);
 		return s > 0.f ? phi_m : phi_p;
 	}
 	else{
 		printf("Error! Shouldn't reach here \n");
 		exit(1);
 	}

}

static float HJ_WENO_Coefficients(float q1, float q2, float q3, float q4, float q5);

void VectorField3D::FindSmoothFactor(const Voxel &voxel, char *needed, char *valid, float *value0, float *S){

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
				if(!valid[INDEX(i-1,j,k)])
					solid_xm = true;
				vel_plus_x = value0[INDEX(i+1, j, k)];
				if(!valid[INDEX(i+1,j,k)])
					solid_xp = true;
				vel_minus_y = value0[INDEX(i, j-1, k)];
				if(!valid[INDEX(i,j-1,k)])
					solid_ym = true;
				vel_plus_y = value0[INDEX(i, j+1, k)];
				if(!valid[INDEX(i,j+1,k)])
					solid_yp = true;
				vel_minus_z = value0[INDEX(i, j, k-1)];
				if(!valid[INDEX(i,j,k-1)])
					solid_zm = true;
				vel_plus_z = value0[INDEX(i, j, k+1)];
				if(!valid[INDEX(i,j,k+1)])
					solid_zp = true;
			}
			else{
				for(int m=0; m<neighbors.size(); m++){
					if(minTentative.LeftNeighbor(neighbors[m])){
						vel_minus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_xm = true;
					}
					if(minTentative.RightNeighbor(neighbors[m])){
						vel_plus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_xp = true;
					}
					if(minTentative.BackNeighbor(neighbors[m])){
						vel_minus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_ym = true;
					}
					if(minTentative.FrontNeighbor(neighbors[m])){
						vel_plus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_yp = true;
					}
					if(minTentative.BottomNeighbor(neighbors[m])){
						vel_minus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_zm = true;
					}
					if(minTentative.TopNeighbor(neighbors[m])){
						vel_plus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_zp = true;
					}
				}
			}

//				float S = value0[pos] / sqrtf(value0[pos] * value0[pos] + delta2);
			float s = value0[pos];

			if( !valid[INDEX(i-3,j,k)] || !valid[INDEX(i-2,j,k)] ||
				!valid[INDEX(i-1,j,k)] || !valid[INDEX(i+1,j,k)] ||
				!valid[INDEX(i+2,j,k)] ){
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

			if( !valid[INDEX(i+3,j,k)] || !valid[INDEX(i+2,j,k)] ||
				!valid[INDEX(i+1,j,k)] || !valid[INDEX(i-1,j,k)] ||
				!valid[INDEX(i-2,j,k)] ){
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

			if( !valid[INDEX(i,j-3,k)] || !valid[INDEX(i,j-2,k)] ||
				!valid[INDEX(i,j-1,k)] || !valid[INDEX(i,j+1,k)] ||
				!valid[INDEX(i,j+2,k)] ){
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
			if( !valid[INDEX(i,j+3,k)] || !valid[INDEX(i,j+2,k)] ||
				!valid[INDEX(i,j+1,k)] || !valid[INDEX(i,j-1,k)] ||
				!valid[INDEX(i,j-2,k)] ){
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

			if( !valid[INDEX(i,j,k-3)] || !valid[INDEX(i,j,k-2)] ||
				!valid[INDEX(i,j,k-1)] || !valid[INDEX(i,j,k+1)] ||
				!valid[INDEX(i,j,k+2)]  ){
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
			if( !valid[INDEX(i,j,k+3)] || !valid[INDEX(i,j,k+2)] ||
				!valid[INDEX(i,j,k+1)] || !valid[INDEX(i,j,k-1)] ||
				!valid[INDEX(i,j,k-2)]  ){
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

void VectorField3D::EulerStep(const Voxel &voxel, char *needed, char *valid, float *S, float *value, float *value0){

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
				if(!valid[INDEX(i-1,j,k)])
					solid_xm = true;
				vel_plus_x = value0[INDEX(i+1, j, k)];
				if(!valid[INDEX(i+1,j,k)])
					solid_xp = true;
				vel_minus_y = value0[INDEX(i, j-1, k)];
				if(!valid[INDEX(i,j-1,k)])
					solid_ym = true;
				vel_plus_y = value0[INDEX(i, j+1, k)];
				if(!valid[INDEX(i,j+1,k)])
					solid_yp = true;
				vel_minus_z = value0[INDEX(i, j, k-1)];
				if(!valid[INDEX(i,j,k-1)])
					solid_zm = true;
				vel_plus_z = value0[INDEX(i, j, k+1)];
				if(!valid[INDEX(i,j,k+1)])
					solid_zp = true;
			}
			else{
				for(int m=0; m<neighbors.size(); m++){
					if(minTentative.LeftNeighbor(neighbors[m])){
						vel_minus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_xm = true;
					}
					if(minTentative.RightNeighbor(neighbors[m])){
						vel_plus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_xp = true;
					}
					if(minTentative.BackNeighbor(neighbors[m])){
						vel_minus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_ym = true;
					}
					if(minTentative.FrontNeighbor(neighbors[m])){
						vel_plus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_yp = true;
					}
					if(minTentative.BottomNeighbor(neighbors[m])){
						vel_minus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_zm = true;
					}
					if(minTentative.TopNeighbor(neighbors[m])){
						vel_plus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[Index(neighbors[m])])
							solid_zp = true;
					}
				}
			}

//				float S = value0[pos] / sqrtf(value0[pos] * value0[pos] + delta2);
			float s = S[pos];

			if( !valid[INDEX(i-3,j,k)] || !valid[INDEX(i-2,j,k)] ||
				!valid[INDEX(i-1,j,k)] || !valid[INDEX(i+1,j,k)] ||
				!valid[INDEX(i+2,j,k)] ){
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

			if( !valid[INDEX(i+3,j,k)] || !valid[INDEX(i+2,j,k)] ||
				!valid[INDEX(i+1,j,k)] || !valid[INDEX(i-1,j,k)] ||
				!valid[INDEX(i-2,j,k)] ){
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

			if( !valid[INDEX(i,j-3,k)] || !valid[INDEX(i,j-2,k)] ||
				!valid[INDEX(i,j-1,k)] || !valid[INDEX(i,j+1,k)] ||
				!valid[INDEX(i,j+2,k)] ){
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
			if( !valid[INDEX(i,j+3,k)] || !valid[INDEX(i,j+2,k)] ||
				!valid[INDEX(i,j+1,k)] || !valid[INDEX(i,j-1,k)] ||
				!valid[INDEX(i,j-2,k)] ){
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

			if( !valid[INDEX(i,j,k-3)] || !valid[INDEX(i,j,k-2)] ||
				!valid[INDEX(i,j,k-1)] || !valid[INDEX(i,j,k+1)] ||
				!valid[INDEX(i,j,k+2)]  ){
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
			if( !valid[INDEX(i,j,k+3)] || !valid[INDEX(i,j,k+2)] ||
				!valid[INDEX(i,j,k+1)] || !valid[INDEX(i,j,k-1)] ||
				!valid[INDEX(i,j,k-2)]  ){
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
	END_FOR
#endif
}

void VectorField3D::ReInitialize(const Voxel &voxel, int vel_index, float *value, float *value0){

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float delta3 = delta2 * delta;
	float inv_delta = 1.f / delta;
	float dtau = delta / 2 ;
	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;


	float eps = 1.e-4f * delta;

	int iterations = int(round(EXTRAPOLATE_VEL_LIMIT / delta * 2));


	char *needed = needRecomputePhi;
	char *valid = new char[DimX*DimY*DimZ];
	//	memset(valid, 0, DimX*DimY*DimZ);
	float *S = w0;
	SetZero(S);

	SetZero(needed);
	SetZero(valid);
	SetValidPoint(voxel, valid, vel_index);

	u_int N = 0;

	SetEqual(value, value0);

	FOR_EACH_CELL
		if(valid[INDEX(i,j,k)]){
			if( value0[INDEX(i,j,k)] >= 0.f && value0[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ){
				needed[INDEX(i,j,k)] = 1;
				++N;
			}
			if( value0[INDEX(i,j,k)] < 0.f && fabs(value0[INDEX(i,j,k)]) < 2 * delta ){
				needed[INDEX(i,j,k)] = 1;
				++N;
			}
		}
	END_FOR

	float old_phi = value[INDEX(I, J, K)];

	SetEqual(value, u0);
	SetEqual(value, v0);

#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif

	for(int n = 0; n < iterations; ++n){

		SWAP(value, value0);
		FindSmoothFactor(voxel, needed, valid, value0, S);
		EulerStep(voxel, needed, valid, S, u0, value0);
		EulerStep(voxel, needed, valid, S, v0, u0);
#ifdef TBB
		AddTwoBody body1(needed, value0, phiNeg, phiPos, 0.75f, 0.25f, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
			body1);
#else
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
				u0[pos] = 0.75f*value0[pos] + 0.25f * v0[pos];
		END_FOR
#endif
		EulerStep(voxel, needed, valid, S, v0, u0);
#ifdef TBB
		AddTwoBody body2(needed, value0, phiNeg, value, one_third, two_third, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
			body2);
#else
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
				value[pos] = one_third * value0[pos] + two_third * v0[pos];
		END_FOR
#endif

//		printf(" (reinit) phi changes sign at iter = %d i=%d, j=%d, k=%d, needed = %d, old = %f, new = %f \n",
//					n, I, J, K, needed[INDEX(I,J,K)], old_phi, value[INDEX(I,J,K)]);


	}

	delete [] valid;


	printf("\nthere are %u points reinitialized \n\n", N);


}

static float WENOCoeffs(int i, int j, int k, int DimX, int DimY, int DimZ,
						float *x0, float inv_delta, char axis, char s){
	float q1, q2, q3, q4, q5;
	switch(axis){
	case 1:
		if(s < 0){
			q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
			q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
			q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
			q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
			q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
		}
		else{
			q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
			q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
			q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
			q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
			q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
		}
		break;
	case 2:
		if(s < 0){
			q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
			q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
			q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
			q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
			q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
		}
		else{
			q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
			q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
			q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
			q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
			q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
		}
		break;
	case 3:
		if(s < 0){
			q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
			q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
			q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
			q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
			q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
		}
		else{
			q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
			q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
			q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
			q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
			q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
		}
		break;
	}
	return HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
}


void VectorField3D::EulerStepNeeded(const Voxel &voxel, const char *needed, float dtau,
			float *value, float *value0, const float *S){


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


 	float q1, q2, q3, q4, q5;
	float *x0 = value0;

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);

		if(needed[pos]){
			float s;
			if(S == NULL)
			  s = value0[pos];
			else
			  s = S[pos];

			if(x0[pos] <= 0.f){   // liquid point
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
					phi_xm = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1, -1);
//					q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
//					q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
//					phi_xm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
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
					phi_xp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1, 1);
//					q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
//					q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
//					q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
//					phi_xp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
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
					phi_ym = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 2, -1);
//					q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
//					q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
//					phi_ym = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
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
					phi_yp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 2, 1);
//					q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
//					q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
//					q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
//					q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
//					phi_yp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
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
					phi_zm = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 3, -1);
//					q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
//					q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
//					phi_zm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
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
					phi_zp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 3, 1);
//					q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
//					q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
//					q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
//					phi_zp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
				}
			}
			else{  // air point
				bool solid_xm = false, solid_ym = false, solid_zm = false;
				bool solid_xp = false, solid_yp = false, solid_zp = false;
				TentativeGridPoint minTentative(value0[pos], i, j, k);
				vector<TentativeGridPoint> neighbors;
				minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
				if(neighbors.size() == 6){
					vel_minus_x = value0[INDEX(i-1, j, k)];
					if(voxel.InSolid(i-1, j, k))
						solid_xm = true;
					vel_plus_x = value0[INDEX(i+1, j, k)];
					if(voxel.InSolid(i+1, j, k))
						solid_xp = true;
					vel_minus_y = value0[INDEX(i, j-1, k)];
					if(voxel.InSolid(i, j-1, k))
						solid_ym = true;
					vel_plus_y = value0[INDEX(i, j+1, k)];
					if(voxel.InSolid(i, j+1, k))
						solid_yp = true;
					vel_minus_z = value0[INDEX(i, j, k-1)];
					if(voxel.InSolid(i, j, k-1))
						solid_zm = true;
					vel_plus_z = value0[INDEX(i, j, k+1)];
					if(voxel.InSolid(i, j, k+1))
						solid_zp = true;
				}
				else{
					for(int m=0; m<neighbors.size(); m++){
						if(minTentative.LeftNeighbor(neighbors[m])){
							vel_minus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
							if(voxel.InSolid(neighbors[m]))
								solid_xm = true;
						}
						if(minTentative.RightNeighbor(neighbors[m])){
							vel_plus_x = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
							if(voxel.InSolid(neighbors[m]))
								solid_xp = true;
						}
						if(minTentative.BackNeighbor(neighbors[m])){
							vel_minus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
							if(voxel.InSolid(neighbors[m]))
								solid_ym = true;
						}
						if(minTentative.FrontNeighbor(neighbors[m])){
							vel_plus_y = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
							if(voxel.InSolid(neighbors[m]))
								solid_yp = true;
						}
						if(minTentative.BottomNeighbor(neighbors[m])){
							vel_minus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
							if(voxel.InSolid(neighbors[m]))
								solid_zm = true;
						}
						if(minTentative.TopNeighbor(neighbors[m])){
							vel_plus_z = value0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
							if(voxel.InSolid(neighbors[m]))
								solid_zp = true;
						}
					}
				}

				if( voxel.InSolid(i-3,j,k) ||
					voxel.InSolid(i-2,j,k) ||
					voxel.InSolid(i-1,j,k) ||
					voxel.InSolid(i+1,j,k) ||
					voxel.InSolid(i+2,j,k)  ){
					if(solid_xm)
						vel_minus_x = value0[pos] + sign(value0[pos]) * delta;
					phi_xm = (value0[pos] - vel_minus_x) * inv_delta ;
				}
				else{
					phi_xm = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1, -1);
//					q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) * inv_delta;
//					q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
//					phi_xm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
				}

				if( voxel.InSolid(i-2,j,k)	||
					voxel.InSolid(i-1,j,k)  ||
					voxel.InSolid(i+3,j,k)  ||
					voxel.InSolid(i+1,j,k)  ||
					voxel.InSolid(i+2,j,k)  ){
					if(solid_xp)
						vel_plus_x = value0[pos] + sign(value0[pos]) * delta;
					phi_xp = (vel_plus_x - value0[pos]) * inv_delta;
				}
				else{
					phi_xp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1, 1);
//					q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) * inv_delta;
//					q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) * inv_delta;
//					q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) * inv_delta;
//					phi_xp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
				}

				if( voxel.InSolid(i,j-3,k)  ||
					voxel.InSolid(i,j-2,k)  ||
					voxel.InSolid(i,j-1,k)  ||
					voxel.InSolid(i,j+1,k)  ||
					voxel.InSolid(i,j+2,k)  ){
					if(solid_ym)
						vel_minus_y = value0[pos] + sign(value0[pos]) * delta;
					phi_ym = (value0[pos] - vel_minus_y) * inv_delta;
				}
				else{
					phi_ym = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 2, -1);
//					q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) * inv_delta;
//					q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
//					phi_ym = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
				}
				if( voxel.InSolid(i,j-2,k)  ||
					voxel.InSolid(i,j-1,k)  ||
					voxel.InSolid(i,j+3,k)  ||
					voxel.InSolid(i,j+1,k)  ||
					voxel.InSolid(i,j+2,k)  ){
					if(solid_yp)
						vel_plus_y = value0[pos] + sign(value0[pos]) * delta;
					phi_yp = (vel_plus_y - value0[pos]) * inv_delta;
				}
				else{
					phi_yp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 2, 1);
//					q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) * inv_delta;
//					q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) * inv_delta;
//					q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) * inv_delta;
//					q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) * inv_delta;
//					phi_yp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
				}

				if( voxel.InSolid(i,j,k-3)  ||
					voxel.InSolid(i,j,k-2)  ||
					voxel.InSolid(i,j,k-1)  ||
					voxel.InSolid(i,j,k+1)  ||
					voxel.InSolid(i,j,k+2)  ){
					if(solid_zm)
						vel_minus_z = value0[pos] + sign(value0[pos]) * delta;
					phi_zm = (value0[pos] - vel_minus_z) * inv_delta;
				}
				else{
					phi_zm = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 3, -1);
//					q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) * inv_delta;
//					q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
//					phi_zm = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
				}
				if( voxel.InSolid(i,j,k-2)  ||
					voxel.InSolid(i,j,k-1)  ||
					voxel.InSolid(i,j,k+3)  ||
					voxel.InSolid(i,j,k+1)  ||
					voxel.InSolid(i,j,k+2)  ){
					if(solid_zp)
						vel_plus_z = value0[pos] + sign(value0[pos]) * delta;
					phi_zp = (vel_plus_z - value0[pos]) * inv_delta;
				}
				else{
					phi_zp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 3, 1);
//					q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) * inv_delta;
//					q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) * inv_delta;
//					q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) * inv_delta;
//					q4 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) * inv_delta;
//					q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) * inv_delta;
//					phi_zp = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
				}
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
			if(S == NULL)
				s = value0[pos] / sqrtf(value0[pos] * value0[pos] +
								(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) * delta2);

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

}

void VectorField3D::ReInitializeNeeded(Voxel &voxel, const char *needed, float *value, float *value0){

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float delta3 = delta2 * delta;
	float inv_delta = 1.f / delta;
	float dtau = delta / 6;
	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;

	float eps = 1.e-4f * delta;

	int iterations = int(round(FASTMARCH_LIMIT / dtau));

	//	memset(valid, 0, DimX*DimY*DimZ);
//	float *S = new float[DimX*DimY*DimZ];
//	SetZero(S);
	float old_phi = value[INDEX(I, J, K)];

#ifdef CUDA
//	iterations = 2 * int(round(FASTMARCH_LIMIT / dtau));
	char *obj = new char[DimX*DimY*DimZ];
	char *mov_obj = new char[DimX*DimY*DimZ];
	SetZero(obj);
	SetZero(mov_obj);
	FOR_EACH_CELL
			if(voxel.InSolid(i,j,k))
				obj[INDEX(i,j,k)] = 1;
			if(voxel.InMovingSolid(i,j,k))
				mov_obj[INDEX(i,j,k)] = 1;
    END_FOR
    int dims[3];
    dims[0] = DimX; dims[1] = DimY; dims[2] = DimZ;

    ReinitializeCUDANeeded(value, obj, mov_obj, needed, iterations, I, J, K, delta, dtau, dims);

	delete [] obj;
	delete [] mov_obj;

#else
	SetEqual(value, u0);
	SetEqual(value, v0);
	float *S = new float[DimX*DimY*DimZ];
	SetEqual(value, S);
	FOR_EACH_CELL
		if(needed[INDEX(i,j,k)])
			S[INDEX(i,j,k)] = value[INDEX(i,j,k)] / sqrtf(value[INDEX(i,j,k)] * value[INDEX(i,j,k)] + delta2);
	END_FOR
	printf("Solving Reinitialization Equation with dtau = %f, iters = %d \n ", dtau, iterations);

#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif

	for(int n = 0; n < iterations; ++n){

		SWAP(value, value0);
//		FindSmoothFactor(voxel, needed, value0, S);
		EulerStepNeeded(voxel, needed, dtau, u0, value0, S);
		EulerStepNeeded(voxel, needed, dtau, v0, u0, S);
#ifdef TBB
		AddTwoBody body1(needed, value0, phiNeg, phiPos, 0.75f, 0.25f, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
			body1);
#else
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
				u0[pos] = 0.75f*value0[pos] + 0.25f * v0[pos];
		END_FOR
#endif
		EulerStepNeeded(voxel, needed, dtau, v0, u0, S);
#ifdef TBB
		AddTwoBody body2(needed, value0, phiNeg, value, one_third, two_third, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
			body2);
#else
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
				value[pos] = one_third * value0[pos] + two_third * v0[pos];
		END_FOR
#endif

//		printf(" (reinit) phi changes sign at iter = %d i=%d, j=%d, k=%d, needed = %d, old = %f, new = %f \n",
//					n, I, J, K, needed[INDEX(I,J,K)], old_phi, value[INDEX(I,J,K)]);


	}
	delete [] S;

#endif

}

#else
void VectorField3D::ReInitialize(const Voxel &voxel){

	printf(" before init: phi_u = %f, phi_u(i-1) = %f, phi_u(j-1) = %f, phi_u(k-1) = %f \n ",
			phi_u[INDEX(I,J,K)], phi_u[INDEX(I-1,J,K)], phi_u[INDEX(I,J-1,K)], phi_u[INDEX(I,J,K-1)]);
	printf(" phi_u = %f, phi_u(i+1) = %f, phi_u(j+1) = %f, phi_u(k+1) = %f \n ",
			phi_u[INDEX(I,J,K)], phi_u[INDEX(I+1,J,K)], phi_u[INDEX(I,J+1,K)], phi_u[INDEX(I,J,K+1)]);

	printf(" before init: phi_v = %f, phi_v(i-1) = %f, phi_v(j-1) = %f, phi_v(k-1) = %f \n ",
			phi_v[INDEX(I,J,K)], phi_v[INDEX(I-1,J,K)], phi_v[INDEX(I,J-1,K)], phi_v[INDEX(I,J,K-1)]);
	printf(" phi_v = %f, phi_v(i+1) = %f, phi_v(j+1) = %f, phi_v(k+1) = %f \n ",
			phi_v[INDEX(I,J,K)], phi_v[INDEX(I+1,J,K)], phi_v[INDEX(I,J+1,K)], phi_v[INDEX(I,J,K+1)]);

	printf(" before init: phi_w = %f, phi_w(i-1) = %f, phi_w(j-1) = %f, phi_w(k-1) = %f \n ",
			phi_w[INDEX(I,J,K)], phi_w[INDEX(I-1,J,K)], phi_w[INDEX(I,J-1,K)], phi_w[INDEX(I,J,K-1)]);
	printf(" phi_w = %f, phi_w(i+1) = %f, phi_w(j+1) = %f, phi_w(k+1) = %f \n ",
			phi_w[INDEX(I,J,K)], phi_w[INDEX(I+1,J,K)], phi_w[INDEX(I,J+1,K)], phi_w[INDEX(I,J,K+1)]);

	char *valid = new char[DimX*DimY*DimZ];
	//memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	map<u_int, KnownPoint> pos_band;
	map<u_int, KnownPoint> neg_band;
	InitializeInterface(voxel, pos_band, neg_band, phi_u, phi_u_obj, valid, 1);
	FastMarching(pos_band, neg_band, voxel, phi_u, phi_u_obj, valid);

	pos_band.clear();
	neg_band.clear();
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	InitializeInterface(voxel, pos_band, neg_band, phi_v, phi_v_obj, valid, 2);
	FastMarching(pos_band, neg_band, voxel, phi_v, phi_v_obj, valid);

	pos_band.clear();
	neg_band.clear();
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	InitializeInterface(voxel, pos_band, neg_band, phi_w, phi_w_obj, valid, 3);
	FastMarching(pos_band, neg_band, voxel, phi_w, phi_w_obj, valid);

	if(valid) delete [] valid;

	printf(" after init: phi_u = %f, phi_u(i-1) = %f, phi_u(j-1) = %f, phi_u(k-1) = %f \n ",
			phi_u[INDEX(I,J,K)], phi_u[INDEX(I-1,J,K)], phi_u[INDEX(I,J-1,K)], phi_u[INDEX(I,J,K-1)]);
	printf(" phi_u = %f, phi_u(i+1) = %f, phi_u(j+1) = %f, phi_u(k+1) = %f \n ",
			phi_u[INDEX(I,J,K)], phi_u[INDEX(I+1,J,K)], phi_u[INDEX(I,J+1,K)], phi_u[INDEX(I,J,K+1)]);
	printf(" after init: phi_v = %f, phi_v(i-1) = %f, phi_v(j-1) = %f, phi_v(k-1) = %f \n ",
			phi_v[INDEX(I,J,K)], phi_v[INDEX(I-1,J,K)], phi_v[INDEX(I,J-1,K)], phi_v[INDEX(I,J,K-1)]);
	printf(" phi_v = %f, phi_v(i+1) = %f, phi_v(j+1) = %f, phi_v(k+1) = %f \n ",
			phi_v[INDEX(I,J,K)], phi_v[INDEX(I+1,J,K)], phi_v[INDEX(I,J+1,K)], phi_v[INDEX(I,J,K+1)]);

	printf(" after init: phi_w = %f, phi_w(i-1) = %f, phi_w(j-1) = %f, phi_w(k-1) = %f \n ",
			phi_w[INDEX(I,J,K)], phi_w[INDEX(I-1,J,K)], phi_w[INDEX(I,J-1,K)], phi_w[INDEX(I,J,K-1)]);
	printf(" phi_w = %f, phi_w(i+1) = %f, phi_w(j+1) = %f, phi_w(k+1) = %f \n ",
			phi_w[INDEX(I,J,K)], phi_w[INDEX(I+1,J,K)], phi_w[INDEX(I,J+1,K)], phi_w[INDEX(I,J,K+1)]);

}
#endif

void VectorField3D::ReInitializeObject(const Voxel &voxel){
//#ifdef WIN32
//	hash_map<u_int, KnownPoint> pos_band;
//	hash_map<u_int, KnownPoint> neg_band;
//#else
//	unordered_map<u_int, KnownPoint> pos_band;
//	unordered_map<u_int, KnownPoint> neg_band;
//#endif


	char *maskforLevelset = new char[DimX*DimY*DimZ];
	memset(maskforLevelset, 0, DimX*DimY*DimZ*sizeof(char));
	InitializeInterface(voxel, maskforLevelset, phi_c_obj, 0);
	FastMarching(maskforLevelset, voxel, phi_c_obj);
//	pos_band.clear();
//	neg_band.clear();

	memset(maskforLevelset, 0, DimX*DimY*DimZ*sizeof(char));
	InitializeInterface(voxel, maskforLevelset, phi_u_obj, 1);
	FastMarching(maskforLevelset, voxel, phi_u_obj);

//	pos_band.clear();
//	neg_band.clear();
	memset(maskforLevelset, 0, DimX*DimY*DimZ*sizeof(char));
	InitializeInterface(voxel, maskforLevelset, phi_v_obj, 2);
	FastMarching(maskforLevelset, voxel, phi_v_obj);

//	pos_band.clear();
//	neg_band.clear();
	memset(maskforLevelset, 0, DimX*DimY*DimZ*sizeof(char));
	InitializeInterface(voxel, maskforLevelset, phi_w_obj, 3);
	FastMarching(maskforLevelset, voxel, phi_w_obj);

	printf("VectorField3D::ReInitializeObject\n");
	printf("phi_u_obj[%d,%d,%d] = %f \n",I, J, K, phi_u_obj[INDEX(I,J,K)]);
	printf("phi_v_obj[%d,%d,%d] = %f \n",I, J, K, phi_v_obj[INDEX(I,J,K)]);
	printf("phi_w_obj[%d,%d,%d] = %f \n",I, J, K, phi_w_obj[INDEX(I,J,K)]);

	delete [] maskforLevelset;

}



void VectorField3D::InitializeInterface( const Voxel &voxel,
								char *mask,	float *phi, int vel_index){


	FOR_EACH_CELL
		u_int index = INDEX(i,j,k);
		KnownPoint tmp(i, j, k);
		if(phi[INDEX(i,j,k)] == 0.f){
//			pos_band.insert(make_pair(INDEX(i,j,k),tmp));
//			neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			mask[index] |= NEGINDONEBAND;
			mask[index] |= POSINDONEBAND;
		}
		else{
			TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
			vector<TentativeGridPoint> neighbors;
			tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
			for(int m=0; m<6; m++){
				if(tp.NeighborExists(m, DimX, DimY, DimZ)){
					neighbors[m].value = phi[POS(neighbors[m])];
					if(tp.value*neighbors[m].value <= 0.f){
						if(phi[INDEX(i,j,k)] > 0.f)
							mask[index] |= POSINDONEBAND;
//							pos_band.insert(make_pair(INDEX(i,j,k),tmp));
						else
							mask[index] |= NEGINDONEBAND;
//							neg_band.insert(make_pair(INDEX(i,j,k),tmp));
					}

				}
			}
		}
	END_FOR


//	map<u_int, KnownPoint>::iterator found;
//	found = pos_band.find(INDEX(I,J,K));
//	if(found == pos_band.end())
//		printf("(%d,%d,%d) is not in the band \n", I,J,K);
//	else
//		printf("(%d,%d,%d) is in the band \n", I,J,K);
//	FOR_EACH_CELL
//		if(i==I && j==J && k==K && vel_index == 1){
//			printf("init phi = %f \n",phi[INDEX(i,j,k)]);
//			printf("init phi-x = %f phi+x = %f\n",phi[INDEX(i-1,j,k)],phi[INDEX(i+1,j,k)]);
//			printf("init phi-y = %f phi+y = %f\n",phi[INDEX(i,j-1,k)],phi[INDEX(i,j+1,k)]);
//			printf("init phi-z = %f phi+z = %f\n",phi[INDEX(i,j,k-1)],phi[INDEX(i,j,k+1)]);
//		}
//	END_FOR
//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);


}


float VectorField3D::RefinedClosestDistance(const Voxel &voxel, float *phi, int i, int j, int k, int vel_index) const{
	float r;
	TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
	vector<TentativeGridPoint> neighbors;
	tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
	Point p = voxel.VelPosition(vel_index, i, j, k);
	float *candidates = new float[6];
	for(int m=0; m<6; m++)
		candidates[m] = INFINITY;
	Point *candidatesP = new Point[6];
	for(int m=0; m<6; m++){
		if(tp.NeighborExists(m, DimX, DimY, DimZ)){
//						Point pn = voxel.VelPosition(vel_index,
//													neighbors[m].ii,
//													neighbors[m].jj,
//													neighbors[m].kk);
			Point pn = voxel.VelClosePosition(m,vel_index,i,j,k);
//						neighbors[m].value = phi[POS(neighbors[m])];
//						if(phi_obj[POS(neighbors[m])] == 0.f && phi_c[INDEX(i,j,k)] < 0.f)
//							neighbors[m].value = 0.f;
			neighbors[m].value = TriInterp(voxel, pn, phi_c);
//						if(phi_obj[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] >= 0.f){
//							if(tp.value <= 0.f){
//								candidates[m] = min(0.5f*delta, fabs(tp.value));
//								candidatesP[m] = pn;
//							}
//							else
//								continue;
//						}
//						else
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
//						if(i == I && j == J && k == K){
//							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
//								i, j, k, tp.value, neighbors[m].value, candidates[m]);
//							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
//									neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
//									neighbors[m].value, phi[POS(neighbors[m])]);
//
//						}
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
		if(v3.Length() == 0.f){
			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
					i,j,k, tp.value);
			for(int n=0; n < neighbors.size(); ++n)
				printf("%f ", candidates[n]);
			printf("\n");
			for(int n=0; n < neighbors.size(); ++n)
				printf("%f ", neighbors[n].value);
			printf("\n");
			printf(" p0 = (%f, %f, %f), p1 = (%f, %f, %f), p2 = (%f, %f, %f)\n",
					phi_min_p[0].x, phi_min_p[0].y, phi_min_p[0].z,
					phi_min_p[1].x, phi_min_p[1].y, phi_min_p[1].z,
					phi_min_p[2].x, phi_min_p[2].y, phi_min_p[2].z);
			phimin = fabs(tp.value);
//						exit(1);
		}
		else{
			v3 = Normalize(v3);
			Vector v4 = p - phi_min_p[0];
			phimin = AbsDot(v3, v4);
		}
	}
	else if( phi_min[0] != INFINITY && phi_min[1] != INFINITY ){
		// 2 points forming a line segment
		phimin= DistanceToLine(phi_min_p[0], phi_min_p[1], p, tp.value);
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
	if(candidates) delete [] candidates;
	if(candidatesP) delete [] candidatesP;
	return phimin;
}

float VectorField3D::ClosestDistance(const Voxel &voxel, float *phi, const char *valid,
		                        int i, int j, int k, int vel_index) const{
	float r;

	float delta = voxel.VoxelDelta();
	TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
	vector<TentativeGridPoint> neighbors;
	tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
	Point p = voxel.VelPosition(vel_index, i, j, k);
	float *candidates = new float[6];
	for(int m=0; m<6; m++)
		candidates[m] = INFINITY;
	Point *candidatesP = new Point[6];
	for(int m=0; m<6; m++){
		if(tp.NeighborExists(m, DimX, DimY, DimZ)){
			Point pn = voxel.VelPosition(vel_index,
										neighbors[m].ii,
										neighbors[m].jj,
										neighbors[m].kk);
//			Point pn = voxel.VelClosePosition(m,vel_index,i,j,k);
			neighbors[m].value = phi[POS(neighbors[m])];
//						if(phi_obj[POS(neighbors[m])] == 0.f && phi_c[INDEX(i,j,k)] < 0.f)
//							neighbors[m].value = 0.f;
//			neighbors[m].value = TriInterp(voxel, pn, phi_c);
//			if(voxel.InSolid(neighbors[m])){
//				if(vel_index == 1){
//					if(voxel.InSolid(tp)){
//						if(tp.LeftNeighbor(neighbors[m]))
//							continue;
//					}
//				}
//				else if(vel_index == 2){
//
//				}
//				else if(vel_index == 3){
//
//				}
//			}
//			if(tp.value*neighbors[m].value == 0.f){
//				candidates[m] = delta;
//				candidatesP[m] = pn;
//			}
			if(valid[INDEX(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk)] &&
					tp.value*neighbors[m].value < 0.f){

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
//						if(i == I && j == J && k == K){
//							printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
//								i, j, k, tp.value, neighbors[m].value, candidates[m]);
//							printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
//									neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
//									neighbors[m].value, phi[POS(neighbors[m])]);
//
//						}
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
//		Vector v1 = phi_min_p[0] - phi_min_p[1];
//		Vector v2 = phi_min_p[0] - phi_min_p[2];
//		Vector v3 =  Cross(v1, v2);
//		if(v3.Length() == 0.f){
//			printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %f \n",
//					i,j,k, tp.value);
//			for(int n=0; n < neighbors.size(); ++n)
//				printf("%f ", candidates[n]);
//			printf("\n");
//			for(int n=0; n < neighbors.size(); ++n)
//				printf("%f ", neighbors[n].value);
//			printf("\n");
//			printf(" p0 = (%f, %f, %f), p1 = (%f, %f, %f), p2 = (%f, %f, %f)\n",
//					phi_min_p[0].x, phi_min_p[0].y, phi_min_p[0].z,
//					phi_min_p[1].x, phi_min_p[1].y, phi_min_p[1].z,
//					phi_min_p[2].x, phi_min_p[2].y, phi_min_p[2].z);
//			phimin = fabs(tp.value);
////						exit(1);
//		}
//		else{
//			v3 = Normalize(v3);
//			Vector v4 = p - phi_min_p[0];
//			phimin = AbsDot(v3, v4);
//		}
		phimin = phi_min[0]*phi_min[1]*phi_min[2] /
			       sqrt( phi_min[0]*phi_min[0]*phi_min[1]*phi_min[1] +
			    		 phi_min[0]*phi_min[0]*phi_min[2]*phi_min[2] +
			    		 phi_min[2]*phi_min[2]*phi_min[1]*phi_min[1]);
	}
	else if( phi_min[0] != INFINITY && phi_min[1] != INFINITY ){
		// 2 points forming a line segment
//		phimin= DistanceToLine(phi_min_p[0], phi_min_p[1], p, tp.value);
		phimin = TwoPointDistance(phi_min[0], phi_min[1]);
	}
	else if( phi_min[1] != INFINITY && phi_min[2] != INFINITY ){
		// 2 points forming a line segment
//		phimin = DistanceToLine(phi_min_p[1], phi_min_p[2], p, tp.value);
		phimin = TwoPointDistance(phi_min[1], phi_min[2]);
	}
	else if( phi_min[0] != INFINITY && phi_min[2] != INFINITY ){
		// 2 points forming a line segment
//		phimin = DistanceToLine(phi_min_p[0], phi_min_p[2], p, tp.value);
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
//	float phimin = FindMinimum(candidates, 6);
	if(candidates) delete [] candidates;
	if(candidatesP) delete [] candidatesP;
	return phimin;
}


void VectorField3D::ProcessNonSolidPoint(const Voxel &voxel, float *phi_tmp,
		                  float *phi, const char *valid,
		                  int i, int j, int k, int vel_index) const{
	if(phi[INDEX(i,j,k)] == 0.f)
		phi_tmp[INDEX(i,j,k)] = 0.f;
	else{
//				float phimin1 = RefinedClosestDistance(voxel, phi, i, j, k, vel_index);
//				float phimin2 = ClosestDistance(voxel, phi, i, j, k, vel_index);
//				float phimin = min(phimin1, phimin2);
		float phimin = ClosestDistance(voxel, phi, valid, i, j, k, vel_index);
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
		if( i == I && j == J && k == K)
			printf("phi = %f, phimin1 = %f, vel_index = %d \n",
					phi[INDEX(i,j,k)], phimin, vel_index);
		if(phi[INDEX(i,j,k)] > 0.f){
			phi_tmp[INDEX(i,j,k)] = phimin;
			//if(phi_min != INFINITY)
				//pos_band.insert(make_pair(INDEX(i,j,k),tmp));
		}
		else{
			phi_tmp[INDEX(i,j,k)] = -phimin;
//			if(phi_min != INFINITY)
//				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
		}

//		printf("after i = %d, j= %d, k = %d, phi = %f \n",
//				i, j, k, phi[INDEX(i,j,k)]);
	}
}

void VectorField3D::SetValidPoint(const Voxel &voxel, char *valid, int vel_index) const{

	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k)){  // air or water cell
//			valid[INDEX(i,j,k)] = 1;
////			if( vel_index == 3 && source && source->IsSourceNeighbor(i,j,k) )
////				valid[INDEX(i,j,k)] = 0;
//		}
//		else{   // solid cell
//			if( vel_index == 1 ){
//				if( i < DimX-1 && !voxel.InSolid(i+1,j,k))
//					valid[INDEX(i,j,k)] = 1;
//			}
//			else if( vel_index == 2 ){
//				if( j < DimY-1 && !voxel.InSolid(i,j+1,k))
//					valid[INDEX(i,j,k)] = 1;
//			}
//			else if( vel_index == 3 ){
//				if( k < DimZ-1 && !voxel.InSolid(i,j,k+1))
//					valid[INDEX(i,j,k)] = 1;
//			}
//		}
		if( vel_index == 1 ){
//			if(!voxel.InSolid(i,j,k) && !voxel.InSolid(i+1,j,k))
			if( phi_u_obj[INDEX(i,j,k)]+E_EPSIL < 0.f)
				valid[INDEX(i,j,k)] = 1;
//			if(i < DimX-1 && voxel.InSource(i+1,j,k))
//				valid[INDEX(i,j,k)] = 0;
			else if(i < DimX-1 && voxel.InSolid(i,j,k) && (voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k)))
				valid[INDEX(i,j,k)] = 1;
			else if(i < DimX-1 && (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InSolid(i+1,j,k))
				valid[INDEX(i,j,k)] = 1;
		}
		else if( vel_index == 2 ){
//			if(!voxel.InSolid(i,j,k) && !voxel.InSolid(i,j+1,k))
			if( phi_v_obj[INDEX(i,j,k)]+E_EPSIL < 0.f)
				valid[INDEX(i,j,k)] = 1;
//			if(j < DimY-1 && voxel.InSource(i,j+1,k))
//				valid[INDEX(i,j,k)] = 0;
			else if(j < DimY-1 && voxel.InSolid(i,j,k) && (voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k)))
				valid[INDEX(i,j,k)] = 1;
			else if(j < DimY-1 && (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InSolid(i,j+1,k))
				valid[INDEX(i,j,k)] = 1;
		}
		else if( vel_index == 3 ){
			if( phi_w_obj[INDEX(i,j,k)]+E_EPSIL < 0.f)
//			if(!voxel.InSolid(i,j,k) && !voxel.InSolid(i,j,k+1))
				valid[INDEX(i,j,k)] = 1;
//			if(k < DimZ-1 && voxel.InSource(i,j,k+1))
//				valid[INDEX(i,j,k)] = 0;
			else if(k < DimZ-1 && voxel.InSolid(i,j,k) && (voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1)))
				valid[INDEX(i,j,k)] = 1;
			else if(k < DimZ-1 && (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InSolid(i,j,k+1))
				valid[INDEX(i,j,k)] = 1;
		}
		else if(vel_index == 0){
			if(!voxel.InSolid(i,j,k))
				valid[INDEX(i,j,k)] = 1;
		}
//		if(voxel.InSource(i,j,k))
//			valid[INDEX(i,j,k)] = 0;
//		if (i == I && j==J && k == K){
//			if(vel_index == 1){
//				printf("vel_index = %d phiu_obj = %f\n", vel_index, phi_u_obj[INDEX(i,j,k)]);
//				printf("phiu_obj_minus_x = %f, phiu_obj_plus_x = %f\n",
//						phi_u_obj[INDEX(i-1,j,k)],phi_u_obj[INDEX(i+1,j,k)]);
//				printf("phiu_obj_minus_y = %f, phiu_obj_plus_y = %f\n",
//						phi_u_obj[INDEX(i,j-1,k)],phi_u_obj[INDEX(i,j+1,k)]);
//				printf("phiu_obj_minus_z = %f, phiu_obj_plus_z = %f\n",
//						phi_u_obj[INDEX(i,j,k-1)],phi_u_obj[INDEX(i,j,k+1)]);
//			}
//			else if(vel_index == 2){
//				printf("vel_index = %d phiu_obj = %f\n", vel_index, phi_v_obj[INDEX(i,j,k)]);
//				printf("phiv_obj_minus_x = %f, phiv_obj_plus_x = %f\n",
//						phi_v_obj[INDEX(i-1,j,k)],phi_v_obj[INDEX(i+1,j,k)]);
//				printf("phiv_obj_minus_y = %f, phiv_obj_plus_y = %f\n",
//						phi_v_obj[INDEX(i,j-1,k)],phi_v_obj[INDEX(i,j+1,k)]);
//				printf("phiv_obj_minus_z = %f, phiv_obj_plus_z = %f\n",
//						phi_v_obj[INDEX(i,j,k-1)],phi_v_obj[INDEX(i,j,k+1)]);
//			}
//			else if(vel_index == 3){
//				printf("vel_index = %d phiw_obj = %f\n", vel_index, phi_w_obj[INDEX(i,j,k)]);
//				printf("phiw_obj_minus_x = %f, phiw_obj_plus_x = %f\n",
//						phi_w_obj[INDEX(i-1,j,k)],phi_w_obj[INDEX(i+1,j,k)]);
//				printf("phiw_obj_minus_y = %f, phiw_obj_plus_y = %f\n",
//						phi_w_obj[INDEX(i,j-1,k)],phi_w_obj[INDEX(i,j+1,k)]);
//				printf("phiw_obj_minus_z = %f, phiw_obj_plus_z = %f\n",
//						phi_w_obj[INDEX(i,j,k-1)],phi_w_obj[INDEX(i,j,k+1)]);
//			}
//
//		}
	END_FOR
}

/*void VectorField3D::InitializeInterface( const Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band,
								float *phi, const float *phi_obj,
								char *valid, int vel_index){

	float *phi_tmp = new float[DimX*DimY*DimZ];
	memset(phi_tmp, 0, DimX*DimY*DimZ);

	SetValidPoint(voxel, valid, vel_index);

	FOR_EACH_CELL
		if(valid[INDEX(i,j,k)])
			ProcessNonSolidPoint(voxel, phi_tmp, phi, valid, i, j, k, vel_index);
		else
			phi_tmp[INDEX(i,j,k)] = INFINITY;
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
		}
	END_FOR
//	map<u_int, KnownPoint>::iterator found;
//	found = pos_band.find(INDEX(I,J,K));
//	if(found == pos_band.end())
//		printf("(%d,%d,%d) is not in the band \n", I,J,K);
//	else
//		printf("(%d,%d,%d) is in the band \n", I,J,K);
//	FOR_EACH_CELL
//		if(i==I && j==J && k==K && vel_index == 1){
//			printf("init phi = %f \n",phi[INDEX(i,j,k)]);
//			printf("init phi-x = %f phi+x = %f\n",phi[INDEX(i-1,j,k)],phi[INDEX(i+1,j,k)]);
//			printf("init phi-y = %f phi+y = %f\n",phi[INDEX(i,j-1,k)],phi[INDEX(i,j+1,k)]);
//			printf("init phi-z = %f phi+z = %f\n",phi[INDEX(i,j,k-1)],phi[INDEX(i,j,k+1)]);
//		}
//	END_FOR
//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);
	if(phi_tmp) delete [] phi_tmp;


}*/

void VectorField3D::InitializeInterface( const Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band,
								float *phi, const float *phi_obj,
								char *valid, int vel_index){

	SetValidPoint(voxel, valid, vel_index);

	FOR_EACH_CELL

//			if(i == I && j == J && k == K)
//				printf("grid point (%d,%d,%d) in band \n", i, j, k);
		if(valid[INDEX(i,j,k)]){
			KnownPoint tmp(i, j, k);
			if(phi[INDEX(i,j,k)] == 0.f){
				pos_band.insert(make_pair(INDEX(i,j,k),tmp));
				neg_band.insert(make_pair(INDEX(i,j,k),tmp));
			}
			else{
				TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
				vector<TentativeGridPoint> neighbors;
				tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
				for(int m=0; m<6; m++){
					if(tp.NeighborExists(m, DimX, DimY, DimZ)){
						neighbors[m].value = phi[POS(neighbors[m])];
						if(valid[POS(neighbors[m])] && tp.value*neighbors[m].value <= 0.f){
							if(phi[INDEX(i,j,k)] > 0.f)
								pos_band.insert(make_pair(INDEX(i,j,k),tmp));
							else
								neg_band.insert(make_pair(INDEX(i,j,k),tmp));
						}
					}
				}
			}
		}
	END_FOR
//	map<u_int, KnownPoint>::iterator found;
//	found = pos_band.find(INDEX(I,J,K));
//	if(found == pos_band.end())
//		printf("(%d,%d,%d) is not in the band \n", I,J,K);
//	else
//		printf("(%d,%d,%d) is in the band \n", I,J,K);
//	FOR_EACH_CELL
//		if(i==I && j==J && k==K && vel_index == 1){
//			printf("init phi = %f \n",phi[INDEX(i,j,k)]);
//			printf("init phi-x = %f phi+x = %f\n",phi[INDEX(i-1,j,k)],phi[INDEX(i+1,j,k)]);
//			printf("init phi-y = %f phi+y = %f\n",phi[INDEX(i,j-1,k)],phi[INDEX(i,j+1,k)]);
//			printf("init phi-z = %f phi+z = %f\n",phi[INDEX(i,j,k-1)],phi[INDEX(i,j,k+1)]);
//		}
//	END_FOR
//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);


}


/*void VectorField3D::InitializeInterface(const Voxel &voxel,
								map<u_int, KnownPoint> &pos_band,
								map<u_int, KnownPoint> &neg_band,
								float *phi, int vel_index){
	float r;
	float *phi_tmp = new float[DimX*DimY*DimZ];
	memset(phi_tmp, 0, DimX*DimY*DimZ);

	FOR_EACH_CELL
		if(phi[INDEX(i,j,k)] == 0.f)
				phi_tmp[INDEX(i,j,k)] = 0.f;
		else{
			TentativeGridPoint tp(phi[INDEX(i,j,k)], i, j, k);
			vector<TentativeGridPoint> neighbors;
			tp.Neigbor(neighbors, DimX, DimY, DimZ);
	//		Point p = voxel.VoxelCenterPosition(i,j,k);
			Point p = voxel.VelPosition(vel_index, i, j, k);
			float *candidates = new float[neighbors.size()];
			for(int m=0; m<neighbors.size(); m++){
				candidates[m] = INFINITY;
	//			Point pn = voxel.VoxelCenterPosition(neighbors[m].ii,
	//												 neighbors[m].jj,
	//												 neighbors[m].kk);
				Point pn = voxel.VelPosition(vel_index,
											neighbors[m].ii,
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
						printf("wrong! phi_min should be > 0 at [%d, %d, %d] with phi = %16.12f \n",
								i,j,k, tp.value);
						printf(" neighbor[m].value = %16.12f  P at (%f,%f,%f)\n",
								neighbors[m].value, p.x, p.y, p.z);
						printf("candidate = %16.12f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
								candidates[m], r, pn.x, pn.y, pn.z);
						// here if candidates[m] < 0.f, it must be due to machine's
						// truncation error. then it is flipped back to positive
						candidates[m] = fabs(candidates[m]);
					}
	//				printf("A. i = %d, j= %d, k = %d, phi = %f, phi_n = %f, cand = %f \n",
	//							i, j, k, tp.value, neighbors[m].value, candidates[m]);
	//				printf("B. i = %d, j= %d, k = %d, phi = %f, phi_n = %f \n",
	//										neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
	//										neighbors[m].value, phi[POS(neighbors[m])]);
				}
	//			if(i==I && j==J && k==K && vel_index == 1){
	//				printf("init phi = %f, neighbor = %f \n", tp.value, neighbors[m].value );
	////				printf("init phi-x = %f phi+x = %f\n",phi[INDEX(i-1,j,k)],phi[INDEX(i+1,j,k)]);
	////				printf("init phi-y = %f phi+y = %f\n",phi[INDEX(i,j-1,k)],phi[INDEX(i,j+1,k)]);
	////				printf("init phi-z = %f phi+z = %f\n",phi[INDEX(i,j,k-1)],phi[INDEX(i,j,k+1)]);
	//			}
			}
			float phi_min = FindMinimumV(candidates, neighbors.size());
	//		if(i==I && j==J && k==K && vel_index == 1){
	//			printf("phi_min = %f \n", phi_min);
	//			if(phi_min != INFINITY)
	//				printf("not inf \n");
	//			else
	//				printf("is inf \n");
	//		}
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
		}
	END_FOR
//	map<u_int, KnownPoint>::iterator found;
//	found = pos_band.find(INDEX(I,J,K));
//	if(found == pos_band.end())
//		printf("(%d,%d,%d) is not in the band \n", I,J,K);
//	else
//		printf("(%d,%d,%d) is in the band \n", I,J,K);
//	FOR_EACH_CELL
//		if(i==I && j==J && k==K && vel_index == 1){
//			printf("init phi = %f \n",phi[INDEX(i,j,k)]);
//			printf("init phi-x = %f phi+x = %f\n",phi[INDEX(i-1,j,k)],phi[INDEX(i+1,j,k)]);
//			printf("init phi-y = %f phi+y = %f\n",phi[INDEX(i,j-1,k)],phi[INDEX(i,j+1,k)]);
//			printf("init phi-z = %f phi+z = %f\n",phi[INDEX(i,j,k-1)],phi[INDEX(i,j,k+1)]);
//		}
//	END_FOR

//	printf("C: phi = %f \n", phi[INDEX(I,J,K)]);
	if(phi_tmp) delete [] phi_tmp;


}*/



void VectorField3D::FastMarching(char *mask, const Voxel &voxel, float *phi){
	// initialize grid points within the band of the interface
	// already initialized?

	float delta = voxel.VoxelDelta();
	int num_points=0;

	MinHeap<TentativeGridPoint, float> heap;


	//Process grid points at phi > 0
//	for(pos=posband.begin(); pos!=posband.end(); ++pos){
//			KnownPoint &p = pos->second;
//			u_int index = pos->first;
//			TentativeGridPoint tp(0.f, p);
//			vector<TentativeGridPoint> neigbor;
//		    tp.Neigbor(neigbor, DimX, DimY, DimZ);
//			for(int i=0;i<neigbor.size();i++){
//				foundpos1 = posband.find(Index(neigbor[i]));
//				if(foundpos1==posband.end()){ // if this neighbor is not in accepted band
//					foundpos = posadj.find(Index(neigbor[i]));
//					if(foundpos==posadj.end()){ // if this neighbor is not in close band
////						if(!(neigbor[i].ii < x1 && neigbor[i].ii > x0 &&
////						     neigbor[i].jj < y1 && neigbor[i].jj > y0 &&
////					   	     neigbor[i].kk < z1 && neigbor[i].kk > z0)){
//						if(phi[POS(neigbor[i])] > 0.f){
//							KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//							posadj.insert(make_pair(Index(neigbor[i]),tmp)); // add it to close band
//						}
//					}
//				}
//			}
////			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
////			i++;
//	}

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mask[pos] & POSINDONEBAND){
			TentativeGridPoint tp(0.f, i, j, k);
			vector<TentativeGridPoint> neigbor;
			tp.Neigbor(neigbor, DimX, DimY, DimZ);
			for(int n=0;n<neigbor.size();n++){
				int ii = neigbor[n].ii;
				int jj = neigbor[n].jj;
				int kk = neigbor[n].kk;
				if(!(mask[INDEX(ii,jj,kk)] & POSINDONEBAND) ){ // if this neighbor is not in accepted band
					if(!(mask[INDEX(ii,jj,kk)] & POSINTRIALBAND)){ // if this neighbor is not in close band
						if(phi[POS(neigbor[n])] > 0.f){
							mask[INDEX(ii,jj,kk)] |= POSINTRIALBAND;
							++num_points;
						}
					}
				}
			}
		}
	END_FOR

	printf(" adj_points = %d \n", num_points);

//	map<u_int, KnownPoint>::iterator foundpos;
//	foundpos = posadj.find(INDEX(I,J,K));
//	if(foundpos == posadj.end())
//		printf("(%d,%d,%d) is not in adjacent band \n", I,J,K);
//	else
//		printf("(%d,%d,%d) is in adjacent band band \n", I,J,K);

//	for(pos=posadj.begin(); pos!=posadj.end(); ++pos){
//			KnownPoint &p = pos->second;
//			u_int index = pos->first;
//			float tentativeValue = UpdateTentativeValue(p, index, posband, voxel, phi);
//			TentativeGridPoint tp(tentativeValue, p);
//			heap.insert(tp, tentativeValue);
//	}

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mask[pos] & POSINTRIALBAND){
			KnownPoint p(i,j,k);
			float tentativeValue = UpdateTentativeValue(p, pos, mask, (char)POSINDONEBAND, voxel, phi);
			TentativeGridPoint tp(tentativeValue, i, j, k);
			heap.insert(tp, tentativeValue);
		}
	END_FOR


	float minValue;
	do{
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		++num_points;
//		printf("i = %d \n", i);
		u_int indexDoneBand = Index(minTentative);
//		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
//		if(minValue > FASTMARCH_LIMIT)
//			break;
//		if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//			printf("old phi = %f\n",phi[POS(t)]);
//		}
		phi[indexDoneBand] = minValue;
//		if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//			printf("new phi = %f\n",phi[POS(t)]);
//		}
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
//	    posband.insert( make_pair(Index(minTentative), t) );
//		posadj.erase(Index(minTentative));
		mask[indexDoneBand] |= POSINDONEBAND;
		mask[indexDoneBand] ^= POSINTRIALBAND;
		//printf("heap size = %d, pos_adj size = %d\n", heap.size(), posadj.size());
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
	    vector<TentativeGridPoint> neigbor;
	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
	    for(int i=0;i<neigbor.size();i++){
	    	u_int indexTrialBand =  Index(neigbor[i]);
	    	if(!(mask[indexTrialBand] & POSINDONEBAND)){ // not currently in done band
//	    		if(neigbor[i].value == INFINITY){
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
//	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    		}
//		    	found = posadj.find(Index(neigbor[i]));
//	    		if(found == posadj.end()){ // not currently in the heap
	    		if(!(mask[indexTrialBand] & POSINTRIALBAND)){ // not currently in trial band
	    			if(phi[indexTrialBand] > 0.f){
//			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//			    		posadj.insert(make_pair(Index(neigbor[i]), tmp));
		//		    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
	    				mask[indexTrialBand] |= POSINTRIALBAND;
	    				neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(),
	    						indexTrialBand, mask, (char)POSINDONEBAND,
	    					   voxel, phi);
			    		heap.insert(neigbor[i], neigbor[i].value);
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
	    		else{	    			// already in heap,
	    			neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(),
	    					indexTrialBand, mask, (char)POSINDONEBAND,
							   voxel, phi);
	    			heap.update(neigbor[i], neigbor[i].value);	// just update values
													    		// according to newly added point
//					printf("update neigbor i = %d, value = %f \n",
//	    							i, neigbor[i].value );
	    		}
	    	}
	    }
	    neigbor.clear();
//	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
	} while(!heap.empty());

//	printf("pos band size = %ld, num_points=%d \n", posband.size(), num_points);
	// in case when marching stops at a certain distance, the heap is not empty
		// make it empty for negative band
		if(!heap.empty())
			heap.clear();

		//Process grid points at phi < 0
//		for(pos=negband.begin(); pos!=negband.end(); ++pos){
//				KnownPoint &p = pos->second;
//	//			u_int index = pos->first;
//				TentativeGridPoint tp(0.f, p);
//				vector<TentativeGridPoint> neigbor;
//			    tp.Neigbor(neigbor, DimX, DimY, DimZ);
//				for(int i=0;i<neigbor.size();i++){
//					foundpos1 = negband.find(Index(neigbor[i]));
//					if(foundpos1==negband.end()){
//	//				if(!v.IsDone(neigbor[i])){
//						foundpos = negadj.find(Index(neigbor[i]));
//						if(foundpos==negadj.end()){
//	//					if(!v.IsClose(neigbor[i])){
//	//						if((neigbor[i].ii < x1 && neigbor[i].ii > x0 &&
//	//						    neigbor[i].jj < y1 && neigbor[i].jj > y0 &&
//	//					   	    neigbor[i].kk < z1 && neigbor[i].kk > z0)){
//							if(phi[POS(neigbor[i])] < 0.f){
//								KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//								negadj.insert(make_pair(Index(neigbor[i]),tmp));
//							//	voxel.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
//							}
//						}
//					}
//				}
//	//			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
//	//			i++;
//		}

	u_int negAdjSize = 0;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mask[pos] & NEGINDONEBAND){
			TentativeGridPoint tp(0.f, i, j, k);
			vector<TentativeGridPoint> neigbor;
			tp.Neigbor(neigbor, DimX, DimY, DimZ);
			for(int n=0;n<neigbor.size();n++){
				int ii = neigbor[n].ii;
				int jj = neigbor[n].jj;
				int kk = neigbor[n].kk;
				if(!(mask[INDEX(ii,jj,kk)] & NEGINDONEBAND) ){ // if this neighbor is not in accepted band
					if(!(mask[INDEX(ii,jj,kk)] & NEGINTRIALBAND)){ // if this neighbor is not in close band
						if(phi[POS(neigbor[n])] < 0.f){
							mask[INDEX(ii,jj,kk)] |= NEGINTRIALBAND;
							++negAdjSize;
						}
					}
				}
			}
		}
	END_FOR
	num_points=0;
	printf("num_points=%d, adjacent points = %d \n",
			num_points, negAdjSize);
	if(negAdjSize == 0){
		return;
	}
//		for(pos=negband.begin(); pos!=negband.end(); ++pos){
//			KnownPoint &p = pos->second;
//			//printf("i = %d, j = %d, k = %d, phi = %f\n", p.ii, p.jj, p.kk, phi[POS(p)]);
//			phi[POS(p)] *= -1;
//		}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mask[pos] & NEGINDONEBAND)
			phi[pos] *= -1;
	END_FOR
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mask[pos] & NEGINTRIALBAND){
			KnownPoint p(i,j,k);
			float tentativeValue = UpdateTentativeValue(p, pos, mask, (char)NEGINDONEBAND, voxel, phi);
			TentativeGridPoint tp(tentativeValue, i, j, k);
			heap.insert(tp, tentativeValue);
		}
	END_FOR
//		for(pos=negadj.begin(); pos!=negadj.end(); ++pos){
//				KnownPoint &p = pos->second;
//				//printf("i = %d, j = %d, k = %d, phi = %f\n", p.ii, p.jj, p.kk, phi[POS(p)]);
//				u_int index = pos->first;
//				float tatitiveValue = UpdateTentativeValue(p, index, negband, voxel, phi);
//	//			float tatitiveValue = UpdateTentativeValue(p, index, v);
//				TentativeGridPoint tp(tatitiveValue, p);
//				heap.insert(tp, tatitiveValue);
//	//			++i;
//	//			printf("i = %d, index = %d \n", i, index);
//		}
//		printf("Tentative band ready\n");
	do{
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		++num_points;
		u_int indexDoneBand = Index(minTentative);
//		printf("i = %d \n", i);
		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
		phi[indexDoneBand] = minValue;
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
//		negband.insert( make_pair(Index(minTentative), t) );
//		negadj.erase(Index(minTentative));
		mask[indexDoneBand] |= NEGINDONEBAND;
		mask[indexDoneBand] ^= NEGINTRIALBAND;
		//if(!v.IsDone(minTentative))
		//	v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, SET, DONE);
		//if(v.IsClose(minTentative))
		//	v.UpdateCellFastMarchType(minTentative.ii, minTentative.jj, minTentative.kk, CLEAR, CLOSE);
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
		vector<TentativeGridPoint> neigbor;
		minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
		for(int i=0;i<neigbor.size();i++){
			u_int indexTrialBand =  Index(neigbor[i]);
			if(!(mask[indexTrialBand] & NEGINDONEBAND)){ // not currently in done band
//			foundinband = negband.find(Index(neigbor[i]));
//			if(foundinband == negband.end()){ // not currently in the band
//	    	if(!v.IsDone(neigbor[i])){
//	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), v);
//	    		if(neigbor[i].value == INFINITY){
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
//	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    		}
//				found = negadj.find(Index(neigbor[i]));
//				if(found == negadj.end()){ // not currently in the heap
				if(!(mask[indexTrialBand] & NEGINTRIALBAND)){ // not currently in trial band
//	    		if(!v.IsClose(neigbor[i])){ // not currently in the heap
					if(phi[indexTrialBand] < 0.f){
//						KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//						negadj.insert(make_pair(Index(neigbor[i]), tmp));
		//	    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
						mask[indexTrialBand] |= NEGINTRIALBAND;
						neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), indexTrialBand,
								mask, (char) NEGINDONEBAND, voxel, phi);
						heap.insert(neigbor[i], neigbor[i].value);
						//v.UpdateCellFastMarchType(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, SET, CLOSE);
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
					neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), indexTrialBand,
									mask, (char) NEGINDONEBAND, voxel, phi);
					heap.update(neigbor[i], neigbor[i].value);	// just update values
																// according to newly added point
				}
			}
		}
		neigbor.clear();
//	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
	} while(!heap.empty());

	printf("num_points=%d \n", num_points);

	// multiply by -1 to get negative distance in phi < 0 area
//	for(pos=negband.begin(); pos!=negband.end(); ++pos){
//		KnownPoint &p = pos->second;
//		phi[POS(p)] *= -1;
//	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mask[pos] & NEGINDONEBAND)
			phi[pos] *= -1;
	END_FOR

}




void VectorField3D::FastMarching(map<u_int, KnownPoint> &posband,
 								 map<u_int, KnownPoint> &negband,
								 const Voxel &voxel, float *phi,
								 const float *phi_obj,
								 const char *valid){
	// initialize grid points within the band of the interface
	// already initialized?

	float delta = voxel.VoxelDelta();
	float delta2 = 3 * delta;
	int num_points=0;
	MinHeap<TentativeGridPoint, float> heap;
	map<u_int, KnownPoint> posadj;
	map<u_int, KnownPoint> negadj;
	map<u_int, KnownPoint>::iterator pos;


	//Process grid points at phi > 0
	for(pos=posband.begin(); pos!=posband.end(); ++pos){
			KnownPoint &p = pos->second;
			u_int index = pos->first;
			TentativeGridPoint tp(0.f, p);
			vector<TentativeGridPoint> neigbor;
		    tp.Neigbor(neigbor, DimX, DimY, DimZ);
	    	map<u_int, KnownPoint>::iterator foundpos;
	    	map<u_int, KnownPoint>::iterator foundpos1;
			for(int i=0;i<neigbor.size();i++){
				foundpos1 = posband.find(Index(neigbor[i]));
				if(foundpos1==posband.end()){ // if this neighbor is not in accepted band
					foundpos = posadj.find(Index(neigbor[i]));
					if(foundpos==posadj.end()){ // if this neighbor is not in close band
//						if(!(neigbor[i].ii < x1 && neigbor[i].ii > x0 &&
//						     neigbor[i].jj < y1 && neigbor[i].jj > y0 &&
//					   	     neigbor[i].kk < z1 && neigbor[i].kk > z0)){
						if(valid[INDEX(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk)] &&
								phi[POS(neigbor[i])] > 0.f){
							KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
							posadj.insert(make_pair(Index(neigbor[i]),tmp)); // add it to close band
						}
					}
				}
			}
//			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
//			i++;
	}
	printf("pos band size = %ld, adj_points = %ld \n", posband.size(), posadj.size());

//	map<u_int, KnownPoint>::iterator foundpos;
//	foundpos = posadj.find(INDEX(I,J,K));
//	if(foundpos == posadj.end())
//		printf("(%d,%d,%d) is not in adjacent band \n", I,J,K);
//	else
//		printf("(%d,%d,%d) is in adjacent band band \n", I,J,K);

	for(pos=posadj.begin(); pos!=posadj.end(); ++pos){
			KnownPoint &p = pos->second;
			u_int index = pos->first;
			float tentativeValue = UpdateTentativeValue(p, index, posband, voxel, phi);
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
//		if(minValue > EXTRAPOLATE_VEL_LIMIT)
//			break;
		if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K){
			printf("old phi = %f\n",phi[POS(t)]);
		}
		phi[POS(t)] = minValue;
//		if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K)
////		if(minValue < 3.f)
//			printf("point(%d, %d, %d) phi = %f \n",
//					minTentative.ii,minTentative.jj,minTentative.kk, phi[POS(t)]);
//		if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//			printf("new phi = %f\n",phi[POS(t)]);
//		}
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
	    posband.insert( make_pair(Index(minTentative), t) );
		posadj.erase(Index(minTentative));
		//printf("heap size = %d, pos_adj size = %d\n", heap.size(), posadj.size());
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
	    vector<TentativeGridPoint> neigbor;
	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
	    for(int i=0;i<neigbor.size();i++){
	    	map<u_int, KnownPoint>::iterator foundinband = posband.find(Index(neigbor[i]));
	    	if(foundinband == posband.end()){ // not currently in the band
	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), posband, voxel, phi);
//	    		if(neigbor[i].value == INFINITY){
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
//	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
//	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
//	    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    		}
		    	map<u_int, KnownPoint>::iterator found = posadj.find(Index(neigbor[i]));
	    		if(found == posadj.end()){ // not currently in the heap
	    			if(valid[INDEX(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk)] &&
	    					phi[Index(neigbor[i])] > 0.f){
			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
			    		posadj.insert(make_pair(Index(neigbor[i]), tmp));
		//		    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
			    		heap.insert(neigbor[i], neigbor[i].value);
			    		if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
			    			printf("\n add close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
			    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
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
	    			if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
		    			printf("\n update close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
		    				neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
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
//	if(!heap.empty())
//		heap.clear();
////	foundpos = posband.find(INDEX(I,J,K));
////	if(foundpos == posband.end())
////		printf("(%d,%d,%d) is not in final band \n", I,J,K);
////	else
////		printf("(%d,%d,%d) is in final band \n", I,J,K);
//	for(pos=negband.begin(); pos!=negband.end(); ++pos){
//			KnownPoint &p = pos->second;
//			u_int index = pos->first;
//			TentativeGridPoint tp(0.f, p);
//			vector<TentativeGridPoint> neigbor;
//		    tp.Neigbor(neigbor, DimX, DimY, DimZ);
//	    	map<u_int, KnownPoint>::iterator foundpos;
//	    	map<u_int, KnownPoint>::iterator foundpos1;
//			for(int i=0;i<neigbor.size();i++){
//				foundpos1 = negband.find(Index(neigbor[i]));
//				if(foundpos1==negband.end()){ // if this neighbor is not in accepted band
//					foundpos = negadj.find(Index(neigbor[i]));
//					if(foundpos==negadj.end()){ // if this neighbor is not in close band
////						if(!(neigbor[i].ii < x1 && neigbor[i].ii > x0 &&
////						     neigbor[i].jj < y1 && neigbor[i].jj > y0 &&
////					   	     neigbor[i].kk < z1 && neigbor[i].kk > z0)){
//						if(valid[INDEX(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk)] &&
//								phi[POS(neigbor[i])] < 0.f){
//							KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//							negadj.insert(make_pair(Index(neigbor[i]),tmp)); // add it to close band
//						}
//					}
//				}
//			}
////			printf("i = %d, phi = %f \n", index, phi[POS(p)]);
////			i++;
//	}
//	printf("pos band size = %ld, adj_points = %ld \n",negband.size(), negadj.size());
//	for(pos=negband.begin(); pos!=negband.end(); ++pos){
//		KnownPoint &p = pos->second;
//		//printf("i = %d, j = %d, k = %d, phi = %f\n", p.ii, p.jj, p.kk, phi[POS(p)]);
//		phi[POS(p)] *= -1;
//	}
//	if(negadj.size() == 0){
//		return;
//	}
//	num_points = 0;
//	for(pos=negadj.begin(); pos!=negadj.end(); ++pos){
//		KnownPoint &p = pos->second;
//		u_int index = pos->first;
//		float tentativeValue = UpdateTentativeValue(p, index, negband, voxel, phi);
////			printf("pos band phi = %f \n", tentativeValue);
//		TentativeGridPoint tp(tentativeValue, p);
//		heap.insert(tp, tentativeValue);
//	}
//	do{
//		TentativeGridPoint minTentative = heap.extract_min(minValue);
//		++num_points;
////		printf("i = %d \n", i);
//		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
//		if(minValue > delta2)
//			break;
////		if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
////			printf("old phi = %f\n",phi[POS(t)]);
////		}
//		phi[POS(t)] = minValue;
////		if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K)
//////		if(minValue < 3.f)
////			printf("point(%d, %d, %d) phi = %f \n",
////					minTentative.ii,minTentative.jj,minTentative.kk, phi[POS(t)]);
////		if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
////			printf("new phi = %f\n",phi[POS(t)]);
////		}
////		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
//	    negband.insert( make_pair(Index(minTentative), t) );
//		negadj.erase(Index(minTentative));
//		//printf("heap size = %d, pos_adj size = %d\n", heap.size(), posadj.size());
//		//printf("# of points = %d, min value = %f\n", num_points, minValue);
//	    vector<TentativeGridPoint> neigbor;
//	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
//	    for(int i=0;i<neigbor.size();i++){
//	    	map<u_int, KnownPoint>::iterator foundinband = negband.find(Index(neigbor[i]));
//	    	if(foundinband == negband.end()){ // not currently in the band
//	    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(), Index(neigbor[i]), negband, voxel, phi);
////	    		if(neigbor[i].value == INFINITY){
////	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, neigbor[i].value,
////	    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
////	    			printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n", i, minTentative.value,
////	    			minTentative.ii, minTentative.jj, minTentative.kk );
////	    		}
//		    	map<u_int, KnownPoint>::iterator found = negadj.find(Index(neigbor[i]));
//	    		if(found == negadj.end()){ // not currently in the heap
//	    			if(valid[INDEX(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk)] &&
//	    					phi[Index(neigbor[i])] < 0.f){
//			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//			    		negadj.insert(make_pair(Index(neigbor[i]), tmp));
//		//		    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
//			    		heap.insert(neigbor[i], neigbor[i].value);
////			    		if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
////			    			printf("\n add close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
////			    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
//	    			}
//	    			else{
////	    				printf("\nfound negative phi in postive band \n");
////	    				printf("neigbor i = %d, value = %f, phi = %16.13f, id = (%d, %d, %d) \n",
////	    						i, neigbor[i].value,  phi[Index(neigbor[i])],
////			    			neigbor[i].ii, neigbor[i].jj, neigbor[i].kk );
////			    		printf("neigbor i = %d, value = %f, id = (%d, %d, %d) \n\n", i, minTentative.value,
////			    			minTentative.ii, minTentative.jj, minTentative.kk );
//	    			}
//	    		}
//	    		else{											// already in heap,
//	    			heap.update(neigbor[i], neigbor[i].value);	// just update values
//													    		// according to newly added point
////	    			if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
////		    			printf("\n update close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
////		    				neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
////					printf("update neigbor i = %d, value = %f \n",
////	    							i, neigbor[i].value );
//	    		}
//	    	}
//	    }
//	    neigbor.clear();
////	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
//	} while(!heap.empty());
//	printf("neg band size = %ld, num_points=%d \n", negband.size(), num_points);
//
//	// multiply by -1 to get negative distance in phi < 0 area
//	for(pos=negband.begin(); pos!=negband.end(); ++pos){
//		KnownPoint &p = pos->second;
//		phi[POS(p)] *= -1;
//	}
//	FOR_EACH_CELL
//		if(voxel.InAir(i,j,k) && phi[INDEX(i,j,k)] > 0.f
//				&& phi[INDEX(i,j,k)] < 3*delta){
//			if( phi[INDEX(i+1,j,k)]-phi[INDEX(i,j,k)] > 1.e-3 &&
//				phi[INDEX(i-1,j,k)]-phi[INDEX(i,j,k)] > 1.e-3 &&
//				phi[INDEX(i,j+1,k)]-phi[INDEX(i,j,k)] > 1.e-3 &&
//				phi[INDEX(i,j-1,k)]-phi[INDEX(i,j,k)] > 1.e-3 &&
//				phi[INDEX(i,j,k+1)]-phi[INDEX(i,j,k)] > 1.e-3 &&
//				phi[INDEX(i,j,k-1)]-phi[INDEX(i,j,k)] > 1.e-3  ){
//					printf("init phi = %f at (%d,%d,%d)\n",phi[INDEX(i,j,k)], i, j, k);
//					printf("init phi-x = %f phi+x = %f\n",phi[INDEX(i-1,j,k)],phi[INDEX(i+1,j,k)]);
//					printf("init phi-y = %f phi+y = %f\n",phi[INDEX(i,j-1,k)],phi[INDEX(i,j+1,k)]);
//					printf("init phi-z = %f phi+z = %f\n",phi[INDEX(i,j,k-1)],phi[INDEX(i,j,k+1)]);
//					exit(1);
//			}
//
//		}
//	END_FOR

}

void VectorField3D::
	CorrectLevelSet(const Voxel &voxel, const Point &pos,
			        float *phiNeg, float* phiPos,
	                const Particle &ps, char sign,
	                unsigned char vel_index) const{

	TentativeGridPoint tp = voxel.ContainsPointInVelCell(pos, vel_index);
	Point center = voxel.VelPosition(vel_index, tp.ii, tp.jj, tp.kk);
	if(sign < 0){
		float phi_center = ps.EvaluatePhiP(center);
		phiNeg[INDEX(tp.ii, tp.jj, tp.kk)] = min(phi_center, phiNeg[INDEX(tp.ii, tp.jj, tp.kk)]);
	}
	else{
		float phi_center = -ps.EvaluatePhiP(center);
		phiPos[INDEX(tp.ii, tp.jj, tp.kk)] = max(phi_center, phiPos[INDEX(tp.ii, tp.jj, tp.kk)]);
	}

}

void VectorField3D::ErrorCorrection(const Voxel &voxel,
							list<Particle> &particles
#ifdef POSITIVE_PARTICLES
							,list<Particle> &posParticles
#endif
							){

	float delta = voxel.VoxelDelta();

	list<Particle>::iterator iter_particle;

//	float *phiuNeg = new float[DimX*DimY*DimZ];
//	float *phiuPos = new float[DimX*DimY*DimZ];
	float *phiuNeg = u0;
	float *phiuPos = u1;
	SetEqual(phi_u, phiuNeg);
	SetEqual(phi_u, phiuPos);
//	float *phivNeg = new float[DimX*DimY*DimZ];
//	float *phivPos = new float[DimX*DimY*DimZ];
	float *phivNeg = v0;
	float *phivPos = v1;
	SetEqual(phi_v, phivNeg);
	SetEqual(phi_v, phivPos);
//	float *phiwNeg = new float[DimX*DimY*DimZ];
//	float *phiwPos = new float[DimX*DimY*DimZ];
	float *phiwNeg = w0;
	float *phiwPos = w1;
	SetEqual(phi_w, phiwNeg);
	SetEqual(phi_w, phiwPos);


	printf("start error correction ... \n");

	//for(int m=0; m<particles.size(); ++m){
	for(iter_particle = particles.begin();
		iter_particle != particles.end();
		++iter_particle){
//		Point pos = particles[m].Position();
//		float radius = particles[m].Radius();
		Particle &ps = *iter_particle;
		Point pos = ps.Position();
		float radius = ps.Radius();

//		TentativeGridPoint tp = voxel.ContainsPoint(pos);
//		if(tp.ii == I && tp.jj == J && tp.kk == K){
//			if(source->IsSourceCell(tp.ii, tp.jj, tp.kk))
//				printf(" (%d,%d,%d) is source cell \n", tp.ii, tp.jj, tp.kk);
//			else
//				printf(" (%d,%d,%d) is not source cell \n", tp.ii, tp.jj, tp.kk);
//		}
//		if(tp.ii == I && tp.jj == J && tp.kk == K){
//			printf(" negative particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
//					pos.x, pos.y, pos.z, phi_i, radius);
//		}
		float phi_i = TriInterp(voxel, pos, phi_u, 1);
		if( phi_i > radius )
			CorrectLevelSet(voxel, pos, phiuNeg, phiuPos, ps, -1, 1);
		phi_i = TriInterp(voxel, pos, phi_v, 2);
		if( phi_i > radius )
			CorrectLevelSet(voxel, pos, phivNeg, phivPos, ps, -1, 2);
		phi_i = TriInterp(voxel, pos, phi_w, 3);
		if( phi_i > radius )
			CorrectLevelSet(voxel, pos, phiwNeg, phiwPos, ps, -1, 3);


	}
//	printf("negative particles error correction finished  \n");

#ifdef POSITIVE_PARTICLES
		for(iter_particle = posParticles.begin();
			iter_particle != posParticles.end();
			++iter_particle){
	//		Point pos = particles[m].Position();
	//		float radius = particles[m].Radius();
			Particle &ps = *iter_particle;
			Point pos = ps.Position();
			float radius = ps.Radius();

//			TentativeGridPoint tp = voxel.ContainsPoint(pos);
//			if(tp.ii == I && tp.jj == J && tp.kk == K){
//				printf(" positive particle at (%f, %f, %f) phi_i = %f, radius = %f \n",
//						pos.x, pos.y, pos.z, phi_i, radius);
//			}
			float phi_i = TriInterp(voxel, pos, phi_u, 1);
			if( -phi_i > radius )
				CorrectLevelSet(voxel, pos, phiuNeg, phiuPos, ps, 1, 1);
			phi_i = TriInterp(voxel, pos, phi_v, 2);
			if( -phi_i > radius )
				CorrectLevelSet(voxel, pos, phivNeg, phivPos, ps, 1, 2);
			phi_i = TriInterp(voxel, pos, phi_w, 3);
			if( -phi_i > radius )
				CorrectLevelSet(voxel, pos, phiwNeg, phiwPos, ps, 1, 3);


		}
//		printf("positive particles error correction finished  \n");
#endif

#ifdef POSITIVE_PARTICLES
		FOR_EACH_CELL
			phi_u[INDEX(i,j,k)] = fabs(phiuPos[INDEX(i,j,k)]) <= fabs(phiuNeg[INDEX(i,j,k)]) ?
								phiuPos[INDEX(i,j,k)] : phiuNeg[INDEX(i,j,k)];
			phi_v[INDEX(i,j,k)] = fabs(phivPos[INDEX(i,j,k)]) <= fabs(phivNeg[INDEX(i,j,k)]) ?
								phivPos[INDEX(i,j,k)] : phivNeg[INDEX(i,j,k)];
			phi_w[INDEX(i,j,k)] = fabs(phiwPos[INDEX(i,j,k)]) <= fabs(phiwNeg[INDEX(i,j,k)]) ?
								phiwPos[INDEX(i,j,k)] : phiwNeg[INDEX(i,j,k)];
		END_FOR
#endif

//	delete [] phiuNeg; delete [] phiuPos;
//	delete [] phivNeg; delete [] phivPos;
//	delete [] phiwNeg; delete [] phiwPos;
}


/*void VectorField3D::UpdateVoxel(Voxel &voxel) const{

//	u_int NP = 0;
//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi_c[INDEX(i,j,k)] < 0.f)
//				++NP;
//	END_FOR
//	printf("Before UpdateVoxel(), there are %u liquid points \n", NP);

	float cell_phi[6];

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
			cell_phi[0] = phi_u[INDEX(i,j,k)];
			cell_phi[1] = phi_u[INDEX(i-1,j,k)];
			cell_phi[2] = phi_v[INDEX(i,j,k)];
			cell_phi[3] = phi_v[INDEX(i,j-1,k)];
			cell_phi[4] = phi_w[INDEX(i,j,k)];
			cell_phi[5] = phi_w[INDEX(i,j,k-1)];
//			cell_phi[1] = 0.5f*(phi[INDEX(i-1,j,k)]+phi[INDEX(i,j,k)]);
//			cell_phi[2] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i+1,j,k)]);
//			cell_phi[3] = 0.5f*(phi[INDEX(i,j-1,k)]+phi[INDEX(i,j,k)]);
//			cell_phi[4] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j+1,k)]);
//			cell_phi[5] = 0.5f*(phi[INDEX(i,j,k-1)]+phi[INDEX(i,j,k)]);
//			cell_phi[6] = 0.5f*(phi[INDEX(i,j,k)]+phi[INDEX(i,j,k+1)]);
			int negative = 0;
			for(int m=0; m<6; ++m){
				if(cell_phi[m] <= 0.f)
					negative++;
			}
			if(negative == 6){ // liquid cell
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
			else{  // surface cell
//				if(cell_phi[0] < 0.f){   // surface cell
//					if(voxel.InAir(i,j,k))
//						voxel.UpdateCellType(i, j, k, CLEAR, EMPTY);
//					if(!voxel.InSurface(i,j,k))
//						voxel.UpdateCellType(i, j, k, SET, SURFACE);
//					if(!voxel.InLiquid(i,j,k))
//						voxel.UpdateCellType(i, j, k, SET, LIQUID);
//				}
//				else{   // air cell
//					if(!voxel.InAir(i,j,k))
//						voxel.UpdateCellType(i, j, k, SET, EMPTY);
//					if(voxel.InSurface(i,j,k))
//						voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
//					if(voxel.InLiquid(i,j,k))
//						voxel.UpdateCellType(i, j, k, CLEAR, LIQUID);
//				}
				if(voxel.InAir(i,j,k))
					voxel.UpdateCellType(i, j, k, CLEAR, EMPTY);
				if(!voxel.InSurface(i,j,k))
					voxel.UpdateCellType(i, j, k, SET, SURFACE);
				if(!voxel.InLiquid(i,j,k))
					voxel.UpdateCellType(i, j, k, SET, LIQUID);
			}

		}

	END_FOR

	u_int NS = 0;

	FOR_EACH_CELL
//		if(voxel.InSurface(i,j,k)){
//			bool realSurfaceCell = false;
//			TentativeGridPoint p(0.f,i,j,k);
//			vector<TentativeGridPoint> neighbors;
//			p.Neigbor(neighbors, DimX, DimY, DimZ);
//			for(int m=0; m<neighbors.size(); m++){
//				if(voxel.InAir(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk))
//					realSurfaceCell = true;
//			}
//			if(!realSurfaceCell)
//				voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
//		}
		if(voxel.InSurface(i,j,k)){
			if(phi_c[INDEX(i,j,k)] <= 0.f)
				voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
		}
		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
				if(phi_c[INDEX(i,j,k)] > 0.f){
					voxel.UpdateCellType(i, j, k, CLEAR, LIQUID);
					voxel.UpdateCellType(i, j, k, SET, EMPTY);
					printf("at (%d, %d, %d) a trapped air cell \n", i,j,k);
					if(voxel.InSurface(i,j,k))
						voxel.UpdateCellType(i, j, k, CLEAR, SURFACE);
				}
		}
		if(voxel.InAir(i,j,k)){
			if(phi_c[INDEX(i,j,k)] < 0.f){
				printf("at (%d, %d, %d) an isolated liquid cell \n", i,j,k);
				voxel.UpdateCellType(i, j, k, CLEAR, EMPTY);
				voxel.UpdateCellType(i, j, k, SET, LIQUID);
			}
		}

//		if(voxel.InSource(i,j,k))
//			++NS;

		if(i==I && j==J && k==K){
			cell_phi[0] = phi_c[INDEX(i,j,k)];
			cell_phi[1] = phi_u[INDEX(i,j,k)];
			cell_phi[2] = phi_v[INDEX(i,j,k)];
			cell_phi[3] = phi_w[INDEX(i,j,k)];
			if(voxel.InLiquid(I,J,K))
			  printf(" \n (%d, %d, %d) in liquid \n\n", I, J, K);
		    if(voxel.InAir(I,J,K))
			  printf(" \n (%d, %d, %d) in air \n\n", I, J, K);
		    if(voxel.InSurface(I,J,K))
			  printf("\n  (%d, %d, %d) in surface \n\n", I, J, K);
		    if(voxel.InSolid(I,J,K))
		    	printf(" \n (%d, %d, %d) in solid \n\n", I, J, K);
		    if(voxel.InSource(I,J,K))
		    	printf(" \n (%d, %d, %d) in source \n\n", I, J, K);
	    	printf(" \n cell phi = (%f, %f, %f, %f) \n\n",
	    		cell_phi[0], cell_phi[1], cell_phi[2], cell_phi[3]);
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
		if(voxel.InSurface(i,j,k)){
			if( phi_c[INDEX(i,j,k)] <= 0.f ){
				printf("wrong! surface cell (%d, %d, %d) has phi (%f) <= 0 \n", i,j,k, phi_c[INDEX(i,j,k)]);
				exit(1);
			}

		}
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
//	printf("There are %u source points \n", NS);

//	u_int NM = 0;
//	NP = 0;
//	FOR_EACH_CELL
//		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
//			++NM;
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi_c[INDEX(i,j,k)] <= 0.f)
//			++NP;
//	END_FOR
//	printf(" After UpdateVoxel(), there are %u and %u liquid points \n", NM, NP);


}*/

void VectorField3D::UpdateVoxel(Voxel &voxel) const{

//	u_int NP = 0;
//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi_c[INDEX(i,j,k)] < 0.f)
//				++NP;
//	END_FOR
//	printf("Before UpdateVoxel(), there are %u liquid points \n", NP);

	float cell_phi[6];

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
			cell_phi[0] = phi_u[INDEX(i,j,k)];
			cell_phi[1] = phi_u[INDEX(i-1,j,k)];
			cell_phi[2] = phi_v[INDEX(i,j,k)];
			cell_phi[3] = phi_v[INDEX(i,j-1,k)];
			cell_phi[4] = phi_w[INDEX(i,j,k)];
			cell_phi[5] = phi_w[INDEX(i,j,k-1)];
			int negative = 0;
			for(int m=0; m<6; ++m){
				if(cell_phi[m] <= 0.f)
					negative++;
			}
			if(phi_c[INDEX(i,j,k)] <= 0.f){ // liquid cell
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
			else{  // surface cell
				if(voxel.InAir(i,j,k))
					voxel.UpdateCellType(i, j, k, CLEAR, EMPTY);
				if(!voxel.InSurface(i,j,k))
					voxel.UpdateCellType(i, j, k, SET, SURFACE);
				if(!voxel.InLiquid(i,j,k))
					voxel.UpdateCellType(i, j, k, SET, LIQUID);
			}

		}

	END_FOR

	u_int NS = 0;

	FOR_EACH_CELL
		if(i==I && j==J && k==K){
			cell_phi[0] = phi_c[INDEX(i,j,k)];
			cell_phi[1] = phi_u[INDEX(i,j,k)];
			cell_phi[2] = phi_v[INDEX(i,j,k)];
			cell_phi[3] = phi_w[INDEX(i,j,k)];
			if(voxel.InLiquid(I,J,K))
			  printf(" \n (%d, %d, %d) in liquid \n\n", I, J, K);
		    if(voxel.InAir(I,J,K))
			  printf(" \n (%d, %d, %d) in air \n\n", I, J, K);
		    if(voxel.InSurface(I,J,K))
			  printf("\n  (%d, %d, %d) in surface \n\n", I, J, K);
		    if(voxel.InSolid(I,J,K))
		    	printf(" \n (%d, %d, %d) in solid \n\n", I, J, K);
		    if(voxel.InSource(I,J,K))
		    	printf(" \n (%d, %d, %d) in source \n\n", I, J, K);
	    	printf(" \n cell phi = (%f, %f, %f, %f) \n\n",
	    		cell_phi[0], cell_phi[1], cell_phi[2], cell_phi[3]);
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
		if(voxel.InSurface(i,j,k)){
			if( phi_c[INDEX(i,j,k)] <= 0.f ){
				printf("wrong! surface cell (%d, %d, %d) has phi (%f) <= 0 \n", i,j,k, phi_c[INDEX(i,j,k)]);
				exit(1);
			}

		}
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
//	printf("There are %u source points \n", NS);

//	u_int NM = 0;
//	NP = 0;
//	FOR_EACH_CELL
//		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
//			++NM;
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && phi_c[INDEX(i,j,k)] <= 0.f)
//			++NP;
//	END_FOR
//	printf(" After UpdateVoxel(), there are %u and %u liquid points \n", NM, NP);


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
//		   printf("minvalue larger than both phi = %f, phi_minus = %f, phi_plus = %f\n",
//				   minvalue, phi_minus, phi_plus);
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

/*void VectorField3D::ExtrapolateOneVelocity(const Voxel &voxel, float *phi_tmp, float *vel, int vel_index, float delta) {

	float minValue;
	MinHeap<TentativeGridPoint, float> heap;

	SetAirBoundary(voxel, vel_index, vel);

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && phi_tmp[INDEX(i,j,k)] > 0.f && phi_tmp[INDEX(i,j,k)] < 3*delta){
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
		}
	END_FOR

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
	    		phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.RightNeighbor(neighbors[m])){
	    		phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.BackNeighbor(neighbors[m])){
	    		phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.FrontNeighbor(neighbors[m])){
	    		phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.BottomNeighbor(neighbors[m])){
	    		phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.TopNeighbor(neighbors[m])){
	    		phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
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
	    	vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
	    }
	    else{
//	    	printf("Error in VectorField3D::ExtractVelocity ! \n ");
//	    	printf(" ii = %d, jj = %d, kk = %d \n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk);
//	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
//	    			minValue, phi_minus_x, phi_plus_x);
//	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
//		    			minValue, phi_minus_y, phi_plus_y);
//	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
//	    		    	minValue, phi_minus_z, phi_plus_z);
	    }
//	    if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K && vel_index == 3){
////	    if(vel_index == 3){
//	    	printf("after interp: at i = %d, j = %d, k = %d vel = %f \n\n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]);
////	    	printf("before interp: A = %f, B = %f, C = %f, a = %f, b = %f, c = %f \n ",
////	    			A, B, C, a, b, c);
//////	    	printf(" ii = %d, jj = %d, kk = %d \n",
//////	    			minTentative.ii, minTentative.jj, minTentative.kk);
////	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
////	    			minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
////	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
////		    			minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
////	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
////	    		    	minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
//	    }

	} while(!heap.empty());

	SetSolidBoundary(voxel, vel_index, vel);
}*/

void VectorField3D::ExtrapolateOneVelocity(const Voxel &voxel, float *phi_tmp, float *phi_obj, float *vel, int vel_index, float delta) {

	float minValue;
	MinHeap<TentativeGridPoint, float> heap;

//	SetAirBoundary(voxel, vel_index, vel);
//	ApplySourceTerm();

	char *valid = new char[DimX*DimY*DimZ];
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	SetValidPoint(voxel, valid, vel_index);


	FOR_EACH_CELL
		if((voxel.InAir(i,j,k) || voxel.InSurface(i,j,k))
		    && phi_c[INDEX(i,j,k)] < EXTRAPOLATE_VEL_LIMIT &&
		   phi_tmp[INDEX(i,j,k)] > 0.f && valid[INDEX(i,j,k)] ){
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
		}
	    if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= EXTRAPOLATE_VEL_LIMIT
	      && valid[INDEX(i,j,k)]){
	    	TentativeGridPoint it(0.f, i, j, k);
	    	vector<TentativeGridPoint> neighbors;
	        it.Neigbor(neighbors, DimX, DimY, DimZ);
	        for(int m=0; m<neighbors.size(); m++){
	        	if(it.RightNeighbor(neighbors[m])
	        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] < EXTRAPOLATE_VEL_LIMIT
	        	  ){
	        		if(vel_index == 1){
	        			KnownPoint tmp(i,j,k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
						if(phi_tmp[INDEX(i,j,k)] < 0.f){
							printf("Error! shouldn't extrapolate this one\n");
							exit(1);
						}
	        		}
        		}
	        	if(it.FrontNeighbor(neighbors[m])
	        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] < EXTRAPOLATE_VEL_LIMIT
	        	  ){
	        		if(vel_index == 2){
	        			KnownPoint tmp(i,j,k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
						if(phi_tmp[INDEX(i,j,k)]< 0.f){
							printf("Error! shouldn't extrapolate this one\n");
							exit(1);
						}
	        		}
        		}
	        	if(it.TopNeighbor(neighbors[m])
	        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] < EXTRAPOLATE_VEL_LIMIT
	        	 ){
	        		if(vel_index == 3){
	        			KnownPoint tmp(i,j,k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
						if(phi_tmp[INDEX(i,j,k)]< 0.f){
							printf("Error! shouldn't extrapolate this one\n");
							exit(1);
						}
	        		}
        		}
	        }

	    }
	END_FOR

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
	    		phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			phi_minus_x = INFINITY;
	    	}
	    	if(minTentative.RightNeighbor(neighbors[m])){
	    		phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			 phi_plus_x = INFINITY;
//	    		if(vel_index == 1 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii+1,neighbors[m].jj,neighbors[m].kk))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_x = INFINITY;
	    	}
	    	if(minTentative.BackNeighbor(neighbors[m])){
	    		phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_minus_y = INFINITY;
	    	}
	    	if(minTentative.FrontNeighbor(neighbors[m])){
	    		phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			 phi_plus_y = INFINITY;
//	    		if(vel_index == 2 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj+1,neighbors[m].kk))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_y = INFINITY;
	    	}
	    	if(minTentative.BottomNeighbor(neighbors[m])){
	    		phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			  phi_minus_z = INFINITY;
	    	}
	    	if(minTentative.TopNeighbor(neighbors[m])){
	    		phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			  phi_plus_z = INFINITY;
//	    		if(vel_index == 3 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk+1))
    			if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_z = INFINITY;
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
	    if(isnan(a) || isinf(a) || isnan(b) || isinf(b) || isnan(c) || isinf(c)){
	    	printf("Error in VectorField3D::ExtractOneVelocity vel_index = %d! \n ",
	    			vel_index);
	    	printf(" ii = %d, jj = %d, kk = %d \n",
	    			minTentative.ii, minTentative.jj, minTentative.kk);
	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
	    			minValue, phi_minus_x, phi_plus_x);
	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
		    			minValue, phi_minus_y, phi_plus_y);
	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
	    		    	minValue, phi_minus_z, phi_plus_z);
	    	printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
	    		    			vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
	    		    			vel_minus_z, vel_plus_z);
	    	exit(1);
	    }
	    if( (A+B+C) > 1.e-10 ){
//	    	if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K && vel_index == 1){
//	    	if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K){
//		    	float gradphi = sqrt(phi_x*phi_x+phi_y*phi_y+phi_z*phi_z);
//		    	if( fabs(gradphi-1.f) > 1.e-3f ){
//		    		printf("ii = %d, jj = %d, kk = %d, vel_index= %d, phi_x = %f, phi_y = %f, phi_z = %f, |grad(phi)| = %f \n ",
//		    			minTentative.ii, minTentative.jj, minTentative.kk,
//		    			vel_index, phi_x, phi_y, phi_z, sqrt(phi_x*phi_x+phi_y*phi_y+phi_z*phi_z));
//		    		printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
//			    			minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
//			    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
//				    			minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
//			    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
//			    		    	minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
//		    	}
//	    	}
	    	vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
	    }
	    else{
//	    	for(int m=0; m<neighbors.size(); m++){
//	    		if(source->IsSourceCell(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)){
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] =
//	    				vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    			break;
//	    		}
//	    	}
//	    	if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//	    	if(minTentative.kk == K){
//		    	printf("Error in VectorField3D::ExtractOneVelocity vel_index = %d\n ",
//		    			vel_index );
//		    	printf(" ii = %d, jj = %d, kk = %d \n",
//		    			minTentative.ii, minTentative.jj, minTentative.kk);
//		    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
//		    			minValue, phi_minus_x, phi_plus_x);
//		    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
//			    			minValue, phi_minus_y, phi_plus_y);
//		    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
//		    		    	minValue, phi_minus_z, phi_plus_z);
//		    	printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
//		    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//		    		    			vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
//		    		    			vel_minus_z, vel_plus_z);
//		    	exit(1);
//	    	}
	    }
	    if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//	    if(vel_index == 3){
//	    if( isnan(vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]) ||
//	    	isinf(vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]) ){
//	    if(phi_minus_x < 0.f && phi_plus_x < 0.f ||
//    	   phi_minus_y < 0.f && phi_plus_y < 0.f ||
//    	   phi_minus_z < 0.f && phi_plus_z < 0.f ) {
//	    if( phi_x == 0.f && A != 0.f || phi_y == 0.f && B != 0.f || phi_z == 0.f && C != 0.f ){
	    	float gradphi = sqrt(phi_x*phi_x+phi_y*phi_y+phi_z*phi_z);
	    	printf("(ExtrapolateOneVelocity) vel_index = %d, |grad(phi)| = %f\n ", vel_index,gradphi);
	    	printf("after interp: at i = %d, j = %d, k = %d vel = %f \n\n",
	    			minTentative.ii, minTentative.jj, minTentative.kk,
	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]);
	    	printf("before interp: A = %f, B = %f, C = %f, a = %f, b = %f, c = %f \n ",
	    			A, B, C, a, b, c);
	    	printf(" phix = %f, phiy = %f, phiz = %f \n",
	    			phi_x, phi_y, phi_z);
	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
	    			minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
		    			minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
	    		    	minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
//	    	exit(1);
	    }

	} while(!heap.empty());

	if(valid) delete [] valid;

//	SetSolidBoundary(voxel, vel_index, vel);
}



#ifdef SPMD
static void VisitNeighbors(const TentativeGridPoint &minTentative,
					float *phi_tmp, float *vel0, int DimX, int DimY, int DimZ,
		            float &phi_minus_x, float &phi_plus_x, float &phi_minus_y,
		            float &phi_plus_y, float &phi_minus_z, float &phi_plus_z,
		            float &vel_minus_x, float &vel_plus_x, float &vel_minus_y,
		            float &vel_plus_y, float &vel_minus_z, float &vel_plus_z){

	vector<TentativeGridPoint> neighbors;
	minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
	u_int pos;
//	int i = minTentative.ii;
//	int j = minTentative.jj;
//	int k = minTentative.kk;
	if(neighbors.size() == 6){
		pos = INDEX(minTentative.ii-1, minTentative.jj, minTentative.kk);
		phi_minus_x = phi_tmp[pos];
		vel_minus_x = vel0[pos];
		pos = INDEX(minTentative.ii+1, minTentative.jj, minTentative.kk);
		phi_plus_x = phi_tmp[pos];
		vel_plus_x  = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj-1, minTentative.kk);
		phi_minus_y = phi_tmp[pos];
		vel_minus_y = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj+1, minTentative.kk);
		phi_plus_y = phi_tmp[pos];
		vel_plus_y  = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj, minTentative.kk-1);
		phi_minus_z = phi_tmp[pos];
		vel_minus_z = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj, minTentative.kk+1);
		phi_plus_z = phi_tmp[pos];
		vel_plus_z  = vel0[pos];
	}
	else{
		for(int m=0; m<neighbors.size(); m++){
			u_int pos = INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
			if(minTentative.LeftNeighbor(neighbors[m])){
				phi_minus_x = phi_tmp[pos];
				vel_minus_x = vel0[pos];
				continue;
			}
			if(minTentative.RightNeighbor(neighbors[m])){
				phi_plus_x = phi_tmp[pos];
				vel_plus_x = vel0[pos];
				continue;
			}
			if(minTentative.BackNeighbor(neighbors[m])){
				phi_minus_y = phi_tmp[pos];
				vel_minus_y = vel0[pos];
				continue;
			}
			if(minTentative.FrontNeighbor(neighbors[m])){
				phi_plus_y = phi_tmp[pos];
				vel_plus_y = vel0[pos];
				continue;
			}
			if(minTentative.BottomNeighbor(neighbors[m])){
				phi_minus_z = phi_tmp[pos];
				vel_minus_z = vel0[pos];
				continue;
			}
			if(minTentative.TopNeighbor(neighbors[m])){
				phi_plus_z = phi_tmp[pos];
				vel_plus_z = vel0[pos];
				continue;
			}
		}
	}

	return;
}

static void MarchingOut(const TentativeGridPoint &minTentative,
						int DimX, int DimY, int DimZ,
						float dtau, float inv_delta2, float *phi_tmp,
						float *vel, float *vel0, int vel_index) {

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
	VisitNeighbors(minTentative, phi_tmp, vel0,
				   DimX, DimY, DimZ,
			      phi_minus_x, phi_plus_x, phi_minus_y,
			      phi_plus_y, phi_minus_z, phi_plus_z,
				  vel_minus_x, vel_plus_x, vel_minus_y,
				  vel_plus_y,  vel_minus_z, vel_plus_z);


	ProcessCoefficents(phi_minus_x, phi_plus_x, minValue,
					vel_minus_x, vel_plus_x, phi_x, A, a);
	ProcessCoefficents(phi_minus_y, phi_plus_y, minValue,
					vel_minus_y, vel_plus_y, phi_y, B, b);
	ProcessCoefficents(phi_minus_z, phi_plus_z, minValue,
					vel_minus_z, vel_plus_z, phi_z, C, c);


	if(phi_x != 0.f || phi_y != 0.f  || phi_z != 0.f){
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
		float mag = sqrtf(inv_delta2 * (phi_x * phi_x + phi_y * phi_y  + phi_z * phi_z));
		vel[pos] = vel0[pos] - dtau * inv_delta2 / mag * (phi_x * vel_x +  phi_y * vel_y + phi_z * vel_z);
	}
//	if(i == 5 && j == 5 && k == 4){
//		printf("(MarchingOut) vel_index = %d phi = %f\n",
//					vel_index, minValue);
//		printf("after interp: at i = %d, j = %d, k = %d, vel0 = %f, vel = %f \n\n",
//				minTentative.ii, minTentative.jj, minTentative.kk,
//				vel0[INDEX(i, j, k)], vel[INDEX(i, j, k)]);
//		printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
//				minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
//		printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
//					minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
//		printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
//					minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
//	}

}

#ifdef TBB
struct ExtrapolateObjectBody {
	float * const t1;
	float * const t2;
	float * const t3;
	map<u_int, KnownPoint> * const knownPoint;
	int vel_index;
	float delta;
	int DimX, DimY, DimZ;
	ExtrapolateObjectBody(map<u_int, KnownPoint> * const kp,
				float *s1, float *s2, float *s3,
			     int vl, float h, int nx, int ny, int nz)
		: t1(s1), t2(s2), t3(s3), knownPoint(kp),
		   vel_index(vl), delta(h), DimX(nx), DimY(ny), DimZ(nz){
	}
    void operator()( const tbb::blocked_range<int>& range ) const {
    	float *phi_tmp = t1;
    	float *vel = t2;
    	float *vel0 = t3;

		float delta2 = delta * delta;
		float inv_delta2 = 1.f / delta2;
		float dtau = delta / 2 ;

		map<u_int, KnownPoint>::iterator foundpos;

        int k_end = range.end();
        for( int k=range.begin(); k!=k_end; ++k ){
        	for(int j=0;j<DimY;j++) {
        		for(int i=0;i<DimX;i++){
        			u_int pos = INDEX(i,j,k);
        			foundpos = knownPoint->find(pos);
					if( foundpos == knownPoint->end() ){
    					if( phi_tmp[pos] >= 0.f && phi_tmp[pos] < LS_ADV_LIMIT ){
						KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[pos], tmp);
						MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp, vel, vel0, vel_index);
					}
        		}
        	}
        }
     }
  }
};
struct ExtrapolateObjectBody1 {
	float * const t1;
	float * const t2;
	float * const t3;
	float * const t4;
	const Voxel * const oxel;
	int vel_index;
	float delta;
	int DimX, DimY, DimZ;
	ExtrapolateObjectBody1(const Voxel * const vox, float *s1, float *s2, float *s3,
			     float *s4, int vl, float h, int nx, int ny, int nz)
		: t1(s1), t2(s2), t3(s3), t4(s4), voxel(vox),
		   vel_index(vl), delta(h), DimX(nx), DimY(ny), DimZ(nz){
	}
    void operator()( const tbb::blocked_range<int>& range ) const {
    	float *phi_tmp = t1;
    	float *vel = t2;
    	float *vel0 = t3;
    	float *phi_c = t4;

		float delta2 = delta * delta;
		float inv_delta2 = 1.f / delta2;
		float dtau = delta / 2 ;

        int k_end = range.end();
        for( int k=range.begin(); k!=k_end; ++k ){
        	for(int j=0;j<DimY;j++) {
        		for(int i=0;i<DimX;i++){
        			u_int pos = INDEX(i,j,k);
					if( phi_tmp[pos] >= 0.f && phi_tmp[pos] < LS_ADV_LIMIT ){
						KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[pos], tmp);
						MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp, vel, vel0, vel_index);
					}
					else{
						if(vel_index == 1){
							if( (!voxel->InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel->InSolid(i+1,j,k)) ||
								(voxel->InSolid(i,j,k) && !voxel->InSolid(i+1,j,k) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
								KnownPoint tmp(i, j, k);
								TentativeGridPoint tp(phi_tmp[pos], tmp);
								MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
											vel, vel0, vel_index);
							}
						}
						if(vel_index == 2){
							if( (!voxel->InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel->InSolid(i,j+1,k)) ||
								(voxel->InSolid(i,j,k) && !voxel->InSolid(i,j+1,k) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
								KnownPoint tmp(i, j, k);
								TentativeGridPoint tp(phi_tmp[pos], tmp);
								MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
											vel, vel0, vel_index);
							}
						}
						if(vel_index == 3){
							if( (!voxel->InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel->InSolid(i,j,k+1)) ||
								(voxel->InSolid(i,j,k) && !voxel->InSolid(i,j,k+1) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
								KnownPoint tmp(i, j, k);
								TentativeGridPoint tp(phi_tmp[pos], tmp);
								MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
										vel, vel0, vel_index);
							}
						}
					}
        	}
        }
     }
  }
};
#endif
#endif

#define EXTRAPOLATE_VEL_INTO_OBJECT_LIMIT 2*delta

#ifdef SPMD
void VectorField3D::ExtrapolateOneVelocityIntoObject(const Voxel &voxel,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoObject() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(EXTRAPOLATE_VEL_INTO_OBJECT_LIMIT / delta * 2));
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
		FOR_EACH_CELL
			if( phi_tmp[INDEX(i,j,k)] >= 0.f && phi_tmp[INDEX(i,j,k)] < EXTRAPOLATE_VEL_INTO_OBJECT_LIMIT){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
							vel, vel0, vel_index);
			}
		END_FOR
	}
    printf("\n End ExtrapolateOneVelocityIntoObject() ..... \n " );
}

#else
void VectorField3D::ExtrapolateOneVelocityIntoObject(const Voxel &voxel,
		float *phi_tmp, float *vel, int vel_index, float delta) {

	MinHeap<TentativeGridPoint, float> heap;
	printf("\n Start ExtrapolateOneVelocityIntoObject() ..... \n " );
	FOR_EACH_CELL
		if( phi_tmp[INDEX(i,j,k)] >= 0.f && phi_tmp[INDEX(i,j,k)] < EXTRAPOLATE_VEL_INTO_OBJECT_LIMIT){
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
		}
//		if( voxel.InSolid(i,j,k) ){
//			KnownPoint tmp(i, j, k);
//			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
//			heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
//		}
//		else{
//			if(vel_index == 1 && voxel.InSolid(i+1,j,k) ||
//			   vel_index == 2 && voxel.InSolid(i,j+1,k)	||
//			   vel_index == 3 && voxel.InSolid(i,j,k+1)	){
//				KnownPoint tmp(i, j, k);
//				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
//				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
//			}
//
//		}
	 END_FOR

//	printf("Extrapolating velocity: b = %d with %ld in heap \n ",
//			    vel_index, heap.size());

	UpdateCells(heap, vel_index, vel, phi_tmp);
	 printf("\n End ExtrapolateOneVelocityIntoObject() ..... \n " );
}
#endif

#ifdef SPMD
void VectorField3D::ExtrapolateOneAirVelocityIntoObject(const Voxel &voxel,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneAirVelocityIntoObject() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(EXTRAPOLATE_VEL_INTO_OBJECT_LIMIT / delta * 2));
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
		FOR_EACH_CELL
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			if(voxel.InSolid(i,j,k) && voxel.CloseToAir(tp)){
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
							vel, vel0, vel_index);
			}
		END_FOR
	}
    printf("\n End ExtrapolateOneAirVelocityIntoObject() ..... \n " );
}

#else
void VectorField3D::ExtrapolateOneAirVelocityIntoObject(const Voxel &voxel,
		float *phi_tmp, float *vel, int vel_index, float delta) {

	MinHeap<TentativeGridPoint, float> heap;
	printf("\n Start ExtrapolateOneAirVelocityIntoObject() ..... \n " );
	FOR_EACH_CELL
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			if(voxel.InSolid(i,j,k) && voxel.CloseToAir(tp)){
				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
			}
//		if( voxel.InSolid(i,j,k) ){
//			KnownPoint tmp(i, j, k);
//			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
//			heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
//		}
//		else{
//			if(vel_index == 1 && voxel.InSolid(i+1,j,k) ||
//			   vel_index == 2 && voxel.InSolid(i,j+1,k)	||
//			   vel_index == 3 && voxel.InSolid(i,j,k+1)	){
//				KnownPoint tmp(i, j, k);
//				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
//				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
//			}
//
//		}
	 END_FOR

//	printf("Extrapolating velocity: b = %d with %ld in heap \n ",
//			    vel_index, heap.size());

	UpdateCells(heap, vel_index, vel, phi_tmp);
	 printf("\n End ExtrapolateOneAirVelocityIntoObject() ..... \n " );
}
#endif

// extrapolate (moving) solid velocities into liquid region
//
// only cells close (with a phi < EXTRAPOLATE_SOLID_VELOCITY_LIMIT)
// to the (moving) solid are updated
//
// these extrapolated solid velocities are used in Projection to
// compute the r.h.s.
//
// the extrapolation is done based upon flip-signed phi_u_obj, phi_v_obj
// and phi_w_obj


#define EXTRAPOLATE_SOLID_VELOCITY_LIMIT 2*delta

#ifdef SPMD
/*void VectorField3D::ExtrapolateOneVelocityIntoLiquid(const Voxel &voxel,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoLiquid() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(EXTRAPOLATE_SOLID_VELOCITY_LIMIT / delta * 2));
	SetZero(needRecomputePhi);
	FOR_EACH_CELL
		TentativeGridPoint tg(0.f, i, j, k);
		int m;
		if( phi_tmp[INDEX(i,j,k)] > 0.f ){
			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tg, m)){
				if( (vel_index == 1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) ||
					(vel_index == 2 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) ||
					(vel_index == 3 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) ){
					needRecomputePhi[INDEX(i,j,k)] = 1;
				}
			}
			if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) ){
				if( (vel_index == 1 && voxel.InSolid(i+1,j,k)) ||
					(vel_index == 2 && voxel.InSolid(i,j+1,k)) ||
					(vel_index == 3 && voxel.InSolid(i,j,k+1)) ){
					needRecomputePhi[INDEX(i,j,k)] = 1;
				}
			}
		}
	END_FOR
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
		FOR_EACH_CELL
			if(needRecomputePhi[INDEX(i,j,k)]){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
							vel, vel0, vel_index);
			}
		END_FOR
	}
	printf("\n End ExtrapolateOneVelocityIntoLiquid() \n " );
}*/

void VectorField3D::ExtrapolateOneVelocityIntoLiquid(const Voxel &voxel,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoLiquid() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(EXTRAPOLATE_SOLID_VELOCITY_LIMIT / delta * 2));
	SetZero(needRecomputePhi);
	FOR_EACH_CELL
		TentativeGridPoint tg(0.f, i, j, k);
		int m;
		if( phi_tmp[INDEX(i,j,k)] > 0.f && phi_tmp[INDEX(i,j,k)] < 3.f*delta){
			needRecomputePhi[INDEX(i,j,k)] = 1;
		}
	END_FOR
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
		FOR_EACH_CELL
			if(needRecomputePhi[INDEX(i,j,k)]){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
							vel, vel0, vel_index);
			}
		END_FOR
	}
	printf("\n End ExtrapolateOneVelocityIntoLiquid() \n " );
}
#else
void VectorField3D::ExtrapolateOneVelocityIntoLiquid(const Voxel &voxel,
		float *phi_tmp, float *vel, int vel_index, float delta) {

	MinHeap<TentativeGridPoint, float> heap;
	printf("\n Start ExtrapolateOneVelocityIntoLiquid() ..... \n " );
	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
		int n;
		if( phi_tmp[INDEX(i,j,k)] > 0.f){
			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tp, n)){
				if( (vel_index == 1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) ||
					(vel_index == 2 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) ||
					(vel_index == 3 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
				}
			}
			if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) ){
				if( (vel_index == 1 && voxel.InSolid(i+1,j,k)) ||
					(vel_index == 2 && voxel.InSolid(i,j+1,k)) ||
					(vel_index == 3 && voxel.InSolid(i,j,k+1)) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
				}
			}
		}
	END_FOR
	if(heap.size() == 0)
		return;
//	printf("Extrapolating velocity: b = %d with %ld in heap \n ",
//			    vel_index, heap.size());
	UpdateCells(heap, vel_index, vel, phi_tmp);
	printf("\n End ExtrapolateOneVelocityIntoLiquid() \n " );

}
#endif

#ifdef SPMD
void VectorField3D::ExtrapolateOneVelocityIntoNonSolid(const Voxel &voxel,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoLiquid() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(EXTRAPOLATE_SOLID_VELOCITY_LIMIT / delta * 2));
	SetZero(needRecomputePhi);
	FOR_EACH_CELL
		TentativeGridPoint tg(0.f, i, j, k);
		int m;
		if( phi_tmp[INDEX(i,j,k)] > 0.f){
			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tg, m)){
				if( (vel_index == 1 && !voxel.InSolid(i+1,j,k)) ||
					(vel_index == 2 && !voxel.InSolid(i,j+1,k)) ||
					(vel_index == 3 && !voxel.InSolid(i,j,k+1)) ){
					needRecomputePhi[INDEX(i,j,k)] = 1;
				}
			}
			if( !voxel.InSolid(i,j,k) ){
				if( (vel_index == 1 && voxel.InSolid(i+1,j,k)) ||
					(vel_index == 2 && voxel.InSolid(i,j+1,k)) ||
					(vel_index == 3 && voxel.InSolid(i,j,k+1)) ){
					needRecomputePhi[INDEX(i,j,k)] = 1;
				}
			}
		}
	END_FOR
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
		FOR_EACH_CELL
			if(needRecomputePhi[INDEX(i,j,k)]){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
							vel, vel0, vel_index);
			}
		END_FOR
	}
	printf("\n End ExtrapolateOneVelocityIntoLiquid() \n " );
}
#else
void VectorField3D::ExtrapolateOneVelocityIntoNonSolid(const Voxel &voxel,
		float *phi_tmp, float *vel, int vel_index, float delta) {

	MinHeap<TentativeGridPoint, float> heap;
	printf("\n Start ExtrapolateOneVelocityIntoLiquid() ..... \n " );
	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
		int n;
		if( phi_tmp[INDEX(i,j,k)] > 0.f){
			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tp, n)){
				if( (vel_index == 1 && !voxel.InSolid(i+1,j,k)) ||
					(vel_index == 2 && !voxel.InSolid(i,j+1,k)) ||
					(vel_index == 3 && !voxel.InSolid(i,j,k+1)) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
				}
			}
			if(!voxel.InSolid(i,j,k)){
				if( (vel_index == 1 && voxel.InSolid(i+1,j,k)) ||
					(vel_index == 2 && voxel.InSolid(i,j+1,k)) ||
					(vel_index == 3 && voxel.InSolid(i,j,k+1)) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
				}
			}
		}
	END_FOR
	if(heap.size() == 0)
		return;
//	printf("Extrapolating velocity: b = %d with %ld in heap \n ",
//			    vel_index, heap.size());
	UpdateCells(heap, vel_index, vel, phi_tmp);
	printf("\n End ExtrapolateOneVelocityIntoLiquid() \n " );

}
#endif

/*
#ifdef SPMD
void VectorField3D::ExtrapolateOneVelocityIntoObjectForAdvection(const Voxel &voxel,
		map<u_int, KnownPoint> &knownPoint,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	map<u_int, KnownPoint>::iterator foundpos;
	printf("\n Start ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(LS_ADV_LIMIT / delta * 2));
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
#ifdef TBB
		ExtrapolateObjectBody body(&knownPoint, phi_tmp, vel, vel0,
					vel_index, delta, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
						  body);
#else
		FOR_EACH_CELL
			foundpos = knownPoint.find(INDEX(i,j,k));
			if( foundpos == knownPoint.end() ){
				if( phi_tmp[INDEX(i,j,k)] >= 0.f && phi_tmp[INDEX(i,j,k)] < LS_ADV_LIMIT ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
								vel, vel0, vel_index);
				}
			}
		END_FOR
#endif
	}
    printf("\n End ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );

}

#else

// extrapolates the fluid(air) velocities into (moving) solid cells
//
// parameter knownPoint contains the cells which should not be updated
// because these cells have been correctly set in during projection step

void VectorField3D::ExtrapolateOneVelocityIntoObjectForAdvection(const Voxel &voxel,
		map<u_int, KnownPoint> &knownPoint,
		float *phi_tmp, float *vel, int vel_index, float delta) {

	MinHeap<TentativeGridPoint, float> heap;
	map<u_int, KnownPoint>::iterator foundpos;
	printf("\n Start ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );
	FOR_EACH_CELL
		foundpos = knownPoint.find(INDEX(i,j,k));
		if( foundpos == knownPoint.end() ){
			if( phi_tmp[INDEX(i,j,k)] >= 0.f ){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
			}
		}
	 END_FOR

//	printf("Extrapolating velocity: b = %d with %ld in heap \n ",
//			    vel_index, heap.size());

	 UpdateCells(heap, vel_index, vel, phi_tmp);
	 printf("\n End ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );

}
#endif
*/


#ifdef SPMD
/*void VectorField3D::ExtrapolateOneVelocityIntoObjectForAdvection(const Voxel &voxel,
		//map<u_int, KnownPoint> &knownPoint,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(LS_ADV_LIMIT / delta * 2));
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
#ifdef TBB
		ExtrapolateObjectBody1 body(&voxel, phi_tmp, vel, vel0, phi_c,
					vel_index, delta, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
						  body);
#else
	FOR_EACH_CELL
		if( phi_tmp[INDEX(i,j,k)] >= 0.f ){
			if(phi_tmp[INDEX(i,j,k)] < LS_ADV_LIMIT ){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
									vel, vel0, vel_index);
			}
		}
		else{
			if(vel_index == 1){
				if( (!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i+1,j,k)) ||
					(voxel.InSolid(i,j,k) && i < DimX-1 && !voxel.InSolid(i+1,j,k) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
								vel, vel0, vel_index);
				}
			}
			if(vel_index == 2){
				if( (!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i,j+1,k)) ||
					(voxel.InSolid(i,j,k) && j < DimY-1 && !voxel.InSolid(i,j+1,k) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
								vel, vel0, vel_index);
				}
			}
			if(vel_index == 3){
				if( (!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i,j,k+1)) ||
					(voxel.InSolid(i,j,k) && k < DimZ-1 && !voxel.InSolid(i,j,k+1) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
							vel, vel0, vel_index);
				}
			}
		}
	 END_FOR
#endif
	}
	 printf("\n End ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );

}*/

void VectorField3D::ExtrapolateOneVelocityIntoObjectForAdvectionHouston2003(const Voxel &voxel,
		//map<u_int, KnownPoint> &knownPoint,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(LS_ADV_LIMIT / delta * 2));
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
#ifdef TBB
		ExtrapolateObjectBody1 body(&voxel, phi_tmp, vel, vel0, phi_c,
					vel_index, delta, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
						  body);
#else
	FOR_EACH_CELL
		if( phi_tmp[INDEX(i,j,k)]+E_EPSIL >= 0.f ){
			if(phi_tmp[INDEX(i,j,k)] < LS_ADV_LIMIT ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
										vel, vel0, vel_index);
			}
		}
		else{
			if( vel_index == 1 && i < DimX-1 ){
				if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i+1,j,k)) ||
				     (voxel.InSolid(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
									vel, vel0, vel_index);
				}
			}
			if( vel_index == 2 && j < DimY-1 ){
				if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j+1,k)) ||
				    (voxel.InSolid(i,j,k)  && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
									vel, vel0, vel_index);
				}
			}
			if( vel_index == 3 && k < DimZ-1 ){
				if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j,k+1)) ||
				     (voxel.InSolid(i,j,k)  && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
									vel, vel0, vel_index);
				}
			}
		}


	 END_FOR
#endif
	}
	 printf("\n End ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );

}

void VectorField3D::ExtrapolateOneVelocityIntoObjectForAdvection(const Voxel &voxel,
		//map<u_int, KnownPoint> &knownPoint,
		float *phi_tmp, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );
	set_time_base();
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(LS_ADV_LIMIT / delta * 2));
	SetEqual(vel, vel0);

	char *needed = new char[DimX*DimY*DimZ];
	SetZero(needed);

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if( phi_tmp[pos]+E_EPSIL >= 0.f ){
			if(phi_tmp[pos] < LS_ADV_LIMIT ){
//			if(voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) || voxel.InSurface(i,j,k) ){
//				if( (vel_index == 1 && i < DimX-1 && voxel.InSolid(i+1,j,k)) ||
//					(vel_index == 2 && j < DimY-1 && voxel.InSolid(i,j+1,k)) ||
//					(vel_index == 3 && k < DimZ-1 && voxel.InSolid(i,j,k+1)) ){
				if(!voxel.InMovingSolid(i,j,k))
					needed[pos] = 1;
			}
//				}
//				if( (vel_index == 1 && i == DimX-1) ||
//				    (vel_index == 2 && j == DimY-1) ||
//					(vel_index == 3 && k == DimZ-1) ){
//					KnownPoint tmp(i, j, k);
//					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
//					MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
//									vel, vel0, vel_index);
//				}
//			}
		}
		else{
			if( vel_index == 1 && i < DimX-1 ){
				if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1,j,k)) ||
				    (voxel.InSolid(i,j,k) && !voxel.InMovingSolid(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)))
					needed[pos] = 1;
			}
			if( vel_index == 2 && j < DimY-1 ){
				if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i,j+1,k)) ||
					(voxel.InSolid(i,j,k) && !voxel.InMovingSolid(i,j,k) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)))
					needed[pos] = 1;
			}
			if( vel_index == 3 && k < DimZ-1 ){
				if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i,j,k+1)) ||
					(voxel.InSolid(i,j,k) && !voxel.InMovingSolid(i,j,k)  && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)))
					needed[pos] = 1;
			}
		}
	 END_FOR

	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
#ifdef TBB
		ExtrapolateObjectBody1 body(&voxel, phi_tmp, vel, vel0, phi_c,
					vel_index, delta, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
						  body);
#else


	FOR_EACH_CELL
		if( needed[INDEX(i,j,k)] == 1 ){
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
 					    vel, vel0, vel_index);
		}
	END_FOR
#endif
	}

	 delete [] needed;
	 printf("End ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );
	 printf("Time spent in ExtrapolateOneVelocityIntoObjectForAdvection() = %f \n\n", get_time_in_seconds());

}

void VectorField3D::ExtrapolateOneVelocityIntoObjectForAdvection1(const Voxel &voxel,
		//map<u_int, KnownPoint> &knownPoint,
		float *phi_tmp, const char *valid, float *vel, float *vel0, int vel_index, float delta) {

	printf("\n Start ExtrapolateOneVelocityIntoObjectForAdvection1() ..... \n " );
	float dtau = delta / 2;
	float inv_delta2 = 1.f / (delta * delta);
	int iterations = int(round(LS_ADV_LIMIT / delta * 2));
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
#ifdef TBB
		ExtrapolateObjectBody1 body(&voxel, phi_tmp, vel, vel0, phi_c,
					vel_index, delta, DimX, DimY, DimZ);
		tbb::parallel_for(tbb::blocked_range<int>(0, DimZ, DimZ/CHUNK_SIZE),
						  body);
#else
	FOR_EACH_CELL
		if( valid[INDEX(i,j,k)] && phi_tmp[INDEX(i,j,k)] < LS_ADV_LIMIT ){
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp,
							vel, vel0, vel_index);
		}
	 END_FOR
#endif
	}
	 printf("\n End ExtrapolateOneVelocityIntoObjectForAdvection1() ..... \n " );
}
#else
  void VectorField3D::ExtrapolateOneVelocityIntoObjectForAdvection(const Voxel &voxel,
		//map<u_int, KnownPoint> &knownPoint,
		float *phi_tmp, float *vel, int vel_index, float delta) {

	MinHeap<TentativeGridPoint, float> heap;
	map<u_int, KnownPoint>::iterator foundpos;
	printf("\n Start ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );
	FOR_EACH_CELL
			if( phi_tmp[INDEX(i,j,k)] >= 0.f && phi_tmp[INDEX(i,j,k)] < LS_ADV_LIMIT ){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
			}
			else{
				if(vel_index == 1){
					if( (!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)] < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i+1,j,k)) ||
						(voxel.InSolid(i,j,k) && !voxel.InSolid(i+1,j,k) && phi_c[INDEX(i+1,j,k)] < EXTRAPOLATE_VEL_LIMIT) ){
						KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
					}
				}
				if(vel_index == 2){
					if( (!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)] < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i,j+1,k)) ||
						(voxel.InSolid(i,j,k) && !voxel.InSolid(i,j+1,k) && phi_c[INDEX(i,j+1,k)] < EXTRAPOLATE_VEL_LIMIT) ){
						KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
					}
				}
				if(vel_index == 3){
					if( (!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)] < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i,j,k+1)) ||
						(voxel.InSolid(i,j,k) && !voxel.InSolid(i,j,k+1) && phi_c[INDEX(i,j,k+1)] < EXTRAPOLATE_VEL_LIMIT) ){
						KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
					}
				}
			}
	 END_FOR

//	printf("Extrapolating velocity: b = %d with %ld in heap \n ",
//			    vel_index, heap.size());

	 UpdateCells(heap, vel_index, vel, phi_tmp);
	 printf("\n End ExtrapolateOneVelocityIntoObjectForAdvection() ..... \n " );

}
#endif





#ifdef SPMD

static void VisitNeighbors(const TentativeGridPoint &minTentative,
					float *phi_tmp, float *vel0, char *valid,
					int DimX, int DimY, int DimZ,
		            float &phi_minus_x, float &phi_plus_x, float &phi_minus_y,
		            float &phi_plus_y, float &phi_minus_z, float &phi_plus_z,
		            float &vel_minus_x, float &vel_plus_x, float &vel_minus_y,
		            float &vel_plus_y, float &vel_minus_z, float &vel_plus_z){

	vector<TentativeGridPoint> neighbors;
	minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
	u_int pos;
	if(neighbors.size() == 6){
		pos = INDEX(minTentative.ii-1, minTentative.jj, minTentative.kk);
		if(valid[pos])
			phi_minus_x = phi_tmp[pos];
		vel_minus_x = vel0[pos];
		pos = INDEX(minTentative.ii+1, minTentative.jj, minTentative.kk);
		if(valid[pos])
			phi_plus_x = phi_tmp[pos];
		vel_plus_x  = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj-1, minTentative.kk);
		if(valid[pos])
			phi_minus_y = phi_tmp[pos];
		vel_minus_y = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj+1, minTentative.kk);
		if(valid[pos])
			phi_plus_y = phi_tmp[pos];
		vel_plus_y  = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj, minTentative.kk-1);
		if(valid[pos])
			phi_minus_z = phi_tmp[pos];
		vel_minus_z = vel0[pos];
		pos = INDEX(minTentative.ii, minTentative.jj, minTentative.kk+1);
		if(valid[pos])
			phi_plus_z = phi_tmp[pos];
		vel_plus_z  = vel0[pos];
	}
	else{
		for(int m=0; m<neighbors.size(); m++){
			u_int pos = INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
			if(minTentative.LeftNeighbor(neighbors[m])){
				phi_minus_x = phi_tmp[pos];
				vel_minus_x = vel0[pos];
				if(!valid[pos])
					 phi_minus_x = INFINITY;
				continue;
			}
			if(minTentative.RightNeighbor(neighbors[m])){
				phi_plus_x = phi_tmp[pos];
				vel_plus_x = vel0[pos];
				if(!valid[pos])
					 phi_plus_x = INFINITY;
				continue;
			}
			if(minTentative.BackNeighbor(neighbors[m])){
				phi_minus_y = phi_tmp[pos];
				vel_minus_y = vel0[pos];
				if(!valid[pos])
					phi_minus_y = INFINITY;
				continue;
			}
			if(minTentative.FrontNeighbor(neighbors[m])){
				phi_plus_y = phi_tmp[pos];
				vel_plus_y = vel0[pos];
				if(!valid[pos])
					 phi_plus_y = INFINITY;
				continue;
			}
			if(minTentative.BottomNeighbor(neighbors[m])){
				phi_minus_z = phi_tmp[pos];
				vel_minus_z = vel0[pos];
				if(!valid[pos])
					  phi_minus_z = INFINITY;
				continue;
			}
			if(minTentative.TopNeighbor(neighbors[m])){
				phi_plus_z = phi_tmp[pos];
				vel_plus_z = vel0[pos];
				if(!valid[pos])
					 phi_plus_z = INFINITY;
				continue;
			}
		}
	}
}


static void MarchingOut(const TentativeGridPoint &minTentative,
						int DimX, int DimY, int DimZ,
						float dtau, float inv_delta2, float *phi_tmp, char *valid,
						float *vel, float *vel0, int vel_index) {

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
	VisitNeighbors(minTentative, phi_tmp, vel0, valid,
			       DimX, DimY, DimZ,
			      phi_minus_x, phi_plus_x, phi_minus_y,
			      phi_plus_y, phi_minus_z, phi_plus_z,
				  vel_minus_x, vel_plus_x, vel_minus_y,
				  vel_plus_y,  vel_minus_z, vel_plus_z);


	ProcessCoefficents(phi_minus_x, phi_plus_x, minValue,
					vel_minus_x, vel_plus_x, phi_x, A, a);
	ProcessCoefficents(phi_minus_y, phi_plus_y, minValue,
					vel_minus_y, vel_plus_y, phi_y, B, b);
	ProcessCoefficents(phi_minus_z, phi_plus_z, minValue,
					vel_minus_z, vel_plus_z, phi_z, C, c);

	if(phi_x == 0.f && phi_y == 0.f && phi_z == 0.f)
		valid[pos] = 0;
	else{
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
		float mag = sqrtf(inv_delta2 * (phi_x * phi_x + phi_y * phi_y  + phi_z * phi_z));
		vel[pos] = vel0[pos] - dtau * inv_delta2 / mag * (phi_x * vel_x +  phi_y * vel_y + phi_z * vel_z);
	}

}


static bool DetectNonValid(const TentativeGridPoint &minTentative,
						int DimX, int DimY, int DimZ,
						float *phi_tmp, char *valid, float *vel0, int vel_index) {

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

	VisitNeighbors(minTentative, phi_tmp, vel0, valid,
					DimX, DimY, DimZ,
				   phi_minus_x, phi_plus_x, phi_minus_y,
				   phi_plus_y, phi_minus_z, phi_plus_z,
				   vel_minus_x, vel_plus_x, vel_minus_y,
				   vel_plus_y,  vel_minus_z, vel_plus_z);


	ProcessCoefficents(phi_minus_x, phi_plus_x, minValue,
					vel_minus_x, vel_plus_x, phi_x, A, a);
	ProcessCoefficents(phi_minus_y, phi_plus_y, minValue,
					vel_minus_y, vel_plus_y, phi_y, B, b);
	ProcessCoefficents(phi_minus_z, phi_plus_z, minValue,
					vel_minus_z, vel_plus_z, phi_z, C, c);

	if(phi_x == 0.f && phi_y == 0.f && phi_z == 0.f){
		valid[INDEX(i,j,k)] = 0;
//	 if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K){
//		printf(" ii = %d, jj = %d, kk = %d vel_index = %d\n",
//				minTentative.ii, minTentative.jj, minTentative.kk, vel_index);
//		printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
//				minValue, phi_minus_x, phi_plus_x);
//		printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
//					minValue, phi_minus_y, phi_plus_y);
//		printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
//					minValue, phi_minus_z, phi_plus_z);
//		printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
//				vel0[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//							vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
//							vel_minus_z, vel_plus_z);
//		printf(" valid = %d, valid_minus_x = %d, valid_plus_x = %d \n ",
//			    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
//			    	valid[INDEX(minTentative.ii-1,minTentative.jj,minTentative.kk)],
//			    	valid[INDEX(minTentative.ii+1,minTentative.jj,minTentative.kk)]);
//		printf(" valid = %d, valid_minus_y = %d, valid_plus_y = %d\n ",
//			    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
//			    	valid[INDEX(minTentative.ii,minTentative.jj-1,minTentative.kk)],
//			    	valid[INDEX(minTentative.ii,minTentative.jj+1,minTentative.kk)]);
//		printf(" valid = %d, valid_minus_z = %d, valid_plus_z = %d\n ",
//			    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
//			    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk-1)],
//			    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk+1)]);
//		 }
		return true;
	}
	else
		return false;

}
#endif

#ifdef SPMD
void VectorField3D::ExtrapolateOneForProjection(const Voxel &voxel,
		                         map<u_int, KnownPoint> &knownPoint,
		                        float *phi_tmp, float *phiobj, float *vel, float *vel0, int vel_index) {

	float delta = voxel.VoxelDelta();
	float inv_delta2 = 1.f / (delta * delta);

	map<u_int, KnownPoint>::iterator pos;
	map<u_int, KnownPoint>::iterator foundpos;

	char *valid = new char[DimX*DimY*DimZ];
	char *needed = new char[DimX*DimY*DimZ];
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	SetZero(needed);
	SetValidPoint(voxel, valid, vel_index);

//	u_int N2 = 0;
//	FOR_EACH_CELL
//		if( (vel_index == 1 && phi_u_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 2 && phi_v_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 3 && phi_w_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ){
//			++N2;
////			printf(" ii = %d, jj = %d, kk = %d  phic = %f\n",
////					i, j, k, phi_c[INDEX(i,j,k)]);
//
//		}
//	END_FOR
//
//	printf("before DetectNonValid there are on %u non-updated points \n ", N2);



	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
		int n;
		foundpos = knownPoint.find(INDEX(i,j,k));
		if( foundpos == knownPoint.end() ){
			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tp, n) && valid[INDEX(i,j,k)] && phi_tmp[INDEX(i,j,k)] > 0.f){
				if( (vel_index == 1 && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 2 && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 3 && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
					needed[INDEX(i,j,k)] = 1;
				}
			}
			if( (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) &&
			     phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT &&
				 phi_tmp[INDEX(i,j,k)] > 0.f && valid[INDEX(i,j,k)]
			  ) {
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
				needed[INDEX(i,j,k)] = 1;
//				if(i == I && j == J && k == K)
//					printf(" wrong! (%d, %d, %d) shouldn't included with phi_c = %f \n",
//							i,j,k, phi_c[INDEX(i,j,k)]);
			}
			if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= EXTRAPOLATE_VEL_LIMIT &&
					valid[INDEX(i,j,k)]){
				TentativeGridPoint it(0.f, i, j, k);
				vector<TentativeGridPoint> neighbors;
				it.Neigbor(neighbors, DimX, DimY, DimZ);
				for(int m=0; m<neighbors.size(); m++){
					if(it.RightNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
					  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
						if(vel_index == 1){
							KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
							needed[INDEX(i,j,k)] = 1;
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
						}
					}
					if(it.FrontNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
					  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
						if(vel_index == 2){
							KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
							needed[INDEX(i,j,k)] = 1;
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
						}
					}
					if(it.TopNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
					  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
						if(vel_index == 3){
							KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
							needed[INDEX(i,j,k)] = 1;
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
						}
					}
				}

			}
		}
	END_FOR

//	u_int N1 = 0;
//	FOR_EACH_CELL
//		if( (vel_index == 1 && phi_u_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 2 && phi_v_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 3 && phi_w_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ){
//			++N1;
//		}
//	END_FOR
//
//	printf("before performing iterations there are on %u non-updated points \n ", N1);

	u_int N = 0;
	FOR_EACH_CELL
		if(needed[INDEX(i,j,k)])
			++N;
	END_FOR
	float dtau = delta / 2;

	int iterations = int(round(EXTRAPOLATE_VEL_LIMIT / delta * 2));
	printf("performing %d iterations on %u points \n ", iterations, N);
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
		FOR_EACH_CELL
			if(needed[INDEX(i,j,k)]){
				KnownPoint tmp(i,j,k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp, valid, vel, vel0, vel_index);
//				if(i == I && j == J && k == K){
//					printf(" ietr = %d, at (%d, %d, %d) vel = %f, with phi = %f, vel_index = %d \n",
//						n, i, j, k, vel[INDEX(i,j,k)], phi_tmp[INDEX(i,j,k)], vel_index);
//					printf(" phi_tmp[i+1] = %f, vel0[i+1] = %f \n", phi_tmp[INDEX(i+1,j,k)], vel0[INDEX(i+1,j,k)] );
//					printf(" phi_tmp[i-1] = %f, vel0[i-1] = %f \n", phi_tmp[INDEX(i-1,j,k)], vel0[INDEX(i-1,j,k)] );
//					printf(" phi_tmp[j+1] = %f, vel0[j+1] = %f \n", phi_tmp[INDEX(i,j+1,k)], vel0[INDEX(i,j+1,k)] );
//					printf(" phi_tmp[j-1] = %f, vel0[j-1] = %f \n", phi_tmp[INDEX(i,j-1,k)], vel0[INDEX(i,j-1,k)] );
//					printf(" phi_tmp[k+1] = %f, vel0[k+1] = %f \n", phi_tmp[INDEX(i,j,k+1)], vel0[INDEX(i,j,k+1)] );
//					printf(" phi_tmp[k-1] = %f, vel0[k-1] = %f \n", phi_tmp[INDEX(i,j,k-1)], vel0[INDEX(i,j,k-1)] );
//				}
			}
		END_FOR
	}
//	FOR_EACH_CELL
//		if(needed[INDEX(i,j,k)] && phi_tmp[INDEX(i,j,k)] < delta){
//			printf(" at (%d, %d, %d) vel = %f, with phi = %f, vel_index = %d \n",
//					i,j,k,vel[INDEX(i,j,k)], phi_tmp[INDEX(i,j,k)], vel_index);
//		}
//	END_FOR
#ifdef WIN32
	hash_map<u_int, TentativeGridPoint> heap;
	hash_map<u_int, TentativeGridPoint>::iterator heap_iter;
#else
	unordered_map<u_int, TentativeGridPoint> heap;
	unordered_map<u_int, TentativeGridPoint>::iterator heap_iter;
#endif
	int iters = 0;
	while(iters < 10 && HasNonUpdatedPoints(voxel, valid, phi_tmp, vel_index, heap) > 0){
		printf("A. at iter %d, there are %u points not updated \n ", iters, heap.size());
		++iters;
		for(heap_iter=heap.begin(); heap_iter!=heap.end(); ++heap_iter){
			TentativeGridPoint &tp = heap_iter->second;
			int i = tp.ii;
			int j = tp.jj;
			int k = tp.kk;
			u_int index = heap_iter->first;
			vector<TentativeGridPoint> neighbors;
		    tp.Neigbor(neighbors, DimX, DimY, DimZ);
			float vel_avg = 0.f;
			int NV = 0;
			for(int m=0; m<neighbors.size(); m++){
				if(valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]){
					++NV;
					vel_avg += vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
				}
			}
			if(NV == 0){

			}
			else{
				vel[index] = vel_avg / NV;
				valid[index] = 1;
			}
//			printf(" NV = %d, at (%d, %d, %d) vel = %f, with phi = %f, vel_index = %d \n",
//				 NV, i, j, k, vel[index], phi_tmp[index], vel_index);
//			printf(" phi_tmp[i+1] = %f, phiobj[i+1] = %f, vel[i+1] = %f valid[i+1] = %d\n",
//					phi_tmp[INDEX(i+1,j,k)], phiobj[INDEX(i+1,j,k)], vel[INDEX(i+1,j,k)], valid[INDEX(i+1,j,k)]);
//			printf(" phi_tmp[i-1] = %f, phiobj[i-1] = %f, vel[i-1] = %f valid[i-1] = %d\n",
//					phi_tmp[INDEX(i-1,j,k)], phiobj[INDEX(i-1,j,k)], vel[INDEX(i-1,j,k)], valid[INDEX(i-1,j,k)]);
//			printf(" phi_tmp[j+1] = %f, phiobj[j+1] = %f, vel[j+1] = %f valid[j+1] = %d\n",
//					phi_tmp[INDEX(i,j+1,k)], phiobj[INDEX(i,j+1,k)], vel[INDEX(i,j+1,k)], valid[INDEX(i,j+1,k)]);
//			printf(" phi_tmp[j-1] = %f, phiobj[j-1] = %f, vel[j-1] = %f valid[j-1] = %d\n",
//					phi_tmp[INDEX(i,j-1,k)], phiobj[INDEX(i,j-1,k)], vel[INDEX(i,j-1,k)], valid[INDEX(i,j-1,k)]);
//			printf(" phi_tmp[k+1] = %f, phiobj[k+1] = %f, vel[k+1] = %f valid[k+1] = %d\n",
//					phi_tmp[INDEX(i,j,k+1)], phiobj[INDEX(i,j,k+1)], vel[INDEX(i,j,k+1)], valid[INDEX(i,j,k+1)]);
//			printf(" phi_tmp[k-1] = %f, phiobj[k-1] = %f, vel[k-1] = %f valid[k-1] = %d\n",
//					phi_tmp[INDEX(i,j,k-1)], phiobj[INDEX(i,j,k-1)], vel[INDEX(i,j,k-1)], valid[INDEX(i,j,k-1)] );
		}
		heap.clear();
	}

	if(valid) delete [] valid;
	if(needed) delete [] needed;

	knownPoint.clear();

}

void VectorField3D::ExtrapolateOneForProjection(const Voxel &voxel, float *phi_tmp,
						float *vel, float *vel0, int vel_index) {

	float delta = voxel.VoxelDelta();
	float inv_delta2 = 1.f / (delta * delta);

	map<u_int, KnownPoint>::iterator pos;
	map<u_int, KnownPoint>::iterator foundpos;

	char *valid = new char[DimX*DimY*DimZ];
	char *needed = new char[DimX*DimY*DimZ];
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	SetZero(needed);
	SetValidPoint(voxel, valid, 0);

//	u_int N2 = 0;
//	FOR_EACH_CELL
//		if( (vel_index == 1 && phi_u_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 2 && phi_v_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 3 && phi_w_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ){
//			++N2;
////			printf(" ii = %d, jj = %d, kk = %d  phic = %f\n",
////					i, j, k, phi_c[INDEX(i,j,k)]);
//
//		}
//	END_FOR
//
//	printf("before DetectNonValid there are on %u non-updated points \n ", N2);



	FOR_EACH_CELL

		if( (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && phi_tmp[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ){
			KnownPoint tmp(i, j, k);
			TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
			DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
			needed[INDEX(i,j,k)] = 1;
		}

		if(voxel.InAir(i,j,k) && phi_tmp[INDEX(i,j,k)]+E_EPSIL >= EXTRAPOLATE_VEL_LIMIT){
		    	TentativeGridPoint it(0.f, i, j, k);
		    	vector<TentativeGridPoint> neighbors;
		        it.Neigbor(neighbors, DimX, DimY, DimZ);
		        for(int m=0; m<neighbors.size(); m++){
		        	if(it.RightNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
		        	  && phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        		KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
						needed[INDEX(i,j,k)] = 1;
						if(phi_tmp[INDEX(i,j,k)]< 0.f){
							printf("Error A ! shouldn't extrapolate this one\n");
							exit(1);
						}

	        		}
		        	if(it.FrontNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
		        	  && phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        		KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
						needed[INDEX(i,j,k)] = 1;
						if(phi_tmp[INDEX(i,j,k)]< 0.f){
							printf("Error B ! shouldn't extrapolate this one\n");
							exit(1);
						}
	        		}
		        	if(it.TopNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
		        	  && phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        		KnownPoint tmp(i, j, k);
						TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
						DetectNonValid(tp, DimX, DimY, DimZ, phi_tmp, valid, vel0, vel_index);
						needed[INDEX(i,j,k)] = 1;
						if(phi_tmp[INDEX(i,j,k)]< 0.f){
							printf("Error C ! shouldn't extrapolate this one\n");
							exit(1);
						}
	        		}
		        }
		}

	END_FOR

//	u_int N1 = 0;
//	FOR_EACH_CELL
//		if( (vel_index == 1 && phi_u_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 2 && phi_v_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
//			(vel_index == 3 && phi_w_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ){
//			++N1;
//		}
//	END_FOR
//
//	printf("before performing iterations there are on %u non-updated points \n ", N1);

	u_int N = 0;
	FOR_EACH_CELL
		if(needed[INDEX(i,j,k)])
			++N;
	END_FOR
	float dtau = delta / 2;

	int iterations = int(round(EXTRAPOLATE_VEL_LIMIT / delta * 2));
	printf("performing %d iterations on %u points \n ", iterations, N);
	SetEqual(vel, vel0);
	for (int n = 0; n < iterations; ++n){
		SWAP(vel, vel0);
		FOR_EACH_CELL
			if(needed[INDEX(i,j,k)]){
				KnownPoint tmp(i,j,k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				MarchingOut(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_tmp, valid, vel, vel0, vel_index);
//				if(i == I && j == J && k == K){
//					printf(" ietr = %d, at (%d, %d, %d) vel = %f, with phi = %f, vel_index = %d \n",
//						n, i, j, k, vel[INDEX(i,j,k)], phi_tmp[INDEX(i,j,k)], vel_index);
//					printf(" phi_tmp[i+1] = %f, vel0[i+1] = %f \n", phi_tmp[INDEX(i+1,j,k)], vel0[INDEX(i+1,j,k)] );
//					printf(" phi_tmp[i-1] = %f, vel0[i-1] = %f \n", phi_tmp[INDEX(i-1,j,k)], vel0[INDEX(i-1,j,k)] );
//					printf(" phi_tmp[j+1] = %f, vel0[j+1] = %f \n", phi_tmp[INDEX(i,j+1,k)], vel0[INDEX(i,j+1,k)] );
//					printf(" phi_tmp[j-1] = %f, vel0[j-1] = %f \n", phi_tmp[INDEX(i,j-1,k)], vel0[INDEX(i,j-1,k)] );
//					printf(" phi_tmp[k+1] = %f, vel0[k+1] = %f \n", phi_tmp[INDEX(i,j,k+1)], vel0[INDEX(i,j,k+1)] );
//					printf(" phi_tmp[k-1] = %f, vel0[k-1] = %f \n", phi_tmp[INDEX(i,j,k-1)], vel0[INDEX(i,j,k-1)] );
//				}
			}
		END_FOR
	}
//	FOR_EACH_CELL
//		if(needed[INDEX(i,j,k)] && phi_tmp[INDEX(i,j,k)] < delta){
//			printf(" at (%d, %d, %d) vel = %f, with phi = %f, vel_index = %d \n",
//					i,j,k,vel[INDEX(i,j,k)], phi_tmp[INDEX(i,j,k)], vel_index);
//		}
//	END_FOR

#ifdef WIN32
	hash_map<u_int, TentativeGridPoint> heap;
	hash_map<u_int, TentativeGridPoint>::iterator heap_iter;
#else
	unordered_map<u_int, TentativeGridPoint> heap;
	unordered_map<u_int, TentativeGridPoint>::iterator heap_iter;
#endif
	int iters = 0;
	while(iters < 10 && HasNonUpdatedPoints(voxel, valid, phi_tmp, 0, heap) > 0){
		printf("B. at iter %d, there are %u points not updated \n ", iters, heap.size());
		++iters;
		for(heap_iter=heap.begin(); heap_iter!=heap.end(); ++heap_iter){
			TentativeGridPoint &tp = heap_iter->second;
			int i = tp.ii;
			int j = tp.jj;
			int k = tp.kk;
			u_int index = heap_iter->first;
			vector<TentativeGridPoint> neighbors;
		    tp.Neigbor(neighbors, DimX, DimY, DimZ);
			float vel_avg = 0.f;
			int NV = 0;
			for(int m=0; m<neighbors.size(); m++){
				if(valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]){
					++NV;
					vel_avg += vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
				}
			}
			if(NV == 0){

			}
			else{
				vel[index] = vel_avg / NV;
				valid[index] = 1;
			}
//			printf(" NV = %d, at (%d, %d, %d) vel = %f, with phi = %f, vel_index = %d \n",
//				 NV, i, j, k, vel[index], phi_tmp[index], vel_index);
//			printf(" phi_tmp[i+1] = %f, phiobj[i+1] = %f, vel[i+1] = %f valid[i+1] = %d\n",
//					phi_tmp[INDEX(i+1,j,k)], phiobj[INDEX(i+1,j,k)], vel[INDEX(i+1,j,k)], valid[INDEX(i+1,j,k)]);
//			printf(" phi_tmp[i-1] = %f, phiobj[i-1] = %f, vel[i-1] = %f valid[i-1] = %d\n",
//					phi_tmp[INDEX(i-1,j,k)], phiobj[INDEX(i-1,j,k)], vel[INDEX(i-1,j,k)], valid[INDEX(i-1,j,k)]);
//			printf(" phi_tmp[j+1] = %f, phiobj[j+1] = %f, vel[j+1] = %f valid[j+1] = %d\n",
//					phi_tmp[INDEX(i,j+1,k)], phiobj[INDEX(i,j+1,k)], vel[INDEX(i,j+1,k)], valid[INDEX(i,j+1,k)]);
//			printf(" phi_tmp[j-1] = %f, phiobj[j-1] = %f, vel[j-1] = %f valid[j-1] = %d\n",
//					phi_tmp[INDEX(i,j-1,k)], phiobj[INDEX(i,j-1,k)], vel[INDEX(i,j-1,k)], valid[INDEX(i,j-1,k)]);
//			printf(" phi_tmp[k+1] = %f, phiobj[k+1] = %f, vel[k+1] = %f valid[k+1] = %d\n",
//					phi_tmp[INDEX(i,j,k+1)], phiobj[INDEX(i,j,k+1)], vel[INDEX(i,j,k+1)], valid[INDEX(i,j,k+1)]);
//			printf(" phi_tmp[k-1] = %f, phiobj[k-1] = %f, vel[k-1] = %f valid[k-1] = %d\n",
//					phi_tmp[INDEX(i,j,k-1)], phiobj[INDEX(i,j,k-1)], vel[INDEX(i,j,k-1)], valid[INDEX(i,j,k-1)] );
		}
		heap.clear();
	}

	if(valid) delete [] valid;
	if(needed) delete [] needed;

}

bool VectorField3D::TestVelPointIsLiquid(const Voxel &voxel, int i, int j, int k, int vel_index,
		                                 const char *valid) const{
	if(valid[INDEX(i,j,k)]){
		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
			return true;
		else{
			if( (vel_index == 1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) ||
				(vel_index == 2 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) ||
				(vel_index == 3 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) )
				return true;
			else
				return false;
		}
	}
	else
		return false;
}

void VectorField3D::DetectIsolatedPoints(const Voxel &voxel,
	                         const char *airVelSet, const float *phi_tmp,
	                         float *vel, int vel_index) {

	u_int N = 0;
	float delta = voxel.VoxelDelta();

	char *valid = new char[DimX*DimY*DimZ];
	SetZero(valid);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(vel_index == 1 && phi_u_obj[pos]+E_EPSIL < 0.f)
			valid[pos] = 1;
		if(vel_index == 2 && phi_v_obj[pos]+E_EPSIL < 0.f)
			valid[pos] = 1;
		if(vel_index == 3 && phi_w_obj[pos]+E_EPSIL < 0.f)
			valid[pos] = 1;
	END_FOR

	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
		u_int pos = INDEX(i,j,k);
		int n;
		//foundpos = knownPoint.find(INDEX(i,j,k));
		if( vel_index == 1 && !(airVelSet[pos] & AIRUVELSET) ||
			vel_index == 2 && !(airVelSet[pos] & AIRVVELSET) ||
			vel_index == 3 && !(airVelSet[pos] & AIRWVELSET) ){
			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tp, n) && valid[pos] && phi_tmp[pos] < E_EPSIL){
				if( (vel_index == 1 && (voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k)) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 2 && (voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k)) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 3 && (voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1)) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ){
//					printf("at (%d,%d,%d), vel = %f ", i,j,k,vel[pos]);
					int NV = 0;
					float vtmp = 0.f;
					if(TestVelPointIsLiquid(voxel,i+1,j,k,vel_index,valid)){
						vtmp += vel[INDEX(i+1,j,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i-1,j,k,vel_index,valid)){
						vtmp += vel[INDEX(i-1,j,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j+1,k,vel_index,valid)){
						vtmp += vel[INDEX(i,j+1,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j-1,k,vel_index,valid)){
						vtmp += vel[INDEX(i,j-1,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j,k+1,vel_index,valid)){
						vtmp += vel[INDEX(i,j,k+1)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j,k-1,vel_index,valid)){
						vtmp += vel[INDEX(i,j,k-1)];
						++NV;
					}
					if(NV > 0)
						vel[pos] = vtmp / NV;
//					printf("after averaging, vel = %f \n ", vel[pos]);
					++N;
				}
			}
			if( (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) &&
				phi_c[pos]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT &&
				phi_tmp[pos] < E_EPSIL && valid[pos] ){
				if( (vel_index == 1 && voxel.InSolid(i+1,j,k)) ||
					(vel_index == 2 && voxel.InSolid(i,j+1,k)) ||
					(vel_index == 3 && voxel.InSolid(i,j,k+1)) ){
//					printf("at (%d,%d,%d), vel = %f ", i,j,k,vel[pos]);
					int NV = 0;
					float vtmp = 0.f;
					if(TestVelPointIsLiquid(voxel,i+1,j,k,vel_index,valid)){
						vtmp += vel[INDEX(i+1,j,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i-1,j,k,vel_index,valid)){
						vtmp += vel[INDEX(i-1,j,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j+1,k,vel_index,valid)){
						vtmp += vel[INDEX(i,j+1,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j-1,k,vel_index,valid)){
						vtmp += vel[INDEX(i,j-1,k)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j,k+1,vel_index,valid)){
						vtmp += vel[INDEX(i,j,k+1)];
						++NV;
					}
					if(TestVelPointIsLiquid(voxel,i,j,k-1,vel_index,valid)){
						vtmp += vel[INDEX(i,j,k-1)];
						++NV;
					}
					if(NV > 0)
						vel[pos] = vtmp / NV;
//					printf("after averaging, vel = %f \n ", vel[pos]);
					++N;
				}
			}
		}
	END_FOR

	printf(" for Vel_index [%d], There are %u isolated vel points \n", vel_index, N);
	delete [] valid;

}


void VectorField3D::ExtrapolateOneUsingFM(const Voxel &voxel,
		                  const char *airVelSet,
		                  float *phi_tmp, float *vel, int vel_index) {


	float delta = voxel.VoxelDelta();

	float minValue;
	MinHeap<TentativeGridPoint, float> heap;

	set_time_base();

	float candidates[6];
	Point candidatesP[6];
	char *valid  = new char[DimX*DimY*DimZ];
	char *needed = new char[DimX*DimY*DimZ];
	char *mask   = new char[DimX*DimY*DimZ];
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	SetZero(needed);
	SetZero(mask);
	SetValidPoint(voxel, valid, vel_index);

//	SetAirBoundary(voxel, vel_index, vel);

	FOR_EACH_CELL
//		map<u_int, KnownPoint>::iterator foundpos;
		TentativeGridPoint tp(0.f, i, j, k);
		u_int pos = INDEX(i,j,k);
		int n;
		if( vel_index == 1 && !(airVelSet[pos] & AIRUVELSET) ||
			vel_index == 2 && !(airVelSet[pos] & AIRVVELSET) ||
			vel_index == 3 && !(airVelSet[pos] & AIRWVELSET) ){
//			if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] < 3*delta &&
//								phi_tmp[INDEX(i,j,k)] > 0.f)
//						printf("phic = %f, phi_tmp = %f at (%d,%d,%d)\n",
//				    		phi_c[INDEX(i,j,k)], phi_tmp[INDEX(i,j,k)], i,j,k);

			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tp, n) && valid[pos] && phi_tmp[pos] > E_EPSIL){
				if( (vel_index == 1 && (voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k)) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 2 && (voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k)) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 3 && (voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1)) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ) ){
					needed[pos] = 1;
				}
			}
			if( (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) &&
				phi_c[pos]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT &&
				phi_tmp[pos] > E_EPSIL && valid[pos] ){
				needed[pos] = 1;
			}
		    if(voxel.InAir(i,j,k) && phi_c[pos]+E_EPSIL >= EXTRAPOLATE_VEL_LIMIT && valid[pos]){
		    	TentativeGridPoint it(0.f, i, j, k);
		    	vector<TentativeGridPoint> neighbors;
		        it.Neigbor(neighbors, DimX, DimY, DimZ);
		        for(int m=0; m<neighbors.size(); m++){
		        	if(it.RightNeighbor(neighbors[m])  && !voxel.InSolid(neighbors[m])
		        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        		if(vel_index == 1){
		        			needed[pos] = 1;
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
		        		}
	        		}
		        	if(it.FrontNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
		        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        		if(vel_index == 2){
		        			needed[pos] = 1;
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
		        		}
	        		}
		        	if(it.TopNeighbor(neighbors[m]) && !voxel.InSolid(neighbors[m])
		        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        		if(vel_index == 3){
		        			needed[pos] = 1;
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
		        		}
	        		}
		        }

		    }
		}
	END_FOR

	u_int N = 0;
	FOR_EACH_CELL
		if(needed[INDEX(i,j,k)])
			++N;
	END_FOR
	printf("there are %u points needing updating velocities \n ", N);
//#ifdef WIN32
//	hash_map<u_int, KnownPoint> pos_band;
//#else
//	unordered_map<u_int, KnownPoint> pos_band;
//#endif
	// add those vel. points with positive phi but already updated in
	// pressure solver into the "done" band
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if( vel_index == 1 && (airVelSet[pos] & AIRUVELSET) ||
			vel_index == 2 && (airVelSet[pos] & AIRVVELSET) ||
			vel_index == 3 && (airVelSet[pos] & AIRWVELSET) ){
//			KnownPoint tmp(i,j,k);
//			pos_band.insert(make_pair(pos,tmp));
			mask[pos] |= POSINDONEBAND;
		}
	END_FOR

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(needed[pos]){
			if(phi_tmp[pos] < 1.e-5f){
//				printf("ExtrapolateOneUsingFM: shouldn't reach here\n");
//				exit(1);
//				Point p = voxel.VelPosition(vel_index, i, j, k);
//				vel[pos] = ProcessTrialPoints(voxel, p, phi_tmp, needed, valid, vel, vel_index);
				TentativeGridPoint tp(phi_tmp[pos], i, j, k);
				vector<TentativeGridPoint> neighbors;
				tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
				int NV = 0;
				float vel_total = 0.f;
				for(int m=0; m<neighbors.size(); m++){
					if(!needed[Index(neighbors[m])] && valid[Index(neighbors[m])]){
						++NV;
						vel_total += vel[Index(neighbors[m])];
					}
				}
				if(NV > 0){
					vel[pos] = vel_total / NV;
//					KnownPoint tmp(i,j,k);
//					pos_band.insert(make_pair(pos,tmp));
					mask[pos] |= POSINDONEBAND;
#ifdef CUDA
					needed[pos] = 0;
#endif
				}
				else{
					printf("something wrong! for vel %d at (%d, %d, %d) phi = %f\n",
						vel_index, i, j, k,	phi_tmp[pos]);
					printf(" no neighbors has valid values \n");
					for(int m=0; m<neighbors.size(); m++)
						printf("neighbor [%d] has phi_tmp = %f, vel = %f, needed = %d, valid = %d \n ",
								m, phi_tmp[Index(neighbors[m])],  vel[Index(neighbors[m])],
								needed[Index(neighbors[m])], valid[Index(neighbors[m])]);
//					exit(1);
				}
			}
			else{
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], i, j, k);
				vector<TentativeGridPoint> neighbors;
				tp.AllNeigbor(neighbors, DimX, DimY, DimZ);
				Point p = voxel.VelPosition(vel_index,i,j,k);
				for(int m=0; m<6; m++){
					candidates[m]  = INFINITY;
				    candidatesP[m] = Point(0.f);
				}
				for(int m=0; m<6; m++){
					if(tp.NeighborExists(m, DimX, DimY, DimZ)){
						Point pn = voxel.VelPosition(vel_index,neighbors[m].ii, neighbors[m].jj, neighbors[m].kk);
						neighbors[m].value = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
							continue;
						else if(tp.value * neighbors[m].value <= 0.f){
							float r;
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
								printf("wrong! phi_min should be > 0 at [%d, %d, %d] with "
										"phi = %f phin[%d] = %f \n",
										i,j,k, tp.value, m, neighbors[m].value);
								printf(" this is %d neighbor in cell (%d, %d, %d) with phi = %f",
										m, neighbors[m].ii, neighbors[m].jj, neighbors[m].kk,
										neighbors[m].value);
								printf("candidate = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
										candidates[m], r, pn.x, pn.y, pn.z);
								printf("phi = %f, r= %f, pn.x = %f, pn.y = %f, pn.z = %f \n",
												tp.value, r, pn.x, pn.y, pn.z);
//								if(m == 3){
//									float r1 =  tp.value * pn.y - neighbors[m].value * p.y;
//									float r2 =  tp.value - neighbors[m].value;
//									printf(" r1 = %f, r 2= %f \n", r1, r2);
//								}
								candidates[m] = INFINITY;
								continue;
//								candidates[m] = fabs(candidates[m]);
							}
						}
					}
				}
				float phi_min[3], phimin = 0.f;
				Point phi_min_p[3];
				Point pi;
				phi_min[0]= min(candidates[0], candidates[1]);
				phi_min_p[0] = (candidates[0] < candidates[1]) ? candidatesP[0]: candidatesP[1];
				phi_min[1]= min(candidates[2], candidates[3]);
				phi_min_p[1] = (candidates[2] < candidates[3]) ? candidatesP[2]: candidatesP[3];
				phi_min[2]= min(candidates[4], candidates[5]);
				phi_min_p[2] = (candidates[4] < candidates[5]) ? candidatesP[4]: candidatesP[5];
				if(fabs(phi_min[0]) < E_EPSIL ||
				   fabs(phi_min[1]) < E_EPSIL ||
				   fabs(phi_min[2]) < E_EPSIL )
					continue;
				if( phi_min[0] != INFINITY &&
					phi_min[1] != INFINITY &&
					phi_min[2] != INFINITY ){  // 3 points forming a plane
//					phi_tmp[pos] = phi_min[0]*phi_min[1]*phi_min[2] /
//						   sqrt( phi_min[0]*phi_min[0]*phi_min[1]*phi_min[1] +
//								 phi_min[0]*phi_min[0]*phi_min[2]*phi_min[2] +
//								 phi_min[2]*phi_min[2]*phi_min[1]*phi_min[1] );
					phimin = phi_tmp[pos];
					float r = 1.f/(phi_min[0]*phi_min[0]);
					float s = 1.f/(phi_min[1]*phi_min[1]);
					float t = 1.f/(phi_min[2]*phi_min[2]);
					pi = phi_min_p[0];
					float rf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(rf)){
						printf("3 rf is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[0].value, neighbors[1].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[0].x, candidatesP[0].y, candidatesP[0].z,
								candidatesP[1].x, candidatesP[1].y, candidatesP[1].z);
						exit(1);
					}
					pi = phi_min_p[1];
					float sf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(sf)){
						printf("3 sf is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[2].value, neighbors[3].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[2].x, candidatesP[2].y, candidatesP[2].z,
								candidatesP[3].x, candidatesP[3].y, candidatesP[3].z);
						exit(1);
					}
					pi = phi_min_p[2];
					float tf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(tf)){
						printf("3 tf is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[4].value, neighbors[5].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[4].x, candidatesP[4].y, candidatesP[4].z,
								candidatesP[5].x, candidatesP[5].y, candidatesP[5].z);
						exit(1);
					}
					if(fabs(r+s+t)<E_EPSIL){
						printf("vel is nan (%d, %d, %d) needed = %d, vel_index = %d \n",
								i,j,k, needed[pos], vel_index);
						printf(" r = %f, s = %f, t =%f \n", r,s,t);
						continue;
					}
					else
						vel[pos] = (r * rf + s * sf + t * tf) / (r + s + t);

				}
				else if( phi_min[0] != INFINITY && phi_min[1] != INFINITY ){
					// 2 points forming a line segment
//					phi_tmp[pos] = TwoPointDistance(phi_min[0], phi_min[1]);
					phimin = phi_tmp[pos];
					float s = 1.f/(phi_min[0]*phi_min[0]);
					float t = 1.f/(phi_min[1]*phi_min[1]);
					pi = phi_min_p[0];
					float sf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(sf)){
						printf("2 [0 1] sf is nan (%d, %d, %d) vel_index = %d \n", i,j, k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[0].value, neighbors[1].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[0].x, candidatesP[0].y, candidatesP[0].z,
								candidatesP[1].x, candidatesP[1].y, candidatesP[1].z);
						exit(1);
					}
					pi = phi_min_p[1];
					float tf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(tf)){
						printf("2 [0 1] tf is nan (%d, %d, %d) vel_index = %d \n", i, j, k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[2].value, neighbors[3].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[2].x, candidatesP[2].y, candidatesP[2].z,
								candidatesP[3].x, candidatesP[3].y, candidatesP[3].z);
						exit(1);
					}
					if(fabs(s+t)<E_EPSIL){
						printf("vel is nan (%d, %d, %d) needed = %d, vel_index = %d \n",
								i,j,k, needed[pos], vel_index);
						printf(" s = %f, t =%f \n", s,t);
						continue;
					}
					else
						vel[pos] = (s * sf + t * tf) / (s + t);
				}
				else if( phi_min[1] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
//					phi_tmp[pos] = TwoPointDistance(phi_min[1], phi_min[2]);
					phimin = phi_tmp[pos];
					float s = 1.f/(phi_min[1]*phi_min[1]);
					float t = 1.f/(phi_min[2]*phi_min[2]);
					pi = phi_min_p[1];
					float sf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(sf)){
						printf("2 [1 2] sf is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[2].value, neighbors[3].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[2].x, candidatesP[2].y, candidatesP[2].z,
								candidatesP[3].x, candidatesP[3].y, candidatesP[3].z);
						exit(1);
					}
					pi = phi_min_p[2];
					float tf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(tf)){
						printf("2 [1 2] tf is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[4].value, neighbors[5].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[4].x, candidatesP[4].y, candidatesP[4].z,
								candidatesP[5].x, candidatesP[5].y, candidatesP[5].z);
						exit(1);
					}
					if(fabs(s+t)<E_EPSIL){
						printf("vel is nan (%d, %d, %d) needed = %d, vel_index = %d \n",
								i,j,k, needed[pos], vel_index);
						printf(" s = %f, t =%f \n", s,t);
						continue;
					}
					else
						vel[pos] = (s * sf + t * tf) / (s + t);
				}
				else if( phi_min[0] != INFINITY && phi_min[2] != INFINITY ){
					// 2 points forming a line segment
//					phi_tmp[pos] = TwoPointDistance(phi_min[0], phi_min[2]);
					phimin = phi_tmp[pos];
					float s = 1.f/(phi_min[0]*phi_min[0]);
					float t = 1.f/(phi_min[2]*phi_min[2]);
					pi = phi_min_p[0];
					float sf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(sf)){
						printf("2 [0 2] sf is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[0].value, neighbors[1].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[0].x, candidatesP[0].y, candidatesP[0].z,
								candidatesP[1].x, candidatesP[1].y, candidatesP[1].z);
						exit(1);
					}
					pi = phi_min_p[2];
					float tf = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(tf)){
						printf("2 [0 2] tf is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[4].value, neighbors[5].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[4].x, candidatesP[4].y, candidatesP[4].z,
								candidatesP[5].x, candidatesP[5].y, candidatesP[5].z);
						exit(1);
					}
					if(fabs(s+t)<E_EPSIL){
						printf("vel is nan (%d, %d, %d) needed = %d, vel_index = %d \n",
								i,j,k, needed[pos], vel_index);
						printf(" s = %f, t =%f \n", s,t);
						continue;
					}
					else
						vel[pos] = (s * sf + t * tf) / (s + t);
				}
				else if(phi_min[0] != INFINITY){
//					phi_tmp[pos] = phi_min[0];
					phimin = phi_tmp[pos];
					pi = phi_min_p[0];
					vel[pos] = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(vel[pos])){
						printf("1 0 vel is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[0].value, neighbors[1].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[0].x, candidatesP[0].y, candidatesP[0].z,
								candidatesP[1].x, candidatesP[1].y, candidatesP[1].z);
						exit(1);
					}
				}
				else if(phi_min[1] != INFINITY){
//					phi_tmp[pos] = phi_min[1];
					phimin = phi_tmp[pos];
					pi = phi_min_p[1];
					vel[pos] = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(vel[pos])){
						printf("1 1 vel is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[2].value, neighbors[3].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[2].x, candidatesP[2].y, candidatesP[2].z,
								candidatesP[3].x, candidatesP[3].y, candidatesP[3].z);
						exit(1);
					}
				}
				else if(phi_min[2] != INFINITY){
//					phi_tmp[pos] = phi_min[2];
					phimin = phi_tmp[pos];
					pi = phi_min_p[2];
					vel[pos] = ProcessTrialPoints(voxel, pi, phi_tmp, needed, valid, vel, vel_index);
					if(isnan(vel[pos])){
						printf("1 2 vel is nan (%d, %d, %d) vel_index = %d \n", i,j,k, vel_index);
						printf(" phi = %f, phi_m = %f, phi_p = %f \n",
								tp.value, neighbors[4].value, neighbors[5].value);
						printf(" p_m = (%f, %f, %f) p_p = (%f, %f, %f) \n",
								candidatesP[4].x, candidatesP[4].y, candidatesP[4].z,
								candidatesP[5].x, candidatesP[5].y, candidatesP[5].z);
						exit(1);
					}
				}
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
				if(phimin != INFINITY){
					if(i == I && j == J && k == K)
						printf("point (%d, %d, %d) in band with vel_index = %d \n", i,j,k,vel_index);
//					KnownPoint tmp(i,j,k);
//					pos_band.insert(make_pair(pos,tmp));
					mask[pos] |= POSINDONEBAND;
#ifdef CUDA
					needed[pos] = 0;
#endif
				}
			}
			if(isnan(vel[pos])){
				printf("vel is nan (%d, %d, %d) needed = %d, vel_index = %d, phi= %f\n",
						i,j,k, needed[pos], vel_index, phi_tmp[pos]);
				printf("X: phi(%d,%d,%d) = %f, phi(%d,%d,%d) = %f \n",
						i-1, j, k, phi_tmp[INDEX(i-1,j,k)], i+1, j, k, phi_tmp[INDEX(i+1,j,k)]);
				printf("Y: phi(%d,%d,%d) = %f, phi(%d,%d,%d) = %f \n",
						i, j-1, k, phi_tmp[INDEX(i,j-1,k)], i, j+1, k, phi_tmp[INDEX(i,j+1,k)]);
				printf("Z: phi(%d,%d,%d) = %f, phi(%d,%d,%d) = %f \n",
						i, j, k-1, phi_tmp[INDEX(i,j,k-1)], i, j, k+1, phi_tmp[INDEX(i,j,k+1)]);
				exit(1);
			}
		}
	END_FOR
	printf("after process trial point (%d, %d, %d) vel = %f, needed = %d, vel_index = %d \n",
				I, J, K, vel[INDEX(I,J,K)], needed[INDEX(I,J,K)], vel_index);
#ifdef CUDA
	int dims[3];
	dims[0] = DimX; dims[1] = DimY; dims[2] = DimZ;
	float dtau = delta / 6;
	ExtrapolateVelCUDA( vel, phi_tmp, needed, valid,
						  I, J, K,
	 				      delta, dtau, INFINITY, FASTMARCH_LIMIT, dims );
#ifdef WIN32
	hash_map<u_int, TentativeGridPoint> nonValid;
	hash_map<u_int, TentativeGridPoint>::iterator nonValid_iter;
#else
	unordered_map<u_int, TentativeGridPoint> nonValid
	unordered_map<u_int, TentativeGridPoint>::iterator nonValid_iter;
#endif
	int iters = 0;
	while(iters < 10 && HasNonUpdatedPoints(voxel, valid, phi_tmp, vel_index, nonValid) > 0){
		printf("C. at iter %d, there are %u points not updated \n ", iters, nonValid.size());
		++iters;
		for(nonValid_iter=nonValid.begin(); nonValid_iter!=nonValid.end(); ++nonValid_iter){
			TentativeGridPoint &tp = nonValid_iter->second;
			int i = tp.ii;
			int j = tp.jj;
			int k = tp.kk;
			u_int index = nonValid_iter->first;
			vector<TentativeGridPoint> neighbors;
		    tp.Neigbor(neighbors, DimX, DimY, DimZ);
			float vel_avg = 0.f;
			int NV = 0;
			for(int m=0; m<neighbors.size(); m++){
				if(valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]){
					++NV;
					vel_avg += vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
				}
			}
			if(NV == 0){

			}
			else{
				vel[index] = vel_avg / NV;
				valid[index] = 1;
			}
//			printf(" NV = %d, at (%d, %d, %d) vel = %f, with phi = %f, vel_index = %d \n",
//				 NV, i, j, k, vel[index], phi_tmp[index], vel_index);
//			printf(" phi_tmp[i+1] = %f, phiobj[i+1] = %f, vel[i+1] = %f valid[i+1] = %d\n",
//					phi_tmp[INDEX(i+1,j,k)], phiobj[INDEX(i+1,j,k)], vel[INDEX(i+1,j,k)], valid[INDEX(i+1,j,k)]);
//			printf(" phi_tmp[i-1] = %f, phiobj[i-1] = %f, vel[i-1] = %f valid[i-1] = %d\n",
//					phi_tmp[INDEX(i-1,j,k)], phiobj[INDEX(i-1,j,k)], vel[INDEX(i-1,j,k)], valid[INDEX(i-1,j,k)]);
//			printf(" phi_tmp[j+1] = %f, phiobj[j+1] = %f, vel[j+1] = %f valid[j+1] = %d\n",
//					phi_tmp[INDEX(i,j+1,k)], phiobj[INDEX(i,j+1,k)], vel[INDEX(i,j+1,k)], valid[INDEX(i,j+1,k)]);
//			printf(" phi_tmp[j-1] = %f, phiobj[j-1] = %f, vel[j-1] = %f valid[j-1] = %d\n",
//					phi_tmp[INDEX(i,j-1,k)], phiobj[INDEX(i,j-1,k)], vel[INDEX(i,j-1,k)], valid[INDEX(i,j-1,k)]);
//			printf(" phi_tmp[k+1] = %f, phiobj[k+1] = %f, vel[k+1] = %f valid[k+1] = %d\n",
//					phi_tmp[INDEX(i,j,k+1)], phiobj[INDEX(i,j,k+1)], vel[INDEX(i,j,k+1)], valid[INDEX(i,j,k+1)]);
//			printf(" phi_tmp[k-1] = %f, phiobj[k-1] = %f, vel[k-1] = %f valid[k-1] = %d\n",
//					phi_tmp[INDEX(i,j,k-1)], phiobj[INDEX(i,j,k-1)], vel[INDEX(i,j,k-1)], valid[INDEX(i,j,k-1)] );
		}
		nonValid.clear();
	}
#else
	ProcessFarAwayPoints(voxel, phi_tmp, valid, needed,	vel_index, mask, vel);
#endif
	printf("after process faraway point (%d, %d, %d) vel = %f, needed = %d, vel_index = %d \n",
				I, J, K, vel[INDEX(I,J,K)], needed[INDEX(I,J,K)], vel_index);
	delete [] valid;
	delete [] needed;
	delete [] mask;
	printf(" Time spent in ExtrapolateOneUsingFM = %f \n", get_time_in_seconds());
}

void VectorField3D::ProcessFarAwayPoints(const Voxel &voxel, float *phi_tmp,
		                         char *valid, char *needed,	 int vel_index,
		                         char *mask, float *vel){
	float delta = voxel.VoxelDelta();

	int num_points = 0;
	//printf(" band has = %d, heap has = %d \n", num_points);
	MinHeap<TentativeGridPoint, float> heap;


//	for(pos=posband.begin(); pos!=posband.end(); ++pos){
//		KnownPoint &p = pos->second;
//		TentativeGridPoint tp(0.f, p);
//		vector<TentativeGridPoint> neigbor;
//	    tp.Neigbor(neigbor, DimX, DimY, DimZ);
//
//		for(int i=0;i<neigbor.size();i++){
//			foundpos1 = posband.find(Index(neigbor[i]));
//			if(foundpos1==posband.end()){ // if this neighbor is not in accepted band
//				foundpos = posadj.find(Index(neigbor[i]));
//				if(foundpos==posadj.end()){ // if this neighbor is not in close band
//					if(needed[Index(neigbor[i])]){
//						KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//						posadj.insert(make_pair(Index(neigbor[i]),tmp)); // add it to close band
//					}
//				}
//			}
//		}
//	}

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//		printf(" at %d, %d, %d \n", i, j, k);
		if(mask[pos] & POSINDONEBAND){
			TentativeGridPoint tp(0.f, i, j, k);
			vector<TentativeGridPoint> neigbor;
			tp.Neigbor(neigbor, DimX, DimY, DimZ);
			for(int n=0;n<neigbor.size();n++){
				int ii = neigbor[n].ii;
				int jj = neigbor[n].jj;
				int kk = neigbor[n].kk;
//				printf(" n = %d for ( %d, %d, %d) \n", n, ii, jj, kk);
				if(!(mask[INDEX(ii,jj,kk)] & POSINDONEBAND)){ // if this neighbor is not in accepted band
					if(!(mask[INDEX(ii,jj,kk)] & POSINTRIALBAND)){ // if this neighbor is not in close band
						if(needed[Index(neigbor[n])]){
							mask[INDEX(ii,jj,kk)] |= POSINTRIALBAND;
							++num_points;
						}
					}
				}
			}
		}
	END_FOR
	printf(" there are %d trial points \n", num_points);
//	for(pos=posadj.begin(); pos!=posadj.end(); ++pos){
//		KnownPoint &p = pos->second;
//		u_int index = pos->first;
//		float tentativeValue = UpdateTentativeValue(p, index, posband, voxel, phi_tmp);
//		TentativeGridPoint tp(tentativeValue, p);
//		heap.insert(tp, tentativeValue);
//	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mask[pos] & POSINTRIALBAND){
			KnownPoint p(i,j,k);
			float tentativeValue = UpdateTentativeValue(p, pos, mask, (char)POSINDONEBAND, voxel, phi_tmp);
			TentativeGridPoint tp(tentativeValue, i, j, k);
			heap.insert(tp, tentativeValue);
		}
	END_FOR
//	printf(" Got here A \n");
	float minValue;
	do{
		TentativeGridPoint minTentative = heap.extract_min(minValue);
		++num_points;
		u_int indexDoneBand = Index(minTentative);
		phi_tmp[indexDoneBand] = minValue;
//		if(minTentative.ii == I && minTentative.jj == J && minTentative.kk == K)
//			printf("\n point(%d, %d, %d) phi = %f \n\n",
//					minTentative.ii,minTentative.jj,minTentative.kk, phi[POS(t)]);
//		printf("point added in band (%d, %d, %d)\n",t.ii, t.jj, t.kk);
//		posband.insert( make_pair(Index(minTentative), kp) );
//		posadj.erase(Index(minTentative));
		mask[indexDoneBand] |= POSINDONEBAND;
		mask[indexDoneBand] ^= POSINTRIALBAND;
//		printf("heap size = %d, pos_adj size = %d\n", heap.size(), num_points);
		//printf("# of points = %d, min value = %f\n", num_points, minValue);
	    vector<TentativeGridPoint> neigbor;
	    minTentative.Neigbor(neigbor, DimX, DimY, DimZ);
//	    assert(neigbor.size() == 6);
	    float s[6] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
	    float t[6] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
	    for(int i=0;i<neigbor.size();i++){
	    	u_int indexTrialBand = Index(neigbor[i]);
			if(!(mask[indexTrialBand] & POSINDONEBAND)){ // not currently in done band
//	    	foundinband = posband.find(Index(neigbor[i]));
//	    	if(foundinband == posband.end()){ // not currently in the band

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

				if(!(mask[indexTrialBand] & POSINTRIALBAND)){ // not currently in trial band
//	    		found = posadj.find(Index(neigbor[i]));
//	    		if(found == posadj.end()){ // not currently in the heap
	    			if(needed[indexTrialBand]){
//			    		KnownPoint tmp(neigbor[i].ii, neigbor[i].jj, neigbor[i].kk);
//			    		posadj.insert(make_pair(Index(neigbor[i]), tmp));
//		    		printf("insert neigbor i = %d, value = %f \n", i, neigbor[i].value );
	    				mask[indexTrialBand] |= POSINTRIALBAND;
			    		neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(),
			    				indexTrialBand, mask, (char)POSINDONEBAND, voxel, phi_tmp);
			    		heap.insert(neigbor[i], neigbor[i].value);
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
	    		else{	// already in heap,
	    			neigbor[i].value = UpdateTentativeValue(neigbor[i].Pos(),
	    					indexTrialBand, mask, (char)POSINDONEBAND, voxel, phi_tmp);
	    			heap.update(neigbor[i], neigbor[i].value);	// just update values
													    		// according to newly added point
//	    			if(neigbor[i].ii == I && neigbor[i].jj == J && neigbor[i].kk == K)
//		    			printf("\n update close neigbor i = %d, value = %f, id = (%d, %d, %d) phi = %f \n\n", i, neigbor[i].value,
//		    					neigbor[i].ii, neigbor[i].jj, neigbor[i].kk, phi[Index(neigbor[i])]);
//					printf("update neigbor i = %d, value = %f \n",
//	    							i, neigbor[i].value );
	    		}
//	    		s[i] = INFINITY;
//	    		t[i] = INFINITY;
	    	}
	    	else{
	    		s[i] = phi_tmp[Index(neigbor[i])];
	    		t[i] = vel[Index(neigbor[i])];
	    	}
	    }
	    float phi_x, A, a;
	    float phi_y, B, b;
	    float phi_z, C, c;
	    ProcessCoefficents(s[0], s[1], minValue, t[0], t[1], phi_x, A, a);
	    ProcessCoefficents(s[2], s[3], minValue, t[2], t[3], phi_y, B, b);
	    ProcessCoefficents(s[4], s[5], minValue, t[4], t[5], phi_z, C, c);

	    if( A != 0.f || B != 0.f || C != 0.f ){
//	    	printf("ii = %d, jj = %d, kk = %d, phi = %f, phi_x = %f, phi_y = %f, phi_z = %f \n ",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			minValue, phi_x, phi_y, phi_z);
	    	vel[Index(minTentative)] = (a+b+c)/(A+B+C);
	    }
	    else{
	    	printf("Error in VectorField3D::ExtractOneVelocityFM vel_index = %d! \n ",
	    			vel_index);
	    	printf(" ii = %d, jj = %d, kk = %d \n",
	    			minTentative.ii, minTentative.jj, minTentative.kk);
	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
	    			minValue, s[0], s[1]);
	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
		    			minValue, s[2], s[3]);
	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
	    		    	minValue, s[4], s[5]);
	    	printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
	    			vel[Index(minTentative)], t[0], t[1], t[2], t[3], t[4], t[5]);
	    	exit(1);
	    }
	    neigbor.clear();
//	    printf("band has = %d, heap has = %d \n", band.size(), heap.size() );
	} while(!heap.empty());

	printf(" heap has = %d \n", num_points);

}



float VectorField3D::ProcessTrialPoints(const Voxel &voxel, const Point &p,
				const float *phi_tmp, char *needed, char *valid, float *vel, int vel_index){

	float delta = voxel.VoxelDelta();

	int i0, j0, k0, i1, j1, k1;
	float s0, t0, s1, t1;
	int x0, y0, z0;
	float r0, r1;
	float xp = p.x/delta;
	float yp = p.y/delta;
	float zp = p.z/delta;


	xp = xp < 0.5f      ? 0.5f      : xp;
	xp = xp > DimX-1.5f ? DimX-1.5f : xp;
	yp = yp < 0.5f      ? 0.5f      : yp;
	yp = yp > DimY-1.5f ? DimY-1.5f : yp;
	zp = zp < 0.5f      ? 0.5f      : zp;
	zp = zp > DimZ-1.5f ? DimZ-1.5f : zp;



	if(vel_index == 1){
		i0=floor(xp)-1; if(i0 < 0) i0 = 0;
		i1 = i0 + 1;
		s1 = xp - 1 - i0;
		s0 = 1 - s1;
		j0 = floor(yp); y0 = j0;
		if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		j1=j0+1;
		if(yp-y0 < 0.5f) t1 = 0.5f+yp-y0;
		else t1 = yp-y0-0.5f;
		t0 = 1-t1;
		k0 = floor(zp); z0 = k0;
		if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
		k1=k0+1;
		if(zp-z0 < 0.5f) r1 = 0.5f+zp-z0;
		else r1 = zp-z0-0.5f;
		r0 = 1-r1;
	}
	else if(vel_index == 2){
		j0=floor(yp)-1; if(j0 < 0) j0 = 0;
		j1 = j0 + 1;
		t1 = yp - 1 - j0;
		t0 = 1 - t1;
		i0 = floor(xp); x0 = i0;
		if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		i1=i0+1;
		if(xp-x0 < 0.5f) s1 = 0.5f+xp-x0;
		else s1 = xp-x0-0.5f;
		s0 = 1-s1;
		k0 = floor(zp); z0 = k0;
		if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
		k1=k0+1;
		if(zp-z0 < 0.5f) r1 = 0.5f+zp-z0;
		else r1 = zp-z0-0.5f;
		r0 = 1-r1;
	}
	else if(vel_index == 3){
		k0=floor(zp)-1; if(k0 < 0) k0 = 0;
		k1 = k0+1;
		r1 = zp - 1 - k0;
		r0 = 1 - r1;
		i0 = floor(xp); x0 = i0;
		if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		i1=i0+1;
		if(xp-x0 < 0.5f) s1 = 0.5f+xp-x0;
		else s1 = xp-x0-0.5f;
		s0 = 1-s1;
		j0 = floor(yp); y0 = j0;
		if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		j1=j0+1;
		if(yp-y0 < 0.5f) t1 = 0.5f+yp-y0;
		else t1 = yp-y0-0.5f;
		t0 = 1-t1;
	}

	float weight[8];
	weight[0] = r0 * s0 * t0;
	weight[1] = r0 * s1 * t0;
	weight[2] = r0 * s0 * t1;
	weight[3] = r0 * s1 * t1;
	weight[4] = r1 * s0 * t0;
	weight[5] = r1 * s1 * t0;
	weight[6] = r1 * s0 * t1;
	weight[7] = r1 * s1 * t1;

	if(!valid[INDEX(i0,j0,k0)] || needed[INDEX(i0,j0,k0)])
		weight[0] = 0.f;
	if(!valid[INDEX(i1,j0,k0)] || needed[INDEX(i1,j0,k0)])
		weight[1] = 0.f;
	if(!valid[INDEX(i0,j1,k0)] || needed[INDEX(i0,j1,k0)])
		weight[2] = 0.f;
	if(!valid[INDEX(i1,j1,k0)] || needed[INDEX(i1,j1,k0)])
		weight[3] = 0.f;
	if(!valid[INDEX(i0,j0,k1)] || needed[INDEX(i0,j0,k1)])
		weight[4] = 0.f;
	if(!valid[INDEX(i1,j0,k1)] || needed[INDEX(i1,j0,k1)])
		weight[5] = 0.f;
	if(!valid[INDEX(i0,j1,k1)] || needed[INDEX(i0,j1,k1)])
		weight[6] = 0.f;
	if(!valid[INDEX(i1,j1,k1)] || needed[INDEX(i1,j1,k1)])
		weight[7] = 0.f;

	float total_weights = 0.f;
	for(int n=0; n<8; ++n)
		total_weights += weight[n];
	if(total_weights == 0.f){
		printf(" Interpolate point is (%f, %f, %f) \n", p.x, p.y, p.z);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f \n", i0, j0, k0,
				valid[INDEX(i0,j0,k0)], i0, j0, k0, needed[INDEX(i0,j0,k0)], weight[0]);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f\n", i1, j0, k0,
				valid[INDEX(i1,j0,k0)], i1, j0, k0, needed[INDEX(i1,j0,k0)], weight[1]);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f\n", i0, j1, k0,
				valid[INDEX(i0,j1,k0)], i0, j1, k0, needed[INDEX(i0,j1,k0)], weight[2]);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f\n", i1, j1, k0,
				valid[INDEX(i1,j1,k0)], i1, j1, k0, needed[INDEX(i1,j1,k0)], weight[3]);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f\n", i0, j0, k1,
				valid[INDEX(i0,j0,k1)], i0, j0, k1, needed[INDEX(i0,j0,k1)], weight[4]);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f\n", i1, j0, k1,
				valid[INDEX(i1,j0,k1)], i1, j0, k1, needed[INDEX(i1,j0,k1)], weight[5]);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f\n", i0, j1, k1,
				valid[INDEX(i0,j1,k1)], i0, j1, k1, needed[INDEX(i0,j1,k1)], weight[6]);
		printf("valid[INDEX(%d,%d,%d)] = %d needed[INDEX(%d,%d,%d)] = %d  weight = %f\n", i1, j1, k1,
				valid[INDEX(i1,j1,k1)], i1, j1, k1, needed[INDEX(i1,j1,k1)], weight[7]);
		printf(" r0 = %f, r1 = %f, s0 = %f, s1 = %f, t0 = %f, t1 = %f \n",
				r0, r1, s0, s1, t0, t1);
		printf(" xp = %f, yp = %f, zp = %f, x0 = %d, y0 = %d, z0 = %d \n",
				xp, yp, zp, x0, y0, z0);

		int n = 0;
		float veln = 0;
		if(valid[INDEX(i0,j0,k0)] && !needed[INDEX(i0,j0,k0)]){
			++n;
			veln += vel[INDEX(i0,j0,k0)];
		}
		if(valid[INDEX(i1,j0,k0)] || !needed[INDEX(i1,j0,k0)]){
			++n;
			veln += vel[INDEX(i1,j0,k0)];
		}
		if(valid[INDEX(i0,j1,k0)] || !needed[INDEX(i0,j1,k0)]){
			++n;
			veln += vel[INDEX(i0,j1,k0)];
		}
		if(valid[INDEX(i1,j1,k0)] || !needed[INDEX(i1,j1,k0)]){
			++n;
			veln += vel[INDEX(i1,j1,k0)];
		}
		if(valid[INDEX(i0,j0,k1)] || !needed[INDEX(i0,j0,k1)]){
			++n;
			veln += vel[INDEX(i0,j0,k1)];
		}
		if(valid[INDEX(i1,j0,k1)] || !needed[INDEX(i1,j0,k1)]){
			++n;
			veln += vel[INDEX(i1,j0,k1)];
		}
		if(valid[INDEX(i0,j1,k1)] || !needed[INDEX(i0,j1,k1)]){
			++n;
			veln += vel[INDEX(i0,j1,k1)];
		}
		if(valid[INDEX(i1,j1,k1)] || !needed[INDEX(i1,j1,k1)]){
			++n;
			veln += vel[INDEX(i1,j1,k1)];
		}
		return veln / n;
	}
	else{
		return ( weight[0]*vel[INDEX(i0,j0,k0)]+weight[2]*vel[INDEX(i0,j1,k0)]+
				 weight[1]*vel[INDEX(i1,j0,k0)]+weight[3]*vel[INDEX(i1,j1,k0)]+
				 weight[4]*vel[INDEX(i0,j0,k1)]+weight[6]*vel[INDEX(i0,j1,k1)]+
				 weight[5]*vel[INDEX(i1,j0,k1)]+weight[7]*vel[INDEX(i1,j1,k1)] )/total_weights;
	}



}

#else
void VectorField3D::UpdateCells(MinHeap<TentativeGridPoint, float> &heap, int vel_index, float *vel, float *phi_tmp){
	float minValue;
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
	    		phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.RightNeighbor(neighbors[m])){
	    		phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.BackNeighbor(neighbors[m])){
	    		phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.FrontNeighbor(neighbors[m])){
	    		phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.BottomNeighbor(neighbors[m])){
	    		phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    	}
	    	if(minTentative.TopNeighbor(neighbors[m])){
	    		phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
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
	    	vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
	    }
	    else{
//	    	printf("Error in VectorField3D::ExtractOneVelocity vel_index = %d! \n ",
//	    			vel_index);
//	    	printf(" ii = %d, jj = %d, kk = %d \n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk);
//	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
//	    			minValue, phi_minus_x, phi_plus_x);
//	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
//		    			minValue, phi_minus_y, phi_plus_y);
//	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
//	    		    	minValue, phi_minus_z, phi_plus_z);
//	    	printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//	    		    			vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
//	    		    			vel_minus_z, vel_plus_z);
//	    	exit(1);
	    }
	    if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//	    if(vel_index == 3){
//	    if(fabs(vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]) > 40.f){
//	    if(phi_minus_x < 0.f && phi_plus_x < 0.f ||
//    	   phi_minus_y < 0.f && phi_plus_y < 0.f ||
//    	   phi_minus_z < 0.f && phi_plus_z < 0.f ) {
//	    if( phi_x == 0.f && A != 0.f || phi_y == 0.f && B != 0.f || phi_z == 0.f && C != 0.f ){
	    	printf("(ExtrapolateOneVelocityIntoObjectForAdvection) vel_index = %d\n", vel_index);
	    	printf("after interp: at i = %d, j = %d, k = %d vel = %f \n\n",
	    			minTentative.ii, minTentative.jj, minTentative.kk,
	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]);
	    	printf("before interp: A = %f, B = %f, C = %f, a = %f, b = %f, c = %f \n ",
	    			A, B, C, a, b, c);
	    	printf(" phix = %f, phiy = %f, phiz = %f \n",
	    			phi_x, phi_y, phi_z);
	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
	    			minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
		    			minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
	    		    	minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
//	    	exit(1);
	    }

	} while(!heap.empty());
}

void VectorField3D::ExtrapolateOneForProjection(const Voxel &voxel,
		                         map<u_int, KnownPoint> &knownPoint,
		                        float *phi_tmp, float *phi_obj, float *vel, int vel_index) {

	float delta = voxel.VoxelDelta();

	float minValue;
	MinHeap<TentativeGridPoint, float> heap;
	map<u_int, KnownPoint>::iterator pos;
	map<u_int, KnownPoint>::iterator foundpos;

	char *valid = new char[DimX*DimY*DimZ];
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	SetValidPoint(voxel, valid, vel_index);

//	SetAirBoundary(voxel, vel_index, vel);

	FOR_EACH_CELL
//		map<u_int, KnownPoint>::iterator foundpos;
		TentativeGridPoint tp(0.f, i, j, k);
		int n;
		foundpos = knownPoint.find(INDEX(i,j,k));
		if( foundpos == knownPoint.end() ){
//			if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] < 3*delta &&
//								phi_tmp[INDEX(i,j,k)] > 0.f)
//						printf("phic = %f, phi_tmp = %f at (%d,%d,%d)\n",
//				    		phi_c[INDEX(i,j,k)], phi_tmp[INDEX(i,j,k)], i,j,k);

			if(voxel.InSolid(i,j,k) && voxel.CloseToNonSolid(tp, n) && valid[INDEX(i,j,k)] && phi_tmp[INDEX(i,j,k)] > 0.f){
				if( (vel_index == 1 && phi_c[INDEX(i+1,j,k)] < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 2 && phi_c[INDEX(i,j+1,k)] < EXTRAPOLATE_VEL_LIMIT ) ||
					(vel_index == 3 && phi_c[INDEX(i,j,k+1)] < EXTRAPOLATE_VEL_LIMIT ) ){
					KnownPoint tmp(i, j, k);
					TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
					heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
				}
			}
			if( (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k))
			   && phi_c[INDEX(i,j,k)] < EXTRAPOLATE_VEL_LIMIT &&
				phi_tmp[INDEX(i,j,k)] > 0.f && valid[INDEX(i,j,k)]
			    ){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
			}
		    if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= EXTRAPOLATE_VEL_LIMIT &&
		    		valid[INDEX(i,j,k)]){
		    	TentativeGridPoint it(0.f, i, j, k);
		    	vector<TentativeGridPoint> neighbors;
		        it.Neigbor(neighbors, DimX, DimY, DimZ);
		        for(int m=0; m<neighbors.size(); m++){
		        	if(it.RightNeighbor(neighbors[m])
		        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] < EXTRAPOLATE_VEL_LIMIT){
		        		if(vel_index == 1){
		        			KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
		        		}
	        		}
		        	if(it.FrontNeighbor(neighbors[m])
		        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] < EXTRAPOLATE_VEL_LIMIT){
		        		if(vel_index == 2){
		        			KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
		        		}
	        		}
		        	if(it.TopNeighbor(neighbors[m])
		        	  && phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] < EXTRAPOLATE_VEL_LIMIT){
		        		if(vel_index == 3){
		        			KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
		        		}
	        		}
		        }

		    }
		}
	END_FOR

	printf("for vel %d, heap size = %d \n", vel_index, heap.size());

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
		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
		knownPoint.insert( make_pair(Index(minTentative), t) );
		vector<TentativeGridPoint> neighbors;
	    minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
	    for(int m=0; m<neighbors.size(); m++){
	    	if(minTentative.LeftNeighbor(neighbors[m])){
	    		phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
    			if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_minus_x = INFINITY;
	    	}
	    	if(minTentative.RightNeighbor(neighbors[m])){
	    		phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			 phi_plus_x = INFINITY;
//	    		if(vel_index == 1 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii+1,neighbors[m].jj,neighbors[m].kk))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_x = INFINITY;
	    	}
	    	if(minTentative.BackNeighbor(neighbors[m])){
	    		phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			phi_minus_y = INFINITY;
	    	}
	    	if(minTentative.FrontNeighbor(neighbors[m])){
	    		phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			 phi_plus_y = INFINITY;
//	    		if(vel_index == 2 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj+1,neighbors[m].kk))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_y = INFINITY;
	    	}
	    	if(minTentative.BottomNeighbor(neighbors[m])){
	    		phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			  phi_minus_z = INFINITY;
	    	}
	    	if(minTentative.TopNeighbor(neighbors[m])){
	    		phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			  phi_plus_z = INFINITY;
//	    		if(vel_index == 3 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk+1))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_z = INFINITY;
	    	}
//	    	if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//	    	    if(voxel.InSolid(neighbors[m]))
//	    	    	printf("neighbor at (%d, %d, %d) in solid \n",
//	    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
//	    	    else
//	    	    	printf("neighbor at (%d, %d, %d) not in solid \n",
//    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
//	    	}
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
	    	vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
	    }
	    else{
//	    	printf("Error in VectorField3D::ExtractOneForProjection vel_index = %d! \n ",
//	    		    			vel_index);
//	    	printf(" ii = %d, jj = %d, kk = %d \n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk);
//	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
//	    			minValue, phi_minus_x, phi_plus_x);
//	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
//		    			minValue, phi_minus_y, phi_plus_y);
//	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
//	    		    	minValue, phi_minus_z, phi_plus_z);
//	    	printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//	    		    			vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
//	    		    			vel_minus_z, vel_plus_z);
//	    	printf(" phic = %f, phi-x = %f, phi+x = %f, phi-y = %f, phi+y = %f, phi-z = %f, phi+z = %f \n ",
//						phi_c[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii-1, minTentative.jj, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii+1, minTentative.jj, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj-1, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj+1, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj, minTentative.kk-1)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj, minTentative.kk+1)]);
	    	valid[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = 0;
//	    	exit(1);
	    }
	    if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
	    	printf("(ExtrapolateOneForProjection) vel_index = %d phi = %f\n",
	    			vel_index, minValue);
	    	printf("after interp: at i = %d, j = %d, k = %d vel = %f \n\n",
	    			minTentative.ii, minTentative.jj, minTentative.kk,
	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]);
	    	printf("before interp: A = %f, B = %f, C = %f, a = %f, b = %f, c = %f \n ",
	    			A, B, C, a, b, c);
	    	printf(" valid = %d, valid_minus_x = %d, valid_plus_x = %d \n ",
	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
	    		    	valid[INDEX(minTentative.ii-1,minTentative.jj,minTentative.kk)],
	    		    	valid[INDEX(minTentative.ii+1,minTentative.jj,minTentative.kk)]);
	    	printf(" valid = %d, valid_minus_y = %d, valid_plus_y = %d\n ",
	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
	    		    	valid[INDEX(minTentative.ii,minTentative.jj-1,minTentative.kk)],
	    		    	valid[INDEX(minTentative.ii,minTentative.jj+1,minTentative.kk)]);
	    	printf(" valid = %d, valid_minus_z = %d, valid_plus_z = %d\n ",
	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk-1)],
	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk+1)]);
	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
	    			minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
		    			minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
	    		    	minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
	    }

	} while(!heap.empty());

	while(HasNonUpdatedPoints(voxel, valid, phi_tmp, vel_index, heap) > 0){

		u_int N0 = heap.size();
		u_int N1 = 0, N2 = 0;
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
				KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
		//		knownPoint.insert( make_pair(Index(minTentative), t) );
				vector<TentativeGridPoint> neighbors;
			    minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
			    for(int m=0; m<neighbors.size(); m++){
			    	if(minTentative.LeftNeighbor(neighbors[m])){
			    		phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			    		vel_minus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		//	    		if(voxel.InSolid(neighbors[m]))
		    			if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
			    			 phi_minus_x = INFINITY;
			    	}
			    	if(minTentative.RightNeighbor(neighbors[m])){
			    		phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			    		vel_plus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		//	    		if(voxel.InSolid(neighbors[m]))
		//	    			 phi_plus_x = INFINITY;
		//	    		if(vel_index == 1 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
		//	    				&& voxel.InSolid(neighbors[m].ii+1,neighbors[m].jj,neighbors[m].kk))
			    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
			    			 phi_plus_x = INFINITY;
			    	}
			    	if(minTentative.BackNeighbor(neighbors[m])){
			    		phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			    		vel_minus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		//	    		if(voxel.InSolid(neighbors[m]))
			    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
			    			phi_minus_y = INFINITY;
			    	}
			    	if(minTentative.FrontNeighbor(neighbors[m])){
			    		phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			    		vel_plus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		//	    		if(voxel.InSolid(neighbors[m]))
		//	    			 phi_plus_y = INFINITY;
		//	    		if(vel_index == 2 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
		//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj+1,neighbors[m].kk))
			    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
			    			 phi_plus_y = INFINITY;
			    	}
			    	if(minTentative.BottomNeighbor(neighbors[m])){
			    		phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			    		vel_minus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		//	    		if(voxel.InSolid(neighbors[m]))
			    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
			    			  phi_minus_z = INFINITY;
			    	}
			    	if(minTentative.TopNeighbor(neighbors[m])){
			    		phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
			    		vel_plus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		//	    		if(voxel.InSolid(neighbors[m]))
		//	    			  phi_plus_z = INFINITY;
		//	    		if(vel_index == 3 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
		//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk+1))
			    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
			    			 phi_plus_z = INFINITY;
			    	}
		//	    	if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
		//	    	    if(voxel.InSolid(neighbors[m]))
		//	    	    	printf("neighbor at (%d, %d, %d) in solid \n",
		//	    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
		//	    	    else
		//	    	    	printf("neighbor at (%d, %d, %d) not in solid \n",
		//    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
		//	    	}
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
			    	vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
			    	valid[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = 1;
			    	++N1;
			    }
			    else{
					float vel_avg = 0.f;
					int NV = 0;
		 			for(int m=0; m<neighbors.size(); m++){
		 				if(valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]){
		 					++NV;
		 					vel_avg += vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		 				}
		 			}
		 			if(NV == 0){
	//	 				printf("serious problem!\n");
	//	 				printf("Error in VectorField3D::ExtractOneForProjection vel_index = %d! \n ",
	//										vel_index);
	//					printf(" ii = %d, jj = %d, kk = %d \n",
	//							minTentative.ii, minTentative.jj, minTentative.kk);
	//					printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
	//							minValue, phi_minus_x, phi_plus_x);
	//					printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
	//								minValue, phi_minus_y, phi_plus_y);
	//					printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
	//								minValue, phi_minus_z, phi_plus_z);
	//					printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
	//							vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
	//										vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
	//										vel_minus_z, vel_plus_z);
		 			}
		 			else{
		 				vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = vel_avg / NV;
		 				valid[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = 1;
		 				++N2;
		 			}
			    }
			}while(!heap.empty());

			printf("\nthere are %u points reupdated, %u through usual manner, %u through one-ring averaging \n\n", N0, N1, N2);
		}
	if(valid) delete [] valid;

//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && phi_tmp[INDEX(i,j,k)] > 0.f){
//			foundpos = knownPoint.find(INDEX(i,j,k));
//			if( foundpos == knownPoint.end() ){
//				vel[INDEX(i,j,k)] = 0.f;
////				printf("Reset point (%d,%d,%d) phi(%f) vel(%d) to zero\n",
////						i, j, k, phi_tmp[INDEX(i,j,k)], vel_index);
//			}
//		}
//	END_FOR
	knownPoint.clear();
//	SetSolidBoundaryNoExtrapolate(voxel, vel_index, vel);
}


#endif

void VectorField3D::ExtrapolateOneForProjection(const Voxel &voxel,
		                        float *phi_tmp, float *vel, int vel_index) {

	float delta = voxel.VoxelDelta();

	float minValue;
	MinHeap<TentativeGridPoint, float> heap;


	char *valid = new char[DimX*DimY*DimZ];
//	memset(valid, 0, DimX*DimY*DimZ);
	SetZero(valid);
	SetValidPoint(voxel, valid, 0);

//	SetAirBoundary(voxel, vel_index, vel);

	FOR_EACH_CELL
//		map<u_int, KnownPoint>::iterator foundpos;


//			if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] < 3*delta &&
//								phi_tmp[INDEX(i,j,k)] > 0.f)
//						printf("phic = %f, phi_tmp = %f at (%d,%d,%d)\n",
//				    		phi_c[INDEX(i,j,k)], phi_tmp[INDEX(i,j,k)], i,j,k);

			if( (voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && phi_tmp[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT ){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
			}
		    if(voxel.InAir(i,j,k) && phi_tmp[INDEX(i,j,k)]+E_EPSIL >= EXTRAPOLATE_VEL_LIMIT){
		    	TentativeGridPoint it(0.f, i, j, k);
		    	vector<TentativeGridPoint> neighbors;
		        it.Neigbor(neighbors, DimX, DimY, DimZ);
		        for(int m=0; m<neighbors.size(); m++){
		        	if(it.RightNeighbor(neighbors[m])
		        	  && phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        			KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
	        		}
		        	if(it.FrontNeighbor(neighbors[m])
		        	  && phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        			KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
	        		}
		        	if(it.TopNeighbor(neighbors[m])
		        	  && phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
		        			KnownPoint tmp(i,j,k);
							TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
							heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
							if(phi_tmp[INDEX(i,j,k)]< 0.f){
								printf("Error! shouldn't extrapolate this one\n");
								exit(1);
							}
	        		}
		        }


		}
	END_FOR

	printf("for vel %d, heap size = %d \n", vel_index, heap.size());

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
		KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
//		knownPoint.insert( make_pair(Index(minTentative), t) );
		vector<TentativeGridPoint> neighbors;
	    minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
	    for(int m=0; m<neighbors.size(); m++){
	    	if(minTentative.LeftNeighbor(neighbors[m])){
	    		phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
    			if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_minus_x = INFINITY;
	    	}
	    	if(minTentative.RightNeighbor(neighbors[m])){
	    		phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			 phi_plus_x = INFINITY;
//	    		if(vel_index == 1 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii+1,neighbors[m].jj,neighbors[m].kk))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_x = INFINITY;
	    	}
	    	if(minTentative.BackNeighbor(neighbors[m])){
	    		phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			phi_minus_y = INFINITY;
	    	}
	    	if(minTentative.FrontNeighbor(neighbors[m])){
	    		phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			 phi_plus_y = INFINITY;
//	    		if(vel_index == 2 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj+1,neighbors[m].kk))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_y = INFINITY;
	    	}
	    	if(minTentative.BottomNeighbor(neighbors[m])){
	    		phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_minus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			  phi_minus_z = INFINITY;
	    	}
	    	if(minTentative.TopNeighbor(neighbors[m])){
	    		phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	    		vel_plus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
//	    		if(voxel.InSolid(neighbors[m]))
//	    			  phi_plus_z = INFINITY;
//	    		if(vel_index == 3 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk+1))
	    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
	    			 phi_plus_z = INFINITY;
	    	}
//	    	if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//	    	    if(voxel.InSolid(neighbors[m]))
//	    	    	printf("neighbor at (%d, %d, %d) in solid \n",
//	    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
//	    	    else
//	    	    	printf("neighbor at (%d, %d, %d) not in solid \n",
//    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
//	    	}
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
	    	vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
	    }
	    else{
//	    	printf("Error in VectorField3D::ExtractOneForProjection vel_index = %d! \n ",
//	    		    			vel_index);
//	    	printf(" ii = %d, jj = %d, kk = %d \n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk);
//	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
//	    			minValue, phi_minus_x, phi_plus_x);
//	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
//		    			minValue, phi_minus_y, phi_plus_y);
//	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
//	    		    	minValue, phi_minus_z, phi_plus_z);
//	    	printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//	    		    			vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
//	    		    			vel_minus_z, vel_plus_z);
//	    	printf(" phic = %f, phi-x = %f, phi+x = %f, phi-y = %f, phi+y = %f, phi-z = %f, phi+z = %f \n ",
//						phi_c[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii-1, minTentative.jj, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii+1, minTentative.jj, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj-1, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj+1, minTentative.kk)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj, minTentative.kk-1)],
//						phi_c[INDEX(minTentative.ii, minTentative.jj, minTentative.kk+1)]);
	    	valid[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = 0;
//	    	exit(1);
	    }
//	    if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
//////	    if(vel_index == 3){
////	    if(fabs(vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]) > 40.f){
////	    if( phi_x == 0.f && A != 0.f || phi_y == 0.f && B != 0.f || phi_z == 0.f && C != 0.f ){
//	    	printf("(ExtrapolateOneForProjection) vel_index = %d phi = %f\n",
//	    			vel_index, minValue);
//	    	printf("after interp: at i = %d, j = %d, k = %d vel = %f \n\n",
//	    			minTentative.ii, minTentative.jj, minTentative.kk,
//	    			vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)]);
//	    	printf("before interp: A = %f, B = %f, C = %f, a = %f, b = %f, c = %f \n ",
//	    			A, B, C, a, b, c);
//	    	printf(" valid = %d, valid_minus_x = %d, valid_plus_x = %d \n ",
//	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
//	    		    	valid[INDEX(minTentative.ii-1,minTentative.jj,minTentative.kk)],
//	    		    	valid[INDEX(minTentative.ii+1,minTentative.jj,minTentative.kk)]);
//	    	printf(" valid = %d, valid_minus_y = %d, valid_plus_y = %d\n ",
//	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
//	    		    	valid[INDEX(minTentative.ii,minTentative.jj-1,minTentative.kk)],
//	    		    	valid[INDEX(minTentative.ii,minTentative.jj+1,minTentative.kk)]);
//	    	printf(" valid = %d, valid_minus_z = %d, valid_plus_z = %d\n ",
//	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk)],
//	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk-1)],
//	    		    	valid[INDEX(minTentative.ii,minTentative.jj,minTentative.kk+1)]);
////	    	printf(" ii = %d, jj = %d, kk = %d \n",
////	    			minTentative.ii, minTentative.jj, minTentative.kk);
//	    	printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f, vel_minus_x = %f, vel_plus_x = %f\n ",
//	    			minValue, phi_minus_x, phi_plus_x, vel_minus_x, vel_plus_x);
//	    	printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f, vel_minus_y = %f, vel_plus_y = %f\n ",
//		    			minValue, phi_minus_y, phi_plus_y, vel_minus_y, vel_plus_y);
//	    	printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f vel_minus_z = %f, vel_plus_z = %f \n\n ",
//	    		    	minValue, phi_minus_z, phi_plus_z,vel_minus_z, vel_plus_z);
////
//////	    	exit(1);
//	    }

	} while(!heap.empty());


	while(HasNonUpdatedPoints(voxel, valid, phi_tmp, 0, heap) > 0){

	u_int N0 = heap.size();
	u_int N1 = 0, N2 = 0;
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
			KnownPoint t(minTentative.ii,minTentative.jj,minTentative.kk);
	//		knownPoint.insert( make_pair(Index(minTentative), t) );
			vector<TentativeGridPoint> neighbors;
		    minTentative.Neigbor(neighbors, DimX, DimY, DimZ);
		    for(int m=0; m<neighbors.size(); m++){
		    	if(minTentative.LeftNeighbor(neighbors[m])){
		    		phi_minus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		    		vel_minus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	//	    		if(voxel.InSolid(neighbors[m]))
	    			if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
		    			 phi_minus_x = INFINITY;
		    	}
		    	if(minTentative.RightNeighbor(neighbors[m])){
		    		phi_plus_x = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		    		vel_plus_x = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	//	    		if(voxel.InSolid(neighbors[m]))
	//	    			 phi_plus_x = INFINITY;
	//	    		if(vel_index == 1 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
	//	    				&& voxel.InSolid(neighbors[m].ii+1,neighbors[m].jj,neighbors[m].kk))
		    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
		    			 phi_plus_x = INFINITY;
		    	}
		    	if(minTentative.BackNeighbor(neighbors[m])){
		    		phi_minus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		    		vel_minus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	//	    		if(voxel.InSolid(neighbors[m]))
		    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
		    			phi_minus_y = INFINITY;
		    	}
		    	if(minTentative.FrontNeighbor(neighbors[m])){
		    		phi_plus_y = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		    		vel_plus_y = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	//	    		if(voxel.InSolid(neighbors[m]))
	//	    			 phi_plus_y = INFINITY;
	//	    		if(vel_index == 2 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
	//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj+1,neighbors[m].kk))
		    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
		    			 phi_plus_y = INFINITY;
		    	}
		    	if(minTentative.BottomNeighbor(neighbors[m])){
		    		phi_minus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		    		vel_minus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	//	    		if(voxel.InSolid(neighbors[m]))
		    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
		    			  phi_minus_z = INFINITY;
		    	}
		    	if(minTentative.TopNeighbor(neighbors[m])){
		    		phi_plus_z = phi_tmp[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
		    		vel_plus_z = vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	//	    		if(voxel.InSolid(neighbors[m]))
	//	    			  phi_plus_z = INFINITY;
	//	    		if(vel_index == 3 && !voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)
	//	    				&& voxel.InSolid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk+1))
		    		if(!valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)])
		    			 phi_plus_z = INFINITY;
		    	}
	//	    	if(minTentative.ii == I & minTentative.jj == J && minTentative.kk == K){
	//	    	    if(voxel.InSolid(neighbors[m]))
	//	    	    	printf("neighbor at (%d, %d, %d) in solid \n",
	//	    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
	//	    	    else
	//	    	    	printf("neighbor at (%d, %d, %d) not in solid \n",
	//    				 neighbors[m].ii,neighbors[m].jj,neighbors[m].kk);
	//	    	}
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
		    	vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = (a+b+c)/(A+B+C);
		    	valid[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = 1;
		    	++N1;
		    }
		    else{
				float vel_avg = 0.f;
				int NV = 0;
	 			for(int m=0; m<neighbors.size(); m++){
	 				if(valid[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)]){
	 					++NV;
	 					vel_avg += vel[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
	 				}
	 			}
	 			if(NV == 0){
//	 				printf("serious problem!\n");
//	 				printf("Error in VectorField3D::ExtractOneForProjection vel_index = %d! \n ",
//										vel_index);
//					printf(" ii = %d, jj = %d, kk = %d \n",
//							minTentative.ii, minTentative.jj, minTentative.kk);
//					printf(" minValue = %f, phi_minus_x = %f, phi_plus_x = %f\n ",
//							minValue, phi_minus_x, phi_plus_x);
//					printf(" minValue = %f, phi_minus_y = %f, phi_plus_y = %f\n ",
//								minValue, phi_minus_y, phi_plus_y);
//					printf(" minValue = %f, phi_minus_z = %f, phi_plus_z = %f\n ",
//								minValue, phi_minus_z, phi_plus_z);
//					printf("vel = %f, vel-x = %f, vel+x = %f ,vel-y = %f, vel+y = %f, vel-z = %f, vel+z = %f \n ",
//							vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)],
//										vel_minus_x, vel_plus_x, vel_minus_y, vel_plus_y,
//										vel_minus_z, vel_plus_z);
	 			}
	 			else{
	 				vel[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = vel_avg / NV;
	 				valid[INDEX(minTentative.ii, minTentative.jj, minTentative.kk)] = 1;
	 				++N2;
	 			}
		    }
		}while(!heap.empty());

		printf("\nthere are %u points reupdated, %u through usual manner, %u through one-ring averaging \n\n", N0, N1, N2);
	}



	if(valid) delete [] valid;

//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && phi_tmp[INDEX(i,j,k)] > 0.f){
//			foundpos = knownPoint.find(INDEX(i,j,k));
//			if( foundpos == knownPoint.end() ){
//				vel[INDEX(i,j,k)] = 0.f;
////				printf("Reset point (%d,%d,%d) phi(%f) vel(%d) to zero\n",
////						i, j, k, phi_tmp[INDEX(i,j,k)], vel_index);
//			}
//		}
//	END_FOR
//	knownPoint.clear();
//	SetSolidBoundaryNoExtrapolate(voxel, vel_index, vel);
}


#ifdef SPMD
u_int VectorField3D::HasNonUpdatedPoints(const Voxel &voxel, char *valid,
				 float *phi_tmp, int vel_index,
#ifdef WIN32
							hash_map<u_int, TentativeGridPoint> &heap) const{
#else
							unordered_map<u_int, TentativeGridPoint> &heap) const{
#endif
	u_int N = 0;
	if( vel_index == 0 ){
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && valid[INDEX(i,j,k)] == 0){
				KnownPoint tmp(i,j,k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(make_pair(INDEX(i,j,k),tp));
				++N;
//			printf(" ii = %d, jj = %d, kk = %d  phic = %f\n",
//					i, j, k, phi_c[INDEX(i,j,k)]);

			}
		END_FOR
	}
	else{
		FOR_EACH_CELL
			if( (vel_index == 1 && phi_u_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
				(vel_index == 2 && phi_v_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
				(vel_index == 3 && phi_w_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ){
				KnownPoint tmp(i,j,k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(make_pair(INDEX(i,j,k),tp));
				++N;
//			printf(" ii = %d, jj = %d, kk = %d  phic = %f\n",
//					i, j, k, phi_c[INDEX(i,j,k)]);

			}
		END_FOR
	}
	return N;
}
#endif
u_int VectorField3D::HasNonUpdatedPoints(const Voxel &voxel, char *valid,
					float *phi_tmp, int vel_index,
							MinHeap<TentativeGridPoint, float> &heap) const{
	u_int N = 0;
	if( vel_index == 0 ){
		FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && valid[INDEX(i,j,k)] == 0){
				KnownPoint tmp(i,j,k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
				++N;
//			printf(" ii = %d, jj = %d, kk = %d  phic = %f\n",
//					i, j, k, phi_c[INDEX(i,j,k)]);

			}
		END_FOR
	}
	else{
		FOR_EACH_CELL
			if( (vel_index == 1 && phi_u_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
				(vel_index == 2 && phi_v_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ||
				(vel_index == 3 && phi_w_obj[INDEX(i,j,k)] < 0.f && valid[INDEX(i,j,k)] == 0 ) ){
				KnownPoint tmp(i,j,k);
				TentativeGridPoint tp(phi_tmp[INDEX(i,j,k)], tmp);
				heap.insert(tp, phi_tmp[INDEX(i,j,k)]);
				++N;
//			printf(" ii = %d, jj = %d, kk = %d  phic = %f\n",
//					i, j, k, phi_c[INDEX(i,j,k)]);

			}
		END_FOR
	}
	return N;
}

void VectorField3D::Extrapolate(const Voxel &voxel, const char *valid) {

	float delta = voxel.VoxelDelta();

	SetEqual(u, u1);
	SetEqual(v, v1);
	SetEqual(w, w1);
	SetEqual(u, u0);
	SetEqual(v, v0);
	SetEqual(w, w0);
//	ExtrapolateOneVelocityIntoObject(voxel, phi_u_obj, u1, 1, delta);
//	ExtrapolateOneVelocityIntoObject(voxel, phi_v_obj, v1, 2, delta);
//	ExtrapolateOneVelocityIntoObject(voxel, phi_w_obj, w1, 3, delta);

//	printf("A. u = %f, u1 = %f, v = %f, v1 = %f \n",
//			u[INDEX(I,J,K)], u[INDEX(I,J,K-1)],
//			v[INDEX(I,J,K)], v[INDEX(I,J,K-1)]);

//	printf("B. u = %f, u1 = %f, v = %f, v1 = %f \n",
//				u[INDEX(I,J,K)], u[INDEX(I,J,K-1)],
//				v[INDEX(I,J,K)], v[INDEX(I,J,K-1)]);
	printf("\n Start Extrapolating velocity ..... \n " );

	FOR_EACH_CELL
		if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
			 voxel.InSource(i,j,k) ){
			u0[INDEX(i,j,k)] = 0.5f * (u1[INDEX(i,j,k)] + u1[INDEX(i-1,j,k)]);
			v0[INDEX(i,j,k)] = 0.5f * (v1[INDEX(i,j,k)] + v1[INDEX(i,j-1,k)]);
			w0[INDEX(i,j,k)] = 0.5f * (w1[INDEX(i,j,k)] + w1[INDEX(i,j,k-1)]);
		}
	END_FOR
//	FOR_EACH_CELL
//	if(i == I && j == J && k == K){
//		printf("1. Extrapolate: at i = %d, j = %d, k = %d, u0 = %f, u0(i-1) = %f, u0(i+1) = %f, u0(j-1) = %f, u0(j+1) = %f,  u0(k-1) = %f, u0(k+1) = %f \n",
//		 			i,j,k, u0[INDEX(i,j,k)], u0[INDEX(i-1,j,k)],u0[INDEX(i+1,j,k)],u0[INDEX(i,j-1,k)],
//		 			u0[INDEX(i,j+1,k)],u0[INDEX(i,j,k-1)],u0[INDEX(i,j,k+1)]);
//		printf("1. Extrapolate: at i = %d, j = %d, k = %d, v0 = %f, v0(i-1) = %f, v0(i+1) = %f, v0(j-1) = %f, v0(j+1) = %f,  v0(k-1) = %f, v0(k+1) = %f \n",
//				 	i,j,k, v0[INDEX(i,j,k)], v0[INDEX(i-1,j,k)],v0[INDEX(i+1,j,k)],v0[INDEX(i,j-1,k)],
//				 	v0[INDEX(i,j+1,k)],v0[INDEX(i,j,k-1)],v0[INDEX(i,j,k+1)]);
//		printf("1. Extrapolate: at i = %d, j = %d, k = %d, w0 = %f, w0(i-1) = %f, w0(i+1) = %f, w0(j-1) = %f, w0(j+1) = %f,  w0(k-1) = %f, w0(k+1) = %f \n",
//				 i,j,k, w0[INDEX(i,j,k)], w0[INDEX(i-1,j,k)],w0[INDEX(i+1,j,k)],w0[INDEX(i,j-1,k)],
//				 w0[INDEX(i,j+1,k)],w0[INDEX(i,j,k-1)],w0[INDEX(i,j,k+1)]);
//	}
//	END_FOR

	// process u component
#ifdef SPMD
	ExtrapolateOneForProjection(voxel, phi_c, u0, temp, 1);
#else
	ExtrapolateOneForProjection(voxel, phi_c, u0, 1);
#endif

	// process v component
#ifdef SPMD
	ExtrapolateOneForProjection(voxel, phi_c, v0, temp, 2);
#else
	ExtrapolateOneForProjection(voxel, phi_c, v0, 2);
#endif

	// process w component
#ifdef SPMD
	ExtrapolateOneForProjection(voxel, phi_c, w0, temp, 3);
#else
	ExtrapolateOneForProjection(voxel, phi_c, w0, 3);
#endif
//	FOR_EACH_CELL
//	if(i == I && j == J && k == K){
//		printf("2. Extrapolate: at i = %d, j = %d, k = %d, u0 = %f, u0(i-1) = %f, u0(i+1) = %f, u0(j-1) = %f, u0(j+1) = %f,  u0(k-1) = %f, u0(k+1) = %f \n",
//					i,j,k, u0[INDEX(i,j,k)], u0[INDEX(i-1,j,k)],u0[INDEX(i+1,j,k)],u0[INDEX(i,j-1,k)],
//					u0[INDEX(i,j+1,k)],u0[INDEX(i,j,k-1)],u0[INDEX(i,j,k+1)]);
//		printf("2. Extrapolate: at i = %d, j = %d, k = %d, v0 = %f, v0(i-1) = %f, v0(i+1) = %f, v0(j-1) = %f, v0(j+1) = %f,  v0(k-1) = %f, v0(k+1) = %f \n",
//					i,j,k, v0[INDEX(i,j,k)], v0[INDEX(i-1,j,k)],v0[INDEX(i+1,j,k)],v0[INDEX(i,j-1,k)],
//					v0[INDEX(i,j+1,k)],v0[INDEX(i,j,k-1)],v0[INDEX(i,j,k+1)]);
//		printf("2. Extrapolate: at i = %d, j = %d, k = %d, w0 = %f, w0(i-1) = %f, w0(i+1) = %f, w0(j-1) = %f, w0(j+1) = %f,  w0(k-1) = %f, w0(k+1) = %f \n",
//				 i,j,k, w0[INDEX(i,j,k)], w0[INDEX(i-1,j,k)],w0[INDEX(i+1,j,k)],w0[INDEX(i,j-1,k)],
//				 w0[INDEX(i,j+1,k)],w0[INDEX(i,j,k-1)],w0[INDEX(i,j,k+1)]);
//	}
//	END_FOR

#ifdef SPMD
	ExtrapolateOneAirVelocityIntoObject(voxel, phi_c_obj, u0, temp, 1, delta);
#else
	ExtrapolateOneAirVelocityIntoObject(voxel, phi_c_obj, u0, 1, delta);
#endif

#ifdef SPMD
	ExtrapolateOneAirVelocityIntoObject(voxel, phi_c_obj, v0, temp, 2, delta);
#else
	ExtrapolateOneAirVelocityIntoObject(voxel, phi_c_obj, v0, 2, delta);
#endif

#ifdef SPMD
	ExtrapolateOneAirVelocityIntoObject(voxel, phi_c_obj, w0, temp, 3, delta);
#else
	ExtrapolateOneAirVelocityIntoObject(voxel, phi_c_obj, w0, 3, delta);
#endif

//	FOR_EACH_CELL
//	if(i == I && j == J && k == K){
//		printf("3. Extrapolate: at i = %d, j = %d, k = %d, u0 = %f, u0(i-1) = %f, u0(i+1) = %f, u0(j-1) = %f, u0(j+1) = %f,  u0(k-1) = %f, u0(k+1) = %f \n",
//					i,j,k, u0[INDEX(i,j,k)], u0[INDEX(i-1,j,k)],u0[INDEX(i+1,j,k)],u0[INDEX(i,j-1,k)],
//					u0[INDEX(i,j+1,k)],u0[INDEX(i,j,k-1)],u0[INDEX(i,j,k+1)]);
//		printf("3. Extrapolate: at i = %d, j = %d, k = %d, v0 = %f, v0(i-1) = %f, v0(i+1) = %f, v0(j-1) = %f, v0(j+1) = %f,  v0(k-1) = %f, v0(k+1) = %f \n",
//					i,j,k, v0[INDEX(i,j,k)], v0[INDEX(i-1,j,k)],v0[INDEX(i+1,j,k)],v0[INDEX(i,j-1,k)],
//					v0[INDEX(i,j+1,k)],v0[INDEX(i,j,k-1)],v0[INDEX(i,j,k+1)]);
//		printf("3. Extrapolate: at i = %d, j = %d, k = %d, w0 = %f, w0(i-1) = %f, w0(i+1) = %f, w0(j-1) = %f, w0(j+1) = %f,  w0(k-1) = %f, w0(k+1) = %f \n",
//				 i,j,k, w0[INDEX(i,j,k)], w0[INDEX(i-1,j,k)],w0[INDEX(i+1,j,k)],w0[INDEX(i,j-1,k)],
//				 w0[INDEX(i,j+1,k)],w0[INDEX(i,j,k-1)],w0[INDEX(i,j,k+1)]);
//	}
//	END_FOR


    if(valid != NULL){
	FOR_EACH_CELL
		if((voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
			if( voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k) )
				u[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
			if( voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k) )
				v[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
			if( voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1) )
				w[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
			if( voxel.InSolid(i+1,j,k) ){
				float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
				CutCellFaceArea(voxel, i, j, k, 1, coeff);
				if(!valid[INDEX(i+1,j,k)] || coeff[1] < E_EPSIL)
					u[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
			}
			if( voxel.InSolid(i,j+1,k) ){
				float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
				CutCellFaceArea(voxel, i, j, k, 3, coeff);
				if(!valid[INDEX(i,j+1,k)] || coeff[3] < E_EPSIL)
					v[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
			}
			if( voxel.InSolid(i,j,k+1) ){
				float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
				CutCellFaceArea(voxel, i, j, k, 5, coeff);
				if(!valid[INDEX(i,j,k+1)] || coeff[5] < E_EPSIL)
					w[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
			}
		}
		if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL >= EXTRAPOLATE_VEL_LIMIT){
			if(voxel.InAir(i+1,j,k) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
				u[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
			if(voxel.InAir(i,j+1,k) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
				v[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
			if(voxel.InAir(i,j,k+1) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
				w[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
		}
		TentativeGridPoint tp(0.f, i, j, k);
		if(voxel.InSolid(i,j,k) && voxel.CloseToAir(tp)){
			if( (voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k)) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
				float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
				CutCellFaceArea(voxel, i, j, k, 1, coeff);
				if(!valid[INDEX(i,j,k)] || coeff[1] < E_EPSIL)
					u[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
			}
			if( (voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k)) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
				float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
				CutCellFaceArea(voxel, i, j, k, 3, coeff);
				if(!valid[INDEX(i,j,k)] || coeff[3] < E_EPSIL)
					v[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
			}
			if( (voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1)) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
				float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
				CutCellFaceArea(voxel, i, j, k, 5, coeff);
				if(!valid[INDEX(i,j,k)] || coeff[5] < E_EPSIL)
					w[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
			}
		}
	END_FOR
    }
    else{
    	FOR_EACH_CELL
			if((voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT){
				if( voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k) ||
					(phi_u_obj[INDEX(i,j,k)] < 0.f && voxel.InSolid(i+1,j,k)) )
					u[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
				if( voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k) ||
					(phi_v_obj[INDEX(i,j,k)] < 0.f && voxel.InSolid(i,j+1,k)) )
					v[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
				if( voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1) ||
					(phi_w_obj[INDEX(i,j,k)] < 0.f && voxel.InSolid(i,j,k+1)) )
					w[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
			}
			if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL >= EXTRAPOLATE_VEL_LIMIT){
				if(voxel.InAir(i+1,j,k) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
					u[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
				if(voxel.InAir(i,j+1,k) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
					v[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
				if(voxel.InAir(i,j,k+1) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
					w[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
			}
			TentativeGridPoint tp(0.f, i, j, k);
			if(voxel.InSolid(i,j,k) && voxel.CloseToAir(tp)){
				if( phi_u_obj[INDEX(i,j,k)] < 0.f &&
					(voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k)) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
					u[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
				if( phi_v_obj[INDEX(i,j,k)] < 0.f &&
					(voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k)) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
					v[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
				if( phi_w_obj[INDEX(i,j,k)] < 0.f &&
					(voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1)) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT)
					w[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
			}
		END_FOR
    }

	printf(" End Extrapolating velocity \n\n " );

//	printf("C. u = %f, u1 = %f, v = %f, v1 = %f \n",
//				u[INDEX(I,J,K)], u[INDEX(I,J,K-1)],
//				v[INDEX(I,J,K)], v[INDEX(I,J,K-1)]);
//	SetSolidBoundary(voxel);
//	printf("D. u = %f, u1 = %f, v = %f, v1 = %f \n",
//				u[INDEX(I,J,K)], u[INDEX(I,J,K-1)],
//				v[INDEX(I,J,K)], v[INDEX(I,J,K-1)]);

//	ProjectionInAirCell(voxel);

}

#ifdef AIR_DIV_FREE

#define EXTRAPOLATE_CELL(i,j,k) (voxel.InSurface((i),(j),(k)) || \
	voxel.InAir((i),(j),(k)) && phi_c[INDEX((i),(j),(k))] < AIR_DIV_FREE_LIMIT)

void VectorField3D::LinearSolverAtAirCell( const Voxel &voxel,
		 						float * x, map<u_int, u_int> &a,
		 						map<u_int, float> &x0 , char *color, u_int n ){

	float delta = voxel.VoxelDelta();
	map<u_int, u_int>::iterator found;
	map<u_int, float>::iterator found1;
	PCGSolver<double> solver;
	int iters = 200;
	if(residual_projection_air > 1.e-2)
		iters = 600;
	else if(residual_projection_air > 1.e-3)
		iters = 500;
	else if(residual_projection_air > 1.e-4)
		iters = 400;
	else if(residual_projection_air > 1.e-6)
		iters = 300;
	solver.set_solver_parameters(1e-13, iters, 0.97, 0.25);
	double residual;
	int iterations;
	SparseMatrixd A(n);
	vector<double> rhs;
	vector<double> p;
	for(u_int i=0; i<n; ++i)
		p.push_back((double)x[i]);
	u_int diag;
	FOR_EACH_CELL
		 u_int index = INDEX(i,j,k);
		 if(EXTRAPOLATE_CELL(i,j,k) && color[INDEX(i,j,k)] != 0){
			 found = a.find(index);
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back((double)found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 int solidNeighbors = 0;
			 for(u_int m=0; m<neighbors.size(); m++){
//				 if( voxel.InSolid(neighbors[m]) || voxel.InSource(neighbors[m]) ||
				 if( voxel.InSource(neighbors[m]) ||
				    (voxel.InLiquid(neighbors[m]) && !voxel.InSurface(neighbors[m])) ){
					 ++solidNeighbors;
				 }
				 else if(EXTRAPOLATE_CELL(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk)
						 && color[INDEX(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk)] != 0){
					 found = a.find(INDEX(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk));
					 A.set_element(diag,found->second,-1.);
				 }
			 }
//			 A.set_element(index,index,(float)(neighbors.size()-solidNeighbors)+total_theta);
			 A.set_element(diag,diag,(double)(neighbors.size()-solidNeighbors));
			 if(solidNeighbors == 6){
				 printf("at (%d, %d, %d) solid neighbors = %d, color = %d \n",
					 i,j,k, solidNeighbors, color[INDEX(i,j,k)]);
				 exit(1);
			 }
		 }
	END_FOR
	//A.write_matlab(cout,"A");
	solver.solve(A, rhs, p, residual, iterations);
	residual_projection_air = residual;
	for(u_int i=0; i<n; ++i)
			x[i] = (float)p[i];
	printf("\n(LinearSolverAtAirCell) PCG iterations = %d, and residual = %16.13f \n\n ",
			iterations, residual);
//	if(residual > 1.f)
//		exit(1);

}

//implementation of the flood fill algorithm
// here we used a queue-based data structure, since a recursive approach,
// although simple to implement, will most likely overflow the stack

void VectorField3D::FloodFill(const Voxel &voxel, char *color,
							   int ii, int jj, int kk,
		                       float delta) const {
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
				if(color[INDEX(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk)] == 0){
					color[INDEX(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk)] = 1;
					TentativeGridPoint tmp(0.f, neighbors[m].ii, neighbors[m].jj, neighbors[m].kk);
					Q.push(tmp);
				}
			}
			Q.pop();
//			printf("there are %lu iterms in queues \n", Q.size());
		}while(!Q.empty());
	}
}


void VectorField3D::SetBoundaryForAirCell(const Voxel &voxel, float *x){
	float delta = voxel.VoxelDelta();
	FOR_EACH_CELL
	  if( voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= AIR_DIV_FREE_LIMIT )
		  x[INDEX(i,j,k)] = 0.f;
//	  if(voxel.InSolid(i,j,k))
//		  x[INDEX(i,j,k)] = 0.f;
	END_FOR
}

// here we project extraploated velocity to be divergence-free.
//
// since we enforce Neumann boundary conditions at the face between
// an air(surface) cell and a liquid cell, large jumps in pressure may occur.
// this is because that pressure must force the air to take on the velocities
// of the liquid to make mass conservation.

void VectorField3D::ProjectionInAirCell(const Voxel &voxel){

	printf("Start ProjectionInAirCell ... \n");
	float delta = voxel.VoxelDelta();
	float residual;

	char *color = new char[DimX*DimY*DimZ];
	SetZero(color);
	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k))
			color[INDEX(i,j,k)] = 2;
		if(voxel.InSource(i,j,k))
			color[INDEX(i,j,k)] = 2;
		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k))
			color[INDEX(i,j,k)] = 2;
	END_FOR
	TentativeGridPoint floodSeed(0.f, 20, 20, 140);
	FloodFill(voxel, color, floodSeed.ii, floodSeed.jj, floodSeed.kk, delta);

	u_int air_pockets = 0, air_points = 0, others = 0;
	FOR_EACH_CELL
		if(color[INDEX(i,j,k)] == 0){
			++air_pockets;
//			printf("at (%d, %d, %d) grid point is in air pocket\n", i,j,k);
		}
		if(color[INDEX(i,j,k)] == 1)
				++air_points;
		if(color[INDEX(i,j,k)] == 2)
				++others;
	END_FOR
	printf("\nthere are %u points in air pocket and %u points in air \n\n", air_pockets, air_points);
	if(air_pockets + air_points + others != DimX*DimY*DimZ){
		printf("wrong! %u points in air pocket and %u points in air and %u others \n",air_pockets, air_points, others);
		exit(1);
	}

	map<u_int, float> rhs;
	map<u_int, u_int> A;
//	map<u_int, TentativeGridPoint> B;
	map<u_int, u_int>::iterator found;

	SetZero(u0);
	SetZero(v0);



	u_int N = 0;
//	CheckVelocity();

	float max_rhs = 0.f;
	int max_rhs_i = 0, max_rhs_j = 0,max_rhs_k = 0;
	FOR_EACH_CELL
	   //if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
	   if( EXTRAPOLATE_CELL(i,j,k) && color[INDEX(i,j,k)] != 0 ){
		   u0[INDEX(i,j,k)] = -(u1[INDEX(i,j,k)]-u1[INDEX(i-1,j,k)]+
				   				v1[INDEX(i,j,k)]-v1[INDEX(i,j-1,k)]+
				   				w1[INDEX(i,j,k)]-w1[INDEX(i,j,k-1)])*delta/dt;
		   if(fabs(max_rhs) < fabs(u0[INDEX(i,j,k)])){
			   max_rhs = u0[INDEX(i,j,k)];
			   max_rhs_i = i; max_rhs_j = j; max_rhs_k = k;
		   }
//		   TentativeGridPoint tp(0.f, i, j, k);
		   rhs.insert( make_pair(INDEX(i,j,k), u0[INDEX(i,j,k)]) );
		   A.insert( make_pair(INDEX(i,j,k), N) );
//		   B.insert( make_pair(INDEX(i,j,k), tp) );
		   ++N;
	   }
	   else
		   u0[INDEX(i,j,k)] = 0.f;
//	   printf("u0 = %f at (%d, %d, %d) \n", u0[INDEX(i,j,k)], i,j,k );
	   if(i == I && j == J && k == K){
//		if(isnan(u0[INDEX(i,j,k)]) || isinf(u0[INDEX(i,j,k)])){
//		   if(voxel.InSurface(i,j,k)) printf("SURFACE CELL\n");
//		   else printf("NOT SURFACE CELL\n");
//		   if(voxel.InAir(i,j,k)) printf(" AIR CELL\n");
//		   	else printf("NOT AIR CELL\n");
//			printf("NOT AIR CELL\n");
		   printf("Before Projection: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
			 			I,J,K, u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
			 			u1[INDEX(i-1,j,k)],v1[INDEX(i,j-1,k)],w1[INDEX(i,j,k-1)]);
		   printf("Before Projection: u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
		   			 			u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
		   			 			u1[INDEX(i+1,j,k)],v1[INDEX(i,j+1,k)],w1[INDEX(i,j,k+1)]);
		   printf("Before Projection: u0 = %f, "
				   		"phi = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
					 		u0[INDEX(i,j,k)], phi_c[INDEX(i,j,k)],
					 		phi_c[INDEX(i,j,k-1)], phi_c[INDEX(i,j,k+1)] );
		   float div = u1[INDEX(i,j,k)]-u1[INDEX(i-1,j,k)]
		   	                      +v1[INDEX(i,j,k)]-v1[INDEX(i,j-1,k)]
		   	                         +w1[INDEX(i,j,k)]-w1[INDEX(i,j,k-1)];
		   printf("divergence = %16.13f \n", div);
//		   printf("Before Projection: "
//   				   		"phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f \n",
//   					 		phi_c[INDEX(i-1,j,k)], phi_c[INDEX(i+1,j,k)],
//   					 		phi_c[INDEX(i,j-1,k)], phi_c[INDEX(i,j+1,k)] );
		}
	   v0[INDEX(i,j,k)] = 0.f;
	END_FOR

	printf("\nmaximum rhs (%f) occurs at (%d, %d, %d) \n", max_rhs,
			max_rhs_i,max_rhs_j,max_rhs_k);
	int i = max_rhs_i, j = max_rhs_j, k = max_rhs_k;
	printf(" at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
				i,j,k, u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
				u1[INDEX(i-1,j,k)],v1[INDEX(i,j-1,k)],w1[INDEX(i,j,k-1)]);
   printf("u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
						u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
						u1[INDEX(i+1,j,k)],v1[INDEX(i,j+1,k)],w1[INDEX(i,j,k+1)]);
   printf("phi = %f, phi[k-1]= %f, phi[k+1] = %f phi[i-1]= %f, phi[i+1] = %f phi[j-1]= %f, phi[j+1] = %f \n\n",
					phi_c[INDEX(i,j,k)], phi_c[INDEX(i,j,k-1)], phi_c[INDEX(i,j,k+1)],
					phi_c[INDEX(i-1,j,k)], phi_c[INDEX(i+1,j,k)],
				    phi_c[INDEX(i,j-1,k)], phi_c[INDEX(i,j+1,k)] );
   printf("phiu = %f, phiu[k-1]= %f, phiu[k+1] = %f phiu[i-1]= %f, phiu[i+1] = %f phiu[j-1]= %f, phiu[j+1] = %f \n\n",
   					phi_u[INDEX(i,j,k)], phi_u[INDEX(i,j,k-1)], phi_u[INDEX(i,j,k+1)],
   					phi_u[INDEX(i-1,j,k)], phi_u[INDEX(i+1,j,k)],
   				    phi_u[INDEX(i,j-1,k)], phi_u[INDEX(i,j+1,k)] );
   printf("phiv = %f, phiv[k-1]= %f, phiv[k+1] = %f phiv[i-1]= %f, phiv[i+1] = %f phiv[j-1]= %f, phiv[j+1] = %f \n\n",
					phi_v[INDEX(i,j,k)], phi_v[INDEX(i,j,k-1)], phi_v[INDEX(i,j,k+1)],
					phi_v[INDEX(i-1,j,k)], phi_v[INDEX(i+1,j,k)],
					phi_v[INDEX(i,j-1,k)], phi_v[INDEX(i,j+1,k)] );
   printf("phiw = %f, phiw[k-1]= %f, phiw[k+1] = %f phiw[i-1]= %f, phiw[i+1] = %f phiw[j-1]= %f, phiw[j+1] = %f \n\n",
				phi_w[INDEX(i,j,k)], phi_w[INDEX(i,j,k-1)], phi_w[INDEX(i,j,k+1)],
				phi_w[INDEX(i-1,j,k)], phi_w[INDEX(i+1,j,k)],
				phi_w[INDEX(i,j-1,k)], phi_w[INDEX(i,j+1,k)] );
   i = max_rhs_i-1; j = max_rhs_j; k = max_rhs_k;
  	printf(" at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
			i,j,k, u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
			u1[INDEX(i-1,j,k)],v1[INDEX(i,j-1,k)],w1[INDEX(i,j,k-1)]);
  printf("u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
					u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
					u1[INDEX(i+1,j,k)],v1[INDEX(i,j+1,k)],w1[INDEX(i,j,k+1)]);
  printf("phiu = %f, phiu[k-1]= %f, phiu[k+1] = %f phiu[i-1]= %f, phiu[i+1] = %f phiu[j-1]= %f, phiu[j+1] = %f \n\n",
					phi_u[INDEX(i,j,k)], phi_u[INDEX(i,j,k-1)], phi_u[INDEX(i,j,k+1)],
					phi_u[INDEX(i-1,j,k)], phi_u[INDEX(i+1,j,k)],
					phi_u[INDEX(i,j-1,k)], phi_u[INDEX(i,j+1,k)] );
   i = max_rhs_i; j = max_rhs_j-1; k = max_rhs_k;
  	printf(" at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
  				i,j,k, u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
  				u1[INDEX(i-1,j,k)],v1[INDEX(i,j-1,k)],w1[INDEX(i,j,k-1)]);
     printf("u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
  						u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
  						u1[INDEX(i+1,j,k)],v1[INDEX(i,j+1,k)],w1[INDEX(i,j,k+1)]);
     printf("phiv = %f, phiv[k-1]= %f, phiv[k+1] = %f phiv[i-1]= %f, phiv[i+1] = %f phiv[j-1]= %f, phiv[j+1] = %f \n\n",
  					phi_v[INDEX(i,j,k)], phi_v[INDEX(i,j,k-1)], phi_v[INDEX(i,j,k+1)],
  					phi_v[INDEX(i-1,j,k)], phi_v[INDEX(i+1,j,k)],
  					phi_v[INDEX(i,j-1,k)], phi_v[INDEX(i,j+1,k)] );
    i = max_rhs_i; j = max_rhs_j; k = max_rhs_k-1;
	printf(" at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
				i,j,k, u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
				u1[INDEX(i-1,j,k)],v1[INDEX(i,j-1,k)],w1[INDEX(i,j,k-1)]);
	printf("u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
						u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)],
						u1[INDEX(i+1,j,k)],v1[INDEX(i,j+1,k)],w1[INDEX(i,j,k+1)]);
	printf("phiw = %f, phiw[k-1]= %f, phiw[k+1] = %f phiw[i-1]= %f, phiw[i+1] = %f phiw[j-1]= %f, phiw[j+1] = %f \n\n",
				phi_w[INDEX(i,j,k)], phi_w[INDEX(i,j,k-1)], phi_w[INDEX(i,j,k+1)],
				phi_w[INDEX(i-1,j,k)], phi_w[INDEX(i+1,j,k)],
				phi_w[INDEX(i,j-1,k)], phi_w[INDEX(i,j+1,k)] );


//	EliminateAirBubble();
//	printf("Start EliminateAirBubble() ...\n");
//	char *visitedpoint = new char[DimX*DimY*DimZ];
//	FOR_EACH_CELL
//	   //if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
//		stack<KnownPoint> visited;
//		memset(visitedpoint, 0, DimX*DimY*DimZ);
//	    if( voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] < 3*delta ){
////		   printf("(%d, %d, %d)\n",i,j,k);
//		   if(!ConnectToAmbientAir(voxel,visited,visitedpoint,i,j,k,delta)){
////			   printf("trapped air cell found at (%d,%d,%d) stacksize = %d \n",
////					   i,j,k, visited.size());
//			   u0[INDEX(i,j,k)] = 0.f;
//		   }
//	   }
//
//	END_FOR
//
//	printf("End EliminateAirBubble()...\n");

	// set velocity divergence boundary condition here
//	SetBoundary(voxel, 0, u0);
//	SetBoundary(voxel, 0, v0);
//	FOR_EACH_CELL
//		if(isnan(u0[INDEX(i,j,k)])){
//			printf("not a number u0 = %f at (%d, %d, %d) \n", u0[INDEX(i,j,k)], i,j,k );
//			exit(1);
//		}
//		if(isinf(u0[INDEX(i,j,k)])){
//			printf("is infinity u0 = %f at (%d, %d, %d) \n", u0[INDEX(i,j,k)], i,j,k );
//			exit(1);
//		}
//	END_FOR
	printf("Before ProjectionInAirCell: there are %u extrapolated air points \n",N);

	float *p = new float[N];
//	memset(p, 0, N);
	SetZero(p, N);
//	for(int i=0;i<N;++i)
//		printf(" p[%d] = %f \n", i, p[i]);
	printf("Start linear solver ...\n");
	LinearSolverAtAirCell(voxel, p, A, rhs, color, N);
	printf("End linear solver ...\n");
	found = A.begin();
	for(; found != A.end(); ++found){
		v0[found->first] = p[found->second];
		if(phi_c[found->first] <= 0.f){
			FOR_EACH_CELL
			  if(INDEX(i,j,k)==found->first)
	                     printf("at i = %d, j = %d, k = %d \n", i,j,k);
			END_FOR
			printf(" wrong! v0 = %f \n",v0[found->first]);
			printf(" wrong! phi_c = %f \n",phi_c[found->first]);
			exit(1);
		}
	}
	delete [] p;



//	FOR_EACH_CELL
//		if(isnan(v0[INDEX(i,j,k)])){
//			printf("not a number v0 = %f at (%d, %d, %d) \n", v0[INDEX(i,j,k)], i,j,k );
//			exit(1);
//		}
//		if(isinf(v0[INDEX(i,j,k)])){
//			printf("is infinity v0 = %f at (%d, %d, %d) \n", v0[INDEX(i,j,k)], i,j,k );
//			exit(1);
//		}
//
//	END_FOR


	// set pressure boundary condition here
//	FOR_EACH_CELL
//	if(i == I & j == J && k == K)
//		printf("before %f, %f, %f, %f\n",v0[INDEX(i-1,j,k)], v0[INDEX(i,j-1,k)],
//						 v0[INDEX(i,j,k-1)], v0[INDEX(i,j,k)]);
//	 END_FOR

	SetBoundaryForAirCell(voxel, v0);
	printf("Air Cell Boundary Conditions Set ...\n");
//	FOR_EACH_CELL
//	if(i == I & j == J && k == K)
//		printf("%f, %f, %f, %f\n",v0[INDEX(i-1,j,k)], v0[INDEX(i,j-1,k)],
//						 v0[INDEX(i,j,k-1)], v0[INDEX(i,j,k)]);
//	 END_FOR


	FOR_EACH_CELL
//		stack<KnownPoint> visited;
//		memset(visitedpoint, 0, DimX*DimY*DimZ);
		if(EXTRAPOLATE_CELL(i,j,k) && color[INDEX(i,j,k)] != 0 ){
				if( !(voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
//				    && !voxel.InSolid(i+1,j,k) &&  !voxel.InSource(i+1,j,k)){
					 &&  !voxel.InSource(i+1,j,k)){
					u1[INDEX(i,j,k)] -= dt*(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)])/delta;
//					u[INDEX(i,j,k)] = u1[INDEX(i,j,k)];
				}
				if( !(voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
//					&& !voxel.InSolid(i,j+1,k) &&  !voxel.InSource(i,j+1,k)){
					 &&  !voxel.InSource(i,j+1,k)){
					v1[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)])/delta;
//					v[INDEX(i,j,k)] = v1[INDEX(i,j,k)];
				}
				if( !(voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
//					&& !voxel.InSolid(i,j,k+1) &&  !voxel.InSource(i,j,k+1)){
					 &&  !voxel.InSource(i,j,k+1)){
					w1[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)])/delta;
//					w[INDEX(i,j,k)] = w1[INDEX(i,j,k)];
				}

		}
		else{
			if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= AIR_DIV_FREE_LIMIT && color[INDEX(i,j,k)] != 0){
				if(phi_u[INDEX(i,j,k)] > 0.f && i < DimX-1 && EXTRAPOLATE_CELL(i+1,j,k) && color[INDEX(i+1,j,k)] != 0){
					u1[INDEX(i,j,k)] -= dt*(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)])/delta;
//					u[INDEX(i,j,k)] = u1[INDEX(i,j,k)];
				}
				if(phi_v[INDEX(i,j,k)] > 0.f && j < DimY-1 && EXTRAPOLATE_CELL(i,j+1,k) && color[INDEX(i,j+1,k)] != 0){
					v1[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)])/delta;
//					v[INDEX(i,j,k)] = v1[INDEX(i,j,k)];
				}
				if(phi_w[INDEX(i,j,k)] > 0.f && k < DimZ-1 && EXTRAPOLATE_CELL(i,j,k+1) && color[INDEX(i,j,k+1)] != 0){
					w1[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)])/delta;
//					w[INDEX(i,j,k)] = w1[INDEX(i,j,k)];
				}
			}
			if(voxel.InSolid(i,j,k)){
				if(i < DimX-1 && EXTRAPOLATE_CELL(i+1,j,k) && color[INDEX(i+1,j,k)] != 0){
					u1[INDEX(i,j,k)] -= dt*(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)])/delta;
					//u[INDEX(i,j,k)] = u1[INDEX(i,j,k)];
				}
				if(j < DimY-1 && EXTRAPOLATE_CELL(i,j+1,k) && color[INDEX(i,j+1,k)] != 0){
					v1[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)])/delta;
					//v[INDEX(i,j,k)] = v1[INDEX(i,j,k)];
				}
				if(k < DimZ-1 && EXTRAPOLATE_CELL(i,j,k+1) && color[INDEX(i,j,k+1)] != 0){
					w1[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)])/delta;
					//w[INDEX(i,j,k)] = w1[INDEX(i,j,k)];
				}
			}
		}

		if(i == I && j == J && k == K){
		float div = u1[INDEX(i,j,k)]-u1[INDEX(i-1,j,k)]
	                      +v1[INDEX(i,j,k)]-v1[INDEX(i,j-1,k)]
	                         +w1[INDEX(i,j,k)]-w1[INDEX(i,j,k-1)];
//		if(voxel.InLiquid(i,j,k) && fabs(div) > 1.e-3){
			printf("(ProjectionInAirCell) \n");
			 printf("After Projection: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
							 			i,j,k, u1[INDEX(i,j,k)],v1[INDEX(i,j,k)],w1[INDEX(i,j,k)]);
			 printf("After Projection: at i = %d, j = %d, k = %d  u1 = %f, v1 = %f, w1 = %f \n",
			 						i,j,k, u1[INDEX(i-1,j,k)],v1[INDEX(i,j-1,k)],w1[INDEX(i,j,k-1)]);
			 printf("After Projection: at i = %d, j = %d, k = %d  phiu = %f, phiv = %f, phiw = %f \n",
			 		 						i,j,k, phi_u[INDEX(i,j,k)],phi_v[INDEX(i,j,k)],phi_w[INDEX(i,j,k)]);
			 printf("After Projection: at i = %d, j = %d, k = %d  phiu1 = %f, phiv1 = %f, phiw1 = %f \n",
			 		 				i,j,k, phi_u[INDEX(i-1,j,k)],phi_v[INDEX(i,j-1,k)],phi_w[INDEX(i,j,k-1)]);
			 printf("%f, %f, %f, %f\n",v0[INDEX(i+1,j,k)], v0[INDEX(i,j+1,k)],
					 v0[INDEX(i,j,k+1)], v0[INDEX(i,j,k)]);
			 printf("%f, %f, %f, %f\n",v0[INDEX(i-1,j,k)], v0[INDEX(i,j-1,k)],
			 				 v0[INDEX(i,j,k-1)], v0[INDEX(i,j,k)]);
			 printf("%d, %d, %d, %d\n",color[INDEX(i+1,j,k)], color[INDEX(i,j+1,k)],
			 					 color[INDEX(i,j,k+1)], color[INDEX(i,j,k)]);
			 printf("%d, %d, %d, %d\n",color[INDEX(i-1,j,k)], color[INDEX(i,j-1,k)],
			 	 				 color[INDEX(i,j,k-1)], color[INDEX(i,j,k)]);
			 printf("divergence = %16.13f \n", div);

//			 exit(1);
		}

	END_FOR


	float max_div = 0.f;
	int i_div = 0, j_div = 0, k_div = 0;
	FOR_EACH_CELL
		if(EXTRAPOLATE_CELL(i,j,k) && color[INDEX(i,j,k)] != 0){
//			if(voxel.InSource(i,j,k)){
//				printf("at C (%d, %d, %d) point mistakenly as in liquid \n", i,j,k);
//				exit(1);
//			}
				float div = u1[INDEX(i,j,k)]-u1[INDEX(i-1,j,k)]
							 +v1[INDEX(i,j,k)]-v1[INDEX(i,j-1,k)]
							 +w1[INDEX(i,j,k)]-w1[INDEX(i,j,k-1)];
				if(max_div < fabs(div)){
					max_div = fabs(div);
					i_div = i;
					j_div = j;
					k_div = k;
				}
//				if(fabs(div) > 1.e-3){
//					printf("(ProjectionInAirCell): (%d,%d,%d) phi = %f, u = %f, v = %f, w = %f, div = %f\n",
//								i, j, k, phi_c[INDEX(i,j,k)],
//						u1[INDEX(i,j,k)], v1[INDEX(i,j,k)], w1[INDEX(i,j,k)], div);
////					exit(1);
//				}
		}
	END_FOR
	printf("maximum divergence occurs at (%d, %d, %d) = %f\n", i_div, j_div, k_div, max_div);
	float max_u = 0.f, max_v = 0.f, max_w = 0.f;
	 int max_w_index = 0, w_ii = 0, w_jj = 0, w_kk = 0;
	 int u_ii = 0, u_jj = 0, u_kk = 0;
	 int v_ii = 0, v_jj = 0, v_kk = 0;
	 int max_u_index = 0, max_v_index = 0;
	 FOR_EACH_CELL
	 //for(int n=0; n<DimX*DimY*DimZ; ++n){
	   u_int n = INDEX(i,j,k);
	 	if( (EXTRAPOLATE_CELL(i,j,k) && color[INDEX(i,j,k)] != 0 && max_u < fabs(u1[n]) ) ||
	 		(voxel.InSolid(i,j,k) && i < DimX-1 && EXTRAPOLATE_CELL(i+1,j,k) &&
	 		 color[INDEX(i+1,j,k)] != 0 && max_u < fabs(u1[n]) ) ||
	 	   (voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= AIR_DIV_FREE_LIMIT && color[INDEX(i,j,k)] != 0 &&
	 	   EXTRAPOLATE_CELL(i+1,j,k) && color[INDEX(i+1,j,k)] != 0 && max_u < fabs(u1[n])) ){
	 		max_u = fabs(u1[n]);
	 		max_u_index = n;
	 		 u_ii = i;
	 	     u_jj = j;
	 	     u_kk = k;
	 	}
	 	if( (EXTRAPOLATE_CELL(i,j,k) && color[INDEX(i,j,k)] != 0 && max_v < fabs(v1[n])) ||
	 		(voxel.InSolid(i,j,k) && j < DimY-1 && EXTRAPOLATE_CELL(i,j+1,k) &&
	 			color[INDEX(i,j+1,k)] != 0 && max_v < fabs(v1[n]) ) ||
 			(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= AIR_DIV_FREE_LIMIT && color[INDEX(i,j,k)] != 0 &&
	 		EXTRAPOLATE_CELL(i,j+1,k) && color[INDEX(i,j+1,k)] != 0 && max_v < fabs(v1[n])) ){
	 		 max_v = fabs(v1[n]);
	 		 max_v_index = n;
	 		 v_ii = i;
	 	     v_jj = j;
	 	     v_kk = k;
	 	}
	 	if( (EXTRAPOLATE_CELL(i,j,k) && color[INDEX(i,j,k)] != 0  && max_w < fabs(w1[n])) ||
	 		(voxel.InSolid(i,j,k) && k < DimZ-1 && EXTRAPOLATE_CELL(i,j,k+1) &&
	 				color[INDEX(i,j,k+1)] != 0 && max_w < fabs(w1[n]) ) ||
 			(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] >= AIR_DIV_FREE_LIMIT && color[INDEX(i,j,k)] != 0 &&
	 		EXTRAPOLATE_CELL(i,j,k+1) && color[INDEX(i,j,k+1)] != 0 && max_w < fabs(w1[n])) ){
	 		max_w = fabs(w1[n]);
	 	    max_w_index = n;
	 	    w_ii = i;
	 	    w_jj = j;
	 	    w_kk = k;
	 	}
	 	END_FOR
//		 printf("max_u at %d, max_v at %d ,max_w at %d \n ",
//				 max_u_index, max_v_index, max_w_index);



	 printf("\n max_u = %f, max_v = %f, max_w = %f \n",
		 u1[max_u_index], v1[max_v_index], w1[max_w_index] );

	 int type;
	 if(phi_c_obj[INDEX(u_ii, u_jj, u_kk)] >= 0.f)
		 type = 0;
	 else{
		 if(phi_c[INDEX(u_ii, u_jj, u_kk)] > 0.f)
			 type = 2;
		 else
			 type = 1;
	 }
	 printf("max_u at (%d,%d,%d) type = %d, phic = %f, phi_u = %f, phi_v = %f, phi_w = %f \n\n",
			 u_ii, u_jj, u_kk, type, phi_c[max_u_index], phi_u[max_u_index],phi_v[max_u_index],phi_w[max_u_index]);
	 printf("max_u at (%d,%d,%d) type = %d, old_u = %f, old_v = %f, old_w = %f \n\n",
	 		 u_ii, u_jj, u_kk, type, u[max_u_index], v[max_u_index], w[max_u_index]);
	 printf("max_u at (%d,%d,%d) type = %d, new_u = %f, new_v = %f, new_w = %f \n\n",
	 	 	u_ii, u_jj, u_kk, type, u1[max_u_index], v1[max_u_index], w1[max_u_index]);
	 printf("max_u at (%d,%d,%d) type = %d, v0 = %f, v0(i+1) = %f phic[i+1] = %f phicobj[i+1] = %f\n\n",
	 			 u_ii, u_jj, u_kk, type, v0[max_u_index], v0[INDEX(u_ii+1,u_jj,u_kk)],
	 			 phi_c[INDEX(u_ii+1,u_jj,u_kk)], phi_c_obj[INDEX(u_ii+1,u_jj,u_kk)]);
	 printf("max_u at (%d,%d,%d) type = %d, u0 = %f, u0[i+1] = %f u0[i-1] = %f u0[j+1] = %f u0[j-1] = %f u0[k+1] = %f, u0[k-1] = %f\n\n",
				 u_ii, u_jj, u_kk, type, u0[max_u_index], u0[INDEX(u_ii+1,u_jj,u_kk)],u0[INDEX(u_ii-1,u_jj,u_kk)],
				 u0[INDEX(u_ii,u_jj+1,u_kk)], u0[INDEX(u_ii,u_jj-1,u_kk)], u0[INDEX(u_ii,u_jj,u_kk+1)], u0[INDEX(u_ii,u_jj,u_kk-1)]);

	 if(phi_c_obj[INDEX(v_ii, v_jj, v_kk)] >= 0.f)
		 type = 0;
	 else{
		 if(phi_c[INDEX(v_ii, v_jj, v_kk)] > 0.f)
			 type = 2;
		 else
			 type = 1;
	 }
	 printf("max_v at (%d,%d,%d) type = %d, phic = %f, phi_u = %f, phi_v = %f, phi_w = %f \n\n",
			 v_ii, v_jj, v_kk, type, phi_c[max_v_index], phi_u[max_v_index],phi_v[max_v_index],phi_w[max_v_index]);
	 printf("max_v at (%d,%d,%d) type = %d, old_u = %f, old_v = %f, old_w = %f \n\n",
			 v_ii, v_jj, v_kk, type, u[max_v_index], v[max_v_index], w[max_v_index]);
	 printf("max_v at (%d,%d,%d) type = %d, new_u = %f, new_v = %f, new_w = %f \n\n",
			v_ii, v_jj, v_kk, type, u1[max_v_index], v1[max_v_index], w1[max_v_index]);
	 printf("max_v at (%d,%d,%d) type = %d, v0 = %f, v0(i+1) = %f phic[j+1] = %f, phicobj[j+1] = %f\n\n",
	 	 	v_ii, v_jj, v_kk, type, v0[max_v_index], v0[INDEX(v_ii,v_jj+1,v_kk)],
	 	 	phi_c[INDEX(v_ii,v_jj+1,v_kk)],phi_c_obj[INDEX(v_ii,v_jj+1,v_kk)] );
	 printf("max_v at (%d,%d,%d) type = %d, u0 = %f, u0[i+1] = %f u0[i-1] = %f u0[j+1] = %f u0[j-1] = %f u0[k+1] = %f, u0[k-1] = %f\n\n",
	 		 v_ii, v_jj, v_kk, type, u0[max_v_index], u0[INDEX(v_ii+1,v_jj,v_kk)],u0[INDEX(v_ii-1,v_jj,v_kk)],
	 		 u0[INDEX(v_ii,v_jj+1,v_kk)], u0[INDEX(v_ii,v_jj-1,v_kk)], u0[INDEX(v_ii,v_jj,v_kk+1)], u0[INDEX(v_ii,v_jj,v_kk-1)]);
	 if(phi_c_obj[INDEX(w_ii, w_jj, w_kk)] >= 0.f)
		 type = 0;
	 else{
		 if(phi_c[INDEX(w_ii, w_jj, w_kk)] > 0.f)
			 type = 2;
		 else
			 type = 1;
	 }
	 printf("max_w at (%d,%d,%d) type = %d, phic = %f, phi_u = %f, phi_v = %f, phi_w = %f \n\n",
			 w_ii, w_jj, w_kk, type, phi_c[max_w_index], phi_u[max_w_index],phi_v[max_w_index],phi_w[max_w_index]);
	 printf("max_w at (%d,%d,%d) type = %d, old_u = %f, old_v = %f, old_w = %f \n\n",
			 w_ii, w_jj, w_kk, type, u[max_w_index], v[max_w_index], w[max_w_index]);
	 printf("max_w at (%d,%d,%d) type = %d, new_u = %f, new_v = %f, new_w = %f \n\n",
			w_ii, w_jj, w_kk, type, u1[max_w_index], v1[max_w_index], w1[max_w_index]);
	 printf("max_w at (%d,%d,%d) type = %d, v0 = %f, v0(i+1) = %f phi_c[k+1] = %f, phicobj[k+1] = %f \n\n",
 			 w_ii, w_jj, w_kk, type, v0[max_w_index], v0[INDEX(w_ii,w_jj,w_kk+1)],
 			 phi_c[INDEX(w_ii,w_jj,w_kk+1)],phi_c_obj[INDEX(w_ii,w_jj,w_kk+1)]);
	 printf("max_w at (%d,%d,%d) type = %d, u0 = %f, u0[i+1] = %f u0[i-1] = %f u0[j+1] = %f u0[j-1] = %f u0[k+1] = %f, u0[k-1] = %f\n\n",
			 w_ii, w_jj, w_kk, type, u0[max_w_index], u0[INDEX(w_ii+1,w_jj,w_kk)],u0[INDEX(w_ii-1,w_jj,w_kk)],
			 u0[INDEX(w_ii,w_jj+1,w_kk)], u0[INDEX(w_ii,w_jj-1,w_kk)], u0[INDEX(w_ii,w_jj,w_kk+1)], u0[INDEX(w_ii,w_jj,w_kk-1)]);

//	 if(fabs(u[max_u_index]) > fabs(v[max_v_index]) &&  fabs(u[max_u_index]) > fabs(w[max_w_index])){
//		 I = u_ii; J = u_jj; K = u_kk;
//	 }
//	 if(fabs(v[max_v_index]) > fabs(u[max_u_index]) &&  fabs(v[max_v_index]) > fabs(w[max_w_index])){
// 		 I = v_ii; J = v_jj; K = v_kk;
// 	 }
//	 if(fabs(w[max_w_index]) > fabs(u[max_u_index]) &&  fabs(w[max_w_index]) > fabs(v[max_v_index])){
// 		 I = w_ii; J = w_jj; K = w_kk;
// 	 }

	 // set velocity boundary condition here

//	SetBoundary(voxel, 1, u);
//	SetBoundary(voxel, 2, v);
//	SetBoundary(voxel, 3, w);
//	SetSolidBoundary(voxel);
	delete [] color;

//	FOR_EACH_CELL
//		if(i == I & j == J && k == K){
//			printf("After Projection: at i = %d, j = %d, k = %d  u1 = %f, v1 = %f, w1 = %f \n",
//			 						I,J,K, u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)]);
//			 printf("divergence = %16.13f \n",
//				 u[INDEX(i,j,k)]-u[INDEX(i-1,j,k)]
//				                   +v[INDEX(i,j,k)]-v[INDEX(i,j-1,k)]
//				                                 +w[INDEX(i,j,k)]-w[INDEX(i,j,k-1)]);
//		}
//	END_FOR

}

#endif

void VectorField3D::LinearSolver( const Voxel &voxel,
					int b, float * x, float * x0, float a, float c )
{
	int m;
	float phi_tmp;
	for ( m=0 ; m<1 ; m++ ) {
		FOR_EACH_CELL
			if(b == 0){
				phi_tmp = phi_c[INDEX(i,j,k)];
			}
			if(b == 1){
				phi_tmp = phi_u[INDEX(i,j,k)];
			}
			if(b == 2){
				phi_tmp = phi_v[INDEX(i,j,k)];
			}
			if(b == 3){
				phi_tmp = phi_w[INDEX(i,j,k)];
			}
			if( b == 0 ){
				if(voxel.InLiquid(i,j,k))
					x[INDEX(i,j,k)] = (x0[INDEX(i,j,k)] + a*( x[INDEX(i-1,j,k)]+x[INDEX(i+1,j,k)]+
	                                x[INDEX(i,j-1,k)]+x[INDEX(i,j+1,k)]+
	                                x[INDEX(i,j,k-1)]+x[INDEX(i,j,k+1)]) ) / c;
//				if(i == I & j == J && k == K && m % 3 == 0){
//					 printf("\nAfter Projection: m = %d at i = %d, j = %d, k = %d  u = %f \n\n",
//							 		m,	I,J,K, x[INDEX(I,J,K)]);
//				}


			}
			else{
//				if(phi_tmp < 0.f)
//					x[INDEX(i,j,k)] = (x0[INDEX(i,j,k)] + a*( x[INDEX(i-1,j,k)]+x[INDEX(i+1,j,k)]+
//			                                x[INDEX(i,j-1,k)]+x[INDEX(i,j+1,k)]+
//			                                x[INDEX(i,j,k-1)]+x[INDEX(i,j,k+1)]) ) / c;
					 x[INDEX(i,j,k)] =  x0[INDEX(i,j,k)];
			}
//			if(i == I && j == J && k == K && b == 3){
//				 printf("\n(LinearSolver): at i = %d, j = %d, k = %d  vel = %f old vel = %f\n",
//						 		i, j, k, x[INDEX(i,j,k)], x0[INDEX(i,j,k)]);
//				 printf(" phi_tmp = %f, a = %f, c = %f \n", phi_tmp, a, c);
//				 printf(" x[i-1,j,k] = %f, x[i+1,j,k] = %f x[i,j-1,k] = %f "
//						 "x[i,j+1,k] = %f x[i,j,k-1] = %f  x[i,j,k+1] = %f \n\n",
//						 x[INDEX(i-1,j,k)],  x[INDEX(i+1,j,k)],  x[INDEX(i,j-1,k)],
//						 x[INDEX(i,j+1,k)],  x[INDEX(i,j,k-1)],  x[INDEX(i,j,k+1)]);
//			}
		END_FOR
//		SetBoundary(voxel, b, x);
		//set_bnd ( N, b, x );
	}
}

void VectorField3D::LinearSolver( const Voxel &voxel,
					int b, float * x, float * x0, u_int n ){

	float delta = voxel.VoxelDelta();

	PCGSolver<float> solver;
	solver.set_solver_parameters(1e-6, 100, 0.97, 0.25);
	float residual;
	int iterations;
	SparseMatrixf A(n);
	vector<float> rhs;
	vector<float> p;
	for(u_int i=0; i<n; ++i)
		p.push_back(x[i]);
	FOR_EACH_CELL
		 u_int index = INDEX(i,j,k);
		 rhs.push_back(x0[INDEX(i,j,k)]);
		 TentativeGridPoint tp(0.f,i,j,k);
		 vector<TentativeGridPoint> neighbors;
		 tp.Neigbor(neighbors, DimX, DimY, DimZ);
		 if(voxel.InLiquid(i,j,k)){
			 int solidNeighbors=0;
			 float total_theta = 0.f;
			 for(u_int m=0; m<neighbors.size(); m++){
				 if(voxel.InSolid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
					 ++solidNeighbors;
					 A.set_element(index,INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk),0.f);
				 }
				 else if(voxel.InAir(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
//					 ++solidNeighbors;
					 A.set_element(index,INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk),0.f);
					 float air_phi = phi_c[INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)];
					 float liquid_phi = phi_c[INDEX(i,j,k)];
					 if(fabs(liquid_phi) < 1.e-6*delta){
						 liquid_phi = -1.e-6*delta;
						 air_phi = delta + liquid_phi;
					 }
					 total_theta += air_phi / (-liquid_phi);
				 }
				 else{
					 A.set_element(index,INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk),-1.f);
				 }
			 }
//			 A.set_element(index,index,(float)(neighbors.size()-solidNeighbors)+total_theta);
			 A.set_element(index,index,(float)(neighbors.size()-solidNeighbors));
//			 printf("non solid neighbors= %d at(%d, %d,%d) \n",
//					 (neighbors.size()-solidNeighbors), i,j,k);
		 }
		 else{
			 for(u_int m=0; m<neighbors.size(); m++){
				 A.set_element(index,INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk),0.f);
			 }
			 A.set_element(index,index,1.f);
		 }
	END_FOR
	//A.write_matlab(cout,"A");
	solver.solve(A, rhs, p, residual, iterations);
	for(u_int i=0; i<n; ++i)
			x[i] = p[i];
	printf("\nPCG iterations = %d, and residual = %16.13f \n\n ",
			iterations, residual);


}


#define VOLUME_FRAC_THRESHOLD 0.01f

// solves the pressure projection equation using the PCG method
//
// first it sets the coefficients of the sparse matrix for pressure
// and then calls the PCG solver to get the solution (pressure)
//
// the key is that the off-diagnal coefficients of the sparse matrix must be symmetric,
// otherwise the solver won't converge

void VectorField3D::LinearSolver( const Voxel &voxel, float *x,
#ifdef WIN32
					 hash_map<u_int, u_int> &a,
#else
					 unordered_map<u_int, u_int> &a,
#endif
					  const float *x0,

					 const char *valid, u_int n ){

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float delta3 = delta * delta2;
#ifdef WIN32
	hash_map<u_int, u_int>::iterator found;
#else
	unordered_map<u_int, u_int>::iterator found;
#endif
	PCGSolver<double> solver;
	int cgIters = 300;
	if(n > 500000){
		cgIters *= int(round((float)n/500000.0f));
	}
	solver.set_solver_parameters(1e-13, cgIters, 0.97, 0.25);
	printf("pressure solver stage 2 \n");
	double residual;
	int iterations;
	SparseMatrixd A(n);
	vector<double> rhs;
	vector<double> p;
	for(u_int i=0; i<n; ++i)
		p.push_back((double)x[i]);
	u_int diag;
	printf("pressure solver stage 3 \n");
	FOR_EACH_CELL
		 u_int index = INDEX(i,j,k);
		 if( valid[index] == 1 ){
			 found = a.find(index);
			 diag = found->second;
			 rhs.push_back((double)x0[index]);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 int solidNeighbors=0;
			 float total_theta = 0.f;
//			 float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
#ifdef USE_MESH
			 Point p0 = voxel.VoxelCenterPosition(i,j,k);
			 double s, t;
#endif
			 for(u_int m=0; m<neighbors.size(); m++){
				 	 int ni = neighbors[m].ii;
			  	 	 int nj = neighbors[m].jj;
				 	 int nk = neighbors[m].kk;
				 	 u_int nindex = INDEX(ni, nj, nk);
#ifdef USE_EMSH
				 	 Point p1 = voxel.VoxelCenterPosition(ni,nj, nk);
#endif
					 if(voxel.InSolid(ni,nj,nk)){
						 if(valid[nindex]){
							 found = a.find(nindex);
							 A.set_element(diag,found->second,-1);
						 }
						 else
							 ++solidNeighbors;
					 }
#ifdef USE_MESH
					 else if(mPhysWorld->SegmentIntersectsMesh(p0, p1, &s, &t))
						 ++solidNeighbors;
#endif
					 else if(voxel.InSource(ni,nj,nk)){
						 ++solidNeighbors;
//						coeff[m] = 0.f;
					 }
					 else if(voxel.InAir(ni,nj,nk) || voxel.InSurface(ni,nj,nk)){
						 float air_phi    = phi_c[nindex];
						 float liquid_phi = phi_c[index];
						 float factor;
						 if(fabs(liquid_phi) < E_EPSIL){
//							 liquid_phi = -1.e-6*delta;
//							 air_phi = delta + liquid_phi;
							 factor = 1.e3;
						 }
						 else
							 factor = - (air_phi / liquid_phi);
//						float factor1 = air_phi/(air_phi - liquid_phi);
//						float factor2 = -liquid_phi/(air_phi - liquid_phi);
//						float factor = factor1 / factor2;
	// 					 printf("at (%d, %d, %d) phi_c = %f, total theta = %f \n", i,j,k,
	// 							 liquid_phi, factor);
	// 					 printf("neightbor at (%d, %d, %d) phi_c = %f, sum = %f \n", neighbors[m].ii,neighbors[m].jj,
	// 							 neighbors[m].kk, air_phi, factor1+factor2);
						 total_theta += factor;
//						coeff[m] = (1.f+factor) * delta3;
					 }
					 else if(voxel.InLiquid(ni,nj,nk) && !voxel.InSurface(ni,nj,nk)){
						 found = a.find(nindex);
//						 A.set_element(diag,found->second,-delta3);
						 A.set_element(diag,found->second,-1);
//						 coeff[m] = delta3;
					 }
			 }
			 if(solidNeighbors == 6){
				printf("at (%d, %d, %d) Newmann b.c. for all directions \n", i,j,k);
				exit(1);
			 }
//			 printf("at (%d, %d, %d) total theta = %f \n", i,j,k,total_theta);
//			 float diag_total = 0.f;
//             for(int m=0; m<6;++m)
//                 diag_total += coeff[m];
//             A.set_element(diag, diag, (double)diag_total);
			 A.set_element(diag,diag,total_theta+(float)(neighbors.size()-solidNeighbors));
//			 printf("non solid neighbors= %d at(%d, %d,%d) \n",
//					 (neighbors.size()-solidNeighbors), i,j,k);
		 }
		 if( valid[index] == 2 ){
			 found = a.find(index);
			 diag = found->second;
			 rhs.push_back((double)x0[index]);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 int solidNeighbors=0;
			 float total_theta = 0.f;
			 float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
			 for(u_int m=0; m<neighbors.size(); m++){
					 if(voxel.InSolid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
						 if(valid[INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)]){
							 CutCellFaceArea(voxel, i, j, k, m, coeff);
//							 if(coeff[m] < E_EPSIL){
//								 printf("(%d, %d, %d) with m = %d full solid face \n", i, j, k, m);
//								 coeff[m] = 1.f;
//							 }
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							 A.set_element(diag,found->second,-coeff[m]);

						 }
						 else
							 coeff[m] = 0.f;
					 }
					 else if(voxel.InSource(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
						 ++solidNeighbors;
						coeff[m] = 0.f;
					 }
	//				 else if(source && source->IsSourceNeighbor(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
	// 					 ++solidNeighbors;
	// 				 }
					 else if(voxel.InAir(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ||
							 voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
	 //					 ++solidNeighbors;
	// 					 A.set_element(index,INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk),0.f);
						 float air_phi = phi_c[INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)];
						 float liquid_phi = phi_c[INDEX(i,j,k)];
						 if(fabs(liquid_phi) < 1.e-6*delta){
							 liquid_phi = -1.e-6*delta;
							 air_phi = delta + liquid_phi;
						 }
	// 					 float factor = air_phi / (-liquid_phi);
						float factor1 = air_phi/(air_phi - liquid_phi);
						float factor2 = -liquid_phi/(air_phi - liquid_phi);
						float factor = factor1 / factor2;
	// 					 printf("at (%d, %d, %d) phi_c = %f, total theta = %f \n", i,j,k,
	// 							 liquid_phi, factor);
	// 					 printf("neightbor at (%d, %d, %d) phi_c = %f, sum = %f \n", neighbors[m].ii,neighbors[m].jj,
	// 							 neighbors[m].kk, air_phi, factor1+factor2);
						 total_theta += factor;
						 CutCellFaceArea(voxel, i, j, k, m, coeff);
						 if(liquid_phi < 0.f)
							 coeff[m] *= 1.f+factor;
						 else
							 coeff[m] *= 1.f;
					 }
	//				 else if(voxel.InAir(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
	//					 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
	//					 A.set_element(diag,found->second,0.f);
	//  			     }
					 else if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) &&
							!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							if(valid[INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)] == 1){
								 coeff[m] = 1.f;
							}
							else{
								 CutCellFaceArea(voxel, i, j, k, m, coeff);
							}
							A.set_element(diag,found->second,-coeff[m]);
					 }
			 }
			 if(solidNeighbors == 6){
				printf("at (%d, %d, %d) Newmann b.c. for all directions \n", i,j,k);
				exit(1);
			 }
 //			 printf("at (%d, %d, %d) total theta = %f \n", i,j,k,total_theta);
			 float diag_total = 0.f;
			  for(int m=0; m<6;++m)
				  diag_total += coeff[m];
              A.set_element(diag, diag, (double)diag_total);
//			 A.set_element(diag,diag,total_theta+(float)(neighbors.size()-solidNeighbors));
 //			 printf("non solid neighbors= %d at(%d, %d,%d) \n",
 //					 (neighbors.size()-solidNeighbors), i,j,k);
		 }


	END_FOR
	printf("pressure solver stage 4 \n");
	//A.write_matlab(cout,"A");
//	printf("\nmatrix A is ready \n\n ");
	A.check_symmetry();
	printf("pressure solver stage 5 \n");
	u_int M = 0;
	if(!A.check_diagnal(M)){
		FOR_EACH_CELL
			u_int index = INDEX(i,j,k);
			found = a.find(index);
			if(found != a.end()){
				diag = found->second;
				if( diag == M ){
					printf(" at (%d, %d, %d) matrix diagnal is too small or negative\n", i, j, k);
					exit(1);
				}
			}
		END_FOR
	}
	printf("pressure solver stage 6 \n");
	solver.solve(A, rhs, p, residual, iterations);
	printf("pressure solver stage 7 \n");
	for(u_int i=0; i<n; ++i)
			x[i] = (float)p[i];
	printf("pressure solver stage 8 \n");
	printf("\nPCG iterations = %d, and residual = %16.13f \n\n ",
			iterations, residual);


}



void VectorField3D::LinearSolverForDiffusion( const Voxel &voxel,
					 float * x, map<u_int, u_int> &a,
					 map<u_int, float> &x0 , u_int n, int vel_index ){

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float coeff = visc * dt / delta2;
	map<u_int, u_int>::iterator found;
	map<u_int, float>::iterator found1;
	PCGSolver<float> solver;
	solver.set_solver_parameters(1e-12, 100, 0.97, 0.25);
	float residual;
	int iterations;
	SparseMatrixf A(n);
	vector<float> rhs;
	vector<float> p;
	for(u_int i=0; i<n; ++i)
		p.push_back(x[i]);
	u_int diag;
	FOR_EACH_CELL
		 u_int index = INDEX(i,j,k);
		 if( voxel.InLiquid(i,j,k)  && !voxel.InSurface(i,j,k) ){
			 found = a.find(index);
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back(found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 for(u_int m=0; m<neighbors.size(); m++){
				 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) &&
					!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
					 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
					 A.set_element(diag,found->second,-1.f * coeff);
				 }
				 else{
					 if(vel_index == 1){
						 if(tp.LeftNeighbor(neighbors[m])){
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							 A.set_element(diag,found->second,-1.f * coeff);
						 }
						 else{
							 if(voxel.InLiquid(neighbors[m].ii+1,neighbors[m].jj, neighbors[m].kk) &&
							 	!voxel.InSurface(neighbors[m].ii+1,neighbors[m].jj, neighbors[m].kk) ){
								 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
								 A.set_element(diag,found->second,-1.f * coeff);
							 }
						 }
					 }
					 if(vel_index == 2){
						 if(tp.BackNeighbor(neighbors[m])){
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							 A.set_element(diag,found->second,-1.f * coeff);
						 }
						 else{
							 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj+1, neighbors[m].kk) &&
								!voxel.InSurface(neighbors[m].ii,neighbors[m].jj+1, neighbors[m].kk) ){
								 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
								 A.set_element(diag,found->second,-1.f * coeff);
							 }
						 }
					 }
					 if(vel_index == 3){
						 if(tp.BottomNeighbor(neighbors[m])){
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							 A.set_element(diag,found->second,-1.f * coeff);
						 }
						 else{
							 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk+1) &&
								!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk+1) ){
								 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
								 A.set_element(diag,found->second,-1.f * coeff);
							 }
						 }
					 }
				 }
			 }
			 A.set_element(diag,diag, 1.f + 6 * coeff);
		 }
		 if(  vel_index == 1 && i < DimX-1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k) ){
			 found = a.find(index);
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back(found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 for(u_int m=0; m<neighbors.size(); m++){
				 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) &&
					!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
					 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
					 A.set_element(diag,found->second,-1.f * coeff);
				 }
				 else{
					 if(voxel.InLiquid(neighbors[m].ii+1,neighbors[m].jj, neighbors[m].kk) &&
						!voxel.InSurface(neighbors[m].ii+1,neighbors[m].jj, neighbors[m].kk) ){
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							 A.set_element(diag,found->second,-1.f * coeff);
					 }
				 }
			 }
			 A.set_element(diag,diag, 1.f + 6 * coeff);
		 }
		 if(  vel_index == 2 && j < DimY-1 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k) ){
			 found = a.find(index);
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back(found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 for(u_int m=0; m<neighbors.size(); m++){
				 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) &&
					!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
					 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
					 A.set_element(diag,found->second,-1.f * coeff);
				 }
				 else{
					 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj+1, neighbors[m].kk) &&
						!voxel.InSurface(neighbors[m].ii,neighbors[m].jj+1, neighbors[m].kk) ){
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							 A.set_element(diag,found->second,-1.f * coeff);
					 }
				 }
			 }
			 A.set_element(diag,diag, 1.f + 6 * coeff);
		 }
		 if(  vel_index == 3 && k < DimZ-1 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1) ){
			 found = a.find(index);
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back(found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 for(u_int m=0; m<neighbors.size(); m++){
				 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) &&
					!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
					 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
					 A.set_element(diag,found->second,-1.f * coeff);
				 }
				 else{
					 if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk+1) &&
						!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk+1) ){
							 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
							 A.set_element(diag,found->second,-1.f * coeff);
					 }
				 }
			 }
			 A.set_element(diag,diag, 1.f + 6 * coeff);
		 }
	END_FOR

	solver.solve(A, rhs, p, residual, iterations);
	for(u_int i=0; i<n; ++i)
			x[i] = p[i];
	printf("\nPCG iterations = %d, and residual = %16.13f \n\n ",
			iterations, residual);


}

void VectorField3D::ExplicitDiffusion(const Voxel &voxel, float *vel, float *vel0, float *phi_tmp, int vel_index){
	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float a = visc * dt / delta2;

	FOR_EACH_CELL
		if(phi_tmp[INDEX(i,j,k)] >= 0.f){
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,vel_index);
			}
			if(vel_index == 1)
				vel0[INDEX(i,j,k)] = obj_vel.x;
			if(vel_index == 2)
				vel0[INDEX(i,j,k)] = obj_vel.y;
			if(vel_index == 3)
				vel0[INDEX(i,j,k)] = obj_vel.z;
		}
	END_FOR

	FOR_EACH_CELL
		if(phi_tmp[INDEX(i,j,k)] < 0.f){
			if( voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) ){
				vel[INDEX(i,j,k)] = vel0[INDEX(i,j,k)] + a * ( -6 * vel0[INDEX(i,j,k)] +
						vel0[INDEX(i-1,j,k)]+ vel0[INDEX(i,j-1,k)] + vel0[INDEX(i,j,k-1)] +
						vel0[INDEX(i+1,j,k)]+ vel0[INDEX(i,j+1,k)] + vel0[INDEX(i,j,k+1)] );

			}
			else{
				if( ( vel_index == 1 && i < DimX-1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k) )  ||
					( vel_index == 2 && j < DimY-1 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k) )  ||
					( vel_index == 3 && k < DimZ-1 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1) ) )
					vel[INDEX(i,j,k)] = vel0[INDEX(i,j,k)] + a * ( -6 * vel0[INDEX(i,j,k)] +
								vel0[INDEX(i-1,j,k)]+ vel0[INDEX(i,j-1,k)] + vel0[INDEX(i,j,k-1)] +
								vel0[INDEX(i+1,j,k)]+ vel0[INDEX(i,j+1,k)] + vel0[INDEX(i,j,k+1)] );
//				else
//					vel[INDEX(i,j,k)] = vel0[INDEX(i,j,k)];
			}
		}
//		else
//			vel[INDEX(i,j,k)] = vel0[INDEX(i,j,k)];
	END_FOR
}

void VectorField3D::ImplicitDiffusion(const Voxel &voxel, float *vel, float *vel0, float *phi_tmp, int vel_index){
	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float a = visc * dt / delta2;
	SetEqual(vel0, vel);

	FOR_EACH_CELL
		if(phi_tmp[INDEX(i,j,k)] >= 0.f){
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,vel_index);
			}
			if(vel_index == 1)
				vel0[INDEX(i,j,k)] = obj_vel.x;
			if(vel_index == 2)
				vel0[INDEX(i,j,k)] = obj_vel.y;
			if(vel_index == 3)
				vel0[INDEX(i,j,k)] = obj_vel.z;
		}
	END_FOR
	map<u_int, float> rhs;
	map<u_int, u_int> A;
	map<u_int, u_int>::iterator found;
	u_int N = 0;
	FOR_EACH_CELL
		if(phi_tmp[INDEX(i,j,k)] < 0.f){
			if( voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) ){
				float tmp = vel0[INDEX(i,j,k)];
				TentativeGridPoint tp(0.f,i,j,k);
				vector<TentativeGridPoint> neighbors;
				tp.Neigbor(neighbors, DimX, DimY, DimZ);
				for(u_int m=0; m<neighbors.size(); m++){
					 if(!(voxel.InLiquid(neighbors[m]) && !voxel.InSurface(neighbors[m]))){
						 if(vel_index == 1){
							 if(!(voxel.InLiquid(neighbors[m].ii+1,neighbors[m].jj,neighbors[m].kk)
							  && !voxel.InSurface(neighbors[m].ii+1, neighbors[m].jj, neighbors[m].kk)))
								 tmp += a * vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						 }
						 if(vel_index == 2){
							 if(!(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj+1,neighbors[m].kk)
							  && !voxel.InSurface(neighbors[m].ii, neighbors[m].jj+1, neighbors[m].kk)))
								 tmp += a * vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						 }
						 if(vel_index == 3){
							 if(!(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk+1)
							  && !voxel.InSurface(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk+1)))
								 tmp += a * vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						 }
					 }
				}
				rhs.insert( make_pair(INDEX(i,j,k), tmp) );
				A.insert( make_pair(INDEX(i,j,k), N) );
				++N;
			}
			else if( ( vel_index == 1 && i < DimX-1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k) )  ||
					 ( vel_index == 2 && j < DimY-1 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k) )  ||
					 ( vel_index == 3 && k < DimZ-1 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1) ) ){
				float tmp = vel0[INDEX(i,j,k)];
				TentativeGridPoint tp(0.f,i,j,k);
				vector<TentativeGridPoint> neighbors;
				tp.Neigbor(neighbors, DimX, DimY, DimZ);
				for(u_int m=0; m<neighbors.size(); m++){
					 if(!(voxel.InLiquid(neighbors[m]) && !voxel.InSurface(neighbors[m]))){
						 if(vel_index == 1){
							 if(!(voxel.InLiquid(neighbors[m].ii+1,neighbors[m].jj,neighbors[m].kk)
							  && !voxel.InSurface(neighbors[m].ii+1, neighbors[m].jj, neighbors[m].kk)))
								 tmp += a * vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						 }
						 if(vel_index == 2){
							 if(!(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj+1,neighbors[m].kk)
							  && !voxel.InSurface(neighbors[m].ii, neighbors[m].jj+1, neighbors[m].kk)))
								 tmp += a * vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						 }
						 if(vel_index == 3){
							 if(!(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk+1)
							  && !voxel.InSurface(neighbors[m].ii, neighbors[m].jj, neighbors[m].kk+1)))
								 tmp += a * vel0[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)];
						 }
					 }
				}
				rhs.insert( make_pair(INDEX(i,j,k), tmp) );
				A.insert( make_pair(INDEX(i,j,k), N) );
				++N;
			}
		}
	END_FOR

	float *p = new float[N];
	SetZero(p, N);
	LinearSolverForDiffusion(voxel, p, A, rhs, N, vel_index);
	found = A.begin();
	for(; found != A.end(); ++found){
		vel[found->first] = p[found->second];
	}
	delete [] p;
}

void VectorField3D::Diffuse(const Voxel &voxel){

	float delta = voxel.VoxelDelta();
	float a = dt*visc/(delta*delta);

//	int i,j,k;
//	FOR_EACH_CELL
//		if(i == I & j == J && k == K){
//			printf("\nbefore diffuse vel is (%f,%f,%f) \n",
//					u0[INDEX(i,j,k)],v0[INDEX(i,j,k)],w0[INDEX(i,j,k)]);
//		}
//	END_FOR
//	SWAP ( u0, u );
//	SWAP ( v0, v );
//	SWAP ( w0, w );
//	SetEqual(u, u0);
//	SetEqual(v, v0);
//	SetEqual(w, w0);
//	SetEqual(u, u1);
//	SetEqual(v, v1);
//	SetEqual(w, w1);
	float visc_cfl_dt = delta * delta /(6 * visc);
	if(visc == 0.f){
		printf("\nperform no diffusion b.c. visc = %f \n\n", visc);
//		SWAP ( u0, u );
//		SWAP ( v0, v );
//		SWAP ( w0, w );
	}
	else if(dt < visc_cfl_dt){
		printf("\nperform explicit diffusion with viscosity = %f \n\n", visc);
		ExplicitDiffusion(voxel, u1, u0, phi_u_obj, 1);
		ExplicitDiffusion(voxel, v1, v0, phi_v_obj, 2);
		ExplicitDiffusion(voxel, w1, w0, phi_w_obj, 3);
	}
	else{
		printf("\nperform implicit diffusion with viscosity = %f \n\n", visc);
		ImplicitDiffusion(voxel, u1, u0, phi_u_obj, 1);
		ImplicitDiffusion(voxel, v1, v0, phi_v_obj, 2);
		ImplicitDiffusion(voxel, w1, w0, phi_w_obj, 3);
	}
//	FOR_EACH_CELL
//		if(i == I & j == J && k == K){
//			printf("after diffuse vel is (%f,%f,%f) \n\n",
//					u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)]);
//		}
//	END_FOR
//	SetBoundary(voxel, 1, u);
//	SetBoundary(voxel, 2, v);
//	SetBoundary(voxel, 3, w);

//	Extrapolate(voxel);
//	SetSolidBoundary(voxel);
//	SetSolidBoundaryForAdvection(voxel);
//	SetSurfaceBoundary(voxel);

//	printf(" VectorField3D::Diffuse() finished! \n");
}

void VectorField3D::ClearData(float *data){
	memset(data, 0, (DimX*DimY*DimZ)*sizeof(float));
}



static float HeavisideFunction(float phi, float delta){
	float eps = 1.5f * delta;
	if(phi < -eps)
		return 0.f;
	else if(phi > eps)
		return 1.f;
	else{
		return 0.5f + 0.5f * phi / eps + 0.5 / M_PI * sin(M_PI * phi / eps);
	}
}

static float DeltaFunction(float phi, float delta){
	float eps = 1.5f * delta ;
	float eps_inv = 1.f / eps;
	if(phi < -eps)
		return 0.f;
	else if(phi > eps)
		return 0.f;
	else{
		return 0.5f * eps_inv + 0.5f * eps_inv * cos(M_PI * phi * eps_inv);
	}
}



static const double Gx[] = {
    1.83434642495649804936e-01,    5.25532409916328985830e-01,
    7.96666477413626739567e-01,    9.60289856497536231661e-01
};

static const double GA[] = {
    3.62683783378361982976e-01,    3.13706645877887287338e-01,
    2.22381034453374470546e-01,    1.01228536290376259154e-01
};

static void Gauss_Legendre_Zeros_8pts( double zeros[] ) {

  zeros[0] = -Gx[3];
  zeros[1] = -Gx[2];
  zeros[2] = -Gx[1];
  zeros[3] = -Gx[0];
  zeros[4] = Gx[0];
  zeros[5] = Gx[1];
  zeros[6] = Gx[2];
  zeros[7] = Gx[3];
}

static void Gauss_Legendre_Coefs_8pts( double coefs[] ) {

   coefs[0] = GA[3];
   coefs[1] = GA[2];
   coefs[2] = GA[1];
   coefs[3] = GA[0];
   coefs[4] = GA[0];
   coefs[5] = GA[1];
   coefs[6] = GA[2];
   coefs[7] = GA[3];
}

static float Gauss_Legendre_Integration_8pts(float xa, float xb,
		float ya, float yb, float za, float zb, float delta, float delta3,
		const Voxel &voxel, const VectorField3D *vel, float *phi)
{
   double integral = 0.0;
   double coeffs[8];
   Gauss_Legendre_Coefs_8pts(coeffs);
   double zeros[8];
   Gauss_Legendre_Zeros_8pts(zeros);
   double a = 0.5*(double)(xb-xa);
   double b = 0.5*(double)(yb-ya);
   double c = 0.5*(double)(zb-za);
   for(int k=0; k<8; ++k){
	   double c8k = coeffs[k];
	   double r8k = zeros[k];
	   for(int j=0; j<8; ++j){
		   double c8j = coeffs[j];
		   double r8j = zeros[j];
		   for(int i=0; i<8; ++i){
			   double c8i = coeffs[i];
			   double r8i = zeros[i];
			   double x = 0.5 * (r8i * (xb-xa) + xa + xb);
			   double y = 0.5 * (r8j * (yb-ya) + ya + yb);
			   double z = 0.5 * (r8k * (zb-za) + za + zb);
			   float f = vel->TriInterp(voxel, Point((float)x, (float)y, (float)z), phi);
			   integral += c8k * c8j * c8i * (1.0 - (double)HeavisideFunction(f, delta));
		   }
	   }
   }
   integral *=  a * b * c;
   return (float)integral;
}



float VectorField3D::EstimateFluidVolume(int ii, int jj, int kk, float delta) const{
	float V = 0.f;
	for(int k=kk-1; k<=kk+1; ++k)
		for(int j=jj-1; j<=jj+1; ++j)
			for(int i=ii-1; i<=ii+1; ++i){
				V += 1.f - HeavisideFunction(phi_c_obj[INDEX(i,j,k)], delta);
			}
	return V;
}

void VectorField3D::CutCellFaceArea(const Voxel &voxel, int i, int j, int k, int m, float *c) const{
		float f[4], g[4];
		container->SetBorderCellPhi(voxel, i, j, k, m, f);
		for(int n = 0; n < movingObjects.size(); ++n){
			MovingObject *obj = movingObjects[n];
			obj->SetBorderCellPhi(voxel, i, j, k, m, g);
		}
		for(int l=0; l < 4; ++l){
			if(f[l] + E_EPSIL < 0.f && g[l] + E_EPSIL < 0.f){
				if(f[l] > g[l])
					g[l] = f[l];
			}
		}
		int N = 0;
		for(int l=0; l < 4; ++l){
//			if(i == I && j == J && k == K)
//				printf(" m = %d, g[%d] = %f \n", m, l, g[l]);
			if(g[l] + E_EPSIL < 0.f){
//				if(i == I && j == J && k == K)
//					printf("g[%d] is less than 0 \n", l);
				++N;
			}
		}
		c[m] = (float)N/4;
	}

float VectorField3D::MarchingCubeArea(const Voxel &voxel, int i, int j, int k) const{

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;

	float val[8];
	Point p[8];

	for(int n=0; n<8; ++n){
		p[n]   = voxel.VoxelCornerPosition(i,j,k,n+1);
		val[n] = TriInterp(voxel, p[n], phi_c_obj);
//		if(i == I && j == J && k == K)
//			printf("val[%d] = %f at point[%d] (%f, %f, %f)\n", n, val[n], n, p[n].x, p[n].y, p[n].z);
	}
	return mc->SurfaceArea(val, p, 0.f, delta) / delta2;


}


// this function implements the variational pressure projection method by Batty et al. (2007)
//
// it deals with both static and moving solids with non-axis-aligned boundary
//
// see Bridson's book 4.5.2 for details


void VectorField3D::FloodFillProjection(int ii, int jj, int kk, char *color) const {
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

void VectorField3D::Projection(const Voxel &voxel){

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float half_delta = delta / 2 ;
	float delta3 = delta * delta2;
	map<u_int, KnownPoint> knownPointU;
	map<u_int, KnownPoint> knownPointV;
	map<u_int, KnownPoint> knownPointW;

#ifdef WIN32
//	hash_map<u_int, float> rhs;
	hash_map<u_int, u_int> A;
	hash_map<u_int, u_int>::iterator found;
#else
//	unordered_map<u_int, float> rhs;
	unordered_map<u_int, u_int> A;
	unordered_map<u_int, u_int>::iterator found;
#endif

	char *maskAirVelSet = new char[DimX*DimY*DimZ];
	SetZero(maskAirVelSet);


//	u_int NL = 0;
//	FOR_EACH_CELL
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//			if(phi_c[INDEX(i,j,k)] <= 0.f)
//				++NL;
//		}
//	END_FOR

	char *valid = new char[DimX*DimY*DimZ];
	SetZero(valid);
	u_int N1 = 0;
	FOR_EACH_CELL
//	 	TentativeGridPoint tp(0.f,i,j,k);
//		if(voxel.InSolid(i,j,k) && voxel.CloseToLiquid(tp)){
//			if(hasMovingObj){
//				bool isBorder1 = false;
//				for(int n = 0; n < movingObjects.size(); ++n){
//					MovingObject *obj = movingObjects[n];
//					if(obj->IsBorderCell(voxel, i, j, k)){
//						isBorder1 = true;
//						break;
//					}
//				}
//				bool isBorder2 = container->IsBorderCell(voxel,i,j,k);
//				if(isBorder1 || isBorder2){
//					valid[INDEX(i,j,k)] = 2;
//					++N1;
//				}
//			}
//			else{
//				if(container->IsBorderCell(voxel,i,j,k)){
//					valid[INDEX(i,j,k)] = 2;
//					++N1;
//				}
//			}
//		}
		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
//			if(hasMovingObj){
//				bool isBorder1 = false;
//				for(int n = 0; n < movingObjects.size(); ++n){
//					MovingObject *obj = movingObjects[n];
//					if(obj->IsBorderCell(voxel, i, j, k)){
//						isBorder1 = true;
//						break;
//					}
//				}
//				bool isBorder2 = container->IsBorderCell(voxel,i,j,k);
//				if(isBorder1 || isBorder2){
//					valid[INDEX(i,j,k)] = 2;
//					++N1;
//				}
//				else
//					valid[INDEX(i,j,k)] = 1;
//			}
//			else{
//				if(container->IsBorderCell(voxel,i,j,k)){
//					valid[INDEX(i,j,k)] = 2;
//					++N1;
//				}
//				else
				valid[INDEX(i,j,k)] = 1;
//			}
		}
	END_FOR
//	printf("Before Projection: there are %u liquid points bordering solids \n",N1);

	SetSolidBoundary(voxel, valid);

	// since complex object boundaries could cause enslaved liquid region
	// use flood fill to get rid of these regions
//	char *color = new char[DimX*DimY*DimZ];
//	SetZero(color);
//	FOR_EACH_CELL
//		if(voxel.InSolid(i,j,k))
//			color[INDEX(i,j,k)] = 2;
//	END_FOR
//	TentativeGridPoint floodSeed(0.f, DimX/2, 390, DimZ-(WALL_THICKNESS+1)-1);
//	FloodFillProjection(floodSeed.ii, floodSeed.jj, floodSeed.kk, color);
//	u_int numEnslavePoints = 0;
//	FOR_EACH_CELL
//		u_int index = INDEX(i,j,k);
//		if(valid[index] == 1 && color[index] == 0){
//			valid[index] = 0;
//			++numEnslavePoints;
//		}
//	END_FOR
//	delete [] color;
//	//if(numEnslavePoints > 0)
//	printf(" There are %u enslaved liquid points \n", numEnslavePoints);

//		SetSolidBoundaryForAdvection(voxel);
	SetZero(u0);
	SetZero(v0);

	u_int N = 0;
	float max_rhs = 0.f;
	int max_rhs_i = 0, max_rhs_j = 0,max_rhs_k = 0;
	float dxdt = delta / dt;
	// compute the r.h.s. of the pressure equation
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//	    float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
//		float coeff[6] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
//	    if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && !valid[INDEX(i,j,k)]){
	 	if(valid[pos]){
//		   float coeff[6];
	       /*bool allfluid = true;
		   TentativeGridPoint tp(0.f,i,j,k);
		   vector<TentativeGridPoint> neighbors;
		   tp.Neigbor(neighbors, DimX, DimY, DimZ);
		   for(u_int m=0; m<neighbors.size(); m++){
			   if(voxel.InSolid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
				 float alpha = phi_c_obj[INDEX(i,j,k)]/(phi_c_obj[INDEX(i,j,k)] - phi_c_obj[POS(neighbors[m])]);
				 if(alpha < VOLUME_FRAC_THRESHOLD) alpha = VOLUME_FRAC_THRESHOLD;
				 coeff[m] = alpha*delta*delta2;
				 if(coeff[m] < 0.f || alpha > 1.f){
					 printf("at (%d %d %d) coeff = %f \n", i,j,k, coeff[m]);
					 exit(1);
				 }

			   }
			   else
				   coeff[m] = delta3;
			   if(coeff[m] != delta3){
   		   			allfluid = false;
   		   		}
		   }*/
//		   for(int m=0; m<6;++m)
//			   coeff[m] = 1.f;
//		   if(valid[INDEX(i,j,k)] == 2){
//			   TentativeGridPoint tp(0.f,i,j,k);
//	   		   vector<TentativeGridPoint> neighbors;
//	   		   tp.Neigbor(neighbors, DimX, DimY, DimZ);
//	   		   for(u_int m=0; m<neighbors.size(); m++){
//	   			   CutCellFaceArea(voxel, i, j, k, m, coeff);
//	   		   }
//	   		   Vector NN = GeometricNormalAt(voxel, i, j, k);
//	   		   Vector obj_vel(0.f,0.f,0.f);
//	   		   if(hasMovingObj)
//	   				SetConstraintVelocity(voxel,obj_vel,i,j,k);
//	   		   float S = delta * DeltaFunction(phi_c_obj[INDEX(i,j,k)], delta);
//	   		   float S1 = MarchingCubeArea(voxel, i, j, k);
////	   		   if(S1 > E_EPSIL)
//	   		   if(i == I && j == J && k == K)
//	   			   printf("at (%d, %d, %d), S = %f, S1 = %f \n", i, j, k, S, S1);
//	   		   u0[INDEX(i,j,k)] += S1 * Dot(obj_vel, NN) * delta / dt;
//	   		   if(coeff[0] < E_EPSIL)
//	   			   u[INDEX(i-1,j,k)] = obj_vel.x;
//	   		   if(coeff[1] < E_EPSIL)
//	   			   u[INDEX(i,j,k)] = obj_vel.x;
//	   		   if(coeff[2] < E_EPSIL)
//	   			   v[INDEX(i,j-1,k)] = obj_vel.y;
//	   		   if(coeff[3] < E_EPSIL)
//	   			   v[INDEX(i,j,k)] = obj_vel.y;
//	   		   if(coeff[4] < E_EPSIL)
//	   			   w[INDEX(i,j,k-1)] = obj_vel.z;
//	   		   if(coeff[5] < E_EPSIL)
//	   			   w[INDEX(i,j,k)] = obj_vel.z;
////	   		   if(voxel.InSolid(i,j,k)){
////	   			  for(u_int m=0; m<neighbors.size(); m++){
////	   				if(voxel.InAir(neighbors[m]) || voxel.InSurface(neighbors[m]))
////	   					coeff[m] = 0.f;
////	   			  }
////	   		   }
//	   		  if(i == I && j == J && k == K){
//	   			  printf("Before Projection, S = %f \n", S);
//	   			  for(int n=0; n<6; ++n)
//	   				  printf("Before Projection coeff[%d] = %f \n", n, coeff[n]);
//	   			  printf("N = (%f, %f, %f) obj_vel = (%f, %f, %f) \n",
//	   					  NN.x, NN.y, NN.z, obj_vel.x, obj_vel.y, obj_vel.z);
//	   			  printf("S*Dot(obj_vel, NN) = %f S1*Dot(obj_vel, NN) = %f \n",S * Dot(obj_vel, NN), S1 * Dot(obj_vel, NN));
//	   		  }
//		   }
		   // first part of the r.h.s. (Bridson's book p.70)
//		   u0[INDEX(i,j,k)] += -(coeff[1]*u[INDEX(i,j,k)]-coeff[0]*u[INDEX(i-1,j,k)]+
//		                         coeff[3]*v[INDEX(i,j,k)]-coeff[2]*v[INDEX(i,j-1,k)]+
//		                         coeff[5]*w[INDEX(i,j,k)]-coeff[4]*w[INDEX(i,j,k-1)])*dxdt;
		   u0[pos] = -(u[pos]-u[INDEX(i-1,j,k)]+
		   		       v[pos]-v[INDEX(i,j-1,k)]+
		   		       w[pos]-w[INDEX(i,j,k-1)])*dxdt;
		  /* float V;
		   if(allfluid)
		   		V = delta3;
		   else{
			   float V1 = delta3 * (1.f - HeavisideFunction(phi_c_obj[INDEX(i,j,k)], delta));
			   float V2 = 0.f;
			   for(int m=0; m<8; ++m){
				   Point subp = voxel.SubVoxelCenterPosition(i,j,k,1,m+1);
				   float phip = TriInterp(voxel, subp, phi_c_obj);
				   float VV = (1.f - HeavisideFunction(phip, delta)) * delta3 / 8;
				   if(phip < 0.f)
					   V2 += 1.f;
				   if(i == I && j == J && k == K)
				   		printf(" m = %d, phip = %f, VV = %f \n", m, phip, VV);
			   }
			   V2 *= delta3 / 8;
	//		   Point pmin = voxel.VoxelCornerPosition(i,j,k,1);
	//		   Point pmax = voxel.VoxelCornerPosition(i,j,k,7);
	//		   float V3 = Gauss_Legendre_Integration_8pts(pmin.x, pmax.x, pmin.y, pmax.y,
	//				   pmin.z, pmax.z, delta, delta3, voxel, this, phi_c_obj);
			   V = V2;
			   if(i == I && j == J && k == K)
				   printf("\nat (%d, %d, %d) V1 = %f, V2 = %f \n\n", i,j,k, V1, V2);
		   }

		   // second part of the r.h.s. (Bridson's book p.70)
		   u0[INDEX(i,j,k)] += ((coeff[1]-V)*u1[INDEX(i,j,k)] - (coeff[0]-V)*u1[INDEX(i-1,j,k)] +
						        (coeff[3]-V)*v1[INDEX(i,j,k)] - (coeff[2]-V)*v1[INDEX(i,j-1,k)] +
								(coeff[5]-V)*w1[INDEX(i,j,k)] - (coeff[4]-V)*w1[INDEX(i,j,k-1)])*delta/dt;
//		   u0[INDEX(i,j,k)] += ((coeff[1])*u1[INDEX(i,j,k)] - (coeff[0])*u1[INDEX(i-1,j,k)] +
//							    (coeff[3])*v1[INDEX(i,j,k)] - (coeff[2])*v1[INDEX(i,j-1,k)] +
//							    (coeff[5])*w1[INDEX(i,j,k)] - (coeff[4])*w1[INDEX(i,j,k-1)])*delta/dt;
  */

		   if(fabs(max_rhs) < fabs(u0[pos])){
			   max_rhs = u0[pos];
			   max_rhs_i = i; max_rhs_j = j; max_rhs_k = k;
		   }
//		   rhs.insert( make_pair(pos, u0[pos]) );
		   A.insert( make_pair(pos, N) );
		   ++N;
//		   valid[INDEX(i,j,k)] = 1;
	   }
		/*if(voxel.InSolid(i,j,k)){
			   TentativeGridPoint tp(0.f,i,j,k);
			   if(voxel.CloseToLiquid(tp)){
				   vector<TentativeGridPoint> neighbors;
				   tp.Neigbor(neighbors, DimX, DimY, DimZ);
				   for(u_int m=0; m<neighbors.size(); m++){
						   if(voxel.InLiquid(neighbors[m]) && !voxel.InSurface(neighbors[m])){
								 float alpha = phi_c_obj[POS(neighbors[m])]/(phi_c_obj[POS(neighbors[m])] - phi_c_obj[INDEX(i,j,k)]);
								 if(alpha < VOLUME_FRAC_THRESHOLD) alpha = VOLUME_FRAC_THRESHOLD;
								 coeff[m] = alpha*delta*delta2;
								 if(coeff[m] < 0.f || alpha > 1.f){
										 printf("at (%d %d %d) coeff = %f \n", i,j,k, coeff[m]);
										 exit(1);
								 }

						   }
						   else
								   coeff[m] = 0.f;
				   }
//                 for(int m=0; m<6;++m)
//                         coeff[m] = 1.f;
				   u0[INDEX(i,j,k)] = -(coeff[1]*u[INDEX(i,j,k)]-coeff[0]*u[INDEX(i-1,j,k)]+
										coeff[3]*v[INDEX(i,j,k)]-coeff[2]*v[INDEX(i,j-1,k)]+
										coeff[5]*w[INDEX(i,j,k)]-coeff[4]*w[INDEX(i,j,k-1)])*delta/dt;

//				   float V = 0.f;
//				   float V = delta3 * (1.f - HeavisideFunction(phi_c_obj[INDEX(i,j,k)], delta));
//				   for(int m=0; m<6;++m)
//						V += coeff[m];
//				   V /= 6;
//				   u0[INDEX(i,j,k)] += ((coeff[1]-V)*u1[INDEX(i,j,k)] - (coeff[0]-V)*u1[INDEX(i-1,j,k)] +
//				  						(coeff[3]-V)*v1[INDEX(i,j,k)] - (coeff[2]-V)*v1[INDEX(i,j-1,k)] +
//				  						(coeff[5]-V)*w1[INDEX(i,j,k)] - (coeff[4]-V)*w1[INDEX(i,j,k-1)])*delta/dt;
				   u0[INDEX(i,j,k)] += ((coeff[1])*u1[INDEX(i,j,k)] - (coeff[0])*u1[INDEX(i-1,j,k)] +
				   					    (coeff[3])*v1[INDEX(i,j,k)] - (coeff[2])*v1[INDEX(i,j-1,k)] +
				   						(coeff[5])*w1[INDEX(i,j,k)] - (coeff[4])*w1[INDEX(i,j,k-1)])*delta/dt;
				   if(fabs(max_rhs) < fabs(u0[INDEX(i,j,k)])){
					   max_rhs = u0[INDEX(i,j,k)];
					   max_rhs_i = i; max_rhs_j = j; max_rhs_k = k;
				   }
				   rhs.insert( make_pair(INDEX(i,j,k), u0[INDEX(i,j,k)]) );
				   A.insert( make_pair(INDEX(i,j,k), N) );
				   ++N;
				   valid[INDEX(i,j,k)] = 1;
			   }
		   }*/

	   if(i == I && j == J && k == K){
//	   if(voxel.InSurface(i,j,k)){
//	   if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
//			   printf("SURFACE CELL\n");
//		   else printf("NOT SURFACE CELL\n");
//		   if(voxel.InAir(i,j,k)) printf(" AIR CELL\n");
//		   	else printf("NOT AIR CELL\n");
		   printf("Before Projection: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
			 			i,j,k, u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)],
			 			u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)]);
		   printf("Before Projection: u = %f, v = %f, w = %f, u_1 = %f, v_1 = %f, w_1 = %f \n",
		   			 			u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)],
		   			 			u[INDEX(i+1,j,k)],v[INDEX(i,j+1,k)],w[INDEX(i,j,k+1)]);
		   printf("Before Projection: u1(%d,%d,%d) = %f, v1(%d,%d,%d) = %f, w1(%d,%d,%d) = %f, "
				   "u1(%d,%d,%d) = %f, v1(%d,%d,%d) = %f, w1(%d,%d,%d) = %f \n",
						i,j,k, u1[INDEX(i,j,k)],i,j,k,v1[INDEX(i,j,k)],i,j,k,w1[INDEX(i,j,k)],
					i-1,j,k, u1[INDEX(i-1,j,k)],i,j-1,k,v1[INDEX(i,j-1,k)],i,j,k-1,w1[INDEX(i,j,k-1)]);
		   printf("Before Projection: u1(%d,%d,%d) = %f, v1(%d,%d,%d) = %f, w1(%d,%d,%d) = %f, "
		   			   "u1(%d,%d,%d) = %f, v1(%d,%d,%d) = %f, w1(%d,%d,%d) = %f \n",
		   				i,j,k, u1[INDEX(i,j,k)],i,j,k,v1[INDEX(i,j,k)],i,j,k,w1[INDEX(i,j,k)],
		   			i+1,j,k, u1[INDEX(i+1,j,k)],i,j+1,k,v1[INDEX(i,j+1,k)],i,j,k+1,w1[INDEX(i,j,k+1)]);
//		   printf("Before Projection: coeff[0] =  %f, coeff[1] = %f, coeff[2] = %f, coeff[3] = %f coeff[4] = %f, coeff[5] = %f \n",
//		   				coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
		   printf("Before Projection: u0 = %f, "
				   		"phi = %f, phi[k-1]= %f, phi[k+1] = %f,  \n",
					 		u0[INDEX(i,j,k)], phi_c[INDEX(i,j,k)],
					 		phi_c[INDEX(i,j,k-1)], phi_c[INDEX(i,j,k+1)] );
		   printf("Before Projection: "
   				   		"phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f \n",
   					 		phi_c[INDEX(i-1,j,k)], phi_c[INDEX(i+1,j,k)],
   					 		phi_c[INDEX(i,j-1,k)], phi_c[INDEX(i,j+1,k)] );
		   printf("valid = %d, valid[i-1] = %d, valid[i+1] = %d, valid[j-1] = %d, valid[j+1] = %d, valid[k-1] = %d, valid[k+1] = %d \n ",
				   valid[INDEX(i,j,k)], valid[INDEX(i-1,j,k)], valid[INDEX(i+1,j,k)], valid[INDEX(i,j-1,k)],
		   	   		valid[INDEX(i,j+1,k)], valid[INDEX(i,j,k-1)], valid[INDEX(i,j,k+1)]);

		}
//	   v0[INDEX(i,j,k)] = 0.f;
	END_FOR
	printf("Before Projection: there are %u liquid points \n",N);
	printf("\nmaximum rhs (%f) occurs at (%d, %d, %d) \n", max_rhs,
				max_rhs_i,max_rhs_j,max_rhs_k);
	if(N == 0) {
//		Extrapolate(voxel);
//		ExtrapolateOneForProjection(voxel, knownPointU, phi_u, phi_u_obj, u, 1);
//		ExtrapolateOneForProjection(voxel, knownPointV, phi_v, phi_v_obj, v, 2);
//		ExtrapolateOneForProjection(voxel, knownPointW, phi_w, phi_w_obj, w, 3);
		DetectIsolatedPoints(voxel, maskAirVelSet, phi_u, u, 1);
		DetectIsolatedPoints(voxel, maskAirVelSet, phi_v, v, 2);
		DetectIsolatedPoints(voxel, maskAirVelSet, phi_w, w, 3);
#ifdef SPMD
//		ExtrapolateOneForProjection(voxel, knownPointUAir, phi_u, phi_u_obj, u, u0, 1);
//		ExtrapolateOneForProjection(voxel, knownPointVAir, phi_v, phi_v_obj, v, v0, 2);
//		ExtrapolateOneForProjection(voxel, knownPointWAir, phi_w, phi_w_obj, w, w0, 3);
//		Extrapolate(voxel, valid);
//#else
		ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_u, u, 1);
		ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_v, v, 2);
		ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_w, w, 3);
//		ExtrapolateOneForProjection(voxel, knownPointUAir, phi_u, phi_u_obj, u, 1);
//		ExtrapolateOneForProjection(voxel, knownPointVAir, phi_v, phi_v_obj, v, 2);
//		ExtrapolateOneForProjection(voxel, knownPointWAir, phi_w, phi_w_obj, w, 3);
#endif
//#ifdef AIR_DIV_FREE
//	ApplySourceTerm(voxel);
//	SetEqual(u, u1);
//	SetEqual(v, v1);
//	SetEqual(w, w1);
//	FOR_EACH_CELL
//		if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
//			u1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,1);
//			}
//			if(obj_vel.x != 0.f)
//				u1[INDEX(i,j,k)] = obj_vel.x;
//		}
//		else if(i < DimX-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i+1,j,k) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i+1,j,k)) )
//			u1[INDEX(i,j,k)] = 0.f;
//
//		if(phi_v_obj[INDEX(i,j,k)] >= 0.f){
//			v1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel = Vector(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,2);
//			}
//			if(obj_vel.y != 0.f)
//				v1[INDEX(i,j,k)] = obj_vel.y;
//		}
//		else if(j < DimY-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i,j+1,k) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i,j+1,k)) )
//			v1[INDEX(i,j,k)] = 0.f;
//
//		if(phi_w_obj[INDEX(i,j,k)] >= 0.f){
//			w1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel = Vector(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,3);
//			}
//			if(obj_vel.z != 0.f)
//				w1[INDEX(i,j,k)] = obj_vel.z;
//		}
//		else if( k < DimZ-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i,j,k+1) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i,j,k+1)) )
//			w1[INDEX(i,j,k)] = 0.f;
//	END_FOR
//
//	ProjectionInAirCell(voxel);
//
//#endif
//		SetSolidBoundaryForAdvection(voxel, knownPointU, knownPointV, knownPointW);
		SetSolidBoundaryForAdvection(voxel);


#ifdef AIR_DIV_FREE
	ApplySourceTerm(voxel);
	SetEqual(u, u1);
	SetEqual(v, v1);
	SetEqual(w, w1);
	ProjectionInAirCell(voxel);
#endif

		ApplySourceTerm(voxel);

//#ifdef AIR_DIV_FREE
		for(int n = 0; n < sources.size();++n)
			sources[n]->SetVelocity(voxel, u1, v1, w1);
//#else
//		SetEqual(u, u1);
//		SetEqual(v, v1);
//		SetEqual(w, w1);
//#endif

		delete [] valid;
		delete [] maskAirVelSet;
		return;
	}
//	if(N != NL){
//		printf("NL(%u) != N(%u) \n ", NL, N);
//		exit(1);
//	}
	// set velocity divergence boundary condition here
//	SetBoundary(voxel, 0, u0);
//	SetBoundary(voxel, 0, v0);

	//LinearSolver(voxel, 0, v0, u0, 1, 6);
//	LinearSolver(voxel, 0, v0, u0, DimX*DimY*DimZ);
	float *p = new float[N];
//	memset(p, 0, N);
	SetZero(p, N);
    printf("pressure solver stage 1 \n");
	LinearSolver(voxel, p, A, u0, valid, N);
	found = A.begin();
	for(; found != A.end(); ++found){
		v0[found->first] = p[found->second];
	}
	delete [] p;

	// set pressure boundary condition here
	SetBoundary(voxel, 0, v0);
	if(source)
		source->SourcePressure(v0);

//	 printf("before update: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
//	   I,J,K, u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)]);


//  make the velocities divergence given the pressure solved

	float dtdx = dt / delta;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
#ifdef USE_MESH
		Point p0 = voxel.VoxelCenterPosition(i,j,k);
		Point p1 = voxel.VoxelCenterPosition(i+1,j,k);
		Point p3 = voxel.VoxelCenterPosition(i,j+1,k);
		Point p5 = voxel.VoxelCenterPosition(i,j,k+1);
		double s, t;
#endif
		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
			if(voxel.InSource(i,j,k)){
				printf("at A (%d, %d, %d) point mistakenly as in liquid \n", i,j,k);
				exit(1);
			}

			if(!voxel.InSolid(i+1,j,k) && !voxel.InSource(i+1,j,k)
#ifdef USE_MESH
				&& !mPhysWorld->SegmentIntersectsMesh(p0, p1, &s, &t)
#endif
				){
//			if(!voxel.InSource(i+1,j,k)){
				if(voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k)){
					float air_phi    = phi_c[INDEX(i+1,j,k)];
					float liquid_phi = phi_c[INDEX(i,j,k)];
					float factor;
					if(fabs(liquid_phi) < E_EPSIL){
						printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
								i, j,k, liquid_phi, i+1, j, k, air_phi, v0[INDEX(i,j,k)]);
//						 liquid_phi = -1.e-6*delta;
//						 air_phi = delta + liquid_phi;
						factor = -1.e3f - 1.f;
					 }
					else
						factor = (air_phi - liquid_phi) / liquid_phi;
					u[INDEX(i,j,k)] -= dtdx * factor * v0[INDEX(i,j,k)];
					if(phi_u[INDEX(i,j,k)] > E_EPSIL){
//						KnownPoint t(i,j,k);
//						knownPointUAir.insert( make_pair(INDEX(i,j,k), t) );
						maskAirVelSet[pos] |= AIRUVELSET;

					}
				}
				else if(valid[INDEX(i+1,j,k)])
					u[INDEX(i,j,k)] -= dtdx *(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)]);
				if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
					KnownPoint t(i,j,k);
					knownPointU.insert( make_pair(INDEX(i,j,k), t) );
				}

			}
			if(!voxel.InSolid(i,j+1,k) && !voxel.InSource(i,j+1,k)
#ifdef USE_MESH
				&& !mPhysWorld->SegmentIntersectsMesh(p0, p3, &s, &t)
#endif
				){
//			if(!voxel.InSource(i,j+1,k)){
				if(voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k)){
					float air_phi    = phi_c[INDEX(i,j+1,k)];
					float liquid_phi = phi_c[INDEX(i,j,k)];
					float factor;
					if(fabs(liquid_phi) < E_EPSIL){
						printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
								i, j, k, liquid_phi, i, j+1, k, air_phi, v0[INDEX(i,j,k)]);
//						 liquid_phi = -1.e-6*delta;
//						 air_phi = delta + liquid_phi;
						 factor = -1.e3f - 1.f;
					 }
					else
						factor = (air_phi - liquid_phi) / liquid_phi;
					v[INDEX(i,j,k)] -= dtdx * factor * v0[INDEX(i,j,k)];
					if(phi_v[INDEX(i,j,k)] > E_EPSIL){
//						KnownPoint t(i,j,k);
//						knownPointVAir.insert( make_pair(INDEX(i,j,k), t) );
						maskAirVelSet[pos] |= AIRVVELSET;
					}
				}
				else if(valid[INDEX(i,j+1,k)])
					v[INDEX(i,j,k)] -= dtdx * (v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)]);
				if(phi_v_obj[INDEX(i,j,k)] >= 0.f){
					KnownPoint t(i,j,k);
					knownPointV.insert( make_pair(INDEX(i,j,k), t) );
				}
			}
			if(!voxel.InSolid(i,j,k+1) && !voxel.InSource(i,j,k+1)
#ifdef USE_MESH
				&& !mPhysWorld->SegmentIntersectsMesh(p0, p5, &s, &t)
#endif
			  ){
//			if(!voxel.InSource(i,j,k+1)){
				if(voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1)){
					float air_phi    = phi_c[INDEX(i,j,k+1)];
					float liquid_phi = phi_c[INDEX(i,j,k)];
					float factor;
					if(fabs(liquid_phi) < E_EPSIL){
						printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
								i, j,  k, liquid_phi, i, j, k+1, air_phi, v0[INDEX(i,j,k)]);
//						 liquid_phi = -1.e-6*delta;
//						 air_phi = delta + liquid_phi;
						 factor = -1.e3f - 1.f;
					 }
					else
						factor = (air_phi - liquid_phi) / liquid_phi;
					w[INDEX(i,j,k)] -= dtdx * factor * v0[INDEX(i,j,k)];
					if(phi_w[INDEX(i,j,k)] > E_EPSIL){
//						KnownPoint t(i,j,k);
//						knownPointWAir.insert( make_pair(INDEX(i,j,k), t) );
						maskAirVelSet[pos] |= AIRWVELSET;
					}
				}
				else if(valid[INDEX(i,j,k+1)])
					w[INDEX(i,j,k)] -= dtdx * (v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)]);
				if(phi_w_obj[INDEX(i,j,k)] >= 0.f){
					KnownPoint t(i,j,k);
					knownPointW.insert( make_pair(INDEX(i,j,k), t) );
				}
			}
		}
		else{
			if(voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)){
//				if(voxel.InSource(i,j,k)){
//					printf("at B (%d, %d, %d) point mistakenly as in air \n", i,j,k);
//					exit(1);
//				}
//					if(!voxel.InSolid(i+1,j,k)){
					if(voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)
#ifdef USE_MESH
					  && !mPhysWorld->SegmentIntersectsMesh(p0, p1, &s, &t)
#endif
					  ){
						float air_phi    = phi_c[INDEX(i,j,k)];
						float liquid_phi = phi_c[INDEX(i+1,j,k)];
						float factor;
						if(fabs(liquid_phi) < E_EPSIL){
							printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
									i+1, j, k,liquid_phi,i, j, k, air_phi, v0[INDEX(i+1,j,k)]);
//							 liquid_phi = -1.e-6*delta;
//							 air_phi = delta + liquid_phi;
							factor = 1.e3 + 1.f;
						 }
						else
							factor = (liquid_phi - air_phi) / liquid_phi;
						u[INDEX(i,j,k)] -= dtdx * factor * v0[INDEX(i+1,j,k)];
						if(phi_u[INDEX(i,j,k)] > E_EPSIL){
//							KnownPoint t(i,j,k);
//							knownPointUAir.insert( make_pair(INDEX(i,j,k), t) );
							maskAirVelSet[pos] |= AIRUVELSET;
						}
					}
//					if(!voxel.InSolid(i,j+1,k))
					if(voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)
#ifdef USE_MESH
					  && !mPhysWorld->SegmentIntersectsMesh(p0, p3, &s, &t)
#endif
					  ){
						float air_phi    = phi_c[INDEX(i,j,k)];
						float liquid_phi = phi_c[INDEX(i,j+1,k)];
						float factor;
						if(fabs(liquid_phi) < E_EPSIL){
							printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
									i, j+1, k, liquid_phi, i, j, k, air_phi, v0[INDEX(i,j+1,k)]);
//							 liquid_phi = -1.e-6*delta;
//							 air_phi = delta + liquid_phi;
							 factor = 1.e3 + 1.f;
						 }
						else
							factor = (liquid_phi - air_phi) / liquid_phi;
						v[INDEX(i,j,k)] -= dtdx * factor * v0[INDEX(i,j+1,k)];
						if(phi_v[INDEX(i,j,k)] > E_EPSIL){
//							KnownPoint t(i,j,k);
//							knownPointVAir.insert( make_pair(INDEX(i,j,k), t) );
							maskAirVelSet[pos] |= AIRVVELSET;
						}
					}
//					if(!voxel.InSolid(i,j,k+1))
					if(voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)
#ifdef USE_MESH
						&& !mPhysWorld->SegmentIntersectsMesh(p0, p5, &s, &t)
#endif
						){
						float air_phi    = phi_c[INDEX(i,j,k)];
						float liquid_phi = phi_c[INDEX(i,j,k+1)];
						float factor;
						if(fabs(liquid_phi) < E_EPSIL){
							printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
									i, j, k+1, liquid_phi, i, j, k, air_phi, v0[INDEX(i,j,k+1)]);
//							 liquid_phi = -1.e-6*delta;
//							 air_phi = delta + liquid_phi;
							 factor = 1.e3 + 1.f;
						 }
						else
							factor = (liquid_phi - air_phi) / liquid_phi;
						w[INDEX(i,j,k)] -= dtdx * factor * v0[INDEX(i,j,k+1)];
						if(phi_w[INDEX(i,j,k)] > E_EPSIL){
//							KnownPoint t(i,j,k);
//							knownPointWAir.insert( make_pair(INDEX(i,j,k), t) );
							maskAirVelSet[pos] |= AIRWVELSET;
						}
					}
			}
					/*if(voxel.InSolid(i+1,j,k) && valid[INDEX(i+1,j,k)]){
						float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
						CutCellFaceArea(voxel, i, j, k, 1, coeff);
						if(coeff[1] > E_EPSIL){
							float air_phi = phi_c[INDEX(i,j,k)];
							float liquid_phi = phi_c[INDEX(i+1,j,k)];
							if(liquid_phi < 0.f){
								if(fabs(liquid_phi) < 1.e-6*delta){
									printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
											i+1, j, k,liquid_phi,i, j, k, air_phi, v0[INDEX(i+1,j,k)]);
									 liquid_phi = -1.e-6*delta;
									 air_phi = delta + liquid_phi;
								 }
								u[INDEX(i,j,k)] -= dt*(1.f+ air_phi / (-liquid_phi))*
												 v0[INDEX(i+1,j,k)]/delta;
							}
							else
								u[INDEX(i,j,k)] -= dt*(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)])/delta;
							if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointU.insert( make_pair(INDEX(i,j,k), t) );
							}
						}
					}
					if(voxel.InSolid(i,j+1,k) && valid[INDEX(i,j+1,k)]){
						float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
						CutCellFaceArea(voxel, i, j, k, 3, coeff);
						if(coeff[3] > E_EPSIL){
							float air_phi = phi_c[INDEX(i,j,k)];
							float liquid_phi = phi_c[INDEX(i,j+1,k)];
							if(liquid_phi < 0.f){
								if(fabs(liquid_phi) < 1.e-6*delta){
									printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
											i+1, j, k,liquid_phi,i, j, k, air_phi, v0[INDEX(i+1,j,k)]);
									 liquid_phi = -1.e-6*delta;
									 air_phi = delta + liquid_phi;
								 }
								v[INDEX(i,j,k)] -= dt*(1.f+ air_phi / (-liquid_phi))*
												 v0[INDEX(i,j+1,k)]/delta;
							}
							else
								v[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)])/delta;
							if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointV.insert( make_pair(INDEX(i,j,k), t) );
							}
						}
					}
					if(voxel.InSolid(i,j,k+1) && valid[INDEX(i,j,k+1)]){
						float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
						CutCellFaceArea(voxel, i, j, k, 5, coeff);
						if(coeff[5] > E_EPSIL){
							float air_phi = phi_c[INDEX(i,j,k)];
							float liquid_phi = phi_c[INDEX(i,j,k+1)];
							if(liquid_phi < 0.f){
								if(fabs(liquid_phi) < 1.e-6*delta){
									printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
											i+1, j, k,liquid_phi,i, j, k, air_phi, v0[INDEX(i+1,j,k)]);
									 liquid_phi = -1.e-6*delta;
									 air_phi = delta + liquid_phi;
								 }
								w[INDEX(i,j,k)] -= dt*(1.f+ air_phi / (-liquid_phi))*
												 v0[INDEX(i,j,k+1)]/delta;
							}
							else
								w[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)])/delta;
							if(phi_w_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointW.insert( make_pair(INDEX(i,j,k), t) );
							}
						}
					}*/
			}
			/*else if(voxel.InSolid(i,j,k) && valid[INDEX(i,j,k)]){
				TentativeGridPoint tp(0.f,i,j,k);
				if(voxel.CloseToLiquid(tp)){
					if(voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)){
							u[INDEX(i,j,k)] -= dt*(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)])/delta;
							if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointU.insert( make_pair(INDEX(i,j,k), t) );
							}
					}
					if(voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)){
							v[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)])/delta;
							if(phi_v_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointV.insert( make_pair(INDEX(i,j,k), t) );
							}
					}
					if(voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)){
							w[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)])/delta;
							if(phi_w_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointW.insert( make_pair(INDEX(i,j,k), t) );
							}
					}
				}
				if(voxel.InAir(i+1,j,k) || voxel.InSurface(i+1,j,k)){
						float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
						CutCellFaceArea(voxel, i, j, k, 1, coeff);
						if(coeff[1] > E_EPSIL){
							float air_phi = phi_c[INDEX(i+1,j,k)];
							float liquid_phi = phi_c[INDEX(i,j,k)];
							if(liquid_phi < 0.f){
								if(fabs(liquid_phi) < 1.e-6*delta){
									printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
											i, j,k, liquid_phi, i+1, j, k, air_phi, v0[INDEX(i,j,k)]);
									 liquid_phi = -1.e-6*delta;
									 air_phi = delta + liquid_phi;
								 }
								u[INDEX(i,j,k)] -= dt*(-air_phi / (-liquid_phi) - 1.f)*
												 v0[INDEX(i,j,k)]/delta;
								}
							else
								u[INDEX(i,j,k)] -= dt*(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)])/delta;
							if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointU.insert( make_pair(INDEX(i,j,k), t) );
							}
						}
				}
				if(voxel.InAir(i,j+1,k) || voxel.InSurface(i,j+1,k)){
						float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
						CutCellFaceArea(voxel, i, j, k, 3, coeff);
						if(coeff[3] > E_EPSIL){
							float air_phi = phi_c[INDEX(i,j+1,k)];
							float liquid_phi = phi_c[INDEX(i,j,k)];
							if(liquid_phi < 0.f){
								if(fabs(liquid_phi) < 1.e-6*delta){
									printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
											i, j,k, liquid_phi, i+1, j, k, air_phi, v0[INDEX(i,j,k)]);
									 liquid_phi = -1.e-6*delta;
									 air_phi = delta + liquid_phi;
								 }
								v[INDEX(i,j,k)] -= dt*(-air_phi / (-liquid_phi) - 1.f)*
												 v0[INDEX(i,j,k)]/delta;
								}
							else
								v[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)])/delta;
							if(phi_v_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointV.insert( make_pair(INDEX(i,j,k), t) );
							}
						}
				}
				if(voxel.InAir(i,j,k+1) || voxel.InSurface(i,j,k+1)){
						float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
						CutCellFaceArea(voxel, i, j, k, 5, coeff);
						if(coeff[5] > E_EPSIL){
							float air_phi = phi_c[INDEX(i,j,k+1)];
							float liquid_phi = phi_c[INDEX(i,j,k)];
							if(liquid_phi < 0.f){
								if(fabs(liquid_phi) < 1.e-6*delta){
									printf("phil(%d, %d, %d) = %16.13f, phia(%d, %d, %d) = %16.13f, v0 = %16.13f  \n ",
											i, j,k, liquid_phi, i+1, j, k, air_phi, v0[INDEX(i,j,k)]);
									 liquid_phi = -1.e-6*delta;
									 air_phi = delta + liquid_phi;
								 }
								w[INDEX(i,j,k)] -= dt*(-air_phi / (-liquid_phi) - 1.f)*
												 v0[INDEX(i,j,k)]/delta;
							}
							else
								w[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)])/delta;
							if(phi_w_obj[INDEX(i,j,k)] >= 0.f){
								KnownPoint t(i,j,k);
								knownPointW.insert( make_pair(INDEX(i,j,k), t) );
							}
						}
				}
				if(voxel.InSolid(i+1,j,k) && valid[INDEX(i+1,j,k)]){
					float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
					CutCellFaceArea(voxel, i, j, k, 1, coeff);
					if(coeff[1] > E_EPSIL)
						u[INDEX(i,j,k)] -= dt*(v0[INDEX(i+1,j,k)]-v0[INDEX(i,j,k)])/delta;
				}
				if(voxel.InSolid(i,j+1,k) && valid[INDEX(i,j+1,k)]){
					float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
					CutCellFaceArea(voxel, i, j, k, 3, coeff);
					if(coeff[3] > E_EPSIL)
						v[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j+1,k)]-v0[INDEX(i,j,k)])/delta;
				}
				if(voxel.InSolid(i,j,k+1) && valid[INDEX(i,j,k+1)]){
					float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
					CutCellFaceArea(voxel, i, j, k, 5, coeff);
					if(coeff[5] > E_EPSIL)
						w[INDEX(i,j,k)] -= dt*(v0[INDEX(i,j,k+1)]-v0[INDEX(i,j,k)])/delta;
				}
			}*/


//		if(voxel.InSurface(i,j,k)){
//		float div = u[INDEX(i,j,k)]-u[INDEX(i-1,j,k)]
//                      +v[INDEX(i,j,k)]-v[INDEX(i,j-1,k)]
//                         +w[INDEX(i,j,k)]-w[INDEX(i,j,k-1)];
//	if(voxel.InLiquid(i,j,k) && fabs(div) > 1.e-3){
		if(i == I && j == J && k == K){
		 printf("After Projection: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
				 			i,j,k, u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)]);
		 printf("After Projection: at i = %d, j = %d, k = %d  u1 = %f, v1 = %f, w1 = %f \n",
		 						i,j,k, u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)]);
		 printf("After Projection: at i = %d, j = %d, k = %d  phiu = %f, phiv = %f, phiw = %f \n",
		 		 						i,j,k, phi_u[INDEX(i,j,k)],phi_v[INDEX(i,j,k)],phi_w[INDEX(i,j,k)]);
		 printf("After Projection: at i = %d, j = %d, k = %d  phiu1 = %f, phiv1 = %f, phiw1 = %f \n",
		 		 				i,j,k, phi_u[INDEX(i-1,j,k)],phi_v[INDEX(i,j-1,k)],phi_w[INDEX(i,j,k-1)]);
		 printf("%f, %f, %f, %f\n",v0[INDEX(i+1,j,k)], v0[INDEX(i,j+1,k)],
				 v0[INDEX(i,j,k+1)], v0[INDEX(i,j,k)]);
		 printf("%f, %f, %f, %f\n",v0[INDEX(i-1,j,k)], v0[INDEX(i,j-1,k)],
		 				 v0[INDEX(i,j,k-1)], v0[INDEX(i,j,k)]);
		 printf("divergence = %16.13f \n",
				 u[INDEX(i,j,k)]-u[INDEX(i-1,j,k)]
				                   +v[INDEX(i,j,k)]-v[INDEX(i,j-1,k)]
				                                 +w[INDEX(i,j,k)]-w[INDEX(i,j,k-1)]);
//		 exit(1);
		}

	END_FOR

	FOR_EACH_CELL
		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
				float div = u[INDEX(i,j,k)]-u[INDEX(i-1,j,k)]
		            		 +v[INDEX(i,j,k)]-v[INDEX(i,j-1,k)]
		            	     +w[INDEX(i,j,k)]-w[INDEX(i,j,k-1)];
				if(fabs(div) > 1.e-3){
					printf("(%d,%d,%d) phi = %f, u = %f, v = %f, w = %f, div = %f\n",
								i, j, k, phi_c[INDEX(i,j,k)],
						u[INDEX(i,j,k)], v[INDEX(i,j,k)], w[INDEX(i,j,k)], div);
//					exit(1);
				}
		}
	END_FOR

	/*FOR_EACH_CELL
		float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
		if(valid[INDEX(i,j,k)]){

	   for(int m=0; m<6;++m)
		   coeff[m] = 1.f;
	   float div = 0.f;
	   float S, S1;
	   Vector NN;
	   Vector obj_vel(0.f, 0.f, 0.f);
	   if(valid[INDEX(i,j,k)] == 2){
		   TentativeGridPoint tp(0.f,i,j,k);
		   vector<TentativeGridPoint> neighbors;
		   tp.Neigbor(neighbors, DimX, DimY, DimZ);
		   for(u_int m=0; m<neighbors.size(); m++){
			   CutCellFaceArea(voxel, i, j, k, m, coeff);
		   }
		   NN = GeometricNormalAt(voxel, i, j, k);
		   if(hasMovingObj)
				SetConstraintVelocity(voxel,obj_vel,i,j,k);
		   S = delta * DeltaFunction(phi_c_obj[INDEX(i,j,k)], delta);
		   S1 = MarchingCubeArea(voxel, i, j, k);
		   div += -S1 * Dot(obj_vel, NN);
	   }
	   // first part of the r.h.s. (Bridson's book p.70)
	   div += (coeff[1]*u[INDEX(i,j,k)]-coeff[0]*u[INDEX(i-1,j,k)]+
			   coeff[3]*v[INDEX(i,j,k)]-coeff[2]*v[INDEX(i,j-1,k)]+
			   coeff[5]*w[INDEX(i,j,k)]-coeff[4]*w[INDEX(i,j,k-1)]);
	   if( fabsf(div) > 1.e-3 ){
//	   if( i == I && j == J && k == K ){
//		   printf(" \nat (%d, %d, %d) div (%f) != 0 \n", i, j, k, div);
		   printf("(%d,%d,%d) phi = %f, div = %f, u = %f, v = %f, w = %f, u1 = %f, v1 = %f, w1 = %f\n",
		   					i, j, k, phi_c[INDEX(i,j,k)], div,
		   				u[INDEX(i,j,k)], v[INDEX(i,j,k)], w[INDEX(i,j,k)],
		   				u[INDEX(i-1,j,k)], v[INDEX(i,j-1,k)], w[INDEX(i,j,k-1)]);
		   printf("Before Projection: coeff[0] =  %f, coeff[1] = %f, coeff[2] = %f, coeff[3] = %f coeff[4] = %f, coeff[5] = %f \n",
		  		   	coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
		   if(valid[INDEX(i,j,k)] == 2){
			  printf("N = (%f, %f, %f) obj_vel = (%f, %f, %f) \n",
					  NN.x, NN.y, NN.z, obj_vel.x, obj_vel.y, obj_vel.z);
			  printf("S = %f, S1 = %f, Dot(obj_vel, NN) = %f \n",S, S1, Dot(obj_vel, NN));
			  printf("-S*Dot(obj_vel, NN) = %f -S1*Dot(obj_vel, NN) = %f\n\n", -S*Dot(obj_vel, NN), -S1*Dot(obj_vel, NN));
		  }
		   exit(1);
	   }
		}

	   END_FOR*/


	// set velocity boundary condition here


//	Extrapolate(voxel);

	DetectIsolatedPoints(voxel, maskAirVelSet, phi_u, u, 1);
	DetectIsolatedPoints(voxel, maskAirVelSet, phi_v, v, 2);
	DetectIsolatedPoints(voxel, maskAirVelSet, phi_w, w, 3);

#ifdef SPMD
//	ExtrapolateOneForProjection(voxel, knownPointUAir, phi_u, phi_u_obj, u, u0, 1);
//	ExtrapolateOneForProjection(voxel, knownPointVAir, phi_v, phi_v_obj, v, v0, 2);
//	ExtrapolateOneForProjection(voxel, knownPointWAir, phi_w, phi_w_obj, w, w0, 3);
//	Extrapolate(voxel, valid);
//#else
	ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_u, u, 1);
	ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_v, v, 2);
	ExtrapolateOneUsingFM(voxel, maskAirVelSet, phi_w, w, 3);
//	ExtrapolateOneForProjection(voxel, knownPointUAir, phi_u, phi_u_obj, u, 1);
//	ExtrapolateOneForProjection(voxel, knownPointVAir, phi_v, phi_v_obj, v, 2);
//	ExtrapolateOneForProjection(voxel, knownPointWAir, phi_w, phi_w_obj, w, 3);
#endif

//#ifdef AIR_DIV_FREE
//	ApplySourceTerm(voxel);
//	SetEqual(u, u1);
//	SetEqual(v, v1);
//	SetEqual(w, w1);
//	FOR_EACH_CELL
//		if(phi_u_obj[INDEX(i,j,k)] >= 0.f){
//			u1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,1);
//			}
//			if(obj_vel.x != 0.f)
//				u1[INDEX(i,j,k)] = obj_vel.x;
//		}
//		else if(i < DimX-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i+1,j,k) ||
//				           voxel.InSolid(i,j,k) && !voxel.InSolid(i+1,j,k)) )
//			u1[INDEX(i,j,k)] = 0.f;
//
//		if(phi_v_obj[INDEX(i,j,k)] >= 0.f){
//			v1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel = Vector(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,2);
//			}
//			if(obj_vel.y != 0.f)
//				v1[INDEX(i,j,k)] = obj_vel.y;
//		}
//		else if(j < DimY-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i,j+1,k) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i,j+1,k)) )
//			v1[INDEX(i,j,k)] = 0.f;
//
//		if(phi_w_obj[INDEX(i,j,k)] >= 0.f){
//			w1[INDEX(i,j,k)] = 0.f;
//			Vector obj_vel = Vector(0.f,0.f,0.f);
//			if(hasMovingObj){
//				SetConstraintVelocity(voxel,obj_vel,i,j,k,3);
//			}
//			if(obj_vel.z != 0.f)
//				w1[INDEX(i,j,k)] = obj_vel.z;
//		}
//		else if( k < DimZ-1 && ( !voxel.InSolid(i,j,k) && voxel.InSolid(i,j,k+1) ||
//						   voxel.InSolid(i,j,k) && !voxel.InSolid(i,j,k+1)) )
//			w1[INDEX(i,j,k)] = 0.f;
//	END_FOR
//	ProjectionInAirCell(voxel);
//
//#endif

//	SetSolidBoundaryForAdvection(voxel, knownPointU, knownPointV, knownPointW);
	SetSolidBoundaryForAdvection(voxel);

#ifdef AIR_DIV_FREE
	ApplySourceTerm(voxel);
	SetEqual(u, u1);
	SetEqual(v, v1);
	SetEqual(w, w1);
	ProjectionInAirCell(voxel);
#endif

	ApplySourceTerm(voxel);

//#ifdef AIR_DIV_FREE
		for(int n = 0; n < sources.size();++n)
			sources[n]->SetVelocity(voxel, u1, v1, w1);
//#else
//		SetEqual(u, u1);
//		SetEqual(v, v1);
//		SetEqual(w, w1);
//		printf("u1, v1 and w1 are set values of u, v, and w respectively\n");
//#endif

	delete [] valid;
	delete [] maskAirVelSet;
//	FOR_EACH_CELL
////		if(i == I & j == J && k == K){
//		if(voxel.InAir(i,j,k) && phi_c[INDEX(i,j,k)] < 3*delta){
//
//			float div = u[INDEX(i,j,k)]-u[INDEX(i-1,j,k)]
//			                   +v[INDEX(i,j,k)]-v[INDEX(i,j-1,k)]
//		                                 +w[INDEX(i,j,k)]-w[INDEX(i,j,k-1)];
//			if(fabs(div) > 1.e-3){
//			 printf("\n (Air) divergence bigger than 1.e-3 = %16.13f and phi = %f \n",
//					 div, phi_c[INDEX(i,j,k)]);
//			 printf("After Projection: at i = %d, j = %d, k = %d  phiu = %f, phiv = %f, phiw = %f \n",
//	 						i,j,k,
//	 						phi_u[INDEX(i,j,k)],phi_v[INDEX(i,j,k)],phi_w[INDEX(i,j,k)]);
//			 printf("After Projection: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
// 						i,j,k,
// 						u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)]);
//			 printf("After Projection: at i = %d, j = %d, k = %d  u1 = %f, v1 = %f, w1 = %f \n\n",
// 						i,j,k,
// 						u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)]);
//			 exit(1);
//			}
//		}
//		if(voxel.InLiquid(i,j,k)){
//
//				float div = u[INDEX(i,j,k)]-u[INDEX(i-1,j,k)]
//				                   +v[INDEX(i,j,k)]-v[INDEX(i,j-1,k)]
//			                                 +w[INDEX(i,j,k)]-w[INDEX(i,j,k-1)];
//				if(fabs(div) > 1.e-3){
//					 printf("\n (Liquid) divergence bigger than 1.e-3 = %16.13f and phi = %f \n",
//							 div, phi_c[INDEX(i,j,k)]);
//					 printf("After Projection: at i = %d, j = %d, k = %d  phiu = %f, phiv = %f, phiw = %f \n",
//			 						i,j,k,
//			 						phi_u[INDEX(i,j,k)],phi_v[INDEX(i,j,k)],phi_w[INDEX(i,j,k)]);
//					 printf("After Projection: at i = %d, j = %d, k = %d  u = %f, v = %f, w = %f \n",
//		 						i,j,k,
//		 						u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)]);
//					 printf("After Projection: at i = %d, j = %d, k = %d  u1 = %f, v1 = %f, w1 = %f \n\n",
//		 						i,j,k,
//		 						u[INDEX(i-1,j,k)],v[INDEX(i,j-1,k)],w[INDEX(i,j,k-1)]);
//					 exit(1);
//				}
//			}
//
//	END_FOR

}


#ifdef COUPLED_FLIP
void VectorField3D::SetGridVelocity(int i, int j, int k, int index, float *vel) const{
	if(index == 1)
		vel[INDEX(i,j,k)] = u[INDEX(i,j,k)];
	if(index == 2)
		vel[INDEX(i,j,k)] = v[INDEX(i,j,k)];
	if(index == 3)
		vel[INDEX(i,j,k)] = w[INDEX(i,j,k)];

}
#endif

void VectorField3D::SetBoundary(const Voxel &voxel, int b, float *x){

	//SetSolidBoundary(voxel, b, x);

	FOR_EACH_CELL

		if(voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)){
			if( b == 0 )
				x[INDEX(i,j,k)] = 0.f;

		}

//		if(voxel.InSurface(i,j,k)){
//			if( b == 0 )
//				x[INDEX(i,j,k)] = 0.f;
//		}

	END_FOR

}




void VectorField3D::SetAirBoundary(const Voxel &voxel, int b, float *x){
	FOR_EACH_CELL

//	if(voxel.InAir(i,j,k)){
//		TentativeGridPoint p(0.f,i,j,k);
//		vector<TentativeGridPoint> neighbors;
//		p.Neigbor(neighbors, DimX, DimY, DimZ);
//		for(int m=0; m<neighbors.size(); m++){
//			if(voxel.InAir(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)){
//				if(p.LeftNeighbor(neighbors[m])){
//					if ( b == 1 ){
//						x[INDEX(i-1,j,k)] = 0.f;
//					}
//				}
//				if(p.RightNeighbor(neighbors[m])){
//					if ( b == 1 ){
//						x[INDEX(i,j,k)] = 0.f;
//					}
//				}
//				if(p.BackNeighbor(neighbors[m])){
//					if ( b == 2 ){
//						x[INDEX(i,j-1,k)] = 0.f;
//					}
//				}
//				if(p.FrontNeighbor(neighbors[m])){
//					if ( b == 2 ){
//						x[INDEX(i,j,k)] = 0.f;
//					}
//				}
//				if(p.BottomNeighbor(neighbors[m])){
//					if ( b == 3 ){
//						x[INDEX(i,j,k-1)] = 0.f;
//					}
//				}
//				if(p.TopNeighbor(neighbors[m])){
//					if ( b == 3 ){
//						x[INDEX(i,j,k)] = 0.f;
//					}
//				}
//			}
//		}
//	}

	if(voxel.InAir(i,j,k) || voxel.InLiquid(i,j,k)){
		if( b == 1 ){
			if( phi_u[INDEX(i,j,k)] > 0.f )
				x[INDEX(i,j,k)] = 0.f;
		}
		else if( b == 2 ){
			if( phi_v[INDEX(i,j,k)] > 0.f )
				x[INDEX(i,j,k)] = 0.f;
		}
		else if( b == 3 ){
			if( phi_w[INDEX(i,j,k)] > 0.f )
				x[INDEX(i,j,k)] = 0.f;
		}
	}
	else{   // in solid
		if( b == 1 ){
			if( phi_u[INDEX(i,j,k)] > 0.f && i < DimX-1 && !voxel.InSolid(i+1,j,k) )
				x[INDEX(i,j,k)] = 0.f;
		}
		else if( b == 2 ){
			if( phi_v[INDEX(i,j,k)] > 0.f && j < DimY-1 && !voxel.InSolid(i,j+1,k) )
				x[INDEX(i,j,k)] = 0.f;
		}
		else if( b == 3 ){
			if( phi_w[INDEX(i,j,k)] > 0.f && k < DimZ-1 && !voxel.InSolid(i,j,k+1) )
				x[INDEX(i,j,k)] = 0.f;
		}

	}
	END_FOR
}

/*void VectorField3D::SetSolidBoundary(const Voxel &voxel, int b, float *x){

	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k))
			x[INDEX(i,j,k)] = 0.f;
	END_FOR

	FOR_EACH_CELL

		if(voxel.InSolid(i,j,k)){
			TentativeGridPoint p(0.f,i,j,k);
			vector<TentativeGridPoint> neighbors;
			p.Neigbor(neighbors, DimX, DimY, DimZ);
			for(int m=0; m<neighbors.size(); m++){
				if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk) ||
				   voxel.InAir(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk) ||
				   voxel.InSurface(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)){
					if(p.LeftNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i-1,j,k)];
						else if ( b == 1 ){
							x[INDEX(i-1,j,k)] = 0.f;
							x[INDEX(i,j,k)] = 0.f;
						}
						else if ( b == 2 || b == 3){
							x[INDEX(i,j,k)] = x[INDEX(i-1,j,k)];
						}
					}
					if(p.RightNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i+1,j,k)];
						else if ( b == 1 ){
							x[INDEX(i,j,k)] = 0.f;
						}
						else if ( b == 2 || b == 3){
							x[INDEX(i,j,k)] = x[INDEX(i+1,j,k)];
						}
					}
					if(p.BackNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j-1,k)];
						else if ( b == 2 ){
							x[INDEX(i,j-1,k)] = 0.f;
							x[INDEX(i,j,k)] = 0.f;
						}
						else if ( b == 1 || b == 3){
							x[INDEX(i,j,k)] = x[INDEX(i,j-1,k)];
						}
					}
					if(p.FrontNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j+1,k)];
						else if ( b == 2 ){
							x[INDEX(i,j,k)] = 0.f;
						}
						else if ( b == 1 || b == 3){
							x[INDEX(i,j,k)] = x[INDEX(i,j+1,k)];
						}
					}
					if(p.BottomNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j,k-1)];
						else if ( b == 3 ){
							x[INDEX(i,j,k-1)] = 0.f;
							x[INDEX(i,j,k)] = 0.f;
						}
						else if ( b == 1 || b == 2){
							x[INDEX(i,j,k)] = x[INDEX(i,j,k-1)];
						}
					}
					if(p.TopNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j,k+1)];
						else if ( b == 3 ){
							x[INDEX(i,j,k)] = 0.f;
						}
						else if ( b == 1 || b == 2){
							x[INDEX(i,j,k)] = x[INDEX(i,j,k+1)];
						}
					}
				}
//				else{
//					x[INDEX(i,j,k)] = 0.f;
//				}

			}
		}

	END_FOR

}*/

#ifndef SPMD

void VectorField3D::SetSolidBoundary(const Voxel &voxel, int b, float *x){

	float delta = voxel.VoxelDelta();

	// first, extrapolate velocity into object using object phi
	if(b == 1)
		ExtrapolateOneVelocityIntoObject(voxel, phi_u_obj, x, 1, delta);
	else if(b == 2)
		ExtrapolateOneVelocityIntoObject(voxel, phi_v_obj, x, 2, delta);
	else if(b == 3)
		ExtrapolateOneVelocityIntoObject(voxel, phi_w_obj, x, 3, delta);

	//apply constraints on extrapolated velocity
	FOR_EACH_CELL

		if(voxel.InSolid(i,j,k)){
			TentativeGridPoint p(0.f,i,j,k);
			vector<TentativeGridPoint> neighbors;
			p.Neigbor(neighbors, DimX, DimY, DimZ);
			for(int m=0; m<neighbors.size(); m++){
				if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk) ||
				   voxel.InAir(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk) ||
				   voxel.InSurface(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)){
					if(p.LeftNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i-1,j,k)];
						else if ( b == 1 ){
//							if(x[INDEX(i,j,k)] > 0.f)
								x[INDEX(i,j,k)] = 0.f;
//							if(x[INDEX(i-1,j,k)] > 0.f)
								x[INDEX(i-1,j,k)] = 0.f;
						}
					}
					if(p.RightNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i+1,j,k)];
						else if ( b == 1 ){
//							if(x[INDEX(i,j,k)] < 0.f)
								x[INDEX(i,j,k)] = 0.f;
						}
					}
					if(p.BackNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j-1,k)];
						else if ( b == 2 ){
//							if(x[INDEX(i,j,k)] > 0.f)
								x[INDEX(i,j,k)] = 0.f;
//							if(x[INDEX(i,j-1,k)] > 0.f)
								x[INDEX(i,j-1,k)] = 0.f;
						}
					}
					if(p.FrontNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j+1,k)];
						else if ( b == 2 ){
//							if(x[INDEX(i,j,k)] < 0.f)
								x[INDEX(i,j,k)] = 0.f;
						}
					}
					if(p.BottomNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j,k-1)];
						else if ( b == 3 ){
//							if(x[INDEX(i,j,k)] > 0.f)
								x[INDEX(i,j,k)] = 0.f;
//							if(x[INDEX(i,j,k-1)] > 0.f)
								x[INDEX(i,j,k-1)] = 0.f;
						}
					}
					if(p.TopNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j,k+1)];
						else if ( b == 3 ){
//							if(x[INDEX(i,j,k)] < 0.f)
								x[INDEX(i,j,k)] = 0.f;
						}
					}
				}
			}
		}

	END_FOR

}

#endif

// set velocities within solid boundary
// this function is called before projection step such that velocities
// bordering a solid cell have proper values.
/*void VectorField3D::SetSolidBoundary(const Voxel &voxel){

	float delta = voxel.VoxelDelta();
	FOR_EACH_CELL
	if(i == I && j == J && k == K){
		printf(" at (%d, %d, %d) \n", i,j,k);
		printf("phiu_obj = %f, phiv_obj = %f, phiw_obj = %f \n",
		   		phi_u_obj[INDEX(i,j,k)], phi_v_obj[INDEX(i,j,k)],
			 		phi_w_obj[INDEX(i,j,k)]);
		printf(" before interpolate into object u = %f, v = %f, w = %f \n",
				u[INDEX(i,j,k)],v[INDEX(i,j,k)],w[INDEX(i,j,k)]);

	}
	END_FOR
	// first, extrapolate velocity into object using object phi
	ExtrapolateOneVelocityIntoObject(voxel, phi_u_obj, u, 1, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_v_obj, v, 2, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_w_obj, w, 3, delta);
//	printf(" after ExtrapolateOneVelocityIntoObject in VectorField3D::SetSolidBoundary()\n");
	//apply constraints on extrapolated velocity
	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k)){
			if(i == I && j == J && k == K){
				printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
//			Vector N1 = ObjectNormalAt(delta, i, j, k);
			Vector N1 = ObjectNormalAt(voxel, i, j, k);
			Vector N = GeometricNormalAt(voxel, i, j, k);
			float u_c, v_c, w_c;
//			if(i == 0)
				u_c = u[INDEX(i,j,k)];
//			else
//				u_c = 0.5f*(u[INDEX(i-1,j,k)]+u[INDEX(i,j,k)]);
//			if(j == 0)
				v_c = v[INDEX(i,j,k)];
//			else
//				v_c = 0.5f*(v[INDEX(i,j-1,k)]+v[INDEX(i,j,k)]);
//			if(k == 0)
				w_c = w[INDEX(i,j,k)];
//			else
//				w_c = 0.5f*(w[INDEX(i,j,k-1)]+w[INDEX(i,j,k)]);
			if(i == I && j == J && k == K){
				printf(" after interpolate into object u = %f, v = %f, w= %f \n",
						u_c, v_c, w_c);
			}
			Vector vel(u_c, v_c, w_c);
			Vector v_para;
			float v_perp;
			v_para = vel - Dot(vel, N) * N;
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj){
				SetConstraintVelocity(obj_vel,i,j,k);
				if(i == I && j == J && k == K){
					printf("object velocity = (%f, %f,  %f) \n",
							obj_vel.x, obj_vel.y, obj_vel.z);
				}
				v_perp = Dot(vel-obj_vel, N);
				if(v_perp < 0.f)
					v_perp = Dot(obj_vel, N);
				else
					v_perp = Dot(vel, N);
			}
			else{
				v_perp = Dot(vel, N);
				if(v_perp != 0.f)
					v_perp = 0.f;
			}
			Vector v_comb = v_para + v_perp * N;
			u_c = Dot(v_comb, Vector(1,0,0));
			v_c = Dot(v_comb, Vector(0,1,0));
			w_c = Dot(v_comb, Vector(0,0,1));
			u[INDEX(i,j,k)] = u_c;
			v[INDEX(i,j,k)] = v_c;
			w[INDEX(i,j,k)] = w_c;
			if(i > 0 && !voxel.InSolid(i-1,j,k))
				u[INDEX(i-1,j,k)] = u_c;
			if(j > 0 && !voxel.InSolid(i,j-1,k) )
				v[INDEX(i,j-1,k)] = v_c;
			if(k > 0 && !voxel.InSolid(i,j,k-1))
				w[INDEX(i,j,k-1)] = w_c;
			if(i == I && j == J && k == K){
				printf("(SetSolidBoundary): u[(%d,%d,%d)] = %f, u[(%d,%d,%d)] = %f\n",
						i, j, k, u[INDEX(i-1,j,k)], i,j,k,u[INDEX(i,j,k)]);
				printf("(SetSolidBoundary): v[(%d,%d,%d)] = %f, v[(%d,%d,%d)] = %f\n",
						i, j, k, v[INDEX(i,j-1,k)], i,j,k,v[INDEX(i,j,k)]);
				printf("(SetSolidBoundary): w[(%d,%d,%d)] = %f, w[(%d,%d,%d)] = %f\n",
						i, j, k, w[INDEX(i,j,k-1)], i,j,k,w[INDEX(i,j,k)]);
				printf("v_comb = (%f, %f, %f) \n",v_comb.x,v_comb.y,v_comb.z);
				printf("v_para = (%f, %f, %f) \n",v_para.x,v_para.y,v_para.z);
				printf("N = (%f, %f, %f) \n",N.x,N.y,N.z);
				printf("N1 = (%f, %f, %f) \n",N1.x,N1.y,N1.z);
				printf("v_perp = %f \n",v_perp);
			}
		}
	END_FOR

}*/

void VectorField3D::SetMovingSolidBoundary(const Voxel &voxel){
	float delta = voxel.VoxelDelta();
	FOR_EACH_CELL
		Vector obj_vel(0.f,0.f,0.f);
		if(hasMovingObj){
			SetConstraintVelocity(voxel,obj_vel,i,j,k,1);
		}
		u1[INDEX(i,j,k)] = obj_vel.x;

		obj_vel = Vector(0.f,0.f,0.f);
		if(hasMovingObj){
			SetConstraintVelocity(voxel,obj_vel,i,j,k,2);
		}
		v1[INDEX(i,j,k)] = obj_vel.y;

		obj_vel = Vector(0.f,0.f,0.f);
		if(hasMovingObj){
			SetConstraintVelocity(voxel,obj_vel,i,j,k,3);
		}
		w1[INDEX(i,j,k)] = obj_vel.z;
	END_FOR

	FOR_EACH_CELL
		phi_u_obj[INDEX(i,j,k)] *= -1;
		phi_v_obj[INDEX(i,j,k)] *= -1;
		phi_w_obj[INDEX(i,j,k)] *= -1;
	END_FOR

#ifdef SPMD
	ExtrapolateOneVelocityIntoLiquid(voxel, phi_u_obj, u1, temp, 1, delta);
	ExtrapolateOneVelocityIntoLiquid(voxel, phi_v_obj, v1, temp, 2, delta);
	ExtrapolateOneVelocityIntoLiquid(voxel, phi_w_obj, w1, temp, 3, delta);
#else
	ExtrapolateOneVelocityIntoLiquid(voxel, phi_u_obj, u1, 1, delta);
	ExtrapolateOneVelocityIntoLiquid(voxel, phi_v_obj, v1, 2, delta);
	ExtrapolateOneVelocityIntoLiquid(voxel, phi_w_obj, w1, 3, delta);
#endif

	FOR_EACH_CELL
		phi_u_obj[INDEX(i,j,k)] *= -1;
		phi_v_obj[INDEX(i,j,k)] *= -1;
		phi_w_obj[INDEX(i,j,k)] *= -1;
	END_FOR
}

void VectorField3D::SetMovingSolidBoundary2(const Voxel &voxel){
	float delta = voxel.VoxelDelta();
	FOR_EACH_CELL
		Vector obj_vel(0.f,0.f,0.f);
		if(hasMovingObj){
			SetConstraintVelocity(voxel,obj_vel,i,j,k,1);
		}
		u1[INDEX(i,j,k)] = obj_vel.x;

		obj_vel = Vector(0.f,0.f,0.f);
		if(hasMovingObj){
			SetConstraintVelocity(voxel,obj_vel,i,j,k,2);
		}
		v1[INDEX(i,j,k)] = obj_vel.y;

		obj_vel = Vector(0.f,0.f,0.f);
		if(hasMovingObj){
			SetConstraintVelocity(voxel,obj_vel,i,j,k,3);
		}
		w1[INDEX(i,j,k)] = obj_vel.z;
	END_FOR

	FOR_EACH_CELL
		phi_u_obj[INDEX(i,j,k)] *= -1;
		phi_v_obj[INDEX(i,j,k)] *= -1;
		phi_w_obj[INDEX(i,j,k)] *= -1;
	END_FOR

#ifdef SPMD
	ExtrapolateOneVelocityIntoNonSolid(voxel, phi_u_obj, u1, temp, 1, delta);
	ExtrapolateOneVelocityIntoNonSolid(voxel, phi_v_obj, v1, temp, 2, delta);
	ExtrapolateOneVelocityIntoNonSolid(voxel, phi_w_obj, w1, temp, 3, delta);
#else
	ExtrapolateOneVelocityIntoNonSolid(voxel, phi_u_obj, u1, 1, delta);
	ExtrapolateOneVelocityIntoNonSolid(voxel, phi_v_obj, v1, 2, delta);
	ExtrapolateOneVelocityIntoNonSolid(voxel, phi_w_obj, w1, 3, delta);
#endif

	FOR_EACH_CELL
		phi_u_obj[INDEX(i,j,k)] *= -1;
		phi_v_obj[INDEX(i,j,k)] *= -1;
		phi_w_obj[INDEX(i,j,k)] *= -1;
	END_FOR
}

// this function performs two tasks:
// 1. it extrapolates fluid velocity into (moving) solid objects (u, v, and w)
// 2. it also extrapolates (moving) solid velocities into fluid region (u1, v1, and w1)
// these velocities are used to compute the r.h.s. of the divergence
// equation in Projection function
//
// the extrapolation is done based on phi_u_obj, phi_v_obj and phi_w_obj (for step 1)
// or their flip-signed version (for step 2)
//
// this function is called before projection step such that velocities
// bordering a solid cell have proper values.
/*void VectorField3D::SetSolidBoundary(const Voxel &voxel){

	float delta = voxel.VoxelDelta();
//	printf(" after ExtrapolateOneVelocityIntoObject in VectorField3D::SetSolidBoundary()\n");
	//apply constraints on extrapolated velocity
//	FOR_EACH_CELL
//		if(voxel.InSolid(i,j,k)){
//
//				Vector obj_vel(0.f,0.f,0.f);
//				if(hasMovingObj){
//					SetConstraintVelocity(voxel,obj_vel,i,j,k);
//					if(i == I && j == J && k == K){
//						printf("object velocity = (%f, %f,  %f) \n",
//								obj_vel.x, obj_vel.y, obj_vel.z);
//					}
//
//				}
//				if(phi_u_obj[INDEX(i,j,k)] >= 0.f)
//					u[INDEX(i,j,k)] = obj_vel.x;
//				if(phi_v_obj[INDEX(i,j,k)] >= 0.f)
//					v[INDEX(i,j,k)] = obj_vel.y;
//				if(phi_w_obj[INDEX(i,j,k)] >= 0.f)
//					w[INDEX(i,j,k)] = obj_vel.z;
//				if(i > 0 && !voxel.InSolid(i-1,j,k) && phi_u_obj[INDEX(i-1,j,k)] >= 0.f)
//					u[INDEX(i-1,j,k)] = obj_vel.x;
//				if(j > 0 && !voxel.InSolid(i,j-1,k) && phi_v_obj[INDEX(i,j-1,k)] >= 0.f)
//					v[INDEX(i,j-1,k)] = obj_vel.y;
//				if(k > 0 && !voxel.InSolid(i,j,k-1) && phi_w_obj[INDEX(i,j,k-1)] >= 0.f)
//					w[INDEX(i,j,k-1)] = obj_vel.z;
//
//		}
//	END_FOR

#ifdef SPMD
	ExtrapolateOneVelocityIntoObject(voxel, phi_u_obj, u, u0, 1, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_v_obj, v, v0, 2, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_w_obj, w, w0, 3, delta);
#else
	ExtrapolateOneVelocityIntoObject(voxel, phi_u_obj, u, 1, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_v_obj, v, 2, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_w_obj, w, 3, delta);
#endif

	// set (moving) solid velocities
	// u1, v1 and w1 get either 0 (no moving solid) or actual
	// solid velocity if the cell is within the moving solid
	SetMovingSolidBoundary(voxel);

}*/

void VectorField3D::SetSolidBoundaryHouston2003(const Voxel &voxel, char *valid){


	float delta = voxel.VoxelDelta();

	printf("VectorField3D::SetSolidBoundaryHouston2003\n\n");
//	FOR_EACH_CELL
//
//		if( k == K && i > 0 && j > 0 ){
////		if( k == K ){
//			if(voxel.InSolid(i,j,k) && !voxel.InSolid(i-1,j,k) && !voxel.InSolid(i,j-1,k))
////			if(voxel.InSolid(i,j,k))
//				printf(" i = %d, j = %d, k = %d, phi_c_obj = %f \n", i,j,k, phi_c_obj[INDEX(i,j,k)]);
//		}
//
//	END_FOR

	// use occlusion velocity when computing divergence on the r.h.s. of pressure equation
	FOR_EACH_CELL
		if( voxel.InSolid(i,j,k) ){
			if(i == I && j == J && k == K){
				printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k);
			}
			u[INDEX(i,j,k)] = obj_vel.x;
			v[INDEX(i,j,k)] = obj_vel.y;
			w[INDEX(i,j,k)] = obj_vel.z;
			if(i > 0 && !voxel.InSolid(i-1,j,k))
				u[INDEX(i-1,j,k)] = obj_vel.x;
			if(j > 0 && !voxel.InSolid(i,j-1,k))
				v[INDEX(i,j-1,k)] = obj_vel.y;
			if(k > 0 && !voxel.InSolid(i,j,k-1))
				w[INDEX(i,j,k-1)] = obj_vel.z;
		}
	END_FOR


}

void VectorField3D::SetSolidBoundary(const Voxel &voxel, char *valid){

	float delta = voxel.VoxelDelta();

//	FOR_EACH_CELL
//		if( k == K && i > 0 && j > 0 ){
////		if( k == K ){
//			if(voxel.InSolid(i,j,k) && !voxel.InSolid(i-1,j,k) && !voxel.InSolid(i,j-1,k))
////			if(voxel.InSolid(i,j,k))
//				printf(" i = %d, j = %d, k = %d, phi_c_obj = %f \n", i,j,k, phi_c_obj[INDEX(i,j,k)]);
//		}
//	END_FOR


	SetEqual(u, u1);
	SetEqual(v, v1);
	SetEqual(w, w1);
#ifdef SPMD
	ExtrapolateOneVelocityIntoObject(voxel, phi_u_obj, u1, temp, 1, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_v_obj, v1, temp, 2, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_w_obj, w1, temp, 3, delta);
#endif

//	memset(needRecomputePhi, 0, DimX*DimY*DimZ);
	SetZero(needRecomputePhi);

//	if(voxel.InSolid(I,J,K))
//		printf("voxel in solid at (%d, %d, %d) \n", I, J, K);
//	else
//		printf("voxel not in solid at (%d, %d, %d) \n", I, J, K);
//	printf(" after ExtrapolateOneVelocityIntoObject in VectorField3D::SetSolidBoundary()\n");
	//apply constraints on extrapolated velocity
	u_int n1 = 0, n2 = 0;
	FOR_EACH_CELL
	    u_int index = INDEX(i,j,k);
		if( voxel.InSolid(i,j,k) ){
			if(i == I && j == J && k == K){
				printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k);
			}
			if(voxel.InMovingSolid(i,j,k)){  // for moving objects, use FF01 approach
				Vector N = GeometricNormalAt(voxel, i, j, k);
				float u_c, v_c, w_c;
				if(i == 0)
					u_c = u1[INDEX(i,j,k)];
				else
					u_c = 0.5f*(u1[INDEX(i-1,j,k)]+u1[INDEX(i,j,k)]);
				if(j == 0)
					v_c = v1[INDEX(i,j,k)];
				else
					v_c = 0.5f*(v1[INDEX(i,j-1,k)]+v1[INDEX(i,j,k)]);
				if(k == 0)
					w_c = w1[INDEX(i,j,k)];
				else
					w_c = 0.5f*(w1[INDEX(i,j,k-1)]+w1[INDEX(i,j,k)]);
				if(i == I && j == J && k == K){
					printf(" after interpolate into object u = %f, v = %f, w= %f \n",
							u_c, v_c, w_c);
				}
				Vector vel(u_c, v_c, w_c);
				Vector v_para;
				float v_perp;
				v_para = vel - Dot(vel, N) * N;

				v_perp = Dot(vel-obj_vel, N);


				float vel_magnitude = (vel-obj_vel).Length();
				if(i == I && j == J && k == K){
					printf("object velocity = (%f, %f,  %f), v_perp = %f, vel_mag = %f\n",
							obj_vel.x, obj_vel.y, obj_vel.z, v_perp, vel_magnitude);
				}
				if(v_perp - 0.1f*vel_magnitude > E_EPSIL){
					needRecomputePhi[INDEX(i,j,k)] = 1;
					if(voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)] < 0.f)
						++n2;
				}
				if(v_perp < E_EPSIL)
					v_perp = Dot(obj_vel, N);
				else{
					TentativeGridPoint tp(0.f, i, j, k);
					if(IsDynamicBorderCell(voxel,i,j,k)){
						v_perp = Dot(vel, N);
	//					valid[INDEX(i,j,k)] = 2;
						++n1;
					}
					else
						v_perp = Dot(obj_vel, N);
	//				++n1;
	//				printf("at (%d, %d, %d) object velocity = (%f, %f,  %f), v_perp = %f, vel_mag = %f\n",
	//						i, j, k, obj_vel.x, obj_vel.y, obj_vel.z, v_perp, vel_magnitude);
	//				v_perp = Dot(vel, N);
	//				valid[INDEX(i,j,k)] = 2;
	//				v_perp = Dot(vel, N);
	//				printf("(at (%d, %d, %d) velocity = (%f, %f,  %f), v_perp = %f, vel_mag = %f\n",
	//						i,j,k, vel.x, vel.y, vel.z, v_perp, vel_magnitude);
				}
	//				v_perp = 0.f;

				Vector v_comb = v_para + v_perp * N;
				u_c = Dot(v_comb, Vector(1,0,0));
				v_c = Dot(v_comb, Vector(0,1,0));
				w_c = Dot(v_comb, Vector(0,0,1));
				u0[INDEX(i,j,k)] = u_c;
				v0[INDEX(i,j,k)] = v_c;
				w0[INDEX(i,j,k)] = w_c;
	//			if(i > 0 && !voxel.InSolid(i-1,j,k))
	//				u[INDEX(i-1,j,k)] = u_c;
	//			if(j > 0 && !voxel.InSolid(i,j-1,k) )
	//				v[INDEX(i,j-1,k)] = v_c;
	//			if(k > 0 && !voxel.InSolid(i,j,k-1))
	//				w[INDEX(i,j,k-1)] = w_c;
				if(i == I && j == J && k == K){
	//				printf("(SetSolidBoundary): u[(%d,%d,%d)] = %f, u[(%d,%d,%d)] = %f\n",
	//						i, j, k+1, u[INDEX(i,j,k+1)], i,j,k,u[INDEX(i,j,k)]);
	//				printf("(SetSolidBoundary): v[(%d,%d,%d)] = %f, v[(%d,%d,%d)] = %f\n",
	//						i, j, k+1, v[INDEX(i,j,k+1)], i,j,k,v[INDEX(i,j,k)]);
	//				printf("(SetSolidBoundary): w[(%d,%d,%d)] = %f, w[(%d,%d,%d)] = %f\n",
	//						i, j, k+1, w[INDEX(i,j,k+1)], i,j,k,w[INDEX(i,j,k)]);
					printf("v_comb = (%f, %f, %f) \n",v_comb.x,v_comb.y,v_comb.z);
					printf("v_para = (%f, %f, %f) \n",v_para.x,v_para.y,v_para.z);
					printf("N = (%f, %f, %f) \n",N.x,N.y,N.z);
					printf("v_perp = %f \n",v_perp);
				}
			}
			else{
				u[INDEX(i,j,k)] = obj_vel.x;
				v[INDEX(i,j,k)] = obj_vel.y;
				w[INDEX(i,j,k)] = obj_vel.z;
				if(i > 0 && !voxel.InSolid(i-1,j,k))
					u[INDEX(i-1,j,k)] = obj_vel.x;
				if(j > 0 && !voxel.InSolid(i,j-1,k))
					v[INDEX(i,j-1,k)] = obj_vel.y;
				if(k > 0 && !voxel.InSolid(i,j,k-1))
					w[INDEX(i,j,k-1)] = obj_vel.z;
			}
		}
#ifdef USE_MESH
		if(valid[index]){
			Point p0 = voxel.VoxelCenterPosition(i,j,k);
			Point p1 = voxel.VoxelCenterPosition(i+1,j,k);
			Point p2 = voxel.VoxelCenterPosition(i-1,j,k);
			Point p3 = voxel.VoxelCenterPosition(i,j+1,k);
			Point p4 = voxel.VoxelCenterPosition(i,j-1,k);
			Point p5 = voxel.VoxelCenterPosition(i,j,k+1);
			Point p6 = voxel.VoxelCenterPosition(i,j,k-1);
			double s, t;
			if(mPhysWorld->SegmentIntersectsMesh(p0, p1, &s, &t)){
				printf(" p0(%f,%f,%f) -> p1(%f,%f,%f) intersects object in "
						"cell(%d,%d,%d) with s = %f t =%f \n",
						p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, i, j, k, s, t);
				Point pint1 = p0 + s*(p1-p0);
				Point pint2 = p1 + t*(p0-p1);
				printf(" p0->p1 intersects at (%f,%f,%f) or (%f,%f,%f) \n ",
						pint1.x, pint1.y, pint1.z, pint2.x, pint2.y, pint2.z);
				u[index] = 0.f;
			}
			if(mPhysWorld->SegmentIntersectsMesh(p0, p2, &s, &t)){
				printf(" p0(%f,%f,%f) -> p2(%f,%f,%f) intersects object in "
				"cell(%d,%d,%d) with s = %f t =%f \n",
				p0.x, p0.y, p0.z, p2.x, p2.y, p2.z, i, j, k, s, t);
				Point pint1 = p0 + s*(p2-p0);
				Point pint2 = p2 + t*(p0-p2);
				printf(" p0->p2 intersects at (%f,%f,%f) or (%f,%f,%f) \n ",
						pint1.x, pint1.y, pint1.z, pint2.x, pint2.y, pint2.z);
				u[INDEX(i-1,j,k)] = 0.f;
			}
			if(mPhysWorld->SegmentIntersectsMesh(p0, p3, &s, &t)){
				printf(" p0(%f,%f,%f) -> p3(%f,%f,%f) intersects object in "
						"cell(%d,%d,%d) with s = %f t =%f \n",
						p0.x, p0.y, p0.z, p3.x, p3.y, p3.z, i, j, k, s, t);
				v[index] = 0.f;
			}
			if(mPhysWorld->SegmentIntersectsMesh(p0, p4, &s, &t)){
				printf(" p0(%f,%f,%f) -> p4(%f,%f,%f) intersects object in "
						"cell(%d,%d,%d) with s = %f t =%f \n",
						p0.x, p0.y, p0.z, p4.x, p4.y, p4.z, i, j, k, s, t);
				v[INDEX(i,j-1,k)] = 0.f;
			}
			if(mPhysWorld->SegmentIntersectsMesh(p0, p5, &s, &t)){
				printf(" p0(%f,%f,%f) -> p5(%f,%f,%f) intersects object in "
						"cell(%d,%d,%d) with s = %f t =%f \n",
						p0.x, p0.y, p0.z, p5.x, p5.y, p5.z, i, j, k, s, t);
				w[index] = 0.f;
			}
			if(mPhysWorld->SegmentIntersectsMesh(p0, p6, &s, &t)){
				printf(" p0(%f,%f,%f) -> p6(%f,%f,%f) intersects object in "
						"cell(%d,%d,%d) with s = %f t =%f \n",
						p0.x, p0.y, p0.z, p6.x, p6.y, p6.z, i, j, k, s, t);
				w[INDEX(i,j,k-1)] = 0.f;
			}
		}
#endif
	END_FOR
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(voxel.InMovingSolid(i,j,k)){  // for moving objects, use FF01 approach
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj)
				SetConstraintVelocity(voxel,obj_vel,i,j,k);
			if(i < DimX-1){
				if(voxel.InSolid(i+1,j,k))
					u[pos] = obj_vel.x;
				else
					u[pos] = u0[pos];
			}
			else
				u[pos] = u0[pos];
			if(j < DimY-1){
				if(voxel.InSolid(i,j+1,k))
					v[pos] = obj_vel.y;
				else
					v[pos] = v0[pos];
			}
			else
				v[pos] = v0[pos];
			if(k < DimZ-1){
				if(voxel.InSolid(i,j,k+1))
					w[pos] = obj_vel.z;
				else
					w[pos] = w0[pos];
			}
			else
				w[pos] = w0[pos];

			if(i > 0 && !voxel.InSolid(i-1,j,k))
				u[INDEX(i-1,j,k)] =  u0[pos];
			if(j > 0 && !voxel.InSolid(i,j-1,k))
				v[INDEX(i,j-1,k)] =  v0[pos];
			if(k > 0 && !voxel.InSolid(i,j,k-1))
				w[INDEX(i,j,k-1)] =  w0[pos];

		}
	END_FOR
//	if(mWaveGenerator)
//		mWaveGenerator->SetBoundaryVelocity(&voxel, DimX, DimY, DimZ,
//				WALL_THICKNESS+1, WALL_THICKNESS+1, WALL_THICKNESS+1,
//				u, v, w);
	printf("Finished Setting velocities at solid boundary cells \n");
}

// set velocities within solid boundary for subsequent levelset and particle advection
//
// this function is called after projection step such that divergence-free velocities
// are propogated


/*void VectorField3D::SetSolidBoundaryForAdvection(const Voxel &voxel,
		map<u_int, KnownPoint> &knownPointU,
		map<u_int, KnownPoint> &knownPointV,
		map<u_int, KnownPoint> &knownPointW){

	float delta = voxel.VoxelDelta();
	map<u_int, KnownPoint>::iterator foundpos;

#ifdef SPMD
#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, u1, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, v1, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, w1, 3, delta);
#else
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, 3, delta);
#endif
	SetEqual(u, u0);
	SetEqual(v, v0);
	SetEqual(w, w0);

//	printf(" got here (1) \n");

	FOR_EACH_CELL

		foundpos = knownPointU.find(INDEX(i,j,k));
		if(phi_u_obj[INDEX(i,j,k)] >= 0.f && foundpos == knownPointU.end()){
			Vector vel;
			vel.x = u0[INDEX(i,j,k)];
			if( j > 0 && i < DimX-1)
				vel.y = 0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i+1,j,k)]+v0[INDEX(i+1,j-1,k)]);
			else
				vel.y = v0[INDEX(i,j,k)];
			if( k > 0 && i < DimX-1)
				vel.z = 0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i+1,j,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i+1,j,k-1)]);
			else
				vel.z = w0[INDEX(i,j,k)];
			Point p = voxel.VelPosition(1, i, j, k);
			Vector N = ObjectNormalAt(voxel, p);
			Vector obj_vel(0.f,0.f,0.f);
			float alpha = 0.f;
			container->Friction(voxel, 1, i, j, k, &alpha);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,1);
				ConstraintFriction(voxel, 1, i, j, k, &alpha);
			}
			float v_perp = Dot(vel-obj_vel, N);
			if(v_perp < 0.f)
				v_perp = Dot(obj_vel, N);
			else
				v_perp = Dot(vel, N);
			Vector v_comb = (1 - alpha) * vel + alpha * obj_vel;
			Vector v_para = v_comb - Dot(v_comb, N) * N;
			Vector v_final = v_para + v_perp * N;
			u[INDEX(i,j,k)] = v_final.x;
		}
//		printf(" at (%d, %d, %d) u finished \n", i,j,k);
		foundpos = knownPointV.find(INDEX(i,j,k));
		if(phi_v_obj[INDEX(i,j,k)] >= 0.f && foundpos == knownPointV.end()){
			Vector vel;
			if( i > 0 && j < DimY-1)
				vel.x = 0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j+1,k)]+u0[INDEX(i-1,j+1,k)]);
			else
				vel.x = u0[INDEX(i,j,k)];
			vel.y = v0[INDEX(i,j,k)];
			if( k > 0 && j < DimY-1)
				vel.z = 0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i,j+1,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i,j+1,k-1)]);
			else
				vel.z = w0[INDEX(i,j,k)];
			Point p = voxel.VelPosition(2, i, j, k);
			Vector N = ObjectNormalAt(voxel, p);
			Vector obj_vel(0.f,0.f,0.f);
			float alpha = 0.f;
			container->Friction(voxel, 2, i, j, k, &alpha);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,2);
				ConstraintFriction(voxel, 2, i, j, k, &alpha);
			}
			float v_perp = Dot(vel-obj_vel, N);
			if(v_perp < 0.f)
				v_perp = Dot(obj_vel, N);
			else
				v_perp = Dot(vel, N);
			Vector v_comb = (1 - alpha) * vel + alpha * obj_vel;
			Vector v_para = v_comb - Dot(v_comb, N) * N;
			Vector v_final = v_para + v_perp * N;
			v[INDEX(i,j,k)] = v_final.y;
		}

//		printf(" at (%d, %d, %d) v finished \n", i,j,k);

		foundpos = knownPointW.find(INDEX(i,j,k));
		if(phi_w_obj[INDEX(i,j,k)] >= 0.f && foundpos == knownPointW.end()){
			Vector vel;
			if( i > 0 && k < DimZ-1)
				vel.x = 0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j,k+1)]+u0[INDEX(i-1,j,k+1)]);
			else
				vel.x = u0[INDEX(i,j,k)];
			if( j > 0 && k < DimZ-1)
				vel.y = 0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i,j,k+1)]+v0[INDEX(i,j-1,k+1)]);
			else
				vel.y = v0[INDEX(i,j,k)];
			vel.z = w0[INDEX(i,j,k)];
			Point p = voxel.VelPosition(3, i, j, k);
			Vector N = ObjectNormalAt(voxel, p);
			Vector obj_vel(0.f,0.f,0.f);
			float alpha = 0.f;
			container->Friction(voxel, 3, i, j, k, &alpha);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,3);
				ConstraintFriction(voxel, 3, i, j, k, &alpha);
			}
			float v_perp = Dot(vel-obj_vel, N);
			if(v_perp < 0.f)
				v_perp = Dot(obj_vel, N);
			else
				v_perp = Dot(vel, N);
			Vector v_comb = (1 - alpha) * vel + alpha * obj_vel;
			Vector v_para = v_comb - Dot(v_comb, N) * N;
			Vector v_final = v_para + v_perp * N;
			w[INDEX(i,j,k)] = v_final.z;
		}

//		printf(" at (%d, %d, %d) w finished \n", i,j,k);

	END_FOR
//	printf(" got here (2) \n");

}*/

void VectorField3D::SetSolidBoundaryForAdvectionHouston2003(const Voxel &voxel, const char *valid){

	float delta = voxel.VoxelDelta();
	printf("VectorField3D::SetSolidBoundaryForAdvectionHouston2003\n\n");

	SetEqual(u, u1);
	SetEqual(v, v1);
	SetEqual(w, w1);

#ifdef SPMD
#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, u0, 1, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, v0, 2, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, w0, 3, delta);
	ExtrapolateOneVelocityIntoObjectForAdvectionHouston2003(voxel, phi_u_obj, u1, u0, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvectionHouston2003(voxel, phi_v_obj, v1, v0, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvectionHouston2003(voxel, phi_w_obj, w1, w0, 3, delta);
#else
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, 1, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, 2, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, 3, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_u_obj, u1, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_v_obj, v1, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_w_obj, w1, 3, delta);
#endif
	SetEqual(u1, u0);
	SetEqual(v1, v0);
	SetEqual(w1, w0);
	printf("before doing adjustment point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
			I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );
//	memset(needRecomputePhi, 0, DimX*DimY*DimZ);
//	SetZero(needRecomputePhi);

//	if(voxel.InSolid(I,J,K))
//		printf("voxel in solid at (%d, %d, %d) \n", I, J, K);
//	else
//		printf("voxel not in solid at (%d, %d, %d) \n", I, J, K);
//	printf(" after ExtrapolateOneVelocityIntoObject in VectorField3D::SetSolidBoundary()\n");
	//apply constraints on extrapolated velocity
	u_int n1 = 0, n2 = 0;
	FOR_EACH_CELL
		if( voxel.InSolid(i,j,k)  ){
//		if( voxel.InSolid(i,j,k) ){
//		if(valid[INDEX(i,j,k)]){
			if(i == I && j == J && k == K){
				printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
//			Vector N = ObjectNormalAt(delta, i, j, k);
//			Vector N = ObjectNormalAt(voxel, i, j, k);
			Vector N = GeometricNormalAt(voxel, i, j, k);
			float u_c, v_c, w_c;
			if(i > 0)
				u_c = 0.5f * ( u1[INDEX(i,j,k)] + u1[INDEX(i-1,j,k)] );
			else
				u_c = u1[INDEX(i,j,k)];
			if(j > 0)
				v_c = 0.5f * ( v1[INDEX(i,j,k)] + v1[INDEX(i,j-1,k)] );
			else
				v_c = v1[INDEX(i,j,k)];
			if(k > 0)
				w_c = 0.5f * ( w1[INDEX(i,j,k)] + w1[INDEX(i,j,k-1)] );
			else
				w_c = w1[INDEX(i,j,k)];
			if(i == I && j == J && k == K){
				printf(" after interpolate into object u = %f, v = %f, w= %f \n",
						u_c, v_c, w_c);
			}
			Vector vel(u_c, v_c, w_c);
			Vector v_para;
			float v_perp;
			v_para = vel - Dot(vel, N) * N;
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj)
				SetConstraintVelocity(voxel,obj_vel,i,j,k);
			v_perp = Dot(vel-obj_vel, N);

			float vel_magnitude = (vel-obj_vel).Length();
			if(i == I && j == J && k == K){
				printf("object velocity = (%f, %f,  %f), v_perp = %f, vel_mag = %f\n",
						obj_vel.x, obj_vel.y, obj_vel.z, v_perp, vel_magnitude);
			}
//			if(v_perp > 0.1f*vel_magnitude){
//				needRecomputePhi[INDEX(i,j,k)] = 1;
//				if(voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)] < 0.f)
//					++n2;
//			}
			if(v_perp < E_EPSIL)
				v_perp = Dot(obj_vel, N);
			else
				v_perp = Dot(vel, N);

			Vector v_comb = v_para + v_perp * N;

			u_c = Dot(v_comb, Vector(1,0,0));
			v_c = Dot(v_comb, Vector(0,1,0));
			w_c = Dot(v_comb, Vector(0,0,1));
			u0[INDEX(i,j,k)] = u_c;
			v0[INDEX(i,j,k)] = v_c;
			w0[INDEX(i,j,k)] = w_c;
//			if(i > 0 && !voxel.InSolid(i-1,j,k))
//				u[INDEX(i-1,j,k)] = u_c;
//			if(j > 0 && !voxel.InSolid(i,j-1,k) )
//				v[INDEX(i,j-1,k)] = v_c;
//			if(k > 0 && !voxel.InSolid(i,j,k-1))
//				w[INDEX(i,j,k-1)] = w_c;
			if(i == I && j == J && k == K){
//				printf("(SetSolidBoundary): u[(%d,%d,%d)] = %f, u[(%d,%d,%d)] = %f\n",
//						i, j, k+1, u[INDEX(i,j,k+1)], i,j,k,u[INDEX(i,j,k)]);
//				printf("(SetSolidBoundary): v[(%d,%d,%d)] = %f, v[(%d,%d,%d)] = %f\n",
//						i, j, k+1, v[INDEX(i,j,k+1)], i,j,k,v[INDEX(i,j,k)]);
//				printf("(SetSolidBoundary): w[(%d,%d,%d)] = %f, w[(%d,%d,%d)] = %f\n",
//						i, j, k+1, w[INDEX(i,j,k+1)], i,j,k,w[INDEX(i,j,k)]);
				printf("v_comb = (%f, %f, %f) \n",v_comb.x,v_comb.y,v_comb.z);
				printf("v_para = (%f, %f, %f) \n",v_para.x,v_para.y,v_para.z);
				printf("N = (%f, %f, %f) \n",N.x,N.y,N.z);
				printf("v_perp = %f \n",v_perp);
			}
		}
	END_FOR
	printf("\n there are %u points with net flow into liquid\n\n",n1);
	printf("\n there are %u solid points need recompute phi \n\n",n2);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if( voxel.InSolid(i,j,k) ){
			if(i < DimX-1){
				if( voxel.InSolid(i+1,j,k) )
					u1[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
				else
					u1[INDEX(i,j,k)] = u0[INDEX(i,j,k)];
			}
			else
				u1[INDEX(i,j,k)] = u0[INDEX(i,j,k)];
			if(j < DimY-1){
				if( voxel.InSolid(i,j+1,k) )
					v1[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
				else
					v1[INDEX(i,j,k)] = v0[INDEX(i,j,k)];
			}
			else
				v1[INDEX(i,j,k)] = v0[INDEX(i,j,k)];
			if(k < DimZ-1){
				if( voxel.InSolid(i,j,k+1) )
					w1[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
				else
					w1[INDEX(i,j,k)] = w0[INDEX(i,j,k)];
			}
			else
				w1[INDEX(i,j,k)] = w0[INDEX(i,j,k)];

			if( i > 0 && !voxel.InSolid(i-1,j,k) )
				u1[INDEX(i-1,j,k)] =  u0[INDEX(i,j,k)];;

			if( j > 0 && !voxel.InSolid(i,j-1,k) )
				v1[INDEX(i,j-1,k)] =  v0[INDEX(i,j,k)];;

			if( k > 0 && !voxel.InSolid(i,j,k-1) )
				w1[INDEX(i,j,k-1)] =  w0[INDEX(i,j,k)];;

		}
	END_FOR
	printf("after doing adjustment point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
			I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );

}

void VectorField3D::SetSolidBoundaryForAdvection(const Voxel &voxel){

	float delta = voxel.VoxelDelta();
	SetEqual(u, u1);
	SetEqual(v, v1);
	SetEqual(w, w1);
#ifdef SPMD
#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, u0, 1, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, v0, 2, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, w0, 3, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_u_obj, u1, u0, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_v_obj, v1, v0, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_w_obj, w1, w0, 3, delta);
#else
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, 1, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, 2, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, 3, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_u_obj, u1, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_v_obj, v1, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_w_obj, w1, 3, delta);
#endif
	SetEqual(u1, u0);
	SetEqual(v1, v0);
	SetEqual(w1, w0);
	printf("before doing adjustment point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
			I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );
//	memset(needRecomputePhi, 0, DimX*DimY*DimZ);
	SetZero(needRecomputePhi);

//	if(voxel.InSolid(I,J,K))
//		printf("voxel in solid at (%d, %d, %d) \n", I, J, K);
//	else
//		printf("voxel not in solid at (%d, %d, %d) \n", I, J, K);
//	printf(" after ExtrapolateOneVelocityIntoObject in VectorField3D::SetSolidBoundary()\n");
	//apply constraints on extrapolated velocity
	u_int n1 = 0, n2 = 0, n3 = 0;
	FOR_EACH_CELL
		if( voxel.InSolid(i,j,k) && !voxel.InMovingSolid(i,j,k) ){
//		if( voxel.InSolid(i,j,k) ){
//		if(valid[INDEX(i,j,k)]){
			if(i == I && j == J && k == K){
				printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
//			Vector N = ObjectNormalAt(delta, i, j, k);
//			Vector N = ObjectNormalAt(voxel, i, j, k);
			Vector N = GeometricNormalAt(voxel, i, j, k);
			float u_c, v_c, w_c;
			if(i > 0)
				u_c = 0.5f * ( u1[INDEX(i,j,k)] + u1[INDEX(i-1,j,k)] );
			else
				u_c = u1[INDEX(i,j,k)];
			if(j > 0)
				v_c = 0.5f * ( v1[INDEX(i,j,k)] + v1[INDEX(i,j-1,k)] );
			else
				v_c = v1[INDEX(i,j,k)];
			if(k > 0)
				w_c = 0.5f * ( w1[INDEX(i,j,k)] + w1[INDEX(i,j,k-1)] );
			else
				w_c = w1[INDEX(i,j,k)];
			if(i == I && j == J && k == K){
				printf(" after interpolate into object u = %f, v = %f, w= %f \n",
						u_c, v_c, w_c);
			}
			Vector vel(u_c, v_c, w_c);
			Vector v_para;
			float v_perp;
			v_para = vel - Dot(vel, N) * N;
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj)
				SetConstraintVelocity(voxel,obj_vel,i,j,k);
			v_perp = Dot(vel-obj_vel, N);

			float vel_magnitude = (vel-obj_vel).Length();
			if(i == I && j == J && k == K){
				printf("object velocity = (%f, %f,  %f), v_perp = %f, vel_mag = %f\n",
						obj_vel.x, obj_vel.y, obj_vel.z, v_perp, vel_magnitude);
			}
			if(v_perp > 0.1f*vel_magnitude){
				needRecomputePhi[INDEX(i,j,k)] = 1;
				++n2;
				if(phi_c[INDEX(i,j,k)] < 0.f)
					++n3;
			}
			if(v_perp < E_EPSIL)
				v_perp = Dot(obj_vel, N);
			else{
//				TentativeGridPoint tp(0.f, i, j, k);
//				if(IsDynamicBorderCell(voxel,i,j,k) || voxel.CloseToAir(tp)){
//					v_perp = Dot(vel, N);
//					++n1;
//				}
//				else
//					v_perp = Dot(obj_vel, N);
				++n1;
				v_perp = Dot(vel, N);
//				printf("at (%d, %d, %d) object velocity = (%f, %f,  %f), v_perp = %f, vel_mag = %f\n",
//						i, j, k, obj_vel.x, obj_vel.y, obj_vel.z, v_perp, vel_magnitude);
//				v_perp = Dot(obj_vel, N);
//				valid[INDEX(i,j,k)] = 2;
//				v_perp = Dot(obj_vel, N);
//				printf("(at (%d, %d, %d) velocity = (%f, %f,  %f), v_perp = %f, vel_mag = %f\n",
//						i,j,k, vel.x, vel.y, vel.z, v_perp, vel_magnitude);
			}

			Vector v_comb = v_para + v_perp * N;
//			Vector v_comb = v_para;

			u_c = Dot(v_comb, Vector(1,0,0));
			v_c = Dot(v_comb, Vector(0,1,0));
			w_c = Dot(v_comb, Vector(0,0,1));
			u0[INDEX(i,j,k)] = u_c;
			v0[INDEX(i,j,k)] = v_c;
			w0[INDEX(i,j,k)] = w_c;
//			if(i > 0 && !voxel.InSolid(i-1,j,k))
//				u[INDEX(i-1,j,k)] = u_c;
//			if(j > 0 && !voxel.InSolid(i,j-1,k) )
//				v[INDEX(i,j-1,k)] = v_c;
//			if(k > 0 && !voxel.InSolid(i,j,k-1))
//				w[INDEX(i,j,k-1)] = w_c;
			if(i == I && j == J && k == K){
//				printf("(SetSolidBoundary): u[(%d,%d,%d)] = %f, u[(%d,%d,%d)] = %f\n",
//						i, j, k+1, u[INDEX(i,j,k+1)], i,j,k,u[INDEX(i,j,k)]);
//				printf("(SetSolidBoundary): v[(%d,%d,%d)] = %f, v[(%d,%d,%d)] = %f\n",
//						i, j, k+1, v[INDEX(i,j,k+1)], i,j,k,v[INDEX(i,j,k)]);
//				printf("(SetSolidBoundary): w[(%d,%d,%d)] = %f, w[(%d,%d,%d)] = %f\n",
//						i, j, k+1, w[INDEX(i,j,k+1)], i,j,k,w[INDEX(i,j,k)]);
				printf("v_comb = (%f, %f, %f) \n",v_comb.x,v_comb.y,v_comb.z);
				printf("v_para = (%f, %f, %f) \n",v_para.x,v_para.y,v_para.z);
				printf("N = (%f, %f, %f) with length = %f \n",N.x, N.y, N.z, N.Length());
				printf("v_perp = %f \n",v_perp);
			}
		}
	END_FOR
	printf("\n there are %u points with net flow into liquid\n\n",n1);
	printf("\n there are %u solid points need recompute phi and %u points with phi < 0 \n\n",n2, n3);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if( voxel.InSolid(i,j,k) && !voxel.InMovingSolid(i,j,k) ){
//		if( voxel.InSolid(i,j,k) ){
			if(i < DimX-1){
				if( voxel.InSolid(i+1,j,k) && !voxel.InMovingSolid(i+1,j,k) )
//				if( voxel.InSolid(i+1,j,k) )
					u1[INDEX(i,j,k)] = 0.5f * (u0[INDEX(i,j,k)] + u0[INDEX(i+1,j,k)]);
				else if(voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
					u1[INDEX(i,j,k)] = u0[INDEX(i,j,k)];
			}
			else
				u1[INDEX(i,j,k)] = u0[INDEX(i,j,k)];
			if(j < DimY-1){
				if(voxel.InSolid(i,j+1,k) && !voxel.InMovingSolid(i,j+1,k))
//				if( voxel.InSolid(i,j+1,k) )
					v1[INDEX(i,j,k)] = 0.5f * (v0[INDEX(i,j,k)] + v0[INDEX(i,j+1,k)]);
				else if(voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
					v1[INDEX(i,j,k)] = v0[INDEX(i,j,k)];
			}
			else
				v1[INDEX(i,j,k)] = v0[INDEX(i,j,k)];
			if(k < DimZ-1){
				if(voxel.InSolid(i,j,k+1) && !voxel.InMovingSolid(i,j,k+1))
//				if( voxel.InSolid(i,j,k+1) )
					w1[INDEX(i,j,k)] = 0.5f * (w0[INDEX(i,j,k)] + w0[INDEX(i,j,k+1)]);
				else if(voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
					w1[INDEX(i,j,k)] = w0[INDEX(i,j,k)];
			}
			else
				w1[INDEX(i,j,k)] = w0[INDEX(i,j,k)];

//			if(i > 0 && !(voxel.InSolid(i-1,j,k) && !voxel.InMovingSolid(i-1,j,k)))
			if(i > 0 && voxel.InLiquid(i-1,j,k) && !voxel.InSurface(i-1,j,k))
				u1[INDEX(i-1,j,k)] =  u0[INDEX(i,j,k)];;
//			if(j > 0 && !(voxel.InSolid(i,j-1,k) && !voxel.InMovingSolid(i,j-1,k)))
			if(j > 0 && voxel.InLiquid(i,j-1,k) && !voxel.InSurface(i,j-1,k))
				v1[INDEX(i,j-1,k)] =  v0[INDEX(i,j,k)];;
//			if(k > 0 && !(voxel.InSolid(i,j,k-1) && !voxel.InMovingSolid(i,j,k-1)))
			if(k > 0 && voxel.InLiquid(i,j,k-1) && !voxel.InSurface(i,j,k-1))
				w1[INDEX(i,j,k-1)] =  w0[INDEX(i,j,k)];;

		}
	END_FOR
	printf("after doing adjustment point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
			I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );

}

void VectorField3D::SetSolidBoundaryForAdvection(const Voxel &voxel,
		map<u_int, KnownPoint> &knownPointU,
		map<u_int, KnownPoint> &knownPointV,
		map<u_int, KnownPoint> &knownPointW){

	float delta = voxel.VoxelDelta();
	map<u_int, KnownPoint>::iterator foundpos;
	SetEqual(u, u1);
	SetEqual(v, v1);
	SetEqual(w, w1);
#ifdef SPMD
#ifdef TBB
	tbb::task_scheduler_init init(THREADS);
#endif
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, u0, 1, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, v0, 2, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, w0, 3, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_u_obj, u1, u0, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_v_obj, v1, v0, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_w_obj, w1, w0, 3, delta);
#else
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointU, phi_u_obj, u, 1, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointV, phi_v_obj, v, 2, delta);
//	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, knownPointW, phi_w_obj, w, 3, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_u_obj, u1, 1, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_v_obj, v1, 2, delta);
	ExtrapolateOneVelocityIntoObjectForAdvection(voxel, phi_w_obj, w1, 3, delta);
#endif
	SetEqual(u1, u0);
	SetEqual(v1, v0);
	SetEqual(w1, w0);

	printf("before doing adjustment point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
				I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );


	FOR_EACH_CELL

//		if( (phi_u_obj[INDEX(i,j,k)] >= 0.f && phi_u_obj[INDEX(i,j,k)] < LS_ADV_LIMIT) ||
//			(!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i+1,j,k)) ||
//			(voxel.InSolid(i,j,k) && i < DimX-1 && !voxel.InSolid(i+1,j,k) && phi_c[INDEX(i+1,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
		bool mod_u = false;
	    // here we only modify velocities on those grid points which border
	    // a solid cell and a liquid cell to make them consistent with the
	    // free-slip boundary condition, as in Houston et al. 03

	    // for those grid points which border a solid cell and an air cell,
	    // we just use extrapolated velocities to advect levelset and particles
	    // doing this way, particles will collide with solid objects and get absorbed
	    // levelsets will also have consistent velocities during advection stage.
		if( phi_u_obj[INDEX(i,j,k)]+E_EPSIL >= 0.f ){
			if(phi_u_obj[INDEX(i,j,k)] < LS_ADV_LIMIT){
//				if(phi_u[INDEX(i,j,k)] <= 0.f){
					if(voxel.InSolid(i,j,k)){
						if(i < DimX-1){
							if(!voxel.InAir(i+1,j,k) && !voxel.InSurface(i+1,j,k))
								mod_u = true;
						}
						else
							mod_u = true;
					}
					else{
						if(!voxel.InAir(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i+1,j,k))
							mod_u = true;
					}
//				}
//				else{
//					if(voxel.InSolid(i,j,k)){
//						if(i < DimX-1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
//							mod_u = true;
//					}
//					else{
//						if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i+1,j,k))
//							mod_u = true;
//					}
//				}
//				if(i < DimX-1){
//					if(voxel.InSolid(i,j,k) && voxel.InSolid(i+1,j,k))
//						mod_u = true;
//				}
//				else{
//					if(voxel.InSolid(i,j,k))
//						mod_u = true;
//				}
			}
		}
		else{
			if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i+1,j,k)) ||
			    (voxel.InSolid(i,j,k) && i < DimX-1 && voxel.InLiquid(i+1,j,k)
			    		&& !voxel.InSurface(i+1,j,k)) )
				mod_u = true;
		}
		if(mod_u){
			Vector vel;
			vel.x = u0[INDEX(i,j,k)];
			if( j > 0 && i < DimX-1)
				vel.y = 0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i+1,j,k)]+v0[INDEX(i+1,j-1,k)]);
			else
				vel.y = v0[INDEX(i,j,k)];
			if( k > 0 && i < DimX-1)
				vel.z = 0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i+1,j,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i+1,j,k-1)]);
			else
				vel.z = w0[INDEX(i,j,k)];
			Point p = voxel.VelPosition(1, i, j, k);
			Vector N = ObjectNormalAt(voxel, p);
			Vector obj_vel(0.f,0.f,0.f);
			float alpha = 0.f;
			container->Friction(voxel, 1, i, j, k, &alpha);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,1);
				ConstraintFriction(voxel, 1, i, j, k, &alpha);
			}
			float v_perp = Dot(vel-obj_vel, N);
			if(v_perp < 0.f)
				v_perp = Dot(obj_vel, N);
			else
//				v_perp = Dot(vel, N);
				v_perp = Dot(obj_vel, N);
			Vector v_comb = (1 - alpha) * vel + alpha * obj_vel;
			Vector v_para = v_comb - Dot(v_comb, N) * N;
			Vector v_final = v_para + v_perp * N;
			u1[INDEX(i,j,k)] = v_final.x;
#ifdef AIR_DIV_FREE
//			u1[INDEX(i,j,k)] = v_final.x;
#endif
		}

//		if( (phi_v_obj[INDEX(i,j,k)] >= 0.f && phi_v_obj[INDEX(i,j,k)] < LS_ADV_LIMIT) ||
//			(!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i,j+1,k)) ||
//			(voxel.InSolid(i,j,k) && j < DimY-1 && !voxel.InSolid(i,j+1,k) && phi_c[INDEX(i,j+1,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
		bool mod_v = false;
		if( phi_v_obj[INDEX(i,j,k)]+E_EPSIL >= 0.f ){
			if(phi_v_obj[INDEX(i,j,k)] < LS_ADV_LIMIT){
//				if(phi_v[INDEX(i,j,k)] <= 0.f){
					if(voxel.InSolid(i,j,k)){
						if(j < DimY-1){
							if(!voxel.InAir(i,j+1,k) && !voxel.InSurface(i,j+1,k))
								mod_v = true;
						}
						else
							mod_v = true;
					}
					else{
						if(!voxel.InAir(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j+1,k))
							mod_v = true;
					}
//				}
//				else{
//					if(voxel.InSolid(i,j,k)){
//						if(j < DimY-1 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
//							mod_v = true;
//					}
//					else{
//						if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j+1,k))
//							mod_v = true;
//					}
//				}
//				if(j < DimY-1){
//					if(voxel.InSolid(i,j,k) && voxel.InSolid(i,j+1,k))
//						mod_v = true;
//				}
//				else{
//					if(voxel.InSolid(i,j,k))
//						mod_v = true;
//				}
			}
		}
		else{
			if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j+1,k)) ||
				(voxel.InSolid(i,j,k) && j < DimY-1 && voxel.InLiquid(i,j+1,k)
						&& !voxel.InSurface(i,j+1,k)) )
				mod_v = true;
		}
		if(mod_v){
			Vector vel;
			if( i > 0 && j < DimY-1)
				vel.x = 0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j+1,k)]+u0[INDEX(i-1,j+1,k)]);
			else
				vel.x = u0[INDEX(i,j,k)];
			vel.y = v0[INDEX(i,j,k)];
			if( k > 0 && j < DimY-1)
				vel.z = 0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i,j+1,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i,j+1,k-1)]);
			else
				vel.z = w0[INDEX(i,j,k)];
			Point p = voxel.VelPosition(2, i, j, k);
			Vector N = ObjectNormalAt(voxel, p);
			Vector obj_vel(0.f,0.f,0.f);
			float alpha = 0.f;
			container->Friction(voxel, 2, i, j, k, &alpha);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,2);
				ConstraintFriction(voxel, 2, i, j, k, &alpha);
			}
			float v_perp = Dot(vel-obj_vel, N);
			if(v_perp < 0.f)
				v_perp = Dot(obj_vel, N);
			else
//				v_perp = Dot(vel, N);
				v_perp = Dot(obj_vel, N);
			Vector v_comb = (1 - alpha) * vel + alpha * obj_vel;
			Vector v_para = v_comb - Dot(v_comb, N) * N;
			Vector v_final = v_para + v_perp * N;
			v1[INDEX(i,j,k)] = v_final.y;
#ifdef AIR_DIV_FREE
//			v1[INDEX(i,j,k)] = v_final.y;
#endif
		}

//		if( (phi_w_obj[INDEX(i,j,k)] >= 0.f && phi_w_obj[INDEX(i,j,k)] < LS_ADV_LIMIT) ||
//			(!voxel.InSolid(i,j,k) && phi_c[INDEX(i,j,k)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT && voxel.InSolid(i,j,k+1)) ||
//			(voxel.InSolid(i,j,k) && k < DimZ-1 && !voxel.InSolid(i,j,k+1) && phi_c[INDEX(i,j,k+1)]+E_EPSIL < EXTRAPOLATE_VEL_LIMIT) ){
		bool mod_w = false;
		if( phi_w_obj[INDEX(i,j,k)]+E_EPSIL >= 0.f ){
			if(phi_w_obj[INDEX(i,j,k)] < LS_ADV_LIMIT){
//				if(phi_w[INDEX(i,j,k)] <= 0.f){
					if(voxel.InSolid(i,j,k)){
						if(k < DimZ-1){
							if(!voxel.InAir(i,j,k+1) && !voxel.InSurface(i,j,k+1))
								mod_w = true;
						}
						else
							mod_w = true;
					}
					else{
						if(!voxel.InAir(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j,k+1))
							mod_w = true;
					}
//				}
//				else{
//					if(voxel.InSolid(i,j,k)){
//						if(k < DimZ-1 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
//							mod_w = true;
//					}
//					else{
//						if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j,k+1))
//							mod_w = true;
//					}
//				}
//				if(k < DimZ-1){
//					if(voxel.InSolid(i,j,k) && voxel.InSolid(i,j,k+1))
//						mod_w = true;
//				}
//				else{
//					if(voxel.InSolid(i,j,k))
//						mod_w = true;
//				}
			}
		}
		else{
			if( (voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) && voxel.InSolid(i,j,k+1)) ||
				(voxel.InSolid(i,j,k) && k < DimZ-1 && voxel.InLiquid(i,j,k+1)
						&& !voxel.InSurface(i,j,k+1)) )
				mod_w = true;
		}
		if(mod_w){
			Vector vel;
			if( i > 0 && k < DimZ-1)
				vel.x = 0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j,k+1)]+u0[INDEX(i-1,j,k+1)]);
			else
				vel.x = u0[INDEX(i,j,k)];
			if( j > 0 && k < DimZ-1)
				vel.y = 0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i,j,k+1)]+v0[INDEX(i,j-1,k+1)]);
			else
				vel.y = v0[INDEX(i,j,k)];
			vel.z = w0[INDEX(i,j,k)];
			if(i == I && j == J && k == K){
				printf("vel = (%f, %f, %f) \n",vel.x, vel.y, vel.z);
				printf("u0 = %f, u0[i-1] = %f, u0[k+1] = %f, u0[i-1,k+1] = %f \n ",
						u0[INDEX(i,j,k)], v0[INDEX(i,j-1,k)], u0[INDEX(i,j,k+1)], u0[INDEX(i-1,j,k+1)]);
				printf("v0 = %f, v0[j-1] = %f, v0[k+1] = %f, v0[j-1,k+1] = %f \n ",
						v0[INDEX(i,j,k)], v0[INDEX(i,j-1,k)], v0[INDEX(i,j,k+1)], v0[INDEX(i,j-1,k+1)]);

			}
			Point p = voxel.VelPosition(3, i, j, k);
			Vector N = ObjectNormalAt(voxel, p);
			Vector obj_vel(0.f,0.f,0.f);
			float alpha = 0.f;
			container->Friction(voxel, 3, i, j, k, &alpha);
			if(hasMovingObj){
				SetConstraintVelocity(voxel,obj_vel,i,j,k,3);
				ConstraintFriction(voxel, 3, i, j, k, &alpha);
			}
			float v_perp = Dot(vel-obj_vel, N);
			if(v_perp < 0.f)
				v_perp = Dot(obj_vel, N);
			else
//				v_perp = Dot(vel, N);
				v_perp = Dot(obj_vel, N);
			Vector v_comb = (1 - alpha) * vel + alpha * obj_vel;
			Vector v_para = v_comb - Dot(v_comb, N) * N;
			Vector v_final = v_para + v_perp * N;
			w1[INDEX(i,j,k)] = v_final.z;
#ifdef AIR_DIV_FREE
//			w1[INDEX(i,j,k)] = v_final.z;
#endif
		}
	END_FOR
	printf("after doing adjustment point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
					I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );

}

/*void VectorField3D::SetSolidBoundaryForAdvection(const Voxel &voxel){

	float delta = voxel.VoxelDelta();

	// first, extrapolate velocity into object using object phi
	ExtrapolateOneVelocityIntoObject(voxel, phi_u_obj, u, 1, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_v_obj, v, 2, delta);
	ExtrapolateOneVelocityIntoObject(voxel, phi_w_obj, w, 3, delta);
//	memset(needRecomputePhi, 0, DimX*DimY*DimZ);
	SetZero(needRecomputePhi);
//	printf(" after ExtrapolateOneVelocityIntoObject in VectorField3D::SetSolidBoundary()\n");
	//apply constraints on extrapolated velocity
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			if(i == I && j == J && k == K){
				printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
			TentativeGridPoint p(0.f, i,j,k);
			if(voxel.CloseToSolid(p)){
	//			Vector N = ObjectNormalAt(delta, i, j, k);
				Vector N = ObjectNormalAt(voxel, i, j, k);
				float u_c, v_c, w_c;
				if(i == 0)
					u_c = u[INDEX(i,j,k)];
				else
					u_c = 0.5f*(u[INDEX(i-1,j,k)]+u[INDEX(i,j,k)]);
				if(j == 0)
					v_c = v[INDEX(i,j,k)];
				else
					v_c = 0.5f*(v[INDEX(i,j-1,k)]+v[INDEX(i,j,k)]);
				if(k == 0)
					w_c = w[INDEX(i,j,k)];
				else
					w_c = 0.5f*(w[INDEX(i,j,k-1)]+w[INDEX(i,j,k)]);
				if(i == I && j == J && k == K){
					printf(" after interpolate into object u = %f, v = %f, w= %f \n",
							u_c, v_c, w_c);
				}
				Vector vel(u_c, v_c, w_c);
				Vector v_para;
				float v_perp;
				v_para = vel - Dot(vel, N) * N;
				Vector obj_vel(0.f,0.f,0.f);
				float vel_magnitude = 0.f;
				if(hasMovingObj){
					SetConstraintVelocity(obj_vel,i,j,k);
					if(i == I && j == J && k == K){
						printf("object velocity = (%f, %f,  %f) \n",
								obj_vel.x, obj_vel.y, obj_vel.z);
					}
					v_perp = Dot(vel-obj_vel, N);

					vel_magnitude = (vel-obj_vel).Length();

//					if(v_perp > 0.1f*vel_magnitude && phi_c[INDEX(i,j,k)] < 0.f){
//						needRecomputePhi[INDEX(i,j,k)] = 1;
//					}
					if(v_perp < 0.f)
						v_perp = Dot(obj_vel, N);
					else
						v_perp = Dot(vel, N);
				}
				else{
					v_perp = Dot(vel, N);
					vel_magnitude = (vel-obj_vel).Length();
//					if(v_perp > 0.1f*vel_magnitude && phi_c[INDEX(i,j,k)] < 0.f){
//							needRecomputePhi[INDEX(i,j,k)] = 1;
//					}
					if(v_perp != 0.f)
						v_perp = 0.f;
				}
				Vector v_comb = v_para + v_perp * N;
				u_c = Dot(v_comb, Vector(1,0,0));
				v_c = Dot(v_comb, Vector(0,1,0));
				w_c = Dot(v_comb, Vector(0,0,1));
				vector<TentativeGridPoint> neighbors;
				p.Neigbor(neighbors, DimX, DimY, DimZ);
				for(int m=0; m<neighbors.size(); ++m){
					if(voxel.InSolid(neighbors[m])){
						if(p.LeftNeighbor(neighbors[m])){
							u[INDEX(i-1,j,k)] = u_c;
						}
						if(p.RightNeighbor(neighbors[m])){
							u[INDEX(i,j,k)] = u_c;
						}
						if(p.FrontNeighbor(neighbors[m])){
							v[INDEX(i,j,k)] = v_c;
						}
						if(p.BackNeighbor(neighbors[m])){
							v[INDEX(i,j-1,k)] = v_c;
						}
						if(p.TopNeighbor(neighbors[m])){
							w[INDEX(i,j,k)] = w_c;
						}
						if(p.BottomNeighbor(neighbors[m])){
							w[INDEX(i,j,k-1)] = w_c;
						}
						if(v_perp > 0.1f*vel_magnitude &&
							phi_c[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] < 0.f){
							needRecomputePhi[INDEX(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)] = 1;
						}
					}
				}

				if(i == I && j == J && k == K){
					printf("(SetSolidBoundary): u[(%d,%d,%d)] = %f, u[(%d,%d,%d)] = %f\n",
							i, j, k+1, u[INDEX(i,j,k+1)], i,j,k,u[INDEX(i,j,k)]);
					printf("(SetSolidBoundary): v[(%d,%d,%d)] = %f, v[(%d,%d,%d)] = %f\n",
							i, j, k+1, v[INDEX(i,j,k+1)], i,j,k,v[INDEX(i,j,k)]);
					printf("(SetSolidBoundary): w[(%d,%d,%d)] = %f, w[(%d,%d,%d)] = %f\n",
							i, j, k+1, w[INDEX(i,j,k+1)], i,j,k,w[INDEX(i,j,k)]);
					printf("v_comb = (%f, %f, %f) \n",v_comb.x,v_comb.y,v_comb.z);
					printf("v_para = (%f, %f, %f) \n",v_para.x,v_para.y,v_para.z);
					printf("N = (%f, %f, %f) \n",N.x,N.y,N.z);
					printf("v_perp = %f \n",v_perp);
				}
			}
		}
		else{
			if(i == I && j == J && k == K){
				printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
//			Vector N = ObjectNormalAt(delta, i, j, k);
			Vector N = ObjectNormalAt(voxel, i, j, k);
			float u_c, v_c, w_c;
//			if(i == 0)
				u_c = u[INDEX(i,j,k)];
//			else
//				u_c = 0.5f*(u[INDEX(i-1,j,k)]+u[INDEX(i,j,k)]);
//			if(j == 0)
				v_c = v[INDEX(i,j,k)];
//			else
//				v_c = 0.5f*(v[INDEX(i,j-1,k)]+v[INDEX(i,j,k)]);
//			if(k == 0)
				w_c = w[INDEX(i,j,k)];
//			else
//				w_c = 0.5f*(w[INDEX(i,j,k-1)]+w[INDEX(i,j,k)]);
			if(i == I && j == J && k == K){
				printf(" after interpolate into object u = %f, v = %f, w= %f \n",
						u_c, v_c, w_c);
			}
			Vector vel(u_c, v_c, w_c);
			Vector v_para;
			float v_perp;
			v_para = vel - Dot(vel, N) * N;
			Vector obj_vel(0.f,0.f,0.f);
			if(hasMovingObj){
				SetConstraintVelocity(obj_vel,i,j,k);
				if(i == I && j == J && k == K){
					printf("object velocity = (%f, %f,  %f) \n",
							obj_vel.x, obj_vel.y, obj_vel.z);
				}
				v_perp = Dot(vel-obj_vel, N);

				float vel_magnitude = (vel-obj_vel).Length();

				if(v_perp > 0.1f*vel_magnitude && phi_c[INDEX(i,j,k)] < 0.f){
					needRecomputePhi[INDEX(i,j,k)] = 1;
				}
				if(v_perp < 0.f)
					v_perp = Dot(obj_vel, N);
				else
					v_perp = Dot(vel, N);

			}
			else{
				v_perp = Dot(vel, N);
				float vel_magnitude = (vel-obj_vel).Length();
				if(v_perp > 0.1f*vel_magnitude && phi_c[INDEX(i,j,k)] < 0.f){
						needRecomputePhi[INDEX(i,j,k)] = 1;
				}
				if(v_perp != 0.f)
					v_perp = 0.f;
			}
			Vector v_comb = v_para + v_perp * N;
			u_c = Dot(v_comb, Vector(1,0,0));
			v_c = Dot(v_comb, Vector(0,1,0));
			w_c = Dot(v_comb, Vector(0,0,1));
			u[INDEX(i,j,k)] = u_c;
			v[INDEX(i,j,k)] = v_c;
			w[INDEX(i,j,k)] = w_c;
			if(i > 0 && !voxel.InSolid(i-1,j,k))
				u[INDEX(i-1,j,k)] = u_c;
			if(j > 0 && !voxel.InSolid(i,j-1,k) )
				v[INDEX(i,j-1,k)] = v_c;
			if(k > 0 && !voxel.InSolid(i,j,k-1))
				w[INDEX(i,j,k-1)] = w_c;
			if(i == I && j == J && k == K){
				printf("(SetSolidBoundary): u[(%d,%d,%d)] = %f, u[(%d,%d,%d)] = %f\n",
						i, j, k+1, u[INDEX(i,j,k+1)], i,j,k,u[INDEX(i,j,k)]);
				printf("(SetSolidBoundary): v[(%d,%d,%d)] = %f, v[(%d,%d,%d)] = %f\n",
						i, j, k+1, v[INDEX(i,j,k+1)], i,j,k,v[INDEX(i,j,k)]);
				printf("(SetSolidBoundary): w[(%d,%d,%d)] = %f, w[(%d,%d,%d)] = %f\n",
						i, j, k+1, w[INDEX(i,j,k+1)], i,j,k,w[INDEX(i,j,k)]);
				printf("v_comb = (%f, %f, %f) \n",v_comb.x,v_comb.y,v_comb.z);
				printf("v_para = (%f, %f, %f) \n",v_para.x,v_para.y,v_para.z);
				printf("N = (%f, %f, %f) \n",N.x,N.y,N.z);
				printf("v_perp = %f \n",v_perp);
			}
		}
	END_FOR

}*/

void VectorField3D::SetSolidBoundaryNoExtrapolate(const Voxel &voxel, int b, float *x){

	float delta = voxel.VoxelDelta();


	//apply constraints on extrapolated velocity
	FOR_EACH_CELL

		if(voxel.InSolid(i,j,k)){
			Vector obj_vel(0.f, 0.f, 0.f);
			if(hasMovingObj)
				SetConstraintVelocity(voxel,obj_vel,i,j,k);
			TentativeGridPoint p(0.f,i,j,k);
			vector<TentativeGridPoint> neighbors;
			p.Neigbor(neighbors, DimX, DimY, DimZ);
			for(int m=0; m<neighbors.size(); m++){
				if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk) ||
				   voxel.InAir(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk) ||
				   voxel.InSurface(neighbors[m].ii,neighbors[m].jj,neighbors[m].kk)){
					if(p.LeftNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i-1,j,k)];
						else if ( b == 1 ){
//							if(x[INDEX(i,j,k)] > 0.f)
								x[INDEX(i,j,k)] = obj_vel.x;
//							if(x[INDEX(i-1,j,k)] > 0.f)
								x[INDEX(i-1,j,k)] = obj_vel.x;
						}
					}
					if(p.RightNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i+1,j,k)];
						else if ( b == 1 ){
//							if(x[INDEX(i,j,k)] < 0.f)
								x[INDEX(i,j,k)] = obj_vel.x;
						}
					}
					if(p.BackNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j-1,k)];
						else if ( b == 2 ){
//							if(x[INDEX(i,j,k)] > 0.f)
								x[INDEX(i,j,k)] = obj_vel.y;
//							if(x[INDEX(i,j-1,k)] > 0.f)
								x[INDEX(i,j-1,k)] = obj_vel.y;
						}
					}
					if(p.FrontNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j+1,k)];
						else if ( b == 2 ){
//							if(x[INDEX(i,j,k)] < 0.f)
								x[INDEX(i,j,k)] = obj_vel.y;
						}
					}
					if(p.BottomNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j,k-1)];
						else if ( b == 3 ){
//							if(x[INDEX(i,j,k)] > 0.f)
								x[INDEX(i,j,k)] = obj_vel.z;
//							if(x[INDEX(i,j,k-1)] > 0.f)
								x[INDEX(i,j,k-1)] = obj_vel.z;
						}
					}
					if(p.TopNeighbor(neighbors[m])){
						if( b == 0 )
							x[INDEX(i,j,k)] = x[INDEX(i,j,k+1)];
						else if ( b == 3 ){
//							if(x[INDEX(i,j,k)] < 0.f)
								x[INDEX(i,j,k)] = obj_vel.z;
						}
					}
				}
			}
		}

	END_FOR

}

void VectorField3D::SetSolidBoundaryNoExtrapolate(const Voxel &voxel){

	float delta = voxel.VoxelDelta();


	//apply constraints on extrapolated velocity
	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k)){
			Vector N = ObjectNormalAt(delta, i, j, k);
			float u_c, v_c, w_c;
			if(i == 0)
				u_c = u[INDEX(i,j,k)];
			else
				u_c = 0.5f*(u[INDEX(i-1,j,k)]+u[INDEX(i,j,k)]);
			if(j == 0)
				v_c = v[INDEX(i,j,k)];
			else
				v_c = 0.5f*(v[INDEX(i,j-1,k)]+v[INDEX(i,j,k)]);
			if(k == 0)
				w_c = w[INDEX(i,j,k)];
			else
				w_c = 0.5f*(w[INDEX(i,j,k-1)]+w[INDEX(i,j,k)]);
			Vector vel(u_c, v_c, w_c);
			Vector v_para;
			float v_perp;
			v_para = vel - Dot(vel, N) * N;
			v_perp = Dot(vel, N);
			if(v_perp != 0.f)
				v_perp = 0.f;
			Vector v_comb = v_para + v_perp*N;
			u_c = Dot(v_comb, Vector(1,0,0));
			v_c = Dot(v_comb, Vector(0,1,0));
			w_c = Dot(v_comb, Vector(0,0,1));
//			u[INDEX(i,j,k)] = u_c;
//			v[INDEX(i,j,k)] = v_c;
//			w[INDEX(i,j,k)] = w_c;
			if(i > 0 && voxel.InAir(i-1,j,k))
				u[INDEX(i-1,j,k)] = u_c;
			if(j > 0 && voxel.InAir(i,j-1,k))
				v[INDEX(i,j-1,k)] = v_c;
			if(k > 0 && voxel.InAir(i,j,k-1))
				w[INDEX(i,j,k-1)] = w_c;
//			if(i == I && j == J && k == K)
//				printf("(SetSolidBoundaryNoExtrapolate): u[(i-1,j,k)] = %f, u[(i,j,k)] = %f\n",
//						u[INDEX(i-1,j,k)], u[INDEX(i,j,k)]);
		}
	END_FOR

}

Vector VectorField3D::ObjectNormalAt(float delta, int i, int j, int k) const{
	float gx, gy, gz;
//	TentativeGridPoint p(0.f,i,j,k);
//	vector<TentativeGridPoint> neighbors;
//	p.Neigbor(neighbors, DimX, DimY, DimZ);
//	if(neighbors.size() == 6){
//		gx = (phi_c_obj[INDEX(i+1,j,k)] - phi_c_obj[INDEX(i-1,j,k)])/(2*delta);
//		gy = (phi_c_obj[INDEX(i,j+1,k)] - phi_c_obj[INDEX(i,j-1,k)])/(2*delta);
//		gz = (phi_c_obj[INDEX(i,j,k+1)] - phi_c_obj[INDEX(i,j,k-1)])/(2*delta);
//	}
//	else{
		if(i == 0)
			gx = (phi_c_obj[INDEX(i+1,j,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
		else if(i == DimX-1)
			gx = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i-1,j,k)])/delta;
		else{
			if(phi_c_obj[INDEX(i-1,j,k)] < phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i+1,j,k)] < phi_c_obj[INDEX(i,j,k)]	){
				if(phi_c_obj[INDEX(i-1,j,k)] < phi_c_obj[INDEX(i+1,j,k)])
					gx = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i-1,j,k)])/delta;
				else
					gx = (phi_c_obj[INDEX(i+1,j,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
			}
			else if( phi_c_obj[INDEX(i-1,j,k)] > phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i+1,j,k)] > phi_c_obj[INDEX(i,j,k)]){
				if(phi_c_obj[INDEX(i-1,j,k)] < phi_c_obj[INDEX(i+1,j,k)])
					gx = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i-1,j,k)])/delta;
				else
					gx = (phi_c_obj[INDEX(i+1,j,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
//				gx = 0.f;
			}
			else if( phi_c_obj[INDEX(i-1,j,k)] == phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i+1,j,k)] == phi_c_obj[INDEX(i,j,k)])
				gx = 0.f;
			else if(phi_c_obj[INDEX(i-1,j,k)] < phi_c_obj[INDEX(i,j,k)])
				gx = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i-1,j,k)])/delta;
			else if(phi_c_obj[INDEX(i+1,j,k)] < phi_c_obj[INDEX(i,j,k)])
				gx = (phi_c_obj[INDEX(i+1,j,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
			else if(phi_c_obj[INDEX(i,j,k)] < phi_c_obj[INDEX(i+1,j,k)])
				gx = (phi_c_obj[INDEX(i+1,j,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
			else if(phi_c_obj[INDEX(i,j,k)] < phi_c_obj[INDEX(i-1,j,k)])
				gx = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i-1,j,k)])/delta;
		}

		if(j == 0)
			gy = (phi_c_obj[INDEX(i,j+1,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
		else if(j == DimY-1)
			gy = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j-1,k)])/delta;
		else{
			if(phi_c_obj[INDEX(i,j-1,k)] < phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i,j+1,k)] < phi_c_obj[INDEX(i,j,k)]	){
				if(phi_c_obj[INDEX(i,j-1,k)] < phi_c_obj[INDEX(i,j+1,k)])
					gy = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j-1,k)])/delta;
				else
					gy = (phi_c_obj[INDEX(i,j+1,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
			}
			else if( phi_c_obj[INDEX(i,j-1,k)] > phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i,j+1,k)] > phi_c_obj[INDEX(i,j,k)]){
				if(phi_c_obj[INDEX(i,j-1,k)] < phi_c_obj[INDEX(i,j+1,k)])
					gy = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j-1,k)])/delta;
				else
					gy = (phi_c_obj[INDEX(i,j+1,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
//				gy = 0.f;
			}
			else if( phi_c_obj[INDEX(i,j-1,k)] == phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i,j+1,k)] == phi_c_obj[INDEX(i,j,k)])
				gy = 0.f;
			else if(phi_c_obj[INDEX(i,j-1,k)] < phi_c_obj[INDEX(i,j,k)])
				gy = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j-1,k)])/delta;
			else if(phi_c_obj[INDEX(i,j+1,k)] < phi_c_obj[INDEX(i,j,k)])
				gy = (phi_c_obj[INDEX(i,j+1,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
			else if(phi_c_obj[INDEX(i,j,k)] < phi_c_obj[INDEX(i,j+1,k)])
				gy = (phi_c_obj[INDEX(i,j+1,k)] - phi_c_obj[INDEX(i,j,k)])/delta;
			else if(phi_c_obj[INDEX(i,j,k)] < phi_c_obj[INDEX(i,j-1,k)])
				gy = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j-1,k)])/delta;
		}
		if(k == 0)
			gz = (phi_c_obj[INDEX(i,j,k+1)] - phi_c_obj[INDEX(i,j,k)])/delta;
		else if(k == DimZ-1)
			gz = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j,k-1)])/delta;
		else{
			if(phi_c_obj[INDEX(i,j,k-1)] < phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i,j,k+1)] < phi_c_obj[INDEX(i,j,k)]	){
				if(phi_c_obj[INDEX(i,j,k-1)] < phi_c_obj[INDEX(i,j,k+1)])
					gz = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j,k-1)])/delta;
				else
					gz = (phi_c_obj[INDEX(i,j,k+1)] - phi_c_obj[INDEX(i,j,k)])/delta;
			}
			else if( phi_c_obj[INDEX(i,j,k-1)] > phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i,j,k+1)] > phi_c_obj[INDEX(i,j,k)]){
				if(phi_c_obj[INDEX(i,j,k-1)] < phi_c_obj[INDEX(i,j,k+1)])
					gz = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j,k-1)])/delta;
				else
					gz = (phi_c_obj[INDEX(i,j,k+1)] - phi_c_obj[INDEX(i,j,k)])/delta;
//				gz = 0.f;
			}
			else if( phi_c_obj[INDEX(i,j,k-1)] == phi_c_obj[INDEX(i,j,k)] &&
			   phi_c_obj[INDEX(i,j,k+1)] == phi_c_obj[INDEX(i,j,k)])
				gz = 0.f;
			else if(phi_c_obj[INDEX(i,j,k-1)] < phi_c_obj[INDEX(i,j,k)])
				gz = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j,k-1)])/delta;
			else if(phi_c_obj[INDEX(i,j,k+1)] < phi_c_obj[INDEX(i,j,k)])
				gz = (phi_c_obj[INDEX(i,j,k+1)] - phi_c_obj[INDEX(i,j,k)])/delta;
			else if(phi_c_obj[INDEX(i,j,k)] < phi_c_obj[INDEX(i,j,k+1)])
				gz = (phi_c_obj[INDEX(i,j,k+1)] - phi_c_obj[INDEX(i,j,k)])/delta;
			else if(phi_c_obj[INDEX(i,j,k)] < phi_c_obj[INDEX(i,j,k-1)])
				gz = (phi_c_obj[INDEX(i,j,k)] - phi_c_obj[INDEX(i,j,k-1)])/delta;
		}
//	}
	Vector N(gx, gy, gz);
	if(N.Length() == 0.f){
		printf("at (%d, %d, %d) object normal undefined!\n", i,j,k);
		printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
		printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
				"phi[k-1] = %f, phi[k+1] = %f \n\n",
		   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
			 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
			 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
	}
	N = Normalize(N);
	return -1*N;  // multiply by -1 because phi at object grid points is positive

}

Vector VectorField3D::ObjectNormalAt(const Voxel &voxel, int i, int j, int k) const{
	float gx, gy, gz;

	float delta = voxel.VoxelDelta();

	if(i == 0 || i == DimX-1 || j == 0 || j == DimY-1 || k == 0 || k == DimZ-1)
		return ObjectNormalAt(delta, i, j, k);


	Point p = voxel.VoxelCenterPosition(i,j,k);

	Point px0(p.x-0.5f*delta, p.y, p.z);
	Point px1(p.x+0.5f*delta, p.y, p.z);
	float phix0 = TriInterp(voxel, px0, phi_c_obj);
	float phix1 = TriInterp(voxel, px1, phi_c_obj);
	gx = (phix1 - phix0) / delta;

	Point py0(p.x, p.y-0.5f*delta, p.z);
	Point py1(p.x, p.y+0.5f*delta, p.z);
	float phiy0 = TriInterp(voxel, py0, phi_c_obj);
	float phiy1 = TriInterp(voxel, py1, phi_c_obj);
	gy = (phiy1 - phiy0) / delta;

	Point pz0(p.x, p.y, p.z-0.5f*delta);
	Point pz1(p.x, p.y, p.z+0.5f*delta);
	float phiz0 = TriInterp(voxel, pz0, phi_c_obj);
	float phiz1 = TriInterp(voxel, pz1, phi_c_obj);
	gz = (phiz1 - phiz0) / delta;

	Vector N(gx, gy, gz);
	if(N.Length() == 0.f){
		printf("at (%d, %d, %d) object normal undefined!\n", i,j,k);
		printf(" at (%d, %d, %d) phi_c = %f \n", i,j,k,phi_c_obj[INDEX(i,j,k)]);
		printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
				"phi[k-1] = %f, phi[k+1] = %f \n\n",
		   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
			 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
			 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
		MTRand mt;
		double rn1 = mt();
		double rn2 = mt();
		double rn3 = mt();
		N.x = (float)rn1;
		N.y = (float)rn2;
		N.z = (float)rn3;
	}
	N = Normalize(N);
	return -1*N;  // multiply by -1 because phi at object grid points is positive

}

Vector VectorField3D::ObjectNormalAt(const Voxel &voxel, const Point &p) const{
	float gx, gy, gz;

	float delta = voxel.VoxelDelta();

	Point px0(p.x-0.5f*delta, p.y, p.z);
	Point px1(p.x+0.5f*delta, p.y, p.z);
	float phix0 = TriInterp(voxel, px0, phi_c_obj);
	float phix1 = TriInterp(voxel, px1, phi_c_obj);
	gx = (phix1 - phix0) / delta;

	Point py0(p.x, p.y-0.5f*delta, p.z);
	Point py1(p.x, p.y+0.5f*delta, p.z);
	float phiy0 = TriInterp(voxel, py0, phi_c_obj);
	float phiy1 = TriInterp(voxel, py1, phi_c_obj);
	gy = (phiy1 - phiy0) / delta;

	Point pz0(p.x, p.y, p.z-0.5f*delta);
	Point pz1(p.x, p.y, p.z+0.5f*delta);
	float phiz0 = TriInterp(voxel, pz0, phi_c_obj);
	float phiz1 = TriInterp(voxel, pz1, phi_c_obj);
	gz = (phiz1 - phiz0) / delta;

	Vector N(gx, gy, gz);
	if(N.Length() == 0.f){
		MTRand mt;
		double rn1 = mt();
		double rn2 = mt();
		double rn3 = mt();
		N.x = (float)rn1;
		N.y = (float)rn2;
		N.z = (float)rn3;
	}
	N = Normalize(N);
	return -1*N;  // multiply by -1 because phi at object grid points is positive

}

Vector VectorField3D::GeometricNormalAt(const Voxel &voxel, int i, int j, int k) const{
//	TentativeGridPoint p(0.f, i,j,k);
//	int nbrs;
//	if(!voxel.CloseToNonSolid(p, nbrs))
//		return ObjectNormalAt(voxel, i, j, k);
//	else{
//		Vector N;
//		container->SetNormal(voxel, N, i, j, k);
//		if(hasMovingObj){
//			for(int n = 0; n < movingObjects.size(); ++n){
//				MovingObject *obj = movingObjects[n];
//				obj->SetNormal(voxel, N, i, j, k);
//			}
//		}
//		return N;
//	}
	Vector N = ObjectNormalAt(voxel, i, j, k);
	if(hasMovingObj){
		for(int n = 0; n < movingObjects.size(); ++n){
			MovingObject *obj = movingObjects[n];
			obj->SetNormal(voxel, N, i, j, k);
		}
	}
	return N;
}

void VectorField3D::RecomputeObjectLevelset(const Voxel &voxel, float *phitmp){

	u_int N1 = 0, N2 = 0, N3 = 0, N4 = 0;
	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k) && phitmp[INDEX(i,j,k)] < 0.f)
			++N1;
		if(voxel.InSolid(i,j,k) && needRecomputePhi[INDEX(i,j,k)])
			++N4;
	END_FOR
	printf("\n before RecomputeObjectLevelset there are %u solid cells needing recompute\n", N4);

	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k)){
			float phi_ori = phitmp[INDEX(i,j,k)];

			if(needRecomputePhi[INDEX(i,j,k)] && phitmp[INDEX(i,j,k)] < 0.f){
//			if(proj_vel > 0.f && phi_c[INDEX(i,j,k)] < 0.f){
				phitmp[INDEX(i,j,k)] = phi_c_obj[INDEX(i,j,k)];
				++N3;
			}

//			else{
//				if(phi_c[INDEX(i,j,k)] < 0.f && proj_vel > 0.f){
//					printf("\nat (%d, %d, %d) phi obj (%f) not reset to phi air (%f) \n",
//							i,j,k, phi_origin, phi_c_obj[INDEX(i,j,k)]);
//					printf("N = (%f, %f, %f) vel = (%f, %f, %f) \n\n",
//							N.x, N.y, N.z, u_c, v_c, w_c);
//				}
//			}
//			phi_c[INDEX(i,j,k)] = phi_c_obj[INDEX(i,j,k)];
			if(i == I && j == J && k == K){
				printf("at (%d, %d, %d) phi obj (%f) reset to phi air (%f) with phi_obj = %f\n",
									i,j,k, phi_ori, phitmp[INDEX(i,j,k)], phi_c_obj[INDEX(i,j,k)]);
//				printf("\nN = (%f, %f, %f) N.length = %f v_perp = %f vel__m = %f, vel = (%f, %f, %f)\n",
//						N.x, N.y, N.z, N.Length(),  v_perp, vel_m, u_c, v_c, w_c);
				printf("phi[i-1] = %f, phi[i+1]= %f, phi[j-1] = %f, phi[j+1] = %f,"
						"phi[k-1] = %f, phi[k+1] = %f \n\n",
   				   		phi_c_obj[INDEX(i-1,j,k)], phi_c_obj[INDEX(i+1,j,k)],
   					 		phi_c_obj[INDEX(i,j-1,k)], phi_c_obj[INDEX(i,j+1,k)],
   					 	phi_c_obj[INDEX(i,j,k-1)], phi_c_obj[INDEX(i,j,k+1)] );
			}
		}
	END_FOR
	FOR_EACH_CELL
		if(voxel.InSolid(i,j,k) && phitmp[INDEX(i,j,k)] < 0.f)
			++N2;
	END_FOR
	printf(" before RecomputeObjectLevelset there are %u solid cells having < 0 phi\n", N1);
	printf(" after RecomputeObjectLevelset there are %u solid cells having < 0 phi\n", N2);
	printf(" after RecomputeObjectLevelset there are %u solid cells resetting phi to > 0 \n\n", N3);
}



void VectorField3D::AdvectOneVelocity(const Voxel &voxel, int b, float *d, const float *d0){


	int i0, j0, k0, i1, j1, k1;
	float s0, t0, s1, t1, dt0;

	float xp, yp, zp, up, vp, wp;

	int x0, y0, z0;

	float r0, r1;

	double p, q;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;

	bool u_active, v_active, w_active;

	dt0 = dt * inv_delta;

#ifdef USE_MESH
	float f[8];
	Point pb(0.f), pi(0.f);
#endif


	FOR_EACH_CELL
		if(b == 1)
			u_active = phi_u_obj[INDEX(i,j,k)] < 0.f && ( phi_u[INDEX(i,j,k)] <= 0.f ||
						(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
						(voxel.InSolid(i,j,k) && i < DimX-1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) ||
						(voxel.InAir(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) ||
						(voxel.InSurface(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) );
		if(b == 2)
			v_active = phi_v_obj[INDEX(i,j,k)] < 0.f && ( phi_v[INDEX(i,j,k)] <= 0.f ||
						(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
						(voxel.InSolid(i,j,k) && j < DimY-1 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) ||
						(voxel.InAir(i,j,k) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) ||
						(voxel.InSurface(i,j,k) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) );
		if(b == 3)
			w_active = phi_w_obj[INDEX(i,j,k)] < 0.f && ( phi_w[INDEX(i,j,k)] <= 0.f ||
						(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
						(voxel.InSolid(i,j,k) && k < DimZ-1 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) ||
						(voxel.InAir(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) ||
						(voxel.InSurface(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) );

		if( (b == 1 && u_active) || (b == 2 && v_active) ||	(b == 3 && w_active) ) {

	        if( b == 1 ){
		        xp = i+1.f-dt0*u0[INDEX(i,j,k)];
			   	yp = j+0.5f-dt0*0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i+1,j,k)]+v0[INDEX(i+1,j-1,k)]);
			   	zp = k+0.5f-dt0*0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i+1,j,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i+1,j,k-1)]);
	        }
	        else if( b == 2 ){
	            xp = i+0.5f-dt0*0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j+1,k)]+u0[INDEX(i-1,j+1,k)]);
	            yp = j+1.f- dt0*v0[INDEX(i,j,k)];
	            zp = k+0.5f-dt0*0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i,j+1,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i,j+1,k-1)]);
	        }
	        else if( b == 3 ){
                xp = i+0.5-dt0*0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j,k+1)]+u0[INDEX(i-1,j,k+1)]);
                yp = j+0.5-dt0*0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i,j,k+1)]+v0[INDEX(i,j-1,k+1)]);
                zp = k+1.f-dt0*w0[INDEX(i,j,k)];
	        }

	        xp = xp < 0.5f      ? 0.5f      : xp;
			xp = xp > DimX-1.5f ? DimX-1.5f : xp;
			yp = yp < 0.5f      ? 0.5f      : yp;
			yp = yp > DimY-1.5f ? DimY-1.5f : yp;
			zp = zp < 0.5f      ? 0.5f      : zp;
			zp = zp > DimZ-1.5f ? DimZ-1.5f : zp;

#ifdef USE_MESH
	        pb = voxel.VelPosition(b, i, j, k);
	        pi = Point(delta*xp, delta*yp, delta*zp);
	        if(mPhysWorld->SegmentIntersectsMesh(pb, pi, &p, &q)){
	        	// clip the semi-lagrangian ray
	        	Vector pib = pi - pb;
	        	float pibl = p * pib.Length();
	        	pib = Normalize(pib);
//	        	Point pin = pb + (p-E_EPSIL)*(pi-pb);
	        	Point pin = pb + (pibl-EPSIL_C) * pib;
	        	xp = pin.x * inv_delta;
	        	yp = pin.y * inv_delta;
	        	zp = pin.z * inv_delta;
	        	pi = pin;
	        }
#endif

//	        if(  i == I && j == J && k == K && b == 1 ){
//				printf(" xp = %f, yp = %f, zp = %f \n", xp, yp, zp);
//			}

	        if( b == 1 ){
	        	i0=floor(xp)-1; if(i0 < 0) i0 = 0;
	        	i1 = i0 + 1;
	            s1 = xp - 1 - i0;
	            s0 = 1 - s1;
	            j0 = floor(yp); y0 = j0;
	            if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
	            j1=j0+1;
	            t1 = yp-y0 < 0.5f ? 0.5f+yp-y0 : yp-y0-0.5f;
	            t0 = 1-t1;
	            k0 = floor(zp); z0 = k0;
	            if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
	            k1=k0+1;
	            r1 = zp-z0 < 0.5f ? 0.5f+zp-z0 : zp-z0-0.5f;
	            r0 = 1-r1;
	        }
	        else if( b == 2 ){
	        	j0=floor(yp)-1; if(j0 < 0) j0 = 0;
	 			j1 = j0 + 1;
	 		    t1 = yp - 1 - j0;
	 		    t0 = 1 - t1;
	 		    i0 = floor(xp); x0 = i0;
	 		    if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
	 		    i1=i0+1;
	 		    s1 = xp-x0 < 0.5f ? 0.5f+xp-x0 : xp-x0-0.5f;
	 		    s0 = 1-s1;
	 		    k0 = floor(zp); z0 = k0;
	 		    if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
	 		    k1=k0+1;
	 		    r1 = zp-z0 < 0.5f ? 0.5f+zp-z0 : zp-z0-0.5f;
	 		    r0 = 1-r1;
	        }
	        else if( b == 3 ){
	        	k0=floor(zp)-1; if(k0 < 0) k0 = 0;
		        k1 = k0+1;
		        r1 = zp - 1 - k0;
		        r0 = 1 - r1;
		        i0 = floor(xp); x0 = i0;
		        if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		        i1=i0+1;
		        s1 = xp-x0 < 0.5f ? 0.5f+xp-x0 : xp-x0-0.5f;
		        s0 = 1-s1;
		        j0 = floor(yp); y0 = j0;
		        if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		        j1=j0+1;
		        t1 = yp-y0 < 0.5f ? 0.5f+yp-y0 : yp-y0-0.5f;
		        t0 = 1-t1;
	        }
#ifdef USE_MESH
	        AssignEightVelocities(i0, j0, k0, i1, j1, k1, d0, f);
			HandleEightRaysForVelocityInterpolation(voxel,
							i0, j0, k0, i1, j1, k1,	b, pi, f);

	        d[INDEX(i,j,k)] = r0*(s0*(t0*f[0]+t1*f[1])+
	        				      s1*(t0*f[2]+t1*f[3])) +
	        			      r1*(s0*(t0*f[4]+t1*f[5])+
	        	     		      s1*(t0*f[6]+t1*f[7]));
#else
			d[INDEX(i,j,k)] = r0*(s0*(t0*d0[INDEX(i0,j0,k0)]+t1*d0[INDEX(i0,j1,k0)])+
				                  s1*(t0*d0[INDEX(i1,j0,k0)]+t1*d0[INDEX(i1,j1,k0)])) +
			                  r1*(s0*(t0*d0[INDEX(i0,j0,k1)]+t1*d0[INDEX(i0,j1,k1)])+
	     		                  s1*(t0*d0[INDEX(i1,j0,k1)]+t1*d0[INDEX(i1,j1,k1)]));
#endif
//			if(  i == I && j == J && k == K && b == 1 ){
//				printf("(AdvectOneVelocity) at (%d, %d, %d) d0 = %f, d = %f \n",
//						i, j, k, d0[INDEX(i,j,k)], d[INDEX(i,j,k)]);
//				printf(" d0[%d,%d,%d] = %f\n", i0, j0, k0, d0[INDEX(i0,j0,k0)]);
//				printf(" d0[%d,%d,%d] = %f\n", i0, j1, k0, d0[INDEX(i0,j1,k0)]);
//				printf(" d0[%d,%d,%d] = %f\n", i1, j0, k0, d0[INDEX(i1,j0,k0)]);
//				printf(" d0[%d,%d,%d] = %f\n", i1, j1, k0, d0[INDEX(i1,j1,k0)]);
//				printf(" d0[%d,%d,%d] = %f\n", i0, j0, k1, d0[INDEX(i0,j0,k1)]);
//				printf(" d0[%d,%d,%d] = %f\n", i0, j1, k1, d0[INDEX(i0,j1,k1)]);
//				printf(" d0[%d,%d,%d] = %f\n", i1, j0, k1, d0[INDEX(i1,j0,k1)]);
//				printf(" d0[%d,%d,%d] = %f\n", i1, j1, k1, d0[INDEX(i1,j1,k1)]);
//				printf(" t0 = %f, t1 = %f, s0 = %f, s1 = %f, r0 = %f, r1 = %f\n",
//						t0, t1, s0, s1, r0, r1);
//			}
			}
		else
			d[INDEX(i,j,k)] = d0[INDEX(i,j,k)];
	END_FOR

}


void VectorField3D::Advect(const Voxel &voxel){

	printf("Before advection point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
			I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );
	printf("Before advection point (%d, %d, %d)  u = %f, v = %f, w = %f \n",
			I, J, K, u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)] );
	printf("\nAdvection of Velocity Starts...\n");
//	SWAP ( u0, u );
//	SWAP ( v0, v );
//	SWAP ( w0, w );
	SetEqual ( u1, u0 );
	SetEqual ( v1, v0 );
	SetEqual ( w1, w0 );
	AdvectOneVelocity(voxel, 1, u, u1);
	AdvectOneVelocity(voxel, 2, v, v1);
	AdvectOneVelocity(voxel, 3, w, w1);

	printf("After advection point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
			I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );
	printf("After advection point (%d, %d, %d)  u = %f, v = %f, w = %f \n",
			I, J, K, u[INDEX(I,J,K)],v[INDEX(I,J,K)],w[INDEX(I,J,K)] );


#ifdef AIR_DIV_FREE
//	SetEqual(u, u_star);
//	SetEqual(v, v_star);
//	SetEqual(w, w_star);
#endif
	printf("\nAdvection of Velocity Ends \n\n");

}

	//  if we use Houston (2003) apporach
void VectorField3D::AdvectHouston2003(const Voxel &voxel){

	printf("Before advection point (%d, %d, %d)  u1 = %f, v1 = %f, w1 = %f \n",
				I, J, K, u1[INDEX(I,J,K)],v1[INDEX(I,J,K)],w1[INDEX(I,J,K)] );

	printf("\nAdvection of Velocity with Houston (2003) approach Starts...\n");
	SetEqual ( u1, u0 );
	SetEqual ( v1, v0 );
	SetEqual ( w1, w0 );
	AdvectOneVelocity(voxel, 1, u, u1);
	AdvectOneVelocity(voxel, 2, v, v1);
	AdvectOneVelocity(voxel, 3, w, w1);

#ifdef AIR_DIV_FREE
//	SetEqual(u, u_star);
//	SetEqual(v, v_star);
//	SetEqual(w, w_star);
#endif
	printf("\nAdvection of Velocity with Houston (2003) approach Ends \n\n");

}

bool VectorField3D::HasNewExtrema(const Voxel &voxel, int i, int j, int k, int b, float d, float *d0) const{


	int i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, s1, t1, dt0;

	float xp, yp, zp, up, vp, wp;

	int x0, y0, z0;

	float r0, r1;

	float delta = voxel.VoxelDelta();

	dt0 = dt/delta;


		bool u_active = phi_u_obj[INDEX(i,j,k)] < 0.f && ( phi_u[INDEX(i,j,k)] <= 0.f ||
						(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
						(voxel.InSolid(i,j,k) && i < DimX-1 && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) ||
						(voxel.InAir(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) ||
						(voxel.InSurface(i,j,k) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)) );
		bool v_active = phi_v_obj[INDEX(i,j,k)] < 0.f && ( phi_v[INDEX(i,j,k)] <= 0.f ||
						(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
						(voxel.InSolid(i,j,k) && j < DimY-1 && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) ||
						(voxel.InAir(i,j,k) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) ||
						(voxel.InSurface(i,j,k) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)) );
		bool w_active = phi_w_obj[INDEX(i,j,k)] < 0.f && ( phi_w[INDEX(i,j,k)] <= 0.f ||
						(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)) ||
						(voxel.InSolid(i,j,k) && k < DimZ-1 && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) ||
						(voxel.InAir(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) ||
						(voxel.InSurface(i,j,k) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)) );

		if( (b == 1 && u_active) || (b == 2 && v_active) ||	(b == 3 && w_active) ) {

	        if( b == 1 ){
		        xp = i+1.f-0.5f*dt0*u0[INDEX(i,j,k)];
			   	yp = j+0.5f-0.5f*dt0*0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i+1,j,k)]+v0[INDEX(i+1,j-1,k)]);
			   	zp = k+0.5f-0.5f*dt0*0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i+1,j,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i+1,j,k-1)]);
	        }
	        else if( b == 2 ){
	            xp = i+0.5f-0.5f*dt0*0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j+1,k)]+u0[INDEX(i-1,j+1,k)]);
	            yp = j+1.f-0.5f*dt0*v0[INDEX(i,j,k)];
	            zp = k+0.5f-0.5f*dt0*0.25f*(w0[INDEX(i,j,k)]+w0[INDEX(i,j+1,k)]+w0[INDEX(i,j,k-1)]+w0[INDEX(i,j+1,k-1)]);
	        }
	        else if( b == 3 ){
                xp = i+0.5-0.5f*dt0*0.25f*(u0[INDEX(i,j,k)]+u0[INDEX(i-1,j,k)]+u0[INDEX(i,j,k+1)]+u0[INDEX(i-1,j,k+1)]);
                yp = j+0.5-0.5f*dt0*0.25f*(v0[INDEX(i,j,k)]+v0[INDEX(i,j-1,k)]+v0[INDEX(i,j,k+1)]+v0[INDEX(i,j-1,k+1)]);
                zp = k+1.f-0.5f*dt0*w0[INDEX(i,j,k)];
	        }
	        if (xp<0.5f) xp=0.5f; if (xp>DimX-1.5f) xp=DimX-1.5f;
	        if (yp<0.5f) yp=0.5f; if (yp>DimY-1.5f) yp=DimY-1.5f;
	        if (zp<0.5f) zp=0.5f; if (zp>DimZ-1.5f) zp=DimZ-1.5f;

//	        if(  i == I && j == J && k == K && b == 1 ){
//				printf(" xp = %f, yp = %f, zp = %f \n", xp, yp, zp);
//			}

	        i0=floor(xp)-1; if(i0 < 0) i0 = 0;
        	i1 = i0 + 1;
            s1 = xp - 1 - i0;
            s0 = 1 - s1;
            j0 = floor(yp); y0 = j0;
            if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
            j1=j0+1;
            if(yp-y0 < 0.5f) t1 = 0.5f+yp-y0;
            else t1 = yp-y0-0.5f;
            t0 = 1-t1;
            k0 = floor(zp); z0 = k0;
            if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
            k1=k0+1;
            if(zp-z0 < 0.5f) r1 = 0.5f+zp-z0;
            else r1 = zp-z0-0.5f;
            r0 = 1-r1;

//			if(i0 < 4){
//				printf("u0 at i0(%d) = %f, u0 at i1(%d) = %f\n",
//						i0, u0[INDEX(i0,j0,k0)], i1, u0[INDEX(i1,j0,k0)]);
//				//exit(1);
//			}
	 		up = r0*(s0*(t0*u0[INDEX(i0,j0,k0)]+t1*u0[INDEX(i0,j1,k0)])+
				  s1*(t0*u0[INDEX(i1,j0,k0)]+t1*u0[INDEX(i1,j1,k0)])) +
		 		 r1*(s0*(t0*u0[INDEX(i0,j0,k1)]+t1*u0[INDEX(i0,j1,k1)])+
					  s1*(t0*u0[INDEX(i1,j0,k1)]+t1*u0[INDEX(i1,j1,k1)]));

	 		j0=floor(yp)-1; if(j0 < 0) j0 = 0;
 			j1 = j0 + 1;
 		    t1 = yp - 1 - j0;
 		    t0 = 1 - t1;
 		    i0 = floor(xp); x0 = i0;
 		    if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
 		    i1=i0+1;
 		    if(xp-x0 < 0.5f) s1 = 0.5f+xp-x0;
 		    else s1 = xp-x0-0.5f;
 		    s0 = 1-s1;
 		    k0 = floor(zp); z0 = k0;
 		    if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
 		    k1=k0+1;
 		    if(zp-z0 < 0.5f) r1 = 0.5f+zp-z0;
 		    else r1 = zp-z0-0.5f;
 		    r0 = 1-r1;

		    vp = r0*(s0*(t0*v0[INDEX(i0,j0,k0)]+t1*v0[INDEX(i0,j1,k0)])+
				    s1*(t0*v0[INDEX(i1,j0,k0)]+t1*v0[INDEX(i1,j1,k0)])) +
				    r1*(s0*(t0*v0[INDEX(i0,j0,k1)]+t1*v0[INDEX(i0,j1,k1)])+
					    s1*(t0*v0[INDEX(i1,j0,k1)]+t1*v0[INDEX(i1,j1,k1)]));

		    k0=floor(zp)-1; if(k0 < 0) k0 = 0;
	        k1 = k0+1;
	        r1 = zp - 1 - k0;
	        r0 = 1 - r1;
	        i0 = floor(xp); x0 = i0;
	        if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
	        i1=i0+1;
	        if(xp-x0 < 0.5f) s1 = 0.5f+xp-x0;
	        else s1 = xp-x0-0.5f;
	        s0 = 1-s1;
	        j0 = floor(yp); y0 = j0;
	        if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
	        j1=j0+1;
	        if(yp-y0 < 0.5f) t1 = 0.5f+yp-y0;
	        else t1 = yp-y0-0.5f;
	        t0 = 1-t1;

		    wp = r0*(s0*(t0*w0[INDEX(i0,j0,k0)]+t1*w0[INDEX(i0,j1,k0)])+
				    	s1*(t0*w0[INDEX(i1,j0,k0)]+t1*w0[INDEX(i1,j1,k0)])) +
			     r1*(s0*(t0*w0[INDEX(i0,j0,k1)]+t1*w0[INDEX(i0,j1,k1)])+
					    s1*(t0*w0[INDEX(i1,j0,k1)]+t1*w0[INDEX(i1,j1,k1)]));

	        if( b == 1 ){
		        x = i+1.f-dt0*up; y = j+0.5f - dt0*vp; z = k+0.5f - dt0*wp;
	        }
	        else if( b == 2 ){
		        x = i+0.5f-dt0*up; y = j+1.f-dt0*vp; z = k+0.5f - dt0*wp;
	        }
	        else if( b == 3 ){
		        x = i+0.5f-dt0*up; y = j+0.5f-dt0*vp; z = k+1.f - dt0*wp;
	        }
	        if (x<0.5f) x=0.5f; if (x>DimX-1.5f) x=DimX-1.5f;
  	        if (y<0.5f) y=0.5f; if (y>DimY-1.5f) y=DimY-1.5f;
  	        if (z<0.5f) z=0.5f; if (z>DimZ-1.5f) z=DimZ-1.5f;

//  	        if(  i == I && j == J && k == K && b == 1 ){
//  	        	printf(" up = %f, vp = %f, wp = %f \n", up, vp, wp);
//  	        	printf(" x = %f, y = %f, z = %f \n", x, y, z);
//			}

	        if( b == 1 ){
	        	i0=floor(x)-1; if(i0 < 0) i0 = 0;
	        	i1 = i0 + 1;
	            j0 = floor(y); y0 = j0;
	            if(y-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
	            j1=j0+1;
	            k0 = floor(z); z0 = k0;
	            if(z-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
	            k1=k0+1;

	        }
	        else if( b == 2 ){
	        	j0=floor(y)-1; if(j0 < 0) j0 = 0;
	 			j1 = j0 + 1;
	 		    i0 = floor(x); x0 = i0;
	 		    if(x-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
	 		    i1=i0+1;
	 		    k0 = floor(z); z0 = k0;
	 		    if(z-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
	 		    k1=k0+1;
	        }
	        else if( b == 3 ){
	        	k0=floor(z)-1; if(k0 < 0) k0 = 0;
		        k1 = k0+1;
		        i0 = floor(x); x0 = i0;
		        if(x-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		        i1=i0+1;
		        j0 = floor(y); y0 = j0;
		        if(y-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		        j1=j0+1;
	        }

			float tm = max(d0[INDEX(i0,j0,k0)], d0[INDEX(i0,j1,k0)]);
			tm = max(tm, d0[INDEX(i1,j0,k0)]);
			tm = max(tm, d0[INDEX(i1,j1,k0)]);
			tm = max(tm, d0[INDEX(i0,j0,k1)]);
			tm = max(tm, d0[INDEX(i0,j1,k1)]);
			tm = max(tm, d0[INDEX(i1,j0,k1)]);
			tm = max(tm, d0[INDEX(i1,j1,k1)]);
			float ts = min(d0[INDEX(i0,j0,k0)], d0[INDEX(i0,j1,k0)]);
			ts = min(ts, d0[INDEX(i1,j0,k0)]);
			ts = min(ts, d0[INDEX(i1,j1,k0)]);
			ts = min(ts, d0[INDEX(i0,j0,k1)]);
			ts = min(ts, d0[INDEX(i0,j1,k1)]);
			ts = min(ts, d0[INDEX(i1,j0,k1)]);
			ts = min(ts, d0[INDEX(i1,j1,k1)]);
			if(d > tm || d < ts)
				return true;
			else
				return false;
		}
		else
			return true;
}

void VectorField3D::AdvectMacCormack(const Voxel &voxel){

	float delta = voxel.VoxelDelta();
	printf("\nAdvection of Velocity Using MacCormack Method Starts...\n");
	float *u2 = new float[DimX*DimY*DimZ];
	float *v2 = new float[DimX*DimY*DimZ];
	float *w2 = new float[DimX*DimY*DimZ];
	SetEqual(u1, u2);
	SetEqual(v1, v2);
	SetEqual(w1, w2);
	SWAP ( u0, u );
	SWAP ( v0, v );
	SWAP ( w0, w );
	// note: we are advecting u1, v1 and w1 using divergence free velocity
	// u, v, w from last time step.
	//
	// note: u1, v1 and w1 within solid cells contains solid velocities
	AdvectOneVelocity(voxel, 1, u, u1);
	AdvectOneVelocity(voxel, 2, v, v1);
	AdvectOneVelocity(voxel, 3, w, w1);
	FOR_EACH_CELL
		u0[INDEX(i,j,k)] = -1.f * u0[INDEX(i,j,k)];
		v0[INDEX(i,j,k)] = -1.f * v0[INDEX(i,j,k)];
		w0[INDEX(i,j,k)] = -1.f * w0[INDEX(i,j,k)];
	END_FOR
	AdvectOneVelocity(voxel, 1, u1, u);
	AdvectOneVelocity(voxel, 2, v1, v);
	AdvectOneVelocity(voxel, 3, w1, w);
	FOR_EACH_CELL
		u0[INDEX(i,j,k)] = -1.f * u0[INDEX(i,j,k)];
		v0[INDEX(i,j,k)] = -1.f * v0[INDEX(i,j,k)];
		w0[INDEX(i,j,k)] = -1.f * w0[INDEX(i,j,k)];
	END_FOR
#define MACCORMACK_LIMIT -5*delta
//	FOR_EACH_CELL
//		if(phi_c[INDEX(i,j,k)] < MACCORMACK_LIMIT && phi_c_obj[INDEX(i,j,k)] < MACCORMACK_LIMIT){
//			u[INDEX(i,j,k)] -= 0.5f*(u1[INDEX(i,j,k)] - u0[INDEX(i,j,k)]);
//			v[INDEX(i,j,k)] -= 0.5f*(v1[INDEX(i,j,k)] - v0[INDEX(i,j,k)]);
//			w[INDEX(i,j,k)] -= 0.5f*(w1[INDEX(i,j,k)] - w0[INDEX(i,j,k)]);
//		}
//	END_FOR
	FOR_EACH_CELL
		if( phi_c_obj[INDEX(i,j,k)] < MACCORMACK_LIMIT ){
			float utmp = u[INDEX(i,j,k)] - 0.5f*(u1[INDEX(i,j,k)] - u2[INDEX(i,j,k)]);
			if(!HasNewExtrema(voxel, i, j, k, 1, utmp, u2))
				u[INDEX(i,j,k)] = utmp;
			float vtmp = v[INDEX(i,j,k)] - 0.5f*(v1[INDEX(i,j,k)] - v2[INDEX(i,j,k)]);
			if(!HasNewExtrema(voxel, i, j, k, 2, vtmp, v2))
				v[INDEX(i,j,k)] = vtmp;
			float wtmp = w[INDEX(i,j,k)] - 0.5f*(w1[INDEX(i,j,k)] - w2[INDEX(i,j,k)]);
			if(!HasNewExtrema(voxel, i, j, k, 3, wtmp, w2))
				w[INDEX(i,j,k)] = wtmp;
		}
	END_FOR

	delete [] u2;
	delete [] v2;
	delete [] w2;
//	SetSolidBoundary(voxel);

//	Extrapolate(voxel);
//	SetSolidBoundary(voxel);
//	SetSolidBoundaryForAdvection(voxel);
//	SetSurfaceBoundary(voxel);


#ifdef AIR_DIV_FREE
//	SetEqual(u, u_star);
//	SetEqual(v, v_star);
//	SetEqual(w, w_star);
#endif
	printf("\nAdvection of Velocity Using MacCormack Method Ends \n\n");

}

// 1st order upwind
/*void VectorField3D::Advect(const Voxel &voxel, float *x, float *x0) {




	float delta = voxel.VoxelDelta();

	float diff_x, diff_y, diff_z;
	float u_c, v_c, w_c;

	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && fabs(x0[INDEX(i,j,k)]) < 3 * delta){
			u_c = 0.5f * (u[INDEX(i,j,k)] + u[INDEX(i-1,j,k)]);
			v_c = 0.5f * (v[INDEX(i,j,k)] + v[INDEX(i,j-1,k)]);
			w_c = 0.5f * (w[INDEX(i,j,k)] + w[INDEX(i,j,k-1)]);

			if(u_c > 0.f)
				diff_x = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) / delta;
			else
				diff_x = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) / delta;
			if(v_c > 0.f)
				diff_y = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) / delta;
			else
				diff_y = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) / delta;
			if(w_c > 0.f)
				diff_z = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) / delta;
			else
				diff_z = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) / delta;

			if(i == I & j == J && k == K){
				 printf("\nBefore Advection: at i = %d, j = %d, k = %d  x0 = %f, u = %f, v = %f, w = %f, diff_x = %f, diff_y = %f, diff_z = %f \n",
				 			I,J,K, x0[INDEX(i,j,k)], u_c ,v_c, w_c, diff_x, diff_y, diff_z);

			}

			x[INDEX(i,j,k)] = x0[INDEX(i,j,k)] -
							dt*(u_c*diff_x + v_c*diff_y + w_c*diff_z);

			if(i == I & j == J && k == K){
				 printf("After Advection: at i = %d, j = %d, k = %d, x = %f,  u = %f, v = %f, w = %f, diff_x = %f, diff_y = %f, diff_z = %f \n\n",
							I,J,K, x[INDEX(i,j,k)], u_c ,v_c, w_c, diff_x, diff_y, diff_z);

			}


		}
		else
			x[INDEX(i,j,k)] = x0[INDEX(i,j,k)];


	END_FOR

}*/

// 3rd order ENO
// see reference in Osher and Fekiw (2002) Chapter 3.3
/*void VectorField3D::Advect(const Voxel &voxel, float *x, float *x0) {


	float delta = voxel.VoxelDelta();

	float diff_x, diff_y, diff_z;
	float u_c, v_c, w_c;

	float q1, q2, q3;

	float c, cstar;
	int kstar;

	FOR_EACH_CELL
			if(!voxel.InSolid(i,j,k) && fabs(x0[INDEX(i,j,k)]) < 3 * delta){
				u_c = 0.5f * (u[INDEX(i,j,k)] + u[INDEX(i-1,j,k)]);
				v_c = 0.5f * (v[INDEX(i,j,k)] + v[INDEX(i,j-1,k)]);
				w_c = 0.5f * (w[INDEX(i,j,k)] + w[INDEX(i,j,k-1)]);
				if(u_c > 0.f){
					q1 = (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) / delta;
					float d2k = (q1 - (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) / delta)
					           / ( 2 * delta );
					float d2kp1 = ((x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) / delta - q1)
					           / ( 2 * delta );
					if(fabs(d2k) <= fabs(d2kp1)){
						c = d2k;
						kstar = i-2;
					}
					else{
						c = d2kp1;
						kstar = i-1;
					}
					float d3kp = ( (x0[INDEX(kstar+2,j,k)] - x0[INDEX(kstar+1,j,k)]) / delta -
								   (x0[INDEX(kstar+1,j,k)] - x0[INDEX(kstar,j,k)]) / delta ) /
								  (2*delta);
					float d3km = ( (x0[INDEX(kstar+1,j,k)] - x0[INDEX(kstar,j,k)]) / delta -
								   (x0[INDEX(kstar,j,k)] - x0[INDEX(kstar-1,j,k)]) / delta ) /
								 (2*delta);
					float d3k = (d3kp - d3km)/ (3*delta);
					d3kp = ( (x0[INDEX(kstar+3,j,k)] - x0[INDEX(kstar+2,j,k)]) / delta -
							   (x0[INDEX(kstar+2,j,k)] - x0[INDEX(kstar+1,j,k)]) / delta ) /
							  (2*delta);
					d3km = ( (x0[INDEX(kstar+2,j,k)] - x0[INDEX(kstar+1,j,k)]) / delta -
							   (x0[INDEX(kstar+1,j,k)] - x0[INDEX(kstar,j,k)]) / delta ) /
							 (2*delta);
					float d3kp1 = (d3kp - d3km)/ (3*delta);
					if(fabs(d3k) <= fabs(d3kp1))
						cstar = d3k;
					else
						cstar = d3kp1;
					q2 = c*delta;
					q3 = cstar*(3*(i-kstar)*(i-kstar)-6*(i-kstar)+2)*(delta*delta);
					diff_x = q1 + q2 + q3;
				}
				else{
					q1 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i,j,k)]) / delta;
					float d2k = (q1 - (x0[INDEX(i,j,k)] - x0[INDEX(i-1,j,k)]) / delta)
			           				/ ( 2 * delta );
					float d2kp1 = ((x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) / delta - q1)
				           				/ ( 2 * delta );
					if(fabs(d2k) <= fabs(d2kp1)){
						c = d2k;
						kstar = i-1;
					}
					else{
						c = d2kp1;
						kstar = i;
					}
					float d3kp = ( (x0[INDEX(kstar+2,j,k)] - x0[INDEX(kstar+1,j,k)]) / delta -
								   (x0[INDEX(kstar+1,j,k)] - x0[INDEX(kstar,j,k)]) / delta ) /
								  (2*delta);
					float d3km = ( (x0[INDEX(kstar+1,j,k)] - x0[INDEX(kstar,j,k)]) / delta -
								   (x0[INDEX(kstar,j,k)] - x0[INDEX(kstar-1,j,k)]) / delta ) /
								 (2*delta);
					float d3k = (d3kp - d3km)/ (3*delta);
					d3kp = ( (x0[INDEX(kstar+3,j,k)] - x0[INDEX(kstar+2,j,k)]) / delta -
							   (x0[INDEX(kstar+2,j,k)] - x0[INDEX(kstar+1,j,k)]) / delta ) /
							  (2*delta);
					d3km = ( (x0[INDEX(kstar+2,j,k)] - x0[INDEX(kstar+1,j,k)]) / delta -
							   (x0[INDEX(kstar+1,j,k)] - x0[INDEX(kstar,j,k)]) / delta ) /
							 (2*delta);
					float d3kp1 = (d3kp - d3km)/ (3*delta);
					if(fabs(d3k) <= fabs(d3kp1))
						cstar = d3k;
					else
						cstar = d3kp1;
					q2 = -c*delta;
					q3 = cstar*(3*(i-kstar)*(i-kstar)-6*(i-kstar)+2)*(delta*delta);
					diff_x = q1 + q2 + q3;
				}
				if(v_c > 0.f){
					q1 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) / delta;
					float d2k = (q1 - (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) / delta)
					           / ( 2 * delta );
					float d2kp1 = ((x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) / delta - q1)
					           / ( 2 * delta );
					if(fabs(d2k) <= fabs(d2kp1)){
						c = d2k;
						kstar = j-2;
					}
					else{
						c = d2kp1;
						kstar = j-1;
					}
					float d3kp = ( (x0[INDEX(i,kstar+2,k)] - x0[INDEX(i,kstar+1,k)]) / delta -
								   (x0[INDEX(i,kstar+1,k)] - x0[INDEX(i,kstar,k)]) / delta ) /
								  (2*delta);
					float d3km = ( (x0[INDEX(i,kstar+1,k)] - x0[INDEX(i,kstar,k)]) / delta -
								   (x0[INDEX(i,kstar,k)] - x0[INDEX(i,kstar-1,k)]) / delta ) /
								 (2*delta);
					float d3k = (d3kp - d3km)/ (3*delta);
					d3kp = ( (x0[INDEX(i,kstar+3,k)] - x0[INDEX(i,kstar+2,k)]) / delta -
							   (x0[INDEX(i,kstar+2,k)] - x0[INDEX(i,kstar+1,k)]) / delta ) /
							  (2*delta);
					d3km = ( (x0[INDEX(i,kstar+2,k)] - x0[INDEX(i,kstar+1,k)]) / delta -
							   (x0[INDEX(i,kstar+1,k)] - x0[INDEX(i,kstar,k)]) / delta ) /
							 (2*delta);
					float d3kp1 = (d3kp - d3km)/ (3*delta);
					if(fabs(d3k) <= fabs(d3kp1))
						cstar = d3k;
					else
						cstar = d3kp1;
					q2 = c*delta;
					q3 = cstar*(3*(j-kstar)*(j-kstar)-6*(j-kstar)+2)*(delta*delta);
					diff_y = q1 + q2 + q3;
				}
				else{
					q1 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j,k)]) / delta;
					float d2k = (q1 - (x0[INDEX(i,j,k)] - x0[INDEX(i,j-1,k)]) / delta)
			           				/ ( 2 * delta );
					float d2kp1 = ((x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) / delta - q1)
				           				/ ( 2 * delta );
					if(fabs(d2k) <= fabs(d2kp1)){
						c = d2k;
						kstar = j-1;
					}
					else{
						c = d2kp1;
						kstar = j;
					}
					float d3kp = ( (x0[INDEX(i,kstar+2,k)] - x0[INDEX(i,kstar+1,k)]) / delta -
								   (x0[INDEX(i,kstar+1,k)] - x0[INDEX(i,kstar,k)]) / delta ) /
								  (2*delta);
					float d3km = ( (x0[INDEX(i,kstar+1,k)] - x0[INDEX(i,kstar,k)]) / delta -
								   (x0[INDEX(i,kstar,k)] - x0[INDEX(i,kstar-1,k)]) / delta ) /
								 (2*delta);
					float d3k = (d3kp - d3km)/ (3*delta);
					d3kp = ( (x0[INDEX(i,kstar+3,k)] - x0[INDEX(i,kstar+2,k)]) / delta -
							   (x0[INDEX(i,kstar+2,k)] - x0[INDEX(i,kstar+1,k)]) / delta ) /
							  (2*delta);
					d3km = ( (x0[INDEX(i,kstar+2,k)] - x0[INDEX(i,kstar+1,k)]) / delta -
							   (x0[INDEX(i,kstar+1,k)] - x0[INDEX(i,kstar,k)]) / delta ) /
							 (2*delta);
					float d3kp1 = (d3kp - d3km)/ (3*delta);
					if(fabs(d3k) <= fabs(d3kp1))
						cstar = d3k;
					else
						cstar = d3kp1;
					q2 = -c*delta;
					q3 = cstar*(3*(j-kstar)*(j-kstar)-6*(j-kstar)+2)*(delta*delta);
					diff_y = q1 + q2 + q3;
				}
				if(w_c > 0.f){
					q1 = (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) / delta;
					float d2k = (q1 - (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) / delta)
					           / ( 2 * delta );
					float d2kp1 = ((x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) / delta - q1)
					           / ( 2 * delta );
					if(fabs(d2k) <= fabs(d2kp1)){
						c = d2k;
						kstar = k-2;
					}
					else{
						c = d2kp1;
						kstar = k-1;
					}
					float d3kp = ( (x0[INDEX(i,j,kstar+2)] - x0[INDEX(i,j,kstar+1)]) / delta -
								   (x0[INDEX(i,j,kstar+1)] - x0[INDEX(i,j,kstar)]) / delta ) /
								  (2*delta);
					float d3km = ( (x0[INDEX(i,j,kstar+1)] - x0[INDEX(i,j,kstar)]) / delta -
								   (x0[INDEX(i,j,kstar)] - x0[INDEX(i,j,kstar-1)]) / delta ) /
								 (2*delta);
					float d3k = (d3kp - d3km)/ (3*delta);
					d3kp = ( (x0[INDEX(i,j,kstar+3)] - x0[INDEX(i,j,kstar+2)]) / delta -
							   (x0[INDEX(i,j,kstar+2)] - x0[INDEX(i,j,kstar+1)]) / delta ) /
							  (2*delta);
					d3km = ( (x0[INDEX(i,j,kstar+2)] - x0[INDEX(i,j,kstar+1)]) / delta -
							   (x0[INDEX(i,j,kstar+1)] - x0[INDEX(i,j,kstar)]) / delta ) /
							 (2*delta);
					float d3kp1 = (d3kp - d3km)/ (3*delta);
					if(fabs(d3k) <= fabs(d3kp1))
						cstar = d3k;
					else
						cstar = d3kp1;
					q2 = c*delta;
					q3 = cstar*(3*(k-kstar)*(k-kstar)-6*(k-kstar)+2)*(delta*delta);
					diff_z = q1 + q2 + q3;
				}
				else{
					q1 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k)]) / delta;
					float d2k = (q1 - (x0[INDEX(i,j,k)] - x0[INDEX(i,j,k-1)]) / delta)
			           				/ ( 2 * delta );
					float d2kp1 = ((x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) / delta - q1)
				           				/ ( 2 * delta );
					if(fabs(d2k) <= fabs(d2kp1)){
						c = d2k;
						kstar = k-1;
					}
					else{
						c = d2kp1;
						kstar = k;
					}
					float d3kp = ( (x0[INDEX(i,j,kstar+2)] - x0[INDEX(i,j,kstar+1)]) / delta -
								   (x0[INDEX(i,j,kstar+1)] - x0[INDEX(i,j,kstar)]) / delta ) /
								  (2*delta);
					float d3km = ( (x0[INDEX(i,j,kstar+1)] - x0[INDEX(i,j,kstar)]) / delta -
								   (x0[INDEX(i,j,kstar)] - x0[INDEX(i,j,kstar-1)]) / delta ) /
								 (2*delta);
					float d3k = (d3kp - d3km)/ (3*delta);
					d3kp = ( (x0[INDEX(i,j,kstar+3)] - x0[INDEX(i,j,kstar+2)]) / delta -
							   (x0[INDEX(i,j,kstar+2)] - x0[INDEX(i,j,kstar+1)]) / delta ) /
							  (2*delta);
					d3km = ( (x0[INDEX(i,j,kstar+2)] - x0[INDEX(i,j,kstar+1)]) / delta -
							   (x0[INDEX(i,j,kstar+1)] - x0[INDEX(i,j,kstar)]) / delta ) /
							 (2*delta);
					float d3kp1 = (d3kp - d3km)/ (3*delta);
					if(fabs(d3k) <= fabs(d3kp1))
						cstar = d3k;
					else
						cstar = d3kp1;
					q2 = -c*delta;
					q3 = cstar*(3*(k-kstar)*(k-kstar)-6*(k-kstar)+2)*(delta*delta);
					diff_z = q1 + q2 + q3;
				}
				if(i == I & j == J && k == K){
					 printf("\nBefore Advection: at i = %d, j = %d, k = %d  x0 = %f, u = %f, v = %f, w = %f, diff_x = %f, diff_y = %f, diff_z = %f \n",
					 			I,J,K, x0[INDEX(i,j,k)], u_c ,v_c, w_c, diff_x, diff_y, diff_z);

				}

				x[INDEX(i,j,k)] = x0[INDEX(i,j,k)] -
								dt*(u_c*diff_x + v_c*diff_y + w_c*diff_z);

				if(i == I & j == J && k == K){
					 printf("After Advection: at i = %d, j = %d, k = %d, x = %f,  u = %f, v = %f, w = %f, diff_x = %f, diff_y = %f, diff_z = %f \n\n",
								I,J,K, x[INDEX(i,j,k)], u_c ,v_c, w_c, diff_x, diff_y, diff_z);

				}
			}
		else
			x[INDEX(i,j,k)] = x0[INDEX(i,j,k)];
	END_FOR

}*/

// 5th order WENO
// see reference in Osher and Fekiw (2002) Chapter 3.3
static float HJ_WENO_Coefficients(float q1, float q2, float q3, float q4, float q5){
	float s1, s2, s3;
	float a1, a2, a3, at;
	float b1, b2, b3;
	float epsil = 1.e-6f;
//	float epsil = 1.e-6*max(q1,max(q2,max(q3,max(q4,q5))))+1.e-99;
	s1 = 13.f/12*(q1-2*q2+q3)*(q1-2*q2+q3) + 0.25f*(q1-4*q2+3*q3)*(q1-4*q2+3*q3);
	s2 = 13.f/12*(q2-2*q3+q4)*(q2-2*q3+q4) + 0.25f*(q2-q4)*(q2-q4);
	s3 = 13.f/12*(q3-2*q4+q5)*(q3-2*q4+q5) + 0.25f*(3*q3-4*q4+q5)*(3*q3-4*q4+q5);
	a1 = 0.1f/((epsil+s1)*(epsil+s1));
	a2 = 0.6f/((epsil+s2)*(epsil+s2));
	a3 = 0.3f/((epsil+s3)*(epsil+s3));
	at = a1 + a2 + a3;
	b1 = a1 / at;
	b2 = a2 / at;
	b3 = a3 / at;
	return b1*(q1/3 - 7.f/6*q2 + 11.f/6*q3) +
	       b2*(-q2/6 + 5.f/6*q3 + q4/3) +
	       b3*(q3/3 + 5.f/6*q4 - q5/6);
}

#ifdef SPMD

static void MarchingOutForPhi(const TentativeGridPoint &minTentative, int DimX, int DimY, int DimZ,
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

void VectorField3D::ExtrapolatePhiIntoObject(const Voxel &voxel,
									         float *phi_tmp,
									         float *phi_tmp0){

	float delta = voxel.VoxelDelta();

	printf("(VectorField3D) Start extrapolating phi into solid object\n ");

	float inv_delta2 = 1.f / (delta * delta);
	float dtau = delta / 2;
	float limit = 2 * delta;
	int iterations = int(round(limit / delta * 2));
	for (int n = 0; n < iterations; ++n){
		SWAP(phi_tmp, phi_tmp0);
		FOR_EACH_CELL
			if(voxel.InSolid(i,j,k) && phi_c_obj[INDEX(i,j,k)] < limit){
				KnownPoint tmp(i, j, k);
				TentativeGridPoint tp(phi_c_obj[INDEX(i,j,k)], tmp);
				MarchingOutForPhi(tp, DimX, DimY, DimZ, dtau, inv_delta2, phi_c_obj, phi_tmp, phi_tmp0);
//					if(i == I && j == J && k == K){
//						printf(" iter = %d, i = %d, j = %d, k = %d, phi = %f, phi[j+1] = %f, phi[i-1] = %f\n",
//								n,i,j,k, phi_tmp[INDEX(i,j,k)], phi_tmp[INDEX(i,j+1,k)], phi_tmp[INDEX(i-1,j,k)]);
//					}
			}
		END_FOR
	}
	FOR_EACH_CELL
	    // if a cell is occupied by a "moving" solid object, we
	    // reset its levelset to phi_c_obj (a positive quantity)
		if(voxel.InMovingSolid(i,j,k)){
			if(phi_tmp[INDEX(i,j,k)] < phi_c_obj[INDEX(i,j,k)])
				phi_tmp[INDEX(i,j,k)] = phi_c_obj[INDEX(i,j,k)];
		}
	END_FOR
	printf("(VectorField3D) End extrapolating phi into solid object\n ");
}

bool VectorField3D::CFLExceedsOne(const Voxel &voxel, float dts, const float *x0, float *f){

	float delta = voxel.VoxelDelta();
	float u_c, v_c, w_c;
	float cfl_max = 0.f;
	int cfl_i = 0, cfl_j = 0, cfl_k = 0;
	float gamma = FASTMARCH_LIMIT;
	float beta = LS_ADV_LIMIT;
	float gb = (gamma-beta)*(gamma-beta)*(gamma-beta);
	SetZero(needRecomputePhi);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && fabs(x0[pos])+E_EPSIL < FASTMARCH_LIMIT ){
			needRecomputePhi[pos] = 1;
			if(fabs(x0[pos]) <= beta){
				temp[pos] = 1.f;
			}
			else{
				temp[pos] = (fabs(x0[pos]) - gamma)*(fabs(x0[pos]) - gamma)
							 * (2 * fabs(x0[pos]) + gamma - 3 * beta) / gb;
			}
		}
	END_FOR
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
			if(needRecomputePhi[pos]){
				float c = temp[pos];
				u_c = 0.5f * c * (u[pos] + u[INDEX(i-1,j,k)]);
				v_c = 0.5f * c * (v[pos] + v[INDEX(i,j-1,k)]);
				w_c = 0.5f * c * (w[pos] + w[INDEX(i,j,k-1)]);

				float cfl = dts / delta * (fabsf(u_c) +fabsf(v_c)+fabsf(w_c));
				if(cfl_max < cfl){
					cfl_max = cfl;
					cfl_i = i;
					cfl_j = j;
					cfl_k = k;
				}
				if(cfl > 1.f){
					printf("\nat (%d, %d, %d) CFL = %f  > 1 \n", i,j,k, cfl);
					printf("In Euler at i = %d, j = %d, k = %d  x0 = %f\n",
								i, j, k, x0[INDEX(i,j,k)]);
					printf(" u_c = %f, v_c = %f, w_c = %f \n", u_c, v_c, w_c);
//					exit(1);
					*f = cfl;
					return true;
				}
			}

	END_FOR
	printf("\nmaximum CFL at this step = %f at cell (%d, %d, %d)\n\n ",
			cfl_max, cfl_i, cfl_j, cfl_k);
	*f = cfl_max;
	return false;
}

float VectorField3D::FindTimeStep(const Voxel &voxel, const float *x0){

	float delta = voxel.VoxelDelta();
	float u_c, v_c, w_c;
	float dt_min = 99999.f;
	int cfl_i = 0, cfl_j = 0, cfl_k = 0;
	float gamma = FASTMARCH_LIMIT;
	float beta = LS_ADV_LIMIT;
	float gb = (gamma-beta)*(gamma-beta)*(gamma-beta);
	char *needed = new char[DimX*DimY*DimZ];
	SetZero(needed);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && fabs(x0[pos])+E_EPSIL < FASTMARCH_LIMIT ){
			needed[pos] = 1;
			if(fabs(x0[pos]) <= beta){
				temp[pos] = 1.f;
			}
			else{
				temp[pos] = (fabs(x0[pos]) - gamma)*(fabs(x0[pos]) - gamma)
							 * (2 * fabs(x0[pos]) + gamma - 3 * beta) / gb;
			}
		}
	END_FOR
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
			if(needed[pos]){
				float c = temp[pos];
				u_c = 0.5f * c * (u[pos] + u[INDEX(i-1,j,k)]);
				v_c = 0.5f * c * (v[pos] + v[INDEX(i,j-1,k)]);
				w_c = 0.5f * c * (w[pos] + w[INDEX(i,j,k-1)]);
				float vel = fabsf(u_c)+ fabsf(v_c)+fabsf(w_c);
				float dts;
				if(vel > E_EPSIL)
					dts = delta / vel;
				else
					dts = 9999.f;
				if(dts < dt_min){
					dt_min = dts;
					cfl_i = i;
					cfl_j = j;
					cfl_k = k;
				}
			}
	END_FOR
	delete [] needed;
	printf("\nminimum time step = %f at cell (%d, %d, %d)\n\n ",
		dt_min, cfl_i, cfl_j, cfl_k);
	return dt_min;
}

void VectorField3D::Euler(const Voxel &voxel, float dts, char *needed, float *x, float *x0, char *cflex1=NULL){

	float delta = voxel.VoxelDelta();
	float u_c, v_c, w_c;
	float q1, q2, q3, q4, q5;
	float diff_x, diff_y, diff_z;
	float cfl_max = 0.f;
	int cfl_i = 0, cfl_j = 0, cfl_k = 0;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
			if(needed[pos]){
				float c = temp[pos];
				u_c = 0.5f * c * (u[pos] + u[INDEX(i-1,j,k)]);
				v_c = 0.5f * c * (v[pos] + v[INDEX(i,j-1,k)]);
				w_c = 0.5f * c * (w[pos] + w[INDEX(i,j,k-1)]);
				if(u_c > 0.f){
					q1 = (x0[INDEX(i-2,j,k)] - x0[INDEX(i-3,j,k)]) / delta;
					q2 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) / delta;
					q3 = (x0[INDEX(i  ,j,k)] - x0[INDEX(i-1,j,k)]) / delta;
					q4 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i  ,j,k)]) / delta;
					q5 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) / delta;
				}
				else{
					q1 = (x0[INDEX(i+3,j,k)] - x0[INDEX(i+2,j,k)]) / delta;
					q2 = (x0[INDEX(i+2,j,k)] - x0[INDEX(i+1,j,k)]) / delta;
					q3 = (x0[INDEX(i+1,j,k)] - x0[INDEX(i  ,j,k)]) / delta;
					q4 = (x0[INDEX(i  ,j,k)] - x0[INDEX(i-1,j,k)]) / delta;
					q5 = (x0[INDEX(i-1,j,k)] - x0[INDEX(i-2,j,k)]) / delta;

				}
				diff_x = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);

				if(v_c > 0.f){
					q1 = (x0[INDEX(i,j-2,k)] - x0[INDEX(i,j-3,k)]) / delta;
					q2 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) / delta;
					q3 = (x0[INDEX(i,j  ,k)] - x0[INDEX(i,j-1,k)]) / delta;
					q4 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j  ,k)]) / delta;
					q5 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) / delta;

				}
				else{
					q1 = (x0[INDEX(i,j+3,k)] - x0[INDEX(i,j+2,k)]) / delta;
					q2 = (x0[INDEX(i,j+2,k)] - x0[INDEX(i,j+1,k)]) / delta;
					q3 = (x0[INDEX(i,j+1,k)] - x0[INDEX(i,j  ,k)]) / delta;
					q4 = (x0[INDEX(i,j  ,k)] - x0[INDEX(i,j-1,k)]) / delta;
					q5 = (x0[INDEX(i,j-1,k)] - x0[INDEX(i,j-2,k)]) / delta;
				}
				diff_y = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);

				if(w_c > 0.f){
					q1 = (x0[INDEX(i,j,k-2)] - x0[INDEX(i,j,k-3)]) / delta;
					q2 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) / delta;
					q3 = (x0[INDEX(i,j,k  )] - x0[INDEX(i,j,k-1)]) / delta;
					q4 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k  )]) / delta;
					q5 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) / delta;
				}
				else{
					q1 = (x0[INDEX(i,j,k+3)] - x0[INDEX(i,j,k+2)]) / delta;
					q2 = (x0[INDEX(i,j,k+2)] - x0[INDEX(i,j,k+1)]) / delta;
					q3 = (x0[INDEX(i,j,k+1)] - x0[INDEX(i,j,k  )]) / delta;
					q4 = (x0[INDEX(i,j,k  )] - x0[INDEX(i,j,k-1)]) / delta;
					q5 = (x0[INDEX(i,j,k-1)] - x0[INDEX(i,j,k-2)]) / delta;
				}
				diff_z = HJ_WENO_Coefficients(q1, q2, q3, q4, q5);

				x[INDEX(i,j,k)] = x0[INDEX(i,j,k)] -
								dts*(u_c*diff_x + v_c*diff_y + w_c*diff_z);
//				if(i == I && j == J && k == K){
//					 printf("\n In Euler at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
//					 			i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
//					 printf(" u_c = %f, v_c = %f, w_c = %f \n", u_c, v_c, w_c);
//					 printf(" diff_x = %f, diff_y = %f, diff_z = %f \n", diff_x, diff_y, diff_z);
//					 printf(" u_c*diff_x = %f, v_c*diff_y = %f, w_c*diff_z = %f \n",
//							 u_c*diff_x, v_c*diff_y, w_c*diff_z);
//				}
				float cfl = dts / delta * (fabsf(u_c) +fabsf(v_c)+fabsf(w_c));
				if(cfl_max < cfl){
					cfl_max = cfl;
					cfl_i = i;
					cfl_j = j;
					cfl_k = k;
				}
				if(cfl > 1.f){
					printf("\n CFL exceeds 1 \n ");
					printf("In Euler at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
								i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
					printf(" u_c = %f, v_c = %f, w_c = %f \n", u_c, v_c, w_c);
					printf(" diff_x = %f, diff_y = %f, diff_z = %f \n", diff_x, diff_y, diff_z);
					printf(" u_c*diff_x = %f, v_c*diff_y = %f, w_c*diff_z = %f \n",
							 u_c*diff_x, v_c*diff_y, w_c*diff_z);
					if( cflex1 != NULL )
						cflex1[pos] = 1;
//					exit(1);
				}
			}
		else
			x[INDEX(i,j,k)] = x0[INDEX(i,j,k)];

	END_FOR
//	printf("\nmaximum CFL at this step = %f at cell (%d, %d, %d)\n\n ",
//		cfl_max, cfl_i, cfl_j, cfl_k);

}

void VectorField3D::Advect_RK3_WENO(const Voxel &voxel, char *needed, float *x, float *x0) {

	float delta = voxel.VoxelDelta();
	int iterations = int(round(LS_ADV_LIMIT / delta));
	float dts  = dt / iterations;
	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;
	printf(" subcycling time step = %f and complete time step = %f\n", dts, dt);

	for(int n = 0; n < iterations; ++n){
		Euler(voxel, dts, needed, u0, x0);
		Euler(voxel, dts, needed, v0, u0);
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
//			if(needed[pos])
			u0[pos] = 0.75f * x0[pos] + 0.25f * v0[pos];
		END_FOR
		Euler(voxel, dts, needed, v0, u0);
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
//			if(needed[pos])
			x[pos] = one_third * x0[pos] + two_third * v0[pos];
		END_FOR
//		printf("\n In Euler before swap iter = %d, at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
//				 		n, I, J, K, x0[INDEX(I,J,K)], x[INDEX(I,J,K)]);
		SetEqual(x, u0);
		ExtrapolatePhiIntoObject(voxel, u0, x0);
		SetEqual(u0, x);
		SetEqual(x, x0);
//		printf("\n In Euler after swap iter = %d, at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
//						n, I, J, K, x0[INDEX(I,J,K)], x[INDEX(I,J,K)]);

	}
}

void VectorField3D::SemiLagrangian(const Voxel &voxel, int i, int j, int k,
				float dt0, float *d, float *d0) {

		int i0, j0, k0, i1, j1, k1;
		float s0, t0, s1, t1;
		float xp, yp, zp;
		int x0, y0, z0;
		float r0, r1;
		u_int pos = INDEX(i,j,k);

		float c = temp[pos];
		xp = i+0.5f-dt0*0.5f*c*(u[INDEX(i,j,k)]+u[INDEX(i-1,j,k)]);
		yp = j+0.5f-dt0*0.5f*c*(v[INDEX(i,j,k)]+v[INDEX(i,j-1,k)]);
		zp = k+0.5f-dt0*0.5f*c*(w[INDEX(i,j,k)]+w[INDEX(i,j,k-1)]);
		if (xp<0.5f) xp=0.5f; if (xp>DimX-1.5f) xp=DimX-1.5f;
		if (yp<0.5f) yp=0.5f; if (yp>DimY-1.5f) yp=DimY-1.5f;
		if (zp<0.5f) zp=0.5f; if (zp>DimZ-1.5f) zp=DimZ-1.5f;

	  i0 = floor(xp); x0 = i0;
	  if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
	  i1=i0+1;
	  if(xp-x0 < 0.5f) s1 = 0.5f+xp-x0;
	  else s1 = xp-x0-0.5f;
	  s0 = 1-s1;
	  j0 = floor(yp); y0 = j0;
	  if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
	  j1=j0+1;
	  if(yp-y0 < 0.5f) t1 = 0.5f+yp-y0;
	  else t1 = yp-y0-0.5f;
	  t0 = 1-t1;
	  k0 = floor(zp); z0 = k0;
	  if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
	  k1=k0+1;
	  if(zp-z0 < 0.5f) r1 = 0.5f+zp-z0;
	  else r1 = zp-z0-0.5f;
	  r0 = 1-r1;


	d[pos] = r0*(s0*(t0*d0[INDEX(i0,j0,k0)]+t1*d0[INDEX(i0,j1,k0)])+
				 s1*(t0*d0[INDEX(i1,j0,k0)]+t1*d0[INDEX(i1,j1,k0)])) +
			 r1*(s0*(t0*d0[INDEX(i0,j0,k1)]+t1*d0[INDEX(i0,j1,k1)])+
				 s1*(t0*d0[INDEX(i1,j0,k1)]+t1*d0[INDEX(i1,j1,k1)]));
}

void VectorField3D::Advect_RK3_WENO_Subcycling(const Voxel &voxel, float dts, char *needed, float *d, float *d0) {

	float delta = voxel.VoxelDelta();

	float dt0 = dts/delta;

	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;
	printf(" subcycling time step = %f and complete time step = %f\n", dts, dt);

	// because advective time step is determined at the beganning of the
	// subcycling, sometimes (although rarely) violation of the CFL condition
	// may occur at some cells starting from the second subcycling step

	// if this does happen, those cells will be marked here by cflex1, and
	// the values of these cells will be updated using a semi-lagrangian scheme,
	// which is guranteed to be stable no matter how big the time step is.

	// since this only happens under some rare conditions, the degraded quality
	// due to semi-lagrangian treatment is not of big concern, especially considering
	// that levelset values will be reinitialized later

	char *cflex1 = new char[DimX*DimY*DimZ];
	SetZero(cflex1);
	Euler(voxel, dts, needed, u0, d0, cflex1);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		// if at any given cell where cfl > 1, use semi-lagrangian advection
	    // scheme to find the value
		if(cflex1[pos])
			SemiLagrangian(voxel, i, j, k, dt0, u0, d0);
	END_FOR
	Euler(voxel, dts, needed, v0, u0);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(cflex1[pos])
			SemiLagrangian(voxel, i, j, k, dt0, v0, u0);
	END_FOR
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//			if(needed[pos])
		u0[pos] = 0.75f * d0[pos] + 0.25f * v0[pos];
	END_FOR
	Euler(voxel, dts, needed, v0, u0);
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(cflex1[pos])
			SemiLagrangian(voxel, i, j, k, dt0, v0, u0);
	END_FOR
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
//			if(needed[pos])
		d[pos] = one_third * d0[pos] + two_third * v0[pos];
	END_FOR
	delete [] cflex1;
//		printf("\n In Euler before swap iter = %d, at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
//				 		n, I, J, K, x0[INDEX(I,J,K)], x[INDEX(I,J,K)]);

//		printf("\n In Euler after swap iter = %d, at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
//						n, I, J, K, x0[INDEX(I,J,K)], x[INDEX(I,J,K)]);


}

void VectorField3D::Advect_FE_WENO(const Voxel &voxel, char *needed, float *x, float *x0) {

	float delta = voxel.VoxelDelta();
	int iterations = int(round(LS_ADV_LIMIT / delta));
	float dts  = dt / iterations;
	printf(" subcycling time step = %f and complete time step = %f\n", dts, dt);

	for(int n = 0; n < iterations; ++n){
		Euler(voxel, dts, needed, x, x0);
//		printf("\n In Euler before swap iter = %d, at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
//				 		n, I, J, K, x0[INDEX(I,J,K)], x[INDEX(I,J,K)]);
		SetEqual(x, u0);
		ExtrapolatePhiIntoObject(voxel, u0, x0);
		SetEqual(u0, x);
		SetEqual(x, x0);
//		printf("\n In Euler after swap iter = %d, at i = %d, j = %d, k = %d  x0 = %f,  x = %f\n",
//						n, I, J, K, x0[INDEX(I,J,K)], x[INDEX(I,J,K)]);

	}
}

void VectorField3D::Advect(const Voxel &voxel, float *x, float *x0) {

	float delta = voxel.VoxelDelta();
	u_int NM = 0, NP = 0;
	SetZero(needRecomputePhi);
	float gamma = FASTMARCH_LIMIT;
	float beta = LS_ADV_LIMIT;
	float gb = (gamma-beta)*(gamma-beta)*(gamma-beta);
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)  && !voxel.InSource(i,j,k)){
			if(x0[INDEX(i,j,k)] <= 0.f)
				++NM;
		}
	END_FOR
	printf("\nbefore advection there are %ld liquid points \n\n", NM);

	SetEqual(x0, w0);
	SetZero(temp);
	u_int N = 0;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && fabs(x0[pos])+E_EPSIL < FASTMARCH_LIMIT ){
			needRecomputePhi[pos] = 1;
			if(fabs(x0[pos]) <= beta){
				temp[pos] = 1.f;
			}
			else{
				temp[pos] = (fabs(x0[pos]) - gamma)*(fabs(x0[pos]) - gamma)
							 * (2 * fabs(x0[pos]) + gamma - 3 * beta) / gb;
			}
			++N;
		}
	END_FOR
	Advect_RK3_WENO(voxel, needRecomputePhi, x, x0);
//	SWAP(x, x0);


	NM = 0;
//	printf("\nthere are %d points liquid->air, %d air->liquid \n", NM, NP);
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			if(x[INDEX(i,j,k)] > 0.f && w0[INDEX(i,j,k)] <= 0.f){
	//			printf("phi changes from liquid to air at (%d, %d, %d) old phi = %f, new phi = %f\n",
	//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
				++NM;
			}
			if(x[INDEX(i,j,k)] <= 0.f && w0[INDEX(i,j,k)] > 0.f){
	//			printf("phi changes from air to liquid at (%d, %d, %d) old phi = %f, new phi = %f\n",
	//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
				++NP;
			}
		}
	END_FOR
	printf("\nafter advection there are %d points liquid->air, %d air->liquid \n\n", NM, NP);
	NM = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && x[INDEX(i,j,k)] <= 0.f){
//			printf(" in Vector3DField::Advect() at (%d, %d, %d) phi = %16.13f \n", i, j, k, x[INDEX(i,j,k)]);
			++NM;
		}
	END_FOR
	printf("\nafter advection there are %ld liquid points \n\n", NM);
	printf("\n there are %d points advected \n\n", N);



}

void VectorField3D::AdvectSubcycling(const Voxel &voxel, float dts, float *x, float *x0) {

	float delta = voxel.VoxelDelta();
	u_int NM = 0, NP = 0;
	char *needed = new char[DimX*DimY*DimZ];
	SetZero(needed);
	float gamma = FASTMARCH_LIMIT;
	float beta = LS_ADV_LIMIT;
	float gb = (gamma-beta)*(gamma-beta)*(gamma-beta);
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)  && !voxel.InSource(i,j,k)){
			if(x0[INDEX(i,j,k)] <= 0.f)
				++NM;
		}
	END_FOR
	printf("\nbefore advection there are %u liquid points \n\n", NM);

	SetEqual(x0, w0);
	SetZero(temp);
	u_int N = 0;
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && fabs(x0[pos])+E_EPSIL < FASTMARCH_LIMIT ){
			needed[pos] = 1;
			if(fabs(x0[pos]) <= beta){
				temp[pos] = 1.f;
			}
			else{
				temp[pos] = (fabs(x0[pos]) - gamma)*(fabs(x0[pos]) - gamma)
							 * (2 * fabs(x0[pos]) + gamma - 3 * beta) / gb;
			}
			++N;
		}
	END_FOR

//	Euler(voxel, dts, needRecomputePhi, x, x0);
	Advect_RK3_WENO_Subcycling(voxel, dts, needed, x, x0);

	NM = 0;
//	printf("\nthere are %d points liquid->air, %d air->liquid \n", NM, NP);
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			if(x[INDEX(i,j,k)] > 0.f && w0[INDEX(i,j,k)] <= 0.f){
	//			printf("phi changes from liquid to air at (%d, %d, %d) old phi = %f, new phi = %f\n",
	//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
				++NM;
			}
			if(x[INDEX(i,j,k)] <= 0.f && w0[INDEX(i,j,k)] > 0.f){
	//			printf("phi changes from air to liquid at (%d, %d, %d) old phi = %f, new phi = %f\n",
	//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
				++NP;
			}
		}
	END_FOR
	printf("\nafter advection there are %d points liquid->air, %d air->liquid \n\n", NM, NP);
	NM = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && x[INDEX(i,j,k)] <= 0.f){
//			printf(" in Vector3DField::Advect() at (%d, %d, %d) phi = %16.13f \n", i, j, k, x[INDEX(i,j,k)]);
			++NM;
		}
	END_FOR
	delete [] needed;
	printf("\nafter advection there are %u liquid points \n\n", NM);
	printf("\n there are %u points advected \n\n", N);
}

//#else
//Semi Lagrangian for Levelset
void VectorField3D::AdvectSemiLagrangian(const Voxel &voxel, float *d, float *d0){


	int i0, j0, k0, i1, j1, k1;
	float s0, t0, s1, t1, dt0;

	float xp, yp, zp;

	int x0, y0, z0;

	float r0, r1;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;
	dt0 = dt/delta;

#ifdef USE_MESH
	float f[8];
	Point pb(0.f), pi(0.f);
#endif

	u_int NM = 0, NP = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)  && !voxel.InSource(i,j,k)){
			if(d0[INDEX(i,j,k)] <= 0.f)
				++NM;
		}
	END_FOR
	printf("\nAdvection of Levelset with Semi-lagrangian Method with dt = %f \n\n", dt);
	printf("\nbefore advection there are %u liquid points \n\n", NM);

	float gamma = FASTMARCH_LIMIT;
	float beta = LS_ADV_LIMIT;
	float inv_gb = 1.f / ((gamma-beta)*(gamma-beta)*(gamma-beta));

	u_int N = 0;
	FOR_EACH_CELL

		if( !voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && fabs(d0[INDEX(i,j,k)])+E_EPSIL < gamma ){
			++N;
			float c = 1.f;
			if(fabs(d0[INDEX(i,j,k)]) > beta)
				c = (fabs(d0[INDEX(i,j,k)]) - gamma)*(fabs(d0[INDEX(i,j,k)]) - gamma)
			         * (2 * fabs(d0[INDEX(i,j,k)]) + gamma - 3 * beta) * inv_gb;

			float dt0c = dt0*0.5f*c;

			xp = i+0.5f-dt0c*(u1[INDEX(i,j,k)]+u1[INDEX(i-1,j,k)]);
			yp = j+0.5f-dt0c*(v1[INDEX(i,j,k)]+v1[INDEX(i,j-1,k)]);
			zp = k+0.5f-dt0c*(w1[INDEX(i,j,k)]+w1[INDEX(i,j,k-1)]);

			xp = xp < 0.5f      ? 0.5f      : xp;
			xp = xp > DimX-1.5f ? DimX-1.5f : xp;
			yp = yp < 0.5f      ? 0.5f      : yp;
			yp = yp > DimY-1.5f ? DimY-1.5f : yp;
			zp = zp < 0.5f      ? 0.5f      : zp;
			zp = zp > DimZ-1.5f ? DimZ-1.5f : zp;

#ifdef USE_MESH
		    pb = voxel.VoxelCenterPosition(i, j, k);
			pi = Point(delta*xp, delta*yp, delta*zp);
			double o, q;
			if(mPhysWorld->SegmentIntersectsMesh(pb, pi, &o, &q)){
				// clip the semi-lagrangian ray
				Vector pib = pi - pb;
				float pibl = o * pib.Length();
				pib = Normalize(pib);
				Point pin = pb + (pibl-EPSIL_C) * pib;
//				Point pin = pb + (o-E_EPSIL)*(pi-pb);
				xp = pin.x * inv_delta;
				yp = pin.y * inv_delta;
				zp = pin.z * inv_delta;
				pi = pin;
			}
#endif

		  i0 = floor(xp); x0 = i0;
		  if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		  i1=i0+1;
		  s1 = xp-x0 < 0.5f ? 0.5f+xp-x0 : xp-x0-0.5f;
		  s0 = 1-s1;
		  j0 = floor(yp); y0 = j0;
		  if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		  j1=j0+1;
		  t1 = yp-y0 < 0.5f ? 0.5f+yp-y0 : yp-y0-0.5f;
		  t0 = 1-t1;
		  k0 = floor(zp); z0 = k0;
		  if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
		  k1=k0+1;
		  r1 = zp-z0 < 0.5f ? 0.5f+zp-z0 : zp-z0-0.5f;
		  r0 = 1-r1;

		  if(i == I && j == J && k == K){
				 printf("\nBefore Advection: at i = %d, j = %d, k = %d  x0 = %f\n",
				 			I,J,K, d0[INDEX(i,j,k)]);


		   }

#ifdef USE_MESH
	 	    AssignEightVelocities(i0, j0, k0, i1, j1, k1, d0, f);

			Point pc = voxel.VoxelCenterPosition(i0, j0, k0);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i0, j0, k0, 0, o, pc, pi, f[0]);
			}
			pc = voxel.VoxelCenterPosition(i0, j1, k0);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i0, j1, k0, 1, o, pc, pi, f[1]);
			}
			pc = voxel.VoxelCenterPosition(i1, j0, k0);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i1, j0, k0, 2, o, pc, pi, f[2]);
			}
			pc = voxel.VoxelCenterPosition(i1, j1, k0);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i1, j1, k0, 3, o, pc, pi, f[3]);
			}
			pc = voxel.VoxelCenterPosition(i0, j0, k1);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i0, j0, k1, 4, o, pc, pi, f[4]);
			}
			pc = voxel.VoxelCenterPosition(i0, j1, k1);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i0, j1, k1, 5, o, pc, pi, f[5]);
			}
			pc = voxel.VoxelCenterPosition(i1, j0, k1);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i1, j0, k1, 6, o, pc, pi, f[6]);
			}
			pc = voxel.VoxelCenterPosition(i1, j1, k1);
			if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &o, &q)){
				ProcessInvisibleNeightborForLevelsetAdvection(voxel, d0,
						i1, j1, k1, 7, o, pc, pi, f[7]);
			}

			d[INDEX(i,j,k)] = r0*(s0*(t0*f[0]+t1*f[1])+
								  s1*(t0*f[2]+t1*f[3])) +
							  r1*(s0*(t0*f[4]+t1*f[5])+
								  s1*(t0*f[6]+t1*f[7]));
#else
			d[INDEX(i,j,k)] = r0*(s0*(t0*d0[INDEX(i0,j0,k0)]+t1*d0[INDEX(i0,j1,k0)])+
				                  s1*(t0*d0[INDEX(i1,j0,k0)]+t1*d0[INDEX(i1,j1,k0)])) +
			                  r1*(s0*(t0*d0[INDEX(i0,j0,k1)]+t1*d0[INDEX(i0,j1,k1)])+
		 		                  s1*(t0*d0[INDEX(i1,j0,k1)]+t1*d0[INDEX(i1,j1,k1)]));
#endif

//			if(d0[INDEX(i,j,k)] < 0.f && d[INDEX(i,j,k)] > 0.f){
//			if(d0[INDEX(i,j,k)] < 0.f && d[INDEX(i,j,k)] > 0.f){
//			if( d0[INDEX(i0,j0,k0)] >= FASTMARCH_LIMIT ||
//				d0[INDEX(i0,j1,k0)] >= FASTMARCH_LIMIT ||
//				d0[INDEX(i1,j0,k0)] >= FASTMARCH_LIMIT ||
//				d0[INDEX(i1,j1,k0)] >= FASTMARCH_LIMIT ||
//				d0[INDEX(i0,j0,k1)] >= FASTMARCH_LIMIT ||
//				d0[INDEX(i0,j1,k1)] >= FASTMARCH_LIMIT ||
//				d0[INDEX(i1,j0,k1)] >= FASTMARCH_LIMIT ||
//				d0[INDEX(i1,j1,k1)] >= FASTMARCH_LIMIT ){
//				float u_c = 0.5f * (u[INDEX(i,j,k)] + u[INDEX(i-1,j,k)]);
//				float v_c = 0.5f * (v[INDEX(i,j,k)] + v[INDEX(i,j-1,k)]);
//				float w_c = 0.5f * (w[INDEX(i,j,k)] + w[INDEX(i,j,k-1)]);
//				printf("d0 = %f, d = %f, u_c = %f, v_c = %f, w_c = %f \n",
//						d0[INDEX(i,j,k)], d[INDEX(i,j,k)], u_c, v_c, w_c);
////				if(voxel.InSolid(i0,j0,k0) && d0[INDEX(i0,j0,k0)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i0, j0, k0, d0[INDEX(i0,j0,k0)]);
////				if(voxel.InSolid(i1,j0,k0) && d0[INDEX(i1,j0,k0)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i1, j0, k0, d0[INDEX(i1,j0,k0)]);
////				if(voxel.InSolid(i0,j1,k0) && d0[INDEX(i0,j1,k0)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i0, j1, k0, d0[INDEX(i0,j1,k0)]);
////				if(voxel.InSolid(i1,j1,k0) && d0[INDEX(i1,j1,k0)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i1, j1, k0, d0[INDEX(i1,j1,k0)]);
////				if(voxel.InSolid(i0,j0,k1) && d0[INDEX(i0,j0,k1)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i0, j0, k1, d0[INDEX(i0,j0,k1)]);
////				if(voxel.InSolid(i1,j0,k1) && d0[INDEX(i1,j0,k1)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i1, j0, k1, d0[INDEX(i1,j0,k1)]);
////				if(voxel.InSolid(i0,j1,k1) && d0[INDEX(i0,j1,k1)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i0, j1, k1, d0[INDEX(i0,j1,k1)]);
////				if(voxel.InSolid(i1,j1,k1) && d0[INDEX(i1,j1,k1)] > 0.f)
//					printf("for (%d, %d, %d) at (%d, %d, %d) phi = %f > 0 \n",
//							i, j, k, i1, j1, k1, d0[INDEX(i1,j1,k1)]);
//					exit(1);
//			}
			 if(i == I && j == J && k == K){
				 printf(" xp = %f, yp = %f, zp = %f \n", xp, yp, zp);
				 printf(" i0 = %d, i1 = %d, j0 = %d, j1 = %d, k0 = %d, k1 = %d \n",
				          i0, i1, j0, j1, k0, k1);
				 printf(" t0 = %f, t1 = %f, s0 = %f, s1 = %f, r0 = %f, r1 = %f\n",
						 t0, t1, s0, s1, r0, r1);
				 printf("After Advection: at i = %d, j = %d, k = %d, x = %f \n\n",
							I,J,K, d[INDEX(i,j,k)]);

			}
		}
		else
			d[INDEX(i,j,k)] = d0[INDEX(i,j,k)];
	END_FOR
	NM = 0;
	NP = 0;
	//	printf("\nthere are %d points liquid->air, %d air->liquid \n", NM, NP);
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k)){
			if(d[INDEX(i,j,k)] > 0.f && d0[INDEX(i,j,k)] <= 0.f){
	//			printf("phi changes from liquid to air at (%d, %d, %d) old phi = %f, new phi = %f\n",
	//					i, j, k, x0[INDEX(i,j,k)], x[INDEX(i,j,k)]);
				++NM;
			}
			if(d[INDEX(i,j,k)] <= 0.f && d0[INDEX(i,j,k)] > 0.f){
//					printf("phi changes from air to liquid at (%d, %d, %d) old phi = %f, new phi = %f\n",
//							i, j, k, d0[INDEX(i,j,k)], d[INDEX(i,j,k)]);
				++NP;
			}
		}
	END_FOR
	printf("\nafter advection there are %d points liquid->air, %d air->liquid \n\n", NM, NP);
	NM = 0;
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
			if(d[INDEX(i,j,k)] < 0.f)
				++NM;
		}
	END_FOR
	printf("\nafter advection there are %u liquid points \n\n", NM);
	printf("\n there are %u points advected \n\n", N);
	//set_bnd ( N, b, d );
}
#endif


void VectorField3D::VortexForce(const Voxel &voxel, VortexParticles *VP){
	SetZero(u0);
	SetZero(v0);
	SetZero(w0);
	VP->VortexForce(voxel, u0, v0, w0, phi_u, phi_v, phi_w,
			phi_u_obj, phi_v_obj, phi_w_obj);

}

void VectorField3D::StretchVortex(const Voxel &voxel, VortexParticles *VP){

	float delta = voxel.VoxelDelta();

	SetZero(u0);
	SetZero(u1);
	SetZero(v0);
	SetZero(v1);
	SetZero(w0);
	SetZero(w1);
	float *u2, *v2, *w2;
	u2 = new float[DimX*DimY*DimZ];
	memset(u2, 0, DimX*DimY*DimZ*sizeof(float));
	v2 = new float[DimX*DimY*DimZ];
	memset(v2, 0, DimX*DimY*DimZ*sizeof(float));
	w2 = new float[DimX*DimY*DimZ];
	memset(w2, 0, DimX*DimY*DimZ*sizeof(float));

	FOR_EACH_CELL
		TentativeGridPoint tp(0.f, i, j, k);
		int n;
		if(!voxel.InSolid(i,j,k) || voxel.CloseToNonSolid(tp, n)){
			u0[INDEX(i,j,k)] = (u[INDEX(i,j,k)] - u[INDEX(i-1,j,k)]) / delta;
			v0[INDEX(i,j,k)] = (u[INDEX(i,j+1,k)] + u[INDEX(i-1,j+1,k)] -
					            u[INDEX(i,j-1,k)] - u[INDEX(i-1,j-1,k)]) / (4*delta);
			w0[INDEX(i,j,k)] = (u[INDEX(i,j,k+1)] + u[INDEX(i-1,j,k+1)] -
								u[INDEX(i,j,k-1)] - u[INDEX(i-1,j,k-1)]) / (4*delta);

			u1[INDEX(i,j,k)] = (v[INDEX(i+1,j,k)] + v[INDEX(i+1,j-1,k)] -
					            v[INDEX(i-1,j,k)] - v[INDEX(i-1,j-1,k)] ) / (4*delta);
			v1[INDEX(i,j,k)] = (v[INDEX(i,j,k)] - v[INDEX(i,j-1,k)]) / delta;
			w1[INDEX(i,j,k)] = (v[INDEX(i,j,k+1)] + v[INDEX(i,j-1,k+1)] -
								v[INDEX(i,j,k-1)] - v[INDEX(i,j-1,k-1)]) / (4*delta);

			u2[INDEX(i,j,k)] = (w[INDEX(i+1,j,k)] + w[INDEX(i+1,j,k-1)] -
					            w[INDEX(i-1,j,k)] - w[INDEX(i-1,j,k-1)] ) / (4*delta);
			v2[INDEX(i,j,k)] = (w[INDEX(i,j+1,k)] + w[INDEX(i,j+1,k-1)] -
			                    w[INDEX(i,j-1,k)] - w[INDEX(i,j-1,k-1)]) / (4*delta);
			w2[INDEX(i,j,k)] = (w[INDEX(i,j,k)] - w[INDEX(i,j,k-1)]) / delta;
		}
	END_FOR

	list<VortexParticle> &vortexParticles = VP->Particles;
	list<VortexParticle>::iterator iter_waterparticle;
	for(iter_waterparticle = vortexParticles.begin();
		iter_waterparticle != vortexParticles.end();
		){
		VortexParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		Vector vor_old = ps.Vorticity();
		float mag_old = ps.Magnitude();
		Vector vor_new;
		float r, s, t;
		r = TriInterp(voxel, pos, u0);
		s = TriInterp(voxel, pos, v0);
		t = TriInterp(voxel, pos, w0);
		vor_new.x = vor_old.x + dt * Dot(vor_old, Vector(r,s,t));
		r = TriInterp(voxel, pos, u1);
		s = TriInterp(voxel, pos, v1);
		t = TriInterp(voxel, pos, w1);
		vor_new.y = vor_old.y + dt * Dot(vor_old, Vector(r,s,t));
		r = TriInterp(voxel, pos, u2);
		s = TriInterp(voxel, pos, v2);
		t = TriInterp(voxel, pos, w2);
		vor_new.z = vor_old.z + dt * Dot(vor_old, Vector(r,s,t));

		Vector d = Normalize(vor_new);
		float mag_new = vor_new.Length();
		if(mag_new > mag_old)
			ps.SetDirection(d);
		else{
			ps.SetMagnitude(mag_new);
			ps.SetDirection(d);
		}
		++iter_waterparticle;
	}

	delete [] u2;
	delete [] v2;
	delete [] w2;
}


void VectorField3D::AdjustOneVelocity(const Voxel &voxel, int vel_index,
							list<WaterParticle> &waterParticles, float *vel){

	float delta = voxel.VoxelDelta();
	float volume = delta * delta * delta;
	list<WaterParticle>::iterator iter_waterparticle;
	int *cellParticles = new int[DimX * DimY *DimZ];
	memset(cellParticles, 0, DimX * DimY *DimZ * sizeof(int));

	typedef list<WaterParticle *> PBUFFER;
	map<u_int, PBUFFER> vCells;
	map<u_int, PBUFFER>::iterator iter_cell;
	list<WaterParticle *>::iterator iter_escape;

	for(iter_waterparticle = waterParticles.begin();
		iter_waterparticle != waterParticles.end();
		){
		WaterParticle &wps = *iter_waterparticle;
		Point pos = wps.Position();
		TentativeGridPoint tp = voxel.ContainsPointInVelCell(pos, vel_index);
		++cellParticles[INDEX(tp.ii, tp.jj, tp.kk)];
		++iter_waterparticle;
	}
	FOR_EACH_CELL
		if(cellParticles[INDEX(i,j,k)] > 0){
			PBUFFER thisheap;
			vCells.insert(make_pair(INDEX(i,j,k), thisheap));
		}
	END_FOR

	for(iter_waterparticle = waterParticles.begin();
		iter_waterparticle != waterParticles.end();
		){
		WaterParticle &wps = *iter_waterparticle;
		Point pos = wps.Position();
		TentativeGridPoint tp = voxel.ContainsPointInVelCell(pos, vel_index);
		iter_cell = vCells.find(INDEX(tp.ii, tp.jj, tp.kk));
		if(iter_cell != vCells.end()){
			PBUFFER &thisheap = iter_cell->second;
			thisheap.push_back(&wps);
		}
		++iter_waterparticle;
	}

	FOR_EACH_CELL
		iter_cell = vCells.find(INDEX(i, j, k));
		if(iter_cell != vCells.end()){
			PBUFFER &thisheap = iter_cell->second;
			// estimating volume loss within the given cell
			// using the Monte Carlo method
			int hit = 0;
			for(int n=0; n < MONTE_CARLO_ITERATIONS; ++n){
				Point p;
				voxel.PlacePoints(i,j,k,p,vel_index);
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
			float total = (volume - V) * vel[INDEX(i,j,k)];
			float total_volume = volume - V;
			for(iter_escape = thisheap.begin();
				iter_escape != thisheap.end();
				++iter_escape){
				WaterParticle *wps = *iter_escape;
				float vp;
				if(vel_index == 1)
					vp = wps->UVelocity();
				else if(vel_index == 2)
					vp = wps->VVelocity();
				else if(vel_index == 3)
					vp = wps->WVelocity();
				float r = wps->Radius();
				float volumep = 4 * M_PI / 3 * r * r * r;
				total += volumep * vp;
				total_volume += volumep;
			}
			vel[INDEX(i,j,k)] = total / total_volume;
		}
	END_FOR

	delete [] cellParticles;
}

void VectorField3D::AdjustVelocity(const Voxel &voxel, list<WaterParticle> &waterParticles){

	AdjustOneVelocity(voxel, 1, waterParticles, u);
	AdjustOneVelocity(voxel, 2, waterParticles, v);
	AdjustOneVelocity(voxel, 3, waterParticles, w);

}

void VectorField3D::Advect(const Voxel &voxel, list<WaterParticle> &waterParticles) const{
	Vector G(0, 0, -g);
	list<WaterParticle>::iterator iter_waterparticle;
	for(iter_waterparticle = waterParticles.begin();
		iter_waterparticle != waterParticles.end();
		){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		float r = ps.Radius();
		float colldist = ps.CollisionDistance();
		Vector particle_vel = ps.Velocity();
		TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos, phi_c_obj);
		particle_vel += dt * G;
		Vector o_normal = ObjectNormalAt(voxel, pos);
		Point new_p = pos + dt * particle_vel;
		TentativeGridPoint new_tp = voxel.ContainsPoint(new_p);
		if(!CheckIndex(new_tp.ii, new_tp.jj, new_tp.kk))
			iter_waterparticle = waterParticles.erase(iter_waterparticle);
		else{
			float phi_i_new = TriInterp(voxel, new_p, phi_c_obj);

//			if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
//				(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
//				(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
			if( phi_i_new >= 0.f || (phi_i_new < 0.f && fabs(phi_i_new) < colldist)){
				float v_normal = Dot(particle_vel, o_normal);
				if(v_normal < 0.f) // particle will drift towards the object
					particle_vel -= v_normal * o_normal;
			}
			pos += dt * particle_vel;
			phi_i_new = TriInterp(voxel, pos, phi_c_obj);
			new_tp = voxel.ContainsPoint(pos);
//			if(phi_i_new > 0.f && voxel.InSolid(new_tp))
			if(phi_i_new > 0.f)
				// if water particles still cannot avoid penetrating the solid object,
				// they are deleted
				iter_waterparticle = waterParticles.erase(iter_waterparticle);
			else{
				ps.SetParticle(pos, r);
				ps.SetVelocity(particle_vel);
				++iter_waterparticle;
			}
		}
	}
	printf("end of Advect water particles\n");
}
#ifdef COUPLED_SPH
void VectorField3D::AdvectSPHParticles(const Voxel *voxel, float dts, list<WaterParticle> &waterParticles) const{

	list<WaterParticle>::iterator iter_waterparticle;
	for(iter_waterparticle = waterParticles.begin();
		iter_waterparticle != waterParticles.end();
		){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		float r = ps.Radius();
		float colldist = ps.CollisionDistance();
		Vector particle_vel = ps.Velocity();
//		TentativeGridPoint o_tp = voxel->ContainsPoint(pos);
//		float phi_i = TriInterp(*voxel, pos, phi_c_obj);
//		Vector o_normal = ObjectNormalAt(*voxel, pos);
		Point new_p = pos + dts * particle_vel;
		TentativeGridPoint new_tp = voxel->ContainsPoint(new_p);
		if(!CheckIndex(new_tp.ii, new_tp.jj, new_tp.kk))
			iter_waterparticle = waterParticles.erase(iter_waterparticle);
		else{
			float phi_i_new = TriInterp(*voxel, new_p, phi_c_obj);

//			if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
//				(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
//				(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
			if( phi_i_new >= 0.f || (phi_i_new < 0.f && fabs(phi_i_new) < colldist)){
				Vector o_normal = ObjectNormalAt(*voxel, new_p);
				float v_normal = Dot(particle_vel, o_normal);
				if(v_normal < 0.f) // particle will drift towards the object
					particle_vel -= v_normal * o_normal;
			}
			pos += dts * particle_vel;
			phi_i_new = TriInterp(*voxel, pos, phi_c_obj);
//			new_tp = voxel.ContainsPoint(pos);
//			if(phi_i_new > 0.f && voxel.InSolid(new_tp))
			if(phi_i_new > 0.f)
				// if water particles still cannot avoid penetrating the solid object,
				// they are deleted
				iter_waterparticle = waterParticles.erase(iter_waterparticle);
			else{
				ps.SetParticle(pos, r);
				ps.SetVelocity(particle_vel);
				++iter_waterparticle;
			}
		}
	}
	printf("end of Advect SPH particles\n");
}
#endif

void VectorField3D::AdvectSubcycling(const Voxel &voxel, float dts,
		WaterParticleList &waterParticles) const{
	Vector G(0, 0, -g);
	WaterParticleList::iterator iter_waterparticle;
	for(iter_waterparticle = waterParticles.begin();
		iter_waterparticle != waterParticles.end();
		){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		float r = ps.Radius();
		float colldist = ps.CollisionDistance();
		Vector particle_vel = ps.Velocity();
		TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos, phi_c_obj);
		particle_vel += dts * G;
		Point new_p = pos + dts * particle_vel;
		TentativeGridPoint new_tp = voxel.ContainsPoint(new_p);
		if(!CheckIndex(new_tp.ii, new_tp.jj, new_tp.kk))
			iter_waterparticle = waterParticles.erase(iter_waterparticle);
		else{
			float phi_i_new = TriInterp(voxel, new_p, phi_c_obj);

//			if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
//				(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
//				(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
			if( phi_i_new >= 0.f || (phi_i_new < 0.f && fabs(phi_i_new) < colldist)){
				Vector o_normal = ObjectNormalAt(voxel, new_p);
				float v_normal = Dot(particle_vel, o_normal);
				if(v_normal < 0.f) // particle will drift towards the object
					particle_vel -= v_normal * o_normal;
				Point projP = new_p + (colldist + phi_i_new) * o_normal;
				new_p = projP;
			}
//			pos += dts * particle_vel;
//			phi_i_new = TriInterp(voxel, pos, phi_c_obj);
//			new_tp = voxel.ContainsPoint(pos);
////			if(phi_i_new > 0.f && voxel.InSolid(new_tp))
//			if(phi_i_new > 0.f)
//				// if water particles still cannot avoid penetrating the solid object,
//				// they are deleted
//				iter_waterparticle = waterParticles.erase(iter_waterparticle);
//			else{
			ps.SetParticle(new_p, r);
			ps.SetVelocity(particle_vel);
			++iter_waterparticle;
//			}
		}
	}
	printf("end of Advect water particles\n");
}

void VectorField3D::Advect(const Voxel &voxel, VortexParticles *VP ) const{

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;

	float xp, yp, zp;
	float old_xp, old_yp, old_zp;
	float up, vp, wp;

	list<VortexParticle> &Particles = VP->Particles;
	list<VortexParticle>::iterator particle;

	for(particle=Particles.begin(); particle != Particles.end(); ){
		VortexParticle &ps = *particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float colldist = ps.CollisionDistance();
		old_xp = xp = pos.x;
		old_yp = yp = pos.y;
		old_zp = zp = pos.z;
		float o_xp = old_xp;
		float o_yp = old_yp;
		float o_zp = old_zp;
		TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos, phi_c_obj);
		TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
		Vector o_normal = ObjectNormalAt(voxel, pos);
		float v_normal;
		Vector v_para;
		xp = old_xp + dt*up; yp = old_yp + dt*vp; zp = old_zp + dt*wp;
		float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
		if( phi_i_new >= 0.f || (phi_i_new < 0.f && fabs(phi_i_new) < colldist)){
			v_normal = Dot(Vector(up, vp, wp), o_normal);
			if(v_normal < 0.f) // particle will drift towards the object
				v_para = Vector(up, vp, wp) - v_normal * o_normal;
			else
				v_para = Vector(up, vp, wp);
			xp = old_xp + dt*v_para.x;
			yp = old_yp + dt*v_para.y;
			zp = old_zp + dt*v_para.z;

		}
		ps.SetParticle(Point(xp,yp,zp), radius);
		++particle;
	}

	printf("end of Advect vortex particles\n");
}


void VectorField3D::AdvectRK3(const Voxel &voxel, list<Particle> &particles, int sign) const{

	float xp, yp, zp;
	float old_xp, old_yp, old_zp;
	float up, vp, wp, up1, vp1, wp1;
	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;
	int N=0;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;
	int iterations = int(round(LS_ADV_LIMIT / delta));
	float dts  = dt / iterations;
	list<Particle>::iterator particle;

	printf("start advecting  %u particles with RK3 : \n",particles.size());
	printf(" subcycling time step = %f and complete time step = %f\n", dts, dt);
	if(sign > 0){ // air particle
		for(particle=particles.begin(); particle != particles.end(); ){
		//		Point pos = particles[m].Postion();
		//		float radius = particles[m].Radius();
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
			float phi_i = TriInterp(voxel, pos, phi_c_obj);
			if(phi_i > 0.f){   // in solid object
//			if(phi_i > 0.f && voxel.InSolid(o_tp)){
				// we only get rid of air particles
				particle = particles.erase(particle);
			}
			else if(voxel.InSource(o_tp))
				particle = particles.erase(particle);
			else{  // particle close to a solid object
	//			float colldist = ps.CollisionDistance();
			bool inside_obj = false;
			for(int n = 0; n < iterations; ++n){
				old_xp = pos.x;
				old_yp = pos.y;
				old_zp = pos.z;
				xp     = pos.x;
				yp     = pos.y;
				zp     = pos.z;
				float o_xp = old_xp;
				float o_yp = old_yp;
				float o_zp = old_zp;
				TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
				xp = old_xp + dts*up; yp = old_yp + dts*vp; zp = old_zp + dts*wp;
				TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
				old_xp = xp + dts*up1;
				old_yp = yp + dts*vp1;
				old_zp = zp + dts*wp1;
				xp = 0.75f * o_xp + 0.25f * old_xp;
				yp = 0.75f * o_yp + 0.25f * old_yp;
				zp = 0.75f * o_zp + 0.25f * old_zp;
				TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
				old_xp = xp + dts*up1;
				old_yp = yp + dts*vp1;
				old_zp = zp + dts*wp1;
				xp = one_third * o_xp + two_third * old_xp;
				yp = one_third * o_yp + two_third * old_yp;
				zp = one_third * o_zp + two_third * old_zp;
				float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
				if(phi_i_new > 0.f){
					particle = particles.erase(particle);
					inside_obj = true;
					break;
				}
				else{
//					ps.SetParticle(Point(xp,yp,zp), radius);
//					++particle;
					pos = Point(xp, yp, zp);
				}
			}
			if(!inside_obj){
				ps.SetParticle(pos, radius);
				++particle;
			}
			}
		}
	}
	else{	// water particles
		for(particle=particles.begin(); particle != particles.end(); ){
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			float colldist = ps.CollisionDistance();
			for(int n = 0; n < iterations; ++n){
				old_xp = pos.x;
				old_yp = pos.y;
				old_zp = pos.z;
				xp     = pos.x;
				yp     = pos.y;
				zp     = pos.z;
				float o_xp = old_xp;
				float o_yp = old_yp;
				float o_zp = old_zp;
				TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
		//		float phi_i = TriInterp(voxel, pos, phi_c_obj);
				TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
		//		    	TentativeGridPoint object = voxel.ClosestObjectCell(o_tp);
		//		    	Vector o_normal = ObjectNormalAt(delta, object.ii, object.jj, object.kk);
		//		    	Vector o_normal = ObjectNormalAt(delta, o_tp.ii, o_tp.jj, o_tp.kk);
		//		    	Vector o_normal = ObjectNormalAt(voxel, o_tp.ii, o_tp.jj, o_tp.kk);
		//		Vector o_normal = ObjectNormalAt(voxel, pos);
		//		    	if(o_tp.ii == I && o_tp.jj == J && o_tp.kk == K){
		//		    		printf("object normal = (%f, %f, %f) \n",
		//		    				o_normal.x, o_normal.y, o_normal.z);
		//		    	}
				xp = old_xp + dts*up; yp = old_yp + dts*vp; zp = old_zp + dts*wp;
				TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
				old_xp = xp + dts * up1;
				old_yp = yp + dts * vp1;
				old_zp = zp + dts * wp1;
				xp = 0.75f * o_xp + 0.25f * old_xp;
				yp = 0.75f * o_yp + 0.25f * old_yp;
				zp = 0.75f * o_zp + 0.25f * old_zp;
				TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
				old_xp = xp + dts*up1;
				old_yp = yp + dts*vp1;
				old_zp = zp + dts*wp1;
				xp = one_third * o_xp + two_third * old_xp;
				yp = one_third * o_yp + two_third * old_yp;
				zp = one_third * o_zp + two_third * old_zp;
				pos = Point(xp, yp, zp);
				float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
		//		if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
		//			(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
		//			(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
				if( phi_i_new >= 0.f ||	fabs(phi_i_new) < colldist ){
					Vector o_normal = ObjectNormalAt(voxel, pos);
					Point projP = pos + (colldist + phi_i_new) * o_normal;
					pos = projP;
				}
			}
			ps.SetParticle(pos, radius);
			++particle;
		}
	}

	printf("end advection particles with RK3: \n");

	printf("There are %d particles in cell(%d,%d,%d) \n ",
			N, I, J, K);


}

void VectorField3D::AdvectRK3NoSubcycling(const Voxel &voxel, list<Particle> &particles, int s) const{

	float xp, yp, zp;
	float old_xp, old_yp, old_zp;
	float up, vp, wp, up1, vp1, wp1;
	float one_third = 1.f / 3.f;
	float two_third = 2.f / 3.f;
	int N=0;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;
	list<Particle>::iterator particle;

	printf("start advecting  %u particles with RK3 no subcycling : \n",particles.size());
	if(s > 0){ // air particle
		for(particle=particles.begin(); particle != particles.end(); ){
		//		Point pos = particles[m].Postion();
		//		float radius = particles[m].Radius();
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
			float phi_i = TriInterp(voxel, pos, phi_c_obj);
			if(phi_i > 0.f){   // in solid object
//			if(phi_i > 0.f && voxel.InSolid(o_tp)){
				// we only get rid of air particles
				particle = particles.erase(particle);
			}
			else if(voxel.InSource(o_tp))
				particle = particles.erase(particle);
			else{  // particle close to a solid object
	//			float colldist = ps.CollisionDistance();
				old_xp = pos.x;
				old_yp = pos.y;
				old_zp = pos.z;
				xp     = pos.x;
				yp     = pos.y;
				zp     = pos.z;
				float o_xp = old_xp;
				float o_yp = old_yp;
				float o_zp = old_zp;
				TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
				xp = old_xp + dt*up; yp = old_yp + dt*vp; zp = old_zp + dt*wp;
				TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
				old_xp = xp + dt*up1;
				old_yp = yp + dt*vp1;
				old_zp = zp + dt*wp1;
				xp = 0.75f * o_xp + 0.25f * old_xp;
				yp = 0.75f * o_yp + 0.25f * old_yp;
				zp = 0.75f * o_zp + 0.25f * old_zp;
				TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
				old_xp = xp + dt*up1;
				old_yp = yp + dt*vp1;
				old_zp = zp + dt*wp1;
				xp = one_third * o_xp + two_third * old_xp;
				yp = one_third * o_yp + two_third * old_yp;
				zp = one_third * o_zp + two_third * old_zp;
				float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
				if(phi_i_new > 0.f){
					particle = particles.erase(particle);
				}
				else{
					pos = Point(xp, yp, zp);
					ps.SetParticle(pos, radius);
					++particle;
				}
			}
		}
	}
	else{	// water particles
		for(particle=particles.begin(); particle != particles.end(); ){
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			float colldist = ps.CollisionDistance();
			old_xp = pos.x;
			old_yp = pos.y;
			old_zp = pos.z;
			xp     = pos.x;
			yp     = pos.y;
			zp     = pos.z;
			float o_xp = old_xp;
			float o_yp = old_yp;
			float o_zp = old_zp;
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
	//		float phi_i = TriInterp(voxel, pos, phi_c_obj);
			TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
	//		    	TentativeGridPoint object = voxel.ClosestObjectCell(o_tp);
	//		    	Vector o_normal = ObjectNormalAt(delta, object.ii, object.jj, object.kk);
	//		    	Vector o_normal = ObjectNormalAt(delta, o_tp.ii, o_tp.jj, o_tp.kk);
	//		    	Vector o_normal = ObjectNormalAt(voxel, o_tp.ii, o_tp.jj, o_tp.kk);
	//		Vector o_normal = ObjectNormalAt(voxel, pos);
	//		    	if(o_tp.ii == I && o_tp.jj == J && o_tp.kk == K){
	//		    		printf("object normal = (%f, %f, %f) \n",
	//		    				o_normal.x, o_normal.y, o_normal.z);
	//		    	}
			xp = old_xp + dt*up; yp = old_yp + dt*vp; zp = old_zp + dt*wp;
			TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
			old_xp = xp + dt * up1;
			old_yp = yp + dt * vp1;
			old_zp = zp + dt * wp1;
			xp = 0.75f * o_xp + 0.25f * old_xp;
			yp = 0.75f * o_yp + 0.25f * old_yp;
			zp = 0.75f * o_zp + 0.25f * old_zp;
			TriInterp(voxel, Point(xp,yp,zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
			old_xp = xp + dt * up1;
			old_yp = yp + dt * vp1;
			old_zp = zp + dt * wp1;
			xp = one_third * o_xp + two_third * old_xp;
			yp = one_third * o_yp + two_third * old_yp;
			zp = one_third * o_zp + two_third * old_zp;
			pos = Point(xp, yp, zp);
			float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
	//		if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
	//			(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
	//			(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
			if( phi_i_new >= 0.f ||	fabs(phi_i_new) < colldist ){
				Vector o_normal = ObjectNormalAt(voxel, pos);
				Point projP = pos + (colldist + phi_i_new) * o_normal;
				pos = projP;
			}
			ps.SetParticle(pos, radius);
			++particle;
		}
	}

	printf("end advection particles with RK3 no subcycling: \n");

	printf("There are %d particles in cell(%d,%d,%d) \n ",
			N, I, J, K);


}

bool VectorField3D::
	ClipSemiLagrangianRay(const Point &pold, Point &pos) const{
	double o, q;
	if(mPhysWorld->SegmentIntersectsMesh(pold, pos, &o, &q)){
		// clip the semi-lagrangian ray
		Vector pib = pos - pold;
		float pibl = o * pib.Length();
		pib = Normalize(pib);
		Point pin = pold + (pibl-EPSIL_C) * pib;
		pos = pin;
		return true;
	}
	else
		return false;
}

void VectorField3D::AdvectRK3Subcycling(const Voxel &voxel,	float dts,
		      ParticleList &particles,
//		      list<Particle> &particles,
		      int s) const{


	const float a1 = 2.f / 9.f;
	const float a2 = 3.f / 9.f;
	const float a3 = 4.f / 9.f;

	int N = 0;


	Vector pvelmax(0.f);
	Point poldmax(0.f);
	Point pmax(0.f);
	TentativeGridPoint tpoldmax(0.f, 0, 0, 0);
	TentativeGridPoint tpmax(0.f, 0, 0, 0);

	Vector k1(0.f), k2(0.f), k3(0.f);

#ifdef USE_MESH
	double r, t;
#endif

	Point pold;

	float dts1 = 0.5f  * dts;
	float dts2 = 0.75f * dts;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;

//	list<Particle>::iterator particle;
	ParticleList::iterator particle;

	if(s < 0)
		printf("\nstart advecting %u NEGATIVE particles with RK3 subcycling : \n", particles.size());
	else
		printf("\nstart advecting %u POSITIVE particles with RK3 subcycling : \n", particles.size());

	if(s > 0){ // air particle
		for( particle  = particles.begin();
		     particle != particles.end();
		     ){
			Particle &ps = *particle;
			Point pos = ps.Position();
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
			if(!CheckIndex(o_tp.ii, o_tp.jj, o_tp.kk))
				particle = particles.erase(particle);
			else{
			float phi_i = TriInterp(voxel, pos, phi_c_obj);
			if(phi_i > 0.f){   // in solid object
				// we only get rid of air particles
				particle = particles.erase(particle);
			}
			else if(voxel.InSource(o_tp))
				particle = particles.erase(particle);
			else{
#ifdef USE_MESH
				pold = pos;
				TriInterp(voxel, pos, k1.x, k1.y, k1.z, inv_delta, u1, v1, w1);
				pos = pold + dts1 * k1;
				if(mPhysWorld->SegmentIntersectsMesh(pold, pos, &r, &t))
					particle = particles.erase(particle);
				else{
					TriInterp(voxel, pos, k2.x, k2.y, k2.z, inv_delta, u1, v1, w1);
					pos = pold + dts2 * k2;
					if(mPhysWorld->SegmentIntersectsMesh(pold, pos, &r, &t))
						particle = particles.erase(particle);
					else{
						TriInterp(voxel, pos, k3.x, k3.y, k3.z, inv_delta, u1, v1, w1);
						pos = pold + dts * (a1 * k1 + a2 * k2 + a3 * k3);
						float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
						TentativeGridPoint tp = voxel.ContainsPoint(pos);
						if(phi_i_new > 0.f || voxel.InSource(tp) ||
						   mPhysWorld->SegmentIntersectsMesh(pold, pos, &r, &t)){
							particle = particles.erase(particle);
						}
						else{
							if(tp.ii == I && tp.jj == J && tp.kk == K)
								++N;
							ps.SetParticle(pos);
							++particle;
						}
					}
				}
#else
				pold = pos;
				TriInterp(voxel, pos, k1.x, k1.y, k1.z, inv_delta, u1, v1, w1);
				pos = pold + dts1 * k1;
				TriInterp(voxel, pos, k2.x, k2.y, k2.z, inv_delta, u1, v1, w1);
				pos = pold + dts2 * k2;
				TriInterp(voxel, pos, k3.x, k3.y, k3.z, inv_delta, u1, v1, w1);
				Vector kave = a1 * k1 + a2 * k2 + a3 * k3;
				pos = pold + dts * kave;
				if(kave.LengthSquared() > pvelmax.LengthSquared()){
					pvelmax = kave;
					tpmax   = voxel.ContainsPoint(pos);
					tpoldmax= voxel.ContainsPoint(pold);
					pmax    = pos;
					poldmax = pold;
				}
				TentativeGridPoint tp = voxel.ContainsPoint(pos);
				if(!CheckIndex(tp.ii, tp.jj, tp.kk))
					particle = particles.erase(particle);
				else{
					float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
					if(phi_i_new > 0.f || voxel.InSource(tp))
						particle = particles.erase(particle);
					else{
						if(tp.ii == I && tp.jj == J && tp.kk == K)
							++N;
						ps.SetParticle(pos);
						++particle;
					}
				}
#endif
			}
		}
		}
	}
	else{	// water particles
		for(particle  = particles.begin();
			particle != particles.end();
			){
			Particle &ps = *particle;
			Point pos = ps.Position();
			float colldist = ps.CollisionDistance();
//			old_xp = pos.x;
//			old_yp = pos.y;
//			old_zp = pos.z;
//			xp     = pos.x;
//			yp     = pos.y;
//			zp     = pos.z;
//			float o_xp = old_xp;
//			float o_yp = old_yp;
//			float o_zp = old_zp;
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
			if(!CheckIndex(o_tp.ii, o_tp.jj, o_tp.kk))
				particle = particles.erase(particle);
			else{
#ifdef USE_MESH
			bool colliding = false;
			pold = pos;
			TriInterp(voxel, pos, k1.x, k1.y, k1.z, inv_delta, u1, v1, w1);
			pos = pold + dts1 * k1;
			colliding = ClipSemiLagrangianRay(pold, pos);
			TriInterp(voxel, pos, k2.x, k2.y, k2.z, inv_delta, u1, v1, w1);
			pos = pold + dts2 * k2;
			colliding = ClipSemiLagrangianRay(pold, pos);
			TriInterp(voxel, pos, k3.x, k3.y, k3.z, inv_delta, u1, v1, w1);
			Vector kave = a1 * k1 + a2 * k2 + a3 * k3;
			Vector objNormal(0.f);
			float distToObj = colldist;
			if(colliding ||
			   mPhysWorld->PointDistToMesh(pold, distToObj, objNormal)){
			   float velNToObj = Dot(kave, objNormal);
			   if(velNToObj < 0.f)
					kave = kave - Dot(kave, objNormal) * objNormal;
			}
			pos = pold + dts * kave;
			distToObj = colldist;
			if(mPhysWorld->PointDistToMesh(pos, distToObj, objNormal)) { // collide particle with objects
				pold = pos;
				pos = pold + (colldist - distToObj) * objNormal;
				double r, t;
				if(mPhysWorld->SegmentIntersectsMesh(pold, pos, &r, &t))
					pos = pold;
			}
#else
			pold = pos;
			TriInterp(voxel, pos, k1.x, k1.y, k1.z, inv_delta, u1, v1, w1);
			pos = pold + dts1 * k1;
			TriInterp(voxel, pos, k2.x, k2.y, k2.z, inv_delta, u1, v1, w1);
			pos = pold + dts2 * k2;
			TriInterp(voxel, pos, k3.x, k3.y, k3.z, inv_delta, u1, v1, w1);
			Vector kave = a1 * k1 + a2 * k2 + a3 * k3;
			pos = pold + dts * kave;
			if(kave.LengthSquared() > pvelmax.LengthSquared()){
				pvelmax = kave;
				tpmax   = voxel.ContainsPoint(pos);
				tpoldmax= voxel.ContainsPoint(pold);
				pmax    = pos;
				poldmax = pold;
			}
#endif
//		float phi_i = TriInterp(voxel, pos, phi_c_obj);
//			TriInterp(up, vp, wp, xp/delta, yp/delta, zp/delta, u1, v1, w1);
//		    	TentativeGridPoint object = voxel.ClosestObjectCell(o_tp);
//		    	Vector o_normal = ObjectNormalAt(delta, object.ii, object.jj, object.kk);
//		    	Vector o_normal = ObjectNormalAt(delta, o_tp.ii, o_tp.jj, o_tp.kk);
//		    	Vector o_normal = ObjectNormalAt(voxel, o_tp.ii, o_tp.jj, o_tp.kk);
//		Vector o_normal = ObjectNormalAt(voxel, pos);
//		    	if(o_tp.ii == I && o_tp.jj == J && o_tp.kk == K){
//		    		printf("object normal = (%f, %f, %f) \n",
//		    				o_normal.x, o_normal.y, o_normal.z);
//		    	}
//			xp = old_xp + dts*up; yp = old_yp + dts*vp; zp = old_zp + dts*wp;
//			TriInterp(up1, vp1, wp1, xp/delta, yp/delta, zp/delta, u1, v1, w1);
//			old_xp = xp + dts * up1;
//			old_yp = yp + dts * vp1;
//			old_zp = zp + dts * wp1;
//			xp = 0.75f * o_xp + 0.25f * old_xp;
//			yp = 0.75f * o_yp + 0.25f * old_yp;
//			zp = 0.75f * o_zp + 0.25f * old_zp;
//			TriInterp(up1, vp1, wp1, xp/delta, yp/delta, zp/delta, u1, v1, w1);
//			old_xp = xp + dts * up1;
//			old_yp = yp + dts * vp1;
//			old_zp = zp + dts * wp1;
//			xp = one_third * o_xp + two_third * old_xp;
//			yp = one_third * o_yp + two_third * old_yp;
//			zp = one_third * o_zp + two_third * old_zp;
//			pos = Point(xp, yp, zp);
			TentativeGridPoint tp = voxel.ContainsPoint(pos);
			if(!CheckIndex(tp.ii, tp.jj, tp.kk))
				particle = particles.erase(particle);
			else{
				float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
//		if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
//			(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
//			(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
				if( phi_i_new >= 0.f ||	fabsf(phi_i_new) < colldist ){
					Vector o_normal = ObjectNormalAt(voxel, pos);
					Point projP = pos + (colldist + phi_i_new) * o_normal;
					pos = projP;
				}
				ps.SetParticle(pos);
				if(tp.ii == I && tp.jj == J && tp.kk == K)
					++N;
				++particle;
			}
			}
		}
	}

	if(s < 0)
		printf("NEGATIVE particle with maximum velocity at old pos (%f,%f,%f),"
				"cell (%d, %d, %d) \n",	poldmax.x, poldmax.y, poldmax.z,
				tpoldmax.ii, tpoldmax.jj, tpoldmax.kk);
	else
		printf("POSITIVE particle with maximum velocity at old pos (%f,%f,%f),"
			"cell (%d, %d, %d) \n",	poldmax.x, poldmax.y, poldmax.z,
			tpoldmax.ii, tpoldmax.jj, tpoldmax.kk);
	printf("new pos (%f, %f, %f) cell (%d, %d, %d) \n",
					pmax.x, pmax.y, pmax.z, tpmax.ii, tpmax.jj, tpmax.kk);
	printf(" particle velocity is (%f, %f, %f) \n",
					pvelmax.x, pvelmax.y, pvelmax.z);

	if(s < 0)
		printf("End advection of NEGATIVE particles with RK3 subcycling: \n");
	else
		printf("End advection of POSITIVE particles with RK3 subcycling: \n");

	printf("There are %d particles in cell(%d,%d,%d) \n ",
			N, I, J, K);


}

void VectorField3D::AdvectForwardEuler(const Voxel &voxel, list<Particle> &particles, int sign) const{

	float xp, yp, zp;
	float old_xp, old_yp, old_zp;
	float up, vp, wp, up1, vp1, wp1;
	int N=0;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;
	int iterations = int(round(LS_ADV_LIMIT / delta));
	float dts  = dt / iterations;
	list<Particle>::iterator particle;

	printf("start advecting  %u particles with Forward Euler : \n",particles.size());
	printf(" subcycling time step = %f and complete time step = %f\n", dts, dt);
	if(sign > 0){ // air particle
		for(particle=particles.begin(); particle != particles.end(); ){
		//		Point pos = particles[m].Postion();
		//		float radius = particles[m].Radius();
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
			float phi_i = TriInterp(voxel, pos, phi_c_obj);
			if(phi_i > 0.f){   // in solid object
//			if(phi_i > 0.f && voxel.InSolid(o_tp)){
				// we only get rid of air particles
				particle = particles.erase(particle);
			}
			else if(voxel.InSource(o_tp))
				particle = particles.erase(particle);
			else{  // particle close to a solid object
	//			float colldist = ps.CollisionDistance();
			bool inside_obj = false;
			for(int n = 0; n < iterations; ++n){
				old_xp = pos.x;
				old_yp = pos.y;
				old_zp = pos.z;
				xp     = pos.x;
				yp     = pos.y;
				zp     = pos.z;
				TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
				xp = old_xp + dts*up; yp = old_yp + dts*vp; zp = old_zp + dts*wp;
				float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
				if(phi_i_new > 0.f){
					particle = particles.erase(particle);
					inside_obj = true;
					break;
				}
				else{
//					ps.SetParticle(Point(xp,yp,zp), radius);
//					++particle;
					pos = Point(xp, yp, zp);
				}
			}
			if(!inside_obj){
				ps.SetParticle(pos, radius);
				++particle;
			}
			}
		}
	}
	else{	// water particles
		for(particle=particles.begin(); particle != particles.end(); ){
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			float colldist = ps.CollisionDistance();
			for(int n = 0; n < iterations; ++n){
				old_xp = pos.x;
				old_yp = pos.y;
				old_zp = pos.z;
				xp     = pos.x;
				yp     = pos.y;
				zp     = pos.z;
				TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
		//		float phi_i = TriInterp(voxel, pos, phi_c_obj);
				TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
		//		    	TentativeGridPoint object = voxel.ClosestObjectCell(o_tp);
		//		    	Vector o_normal = ObjectNormalAt(delta, object.ii, object.jj, object.kk);
		//		    	Vector o_normal = ObjectNormalAt(delta, o_tp.ii, o_tp.jj, o_tp.kk);
		//		    	Vector o_normal = ObjectNormalAt(voxel, o_tp.ii, o_tp.jj, o_tp.kk);
		//		Vector o_normal = ObjectNormalAt(voxel, pos);
		//		    	if(o_tp.ii == I && o_tp.jj == J && o_tp.kk == K){
		//		    		printf("object normal = (%f, %f, %f) \n",
		//		    				o_normal.x, o_normal.y, o_normal.z);
		//		    	}
				xp = old_xp + dts*up; yp = old_yp + dts*vp; zp = old_zp + dts*wp;
				pos = Point(xp, yp, zp);
				float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
		//		if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
		//			(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
		//			(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
				if( phi_i_new >= 0.f ||	fabs(phi_i_new) < colldist ){
					Vector o_normal = ObjectNormalAt(voxel, pos);
					Point projP = pos + (colldist + phi_i_new) * o_normal;
					pos = projP;
				}
			}
			ps.SetParticle(pos, radius);
			++particle;
		}
	}

	printf("end advection particles with Forward Euler : \n");

	printf("There are %d particles in cell(%d,%d,%d) \n ",
			N, I, J, K);


}

void VectorField3D::AdvectForwardEulerSubcycling(const Voxel &voxel, float dts, list<Particle> &particles, int sign) const{

	float xp, yp, zp;
	float old_xp, old_yp, old_zp;
	float up, vp, wp, up1, vp1, wp1;
	int N=0;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;
	list<Particle>::iterator particle;

	printf("start advecting  %u particles with Forward Euler Subcycling : \n",particles.size());
	printf(" subcycling time step = %f and complete time step = %f\n", dts, dt);
	if(sign > 0){ // air particle
		for(particle=particles.begin(); particle != particles.end(); ){
		//		Point pos = particles[m].Postion();
		//		float radius = particles[m].Radius();
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
			float phi_i = TriInterp(voxel, pos, phi_c_obj);
			if(phi_i > 0.f){   // in solid object
//			if(phi_i > 0.f && voxel.InSolid(o_tp)){
				// we only get rid of air particles
				particle = particles.erase(particle);
			}
			else if(voxel.InSource(o_tp))
				particle = particles.erase(particle);
			else{  // particle close to a solid object
	//			float colldist = ps.CollisionDistance();
				old_xp = pos.x;
				old_yp = pos.y;
				old_zp = pos.z;
				xp     = pos.x;
				yp     = pos.y;
				zp     = pos.z;
				TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
				xp = old_xp + dts*up; yp = old_yp + dts*vp; zp = old_zp + dts*wp;
				float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
				if(phi_i_new > 0.f){
					particle = particles.erase(particle);
				}
				else{
					pos = Point(xp, yp, zp);
					ps.SetParticle(pos, radius);
					++particle;
				}
			}
		}
	}
	else{	// water particles
		for(particle=particles.begin(); particle != particles.end(); ){
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
			float colldist = ps.CollisionDistance();
			old_xp = pos.x;
			old_yp = pos.y;
			old_zp = pos.z;
			xp     = pos.x;
			yp     = pos.y;
			zp     = pos.z;
//			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
	//		float phi_i = TriInterp(voxel, pos, phi_c_obj);
			TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
	//		    	TentativeGridPoint object = voxel.ClosestObjectCell(o_tp);
	//		    	Vector o_normal = ObjectNormalAt(delta, object.ii, object.jj, object.kk);
	//		    	Vector o_normal = ObjectNormalAt(delta, o_tp.ii, o_tp.jj, o_tp.kk);
	//		    	Vector o_normal = ObjectNormalAt(voxel, o_tp.ii, o_tp.jj, o_tp.kk);
	//		Vector o_normal = ObjectNormalAt(voxel, pos);
	//		    	if(o_tp.ii == I && o_tp.jj == J && o_tp.kk == K){
	//		    		printf("object normal = (%f, %f, %f) \n",
	//		    				o_normal.x, o_normal.y, o_normal.z);
	//		    	}
			xp = old_xp + dts*up; yp = old_yp + dts*vp; zp = old_zp + dts*wp;
			pos = Point(xp, yp, zp);
			float phi_i_new = TriInterp(voxel, pos, phi_c_obj);
	//		if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
	//			(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
	//			(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
			if( phi_i_new >= 0.f ||	fabs(phi_i_new) < colldist ){
				Vector o_normal;
				TentativeGridPoint tp = voxel.ContainsPoint(pos);
				if(voxel.InSolid(tp))
					o_normal = GeometricNormalAt(voxel, tp.ii, tp.jj, tp.kk);
				else
					o_normal = ObjectNormalAt(voxel, pos);
//				Vector o_normal = ObjectNormalAt(voxel, pos);
				Point projP = pos + (colldist + phi_i_new) * o_normal;
				pos = projP;
			}
			ps.SetParticle(pos, radius);
			++particle;
		}
	}

	printf("end advection particles with Forward Euler Subcycling: \n");

	printf("There are %d particles in cell(%d,%d,%d) \n ",
			N, I, J, K);


}

void VectorField3D::AdvectRK2(const Voxel &voxel, list<Particle> &particles, int sign) const{


	float xp, yp, zp;
	float old_xp, old_yp, old_zp;
	float up, vp, wp, up1, vp1, wp1;
	int N=0;

	float delta = voxel.VoxelDelta();
	float inv_delta = 1.f / delta;

	list<Particle>::iterator particle;

	printf("start advecting  %ld particles with RK2: \n",particles.size());

	if(sign > 0){ // air particle
		for(particle=particles.begin(); particle != particles.end(); ){
	//		Point pos = particles[m].Postion();
	//		float radius = particles[m].Radius();
			Particle &ps = *particle;
			Point pos = ps.Position();
			float radius = ps.Radius();
//			float colldist = ps.CollisionDistance();
			old_xp = xp = pos.x;
			old_yp = yp = pos.y;
			old_zp = zp = pos.z;
			float o_xp = old_xp;
			float o_yp = old_yp;
			float o_zp = old_zp;
			TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
			float phi_i = TriInterp(voxel, pos, phi_c_obj);
			if(phi_i > 0.f){   // in solid object
//			if(phi_i > 0.f && voxel.InSolid(o_tp)){
				// we only get rid of air particles
				particle = particles.erase(particle);
			}
			else if(voxel.InSource(o_tp))
				particle = particles.erase(particle);
			else{  // particle close to a solid object
				TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
//				xp = old_xp + dt*up; yp = old_yp + dt*vp; zp = old_zp + dt*wp;
				xp = old_xp + 0.5f*dt*up; yp = old_yp + 0.5f*dt*vp; zp = old_zp + 0.5f*dt*wp;
//				TentativeGridPoint new_tp = voxel.ContainsPoint(Point(xp,yp,zp));
				float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
				if(phi_i_new > 0.f)
//				if(phi_i_new > 0.f && voxel.InSolid(new_tp))
//		    	if(phi_i > 0.f) //air particle collides with solid objects
					particle = particles.erase(particle);
				else{
					TriInterp(voxel, Point(xp, yp, zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
//					xp = old_xp + 0.5f*dt*(up+up1);
//					yp = old_yp + 0.5f*dt*(vp+vp1);
//					zp = old_zp + 0.5f*dt*(wp+wp1);
					xp = old_xp + dt*up1;
					yp = old_yp + dt*vp1;
					zp = old_zp + dt*wp1;
//					new_tp = voxel.ContainsPoint(Point(xp,yp,zp));
					float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
				    if(phi_i_new > 0.f)
//					if(phi_i_new > 0.f && voxel.InSolid(new_tp))
						particle = particles.erase(particle);
					else{
						ps.SetParticle(Point(xp,yp,zp), radius);
						++particle;
					}
				}
			}

		}

	}
	else{	// water particles
		for(particle=particles.begin(); particle != particles.end(); ){
		Particle &ps = *particle;
		Point pos = ps.Position();
		float radius = ps.Radius();
		float colldist = ps.CollisionDistance();
		old_xp = xp = pos.x;
		old_yp = yp = pos.y;
		old_zp = zp = pos.z;
		float o_xp = old_xp;
		float o_yp = old_yp;
		float o_zp = old_zp;
		TentativeGridPoint o_tp = voxel.ContainsPoint(pos);
		float phi_i = TriInterp(voxel, pos, phi_c_obj);
		TriInterp(voxel, pos, up, vp, wp, inv_delta, u1, v1, w1);
//		    	TentativeGridPoint object = voxel.ClosestObjectCell(o_tp);
//		    	Vector o_normal = ObjectNormalAt(delta, object.ii, object.jj, object.kk);
//		    	Vector o_normal = ObjectNormalAt(delta, o_tp.ii, o_tp.jj, o_tp.kk);
//		    	Vector o_normal = ObjectNormalAt(voxel, o_tp.ii, o_tp.jj, o_tp.kk);
		Vector o_normal = ObjectNormalAt(voxel, pos);
//		    	if(o_tp.ii == I && o_tp.jj == J && o_tp.kk == K){
//		    		printf("object normal = (%f, %f, %f) \n",
//		    				o_normal.x, o_normal.y, o_normal.z);
//		    	}
		float v_normal;
		Vector v_para;
//		xp = old_xp + dt*up; yp = old_yp + dt*vp; zp = old_zp + dt*wp;
		xp = old_xp + 0.5f*dt*up; yp = old_yp + 0.5f*dt*vp; zp = old_zp + 0.5f*dt*wp;
		float phi_i_new = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
		if( phi_i_new >= 0.f || (phi_i_new < 0.f && fabs(phi_i_new) < colldist)){
			v_normal = Dot(Vector(up, vp, wp), o_normal);
			if(v_normal < 0.f) // particle will drift towards the object
				v_para = Vector(up, vp, wp) - v_normal * o_normal;
			else
				v_para = Vector(up, vp, wp);
			xp = old_xp + dt*v_para.x;
			yp = old_yp + dt*v_para.y;
			zp = old_zp + dt*v_para.z;
//			float phip = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
//			if(phip >= 0.f || (phip < 0.f && fabs(phip) < colldist)){
//				printf("\nparticle original pos at (%f, %f, %f) with normal (%f, %f, %f) \n",
//						old_xp, old_yp, old_zp, o_normal.x, o_normal.y, o_normal.z);
//				printf("particle now at (%f, %f, %f), phi_i_new = %f, phip = %f  colldist = %f \n",
//						xp, yp, zp, phi_i_new, phip, colldist);
//				printf(" particle in cell (%d, %d, %d) with phic = %f and phic_obj = %f\n",
//						o_tp.ii, o_tp.jj, o_tp.kk, phi_c[INDEX(o_tp.ii, o_tp.jj, o_tp.kk)], phi_c_obj[INDEX(o_tp.ii, o_tp.jj, o_tp.kk)]);
//				printf("v_para = (%f, %f, %f) up = %f, vp = %f, wp = %f, v_normal = %f\n\n",
//						v_para.x, v_para.y, v_para.z, up, vp, wp, v_normal);
//			}
		}

//		TentativeGridPoint new_tp = voxel.ContainsPoint(Point(xp,yp,zp));
		TriInterp(voxel, Point(xp, yp, zp), up1, vp1, wp1, inv_delta, u1, v1, w1);
//		float u_avg = 0.5f * (up + up1);
//		float v_avg = 0.5f * (vp + vp1);
//		float w_avg = 0.5f * (wp + wp1);
		float u_avg =  up1;
		float v_avg =  vp1;
		float w_avg =  wp1;
		xp = old_xp + dt * u_avg;
		yp = old_yp + dt * v_avg;
		zp = old_zp + dt * w_avg;
		float phi_i_newp = TriInterp(voxel, Point(xp,yp,zp), phi_c_obj);
//		if( (phi_i_new > 0.f && voxel.InSolid(new_tp)) ||  // truly in solid
//			(phi_i < 0.f && fabs(phi_i) < colldist) ||   // truly not in solid
//			(phi_i > 0.f && !voxel.InSolid(o_tp)) ){     // truly not in solid, but very close to solid object
		if( phi_i_newp >= 0.f ||  // truly in solid
			(phi_i_newp < 0.f && fabs(phi_i_newp) < colldist)){
			v_normal = Dot(Vector(u_avg,v_avg,w_avg), o_normal);
			if(v_normal < 0.f) // particle will drift towards the object
				v_para = Vector(u_avg,v_avg,w_avg) - v_normal * o_normal;
			else
				v_para = Vector(u_avg,v_avg,w_avg);
			u_avg = v_para.x;
			v_avg = v_para.y;
			w_avg = v_para.z;
			xp = old_xp + dt * u_avg;
			yp = old_yp + dt * v_avg;
			zp = old_zp + dt * w_avg;
		}
		ps.SetParticle(Point(xp,yp,zp), radius);
		++particle;
		}

	}


	printf("end advection particles with RK2: \n");

	printf("There are %d particles in cell(%d,%d,%d) \n ",
			N, I, J, K);

//	particle=particles.begin();
//	Particle &ps = *particle;
//	Point pos = ps.Position();
//	float radius = ps.Radius();
//	TentativeGridPoint tpp = voxel.ContainsPoint(pos);
//	printf("\n particle at (%f, %f, %f) with radius = %f and cell phi = %f \n\n ",
//			pos.x, pos.y, pos.z, radius, phi_c[INDEX(tpp.ii, tpp.jj, tpp.kk)]);

}


bool VectorField3D::CheckIndex(int i, int j, int k) const {
	if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
//		printf("Error in (VectorField3D::CheckIndex())! index out of range i=%d, j=%d, k=%d \n",
//				i,j,k);
		//exit(1);
		return false;
	}
	else
		return true;
}

void VectorField3D::GetVelocityAtPos(const Voxel &voxel, const Point &p,
		float delta, float &uu, float &vv, float &ww) const{
	float inv_delta = 1.f / delta;
	TriInterp(voxel, p, uu, vv, ww, inv_delta, u, v, w);
}

#if defined(COUPLED_SPH) || defined(COUPLED_FLIP)
void VectorField3D::MomentumConservation(const Voxel &voxel,
//		list<WaterParticle> &absorbed
		WaterParticleList &absorbed
		){

	float delta = voxel.VoxelDelta();
	float cellVolume = delta * delta * delta;

//	list<WaterParticle>::iterator iter_waterparticle;
//	typedef list<WaterParticle> PLIST;

	WaterParticleList::iterator iter_waterparticle;
	typedef WaterParticleList PLIST;

#ifdef WIN32
	hash_map<u_int, PLIST> forDelete;
	hash_map<u_int, PLIST>::iterator iter_cell;
#else
	unordered_map<u_int, PLIST> forDelete;
	unordered_map<u_int, PLIST>::iterator iter_cell;
#endif

	for(iter_waterparticle  = absorbed.begin();
		iter_waterparticle != absorbed.end();
		++iter_waterparticle){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPointInVelCell(pos, 1);
		u_int index = INDEX(tp.ii,tp.jj,tp.kk);
		iter_cell = forDelete.find(index);
		if(iter_cell == forDelete.end()){
			PLIST plist;
			plist.push_back(ps);
			forDelete.insert(make_pair(index, plist));
		}
		else{
			PLIST &plist = iter_cell->second;
			plist.push_back(ps);
		}
	}
	for(iter_cell  = forDelete.begin();
	    iter_cell != forDelete.end();
	    ++iter_cell){
		PLIST &plist = iter_cell->second;
		u_int index  = iter_cell->first;
		float gridu  = u[index];
		float utot   = gridu * cellVolume;
		float voltot = cellVolume;
		for(iter_waterparticle  = plist.begin();
			iter_waterparticle != plist.end();
			++iter_waterparticle){
			WaterParticle &ps = *iter_waterparticle;
			float uvel = ps.UVelocity();
			float rad  = ps.Radius();
			float pvol = 4*M_PI/3*rad*rad*rad;
			utot   += uvel * pvol;
			voltot +=        pvol;
		}
		plist.clear();
		u[index] = utot / voltot;
	}

	forDelete.clear();

	for(iter_waterparticle  = absorbed.begin();
		iter_waterparticle != absorbed.end();
		++iter_waterparticle){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPointInVelCell(pos, 2);
		u_int index = INDEX(tp.ii,tp.jj,tp.kk);
		iter_cell = forDelete.find(index);
		if(iter_cell == forDelete.end()){
			PLIST plist;
			plist.push_back(ps);
			forDelete.insert(make_pair(index, plist));
		}
		else{
			PLIST &plist = iter_cell->second;
			plist.push_back(ps);
		}
	}
	for(iter_cell  = forDelete.begin();
	    iter_cell != forDelete.end();
	    ++iter_cell){
		PLIST &plist = iter_cell->second;
		u_int index  = iter_cell->first;
		float gridv  = v[index];
		float vtot   = gridv * cellVolume;
		float voltot = cellVolume;
		for(iter_waterparticle  = plist.begin();
			iter_waterparticle != plist.end();
			++iter_waterparticle){
			WaterParticle &ps = *iter_waterparticle;
			float vvel = ps.VVelocity();
			float rad  = ps.Radius();
			float pvol = 4*M_PI/3*rad*rad*rad;
			vtot   += vvel * pvol;
			voltot +=        pvol;
		}
		plist.clear();
		v[index] = vtot / voltot;
	}

	forDelete.clear();

	for(iter_waterparticle  = absorbed.begin();
		iter_waterparticle != absorbed.end();
		++iter_waterparticle){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPointInVelCell(pos, 3);
		u_int index = INDEX(tp.ii,tp.jj,tp.kk);
		iter_cell = forDelete.find(index);
		if(iter_cell == forDelete.end()){
			PLIST plist;
			plist.push_back(ps);
			forDelete.insert(make_pair(index, plist));
		}
		else{
			PLIST &plist = iter_cell->second;
			plist.push_back(ps);
		}
	}
	for(iter_cell  = forDelete.begin();
	    iter_cell != forDelete.end();
	    ++iter_cell){
		PLIST &plist = iter_cell->second;
		u_int index  = iter_cell->first;
		float gridw  = w[index];
		float wtot   = gridw * cellVolume;
		float voltot = cellVolume;
		for(iter_waterparticle  = plist.begin();
			iter_waterparticle != plist.end();
			++iter_waterparticle){
			WaterParticle &ps = *iter_waterparticle;
			float wvel = ps.WVelocity();
			float rad  = ps.Radius();
			float pvol = 4*M_PI/3*rad*rad*rad;
			wtot   += wvel * pvol;
			voltot +=        pvol;
		}
		plist.clear();
		w[index] = wtot / voltot;
	}

	forDelete.clear();
}
#endif
#ifdef COUPLED_FLIP

void VectorField3D::AccumulateValues(const Point &pos, float invh,
				float vel, int index,
				float *mWeights, float *vc){
	int i0, j0, k0, i1, j1, k1;
	int x0, y0, z0;

	float s0, t0, s1, t1;
	float r0, r1;
	float weight;

	float xp = pos.x * invh;
	float yp = pos.y * invh;
	float zp = pos.z * invh;

	if (xp<0.5f) xp=0.5f; if (xp>DimX-1.5f) xp=DimX-1.5f;
	if (yp<0.5f) yp=0.5f; if (yp>DimY-1.5f) yp=DimY-1.5f;
	if (zp<0.5f) zp=0.5f; if (zp>DimZ-1.5f) zp=DimZ-1.5f;

	if(index == 1){

		i0=floor(xp)-1; if(i0 < 0) i0 = 0;
		i1 = i0 + 1;
		s1 = xp - 1 - i0;
		s0 = 1 - s1;
		j0 = floor(yp); y0 = j0;
		if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		j1=j0+1;
		if(yp-y0 < 0.5f) t1 = 0.5f+yp-y0;
		else t1 = yp-y0-0.5f;
		t0 = 1-t1;
		k0 = floor(zp); z0 = k0;
		if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
		k1=k0+1;
		if(zp-z0 < 0.5f) r1 = 0.5f+zp-z0;
		else r1 = zp-z0-0.5f;
		r0 = 1-r1;
	}
	if(index == 2){

		j0=floor(yp)-1; if(j0 < 0) j0 = 0;
		j1 = j0 + 1;
		t1 = yp - 1 - j0;
		t0 = 1 - t1;
		i0 = floor(xp); x0 = i0;
		if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		i1=i0+1;
		if(xp-x0 < 0.5f) s1 = 0.5f+xp-x0;
		else s1 = xp-x0-0.5f;
		s0 = 1-s1;
		k0 = floor(zp); z0 = k0;
		if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
		k1=k0+1;
		if(zp-z0 < 0.5f) r1 = 0.5f+zp-z0;
		else r1 = zp-z0-0.5f;
		r0 = 1-r1;

	}
	if(index == 3){

		k0=floor(zp)-1; if(k0 < 0) k0 = 0;
		k1 = k0+1;
		r1 = zp - 1 - k0;
		r0 = 1 - r1;
		i0 = floor(xp); x0 = i0;
		if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
		i1=i0+1;
		if(xp-x0 < 0.5f) s1 = 0.5f+xp-x0;
		else s1 = xp-x0-0.5f;
		s0 = 1-s1;
		j0 = floor(yp); y0 = j0;
		if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
		j1=j0+1;
		if(yp-y0 < 0.5f) t1 = 0.5f+yp-y0;
		else t1 = yp-y0-0.5f;
		t0 = 1-t1;

	}

	weight = r0*s0*t0;
	vc[INDEX(i0,j0,k0)] += weight * vel;
	mWeights[INDEX(i0,j0,k0)] += weight;
	weight = r0*s0*t1;
	vc[INDEX(i0,j1,k0)] += weight * vel;
	mWeights[INDEX(i0,j1,k0)] += weight;
	weight = r0*s1*t0;
	vc[INDEX(i1,j0,k0)] += weight * vel;
	mWeights[INDEX(i1,j0,k0)] += weight;
	weight = r0*s1*t1;
	vc[INDEX(i1,j1,k0)] += weight * vel;
	mWeights[INDEX(i1,j1,k0)] += weight;
	weight = r1*s0*t0;
	vc[INDEX(i0,j0,k1)] += weight * vel;
	mWeights[INDEX(i0,j0,k1)] += weight;
	weight = r1*s0*t1;
	vc[INDEX(i0,j1,k1)] += weight * vel;
	mWeights[INDEX(i0,j1,k1)] += weight;
	weight = r1*s1*t0;
	vc[INDEX(i1,j0,k1)] += weight * vel;
	mWeights[INDEX(i1,j0,k1)] += weight;
	weight = r1*s1*t1;
	vc[INDEX(i1,j1,k1)] += weight * vel;
	mWeights[INDEX(i1,j1,k1)] += weight;

}

#define RADIUS_MAX 0.5f*delta
#define RADIUS_MIN 0.1f*delta

void VectorField3D::
	 ConvertParticlesToLevelSet(const Voxel &voxel,
		WaterParticleList &wpl, ParticleList &pl,
		int nPerCell){

	float delta = voxel.VoxelDelta();
	float invh = 1.f / delta;

	WaterParticleList::iterator iter_waterparticle;
	WaterParticleList tempList;
	int *marker = new int[DimX*DimY*DimZ];
	memset(marker, 0, DimX * DimY * DimZ * sizeof(int));
	SetEqual(phi_c, temp);

	SetZero(u0);
	SetZero(v0);
	SetZero(w0);

	for(iter_waterparticle  = wpl.begin();
		iter_waterparticle != wpl.end();
		){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		++marker[INDEX(tp.ii, tp.jj, tp.kk)];
		++iter_waterparticle;
	}

	for(iter_waterparticle  = wpl.begin();
		iter_waterparticle != wpl.end();
		){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		float rad = ps.Radius();
		TentativeGridPoint tp = voxel.ContainsPoint(pos);
		int i = tp.ii, j = tp.jj, k = tp.kk;
		u_int index = INDEX(i,j,k);
		if(marker[index] >= nPerCell){
			Point center = voxel.VoxelCenterPosition(i, j, k);
			float phi_center = ps.EvaluatePhiP(center);
			temp[index] = min(phi_center, temp[index]);
			tempList.push_back(ps);
			iter_waterparticle = wpl.erase(iter_waterparticle);
		}
		else
			++iter_waterparticle;
	}

	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		phi_c[pos] = fabs(phi_c[pos]) <= fabs(temp[pos]) ? phi_c[pos] : temp[pos];
	END_FOR


	float *mWeights = temp;

	SetZero(mWeights);
	for(iter_waterparticle  = tempList.begin();
	    iter_waterparticle != tempList.end();){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		Vector vel = ps.Velocity();
		AccumulateValues(pos, invh, vel.x, 1, mWeights, u0);
		++iter_waterparticle;
	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mWeights[pos] != 0.f)
			u0[pos] /= mWeights[pos];
	END_FOR
	FOR_EACH_CELL
		u_int index = INDEX(i,j,k);
		if(i < DimX-1)
			if(mWeights[index] != 0.f && voxel.InAir(i,j,k) && voxel.InAir(i+1,j,k))
				u[index] = u0[index];
	END_FOR

	SetZero(mWeights);
	for(iter_waterparticle  = tempList.begin();
	    iter_waterparticle != tempList.end();){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		Vector vel = ps.Velocity();
		AccumulateValues(pos, invh, vel.y, 2, mWeights, v0);
		++iter_waterparticle;
	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mWeights[pos] != 0.f)
			v0[pos] /= mWeights[pos];
	END_FOR
	FOR_EACH_CELL
		u_int index = INDEX(i,j,k);
		if(j < DimY-1)
			if(mWeights[index] != 0.f && voxel.InAir(i,j,k) && voxel.InAir(i,j+1,k))
				v[index] = v0[index];
	END_FOR

	SetZero(mWeights);
	for(iter_waterparticle  = tempList.begin();
	    iter_waterparticle != tempList.end();){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		Vector vel = ps.Velocity();
		AccumulateValues(pos, invh, vel.z, 3, mWeights, w0);
		++iter_waterparticle;
	}
	FOR_EACH_CELL
		u_int pos = INDEX(i,j,k);
		if(mWeights[pos] != 0.f)
			w0[pos] /= mWeights[pos];
	END_FOR
	FOR_EACH_CELL
		u_int index = INDEX(i,j,k);
		if(k < DimZ-1)
			if(mWeights[index] != 0.f && voxel.InAir(i,j,k) && voxel.InAir(i,j,k+1))
				w[index] = w0[index];
	END_FOR




	for(iter_waterparticle  = tempList.begin();
		iter_waterparticle != tempList.end();
		){
		WaterParticle &ps = *iter_waterparticle;
		Point pos = ps.Position();
		float rad = ps.Radius();
		float phi_i = TriInterp(voxel, pos, phi_c);
		if(phi_i < rad){ // accepted as a negative marker particle
			if(phi_i <= 0.f){
				if( -phi_i > RADIUS_MAX ) rad = RADIUS_MAX;
				else if( -phi_i < RADIUS_MIN ) rad = RADIUS_MIN;
				else rad = -phi_i;
			}
			else
				rad = RADIUS_MIN;
			ps.SetParticle(pos, rad);
			MTRand mt;
			double rn = mt();
			Particle p(pos, rad, -1, (0.1f+(float)rn*0.9f)*delta);
			pl.push_back(p);
		}
		++iter_waterparticle;
	}

	delete [] marker;
}
#endif


void VectorField3D::AssignEightVelocities(int i0, int j0, int k0,
		int i1, int j1, int k1,
		const float *vel, float f[8]) const{

	f[0] =  vel[INDEX(i0,j0,k0)];
	f[1] =  vel[INDEX(i0,j1,k0)];
	f[2] =  vel[INDEX(i1,j0,k0)];
	f[3] =  vel[INDEX(i1,j1,k0)];
	f[4] =  vel[INDEX(i0,j0,k1)];
	f[5] =  vel[INDEX(i0,j1,k1)];
	f[6] =  vel[INDEX(i1,j0,k1)];
	f[7] =  vel[INDEX(i1,j1,k1)];

}

void VectorField3D::
	HandleEightRaysForVelocityInterpolation(const Voxel &voxel,
			int i0, int j0, int k0,
			int i1, int j1, int k1,
			int vel_index, const Point &pi,
		    float f[8]) const{

	double p, q;

	Point pc = voxel.VelPosition(vel_index, i0, j0, k0);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[0] = 0.f;
	pc = voxel.VelPosition(vel_index, i0, j1, k0);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[1] = 0.f;
	pc = voxel.VelPosition(vel_index, i1, j0, k0);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[2] = 0.f;
	pc = voxel.VelPosition(vel_index, i1, j1, k0);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[3] = 0.f;
	pc = voxel.VelPosition(vel_index, i0, j0, k1);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[4] = 0.f;
	pc = voxel.VelPosition(vel_index, i0, j1, k1);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[5] = 0.f;
	pc = voxel.VelPosition(vel_index, i1, j0, k1);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[6] = 0.f;
	pc = voxel.VelPosition(vel_index, i1, j1, k1);
	if(mPhysWorld->SegmentIntersectsMesh(pi, pc, &p, &q))
		f[7] = 0.f;
}

void VectorField3D::TriInterp(const Voxel &voxel, const Point &pi,
		float &ui, float &vi, float &wi,
		float inv_delta, float *uvel, float *vvel, float *wvel) const{

	int i0, j0, k0, i1, j1, k1;
	float s0, t0, s1, t1;
	int x0, y0, z0;

	float r0, r1;
	float xp0, yp0, zp0;

#ifdef USE_MESH
	double p,q;
	float f[8];
	int vel_index;
#endif

	float xp = pi.x * inv_delta;
	float yp = pi.y * inv_delta;
	float zp = pi.z * inv_delta;


	xp = xp < 0.5f      ? 0.5f      : xp;
	xp = xp > DimX-1.5f ? DimX-1.5f : xp;
	yp = yp < 0.5f      ? 0.5f      : yp;
	yp = yp > DimY-1.5f ? DimY-1.5f : yp;
	zp = zp < 0.5f      ? 0.5f      : zp;
	zp = zp > DimZ-1.5f ? DimZ-1.5f : zp;

    i0=floor(xp)-1; if(i0 < 0) i0 = 0;
	i1 = i0 + 1;
    s1 = xp - 1 - i0;
    s0 = 1 - s1;
    j0 = floor(yp); y0 = j0;
    if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
    j1=j0+1;
    yp0 = yp-y0;
    t1 = yp0 < 0.5f ? 0.5f+yp0 : yp0-0.5f;
    t0 = 1-t1;
    k0 = floor(zp); z0 = k0;
    if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
    k1=k0+1;
    zp0 = zp-z0;
    r1 = zp0 < 0.5f ? 0.5f+zp0 : zp0-0.5f;
    r0 = 1-r1;

//	CheckIndex(i0,j0,k0);
//	CheckIndex(i0,j1,k0);
//	CheckIndex(i1,j0,k0);
//	CheckIndex(i1,j1,k0);
//	CheckIndex(i0,j0,k1);
//	CheckIndex(i0,j1,k1);
//	CheckIndex(i1,j0,k1);
//	CheckIndex(i1,j1,k1);
#ifdef USE_MESH
    AssignEightVelocities(i0, j0, k0, i1, j1, k1, uvel, f);
	vel_index = 1;
	HandleEightRaysForVelocityInterpolation(voxel,
				i0, j0, k0, i1, j1, k1,
				vel_index, pi, f);

	ui =  r0*(s0*(t0*f[0]+t1*f[1])+
			  s1*(t0*f[2]+t1*f[3])) +
		  r1*(s0*(t0*f[4]+t1*f[5])+
			  s1*(t0*f[6]+t1*f[7]));
#else
	ui = r0*(s0*(t0*uvel[INDEX(i0,j0,k0)]+t1*uvel[INDEX(i0,j1,k0)]) +
		     s1*(t0*uvel[INDEX(i1,j0,k0)]+t1*uvel[INDEX(i1,j1,k0)]))+
 		 r1*(s0*(t0*uvel[INDEX(i0,j0,k1)]+t1*uvel[INDEX(i0,j1,k1)]) +
			 s1*(t0*uvel[INDEX(i1,j0,k1)]+t1*uvel[INDEX(i1,j1,k1)]));
#endif

	j0=floor(yp)-1; if(j0 < 0) j0 = 0;
	j1 = j0 + 1;
    t1 = yp - 1 - j0;
    t0 = 1 - t1;
    i0 = floor(xp); x0 = i0;
    if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
    i1=i0+1;
    xp0 = xp - x0;
    s1 = xp0 < 0.5f ? 0.5f+xp0 : xp0-0.5f;
    s0 = 1-s1;
    k0 = floor(zp); z0 = k0;
    if(zp-k0 < 0.5f) k0--; if(k0 < 0) k0=0;
    k1=k0+1;
    zp0 = zp-z0;
    r1 = zp0 < 0.5f ? 0.5f+zp0 : zp0-0.5f;
    r0 = 1-r1;

//	CheckIndex(i0,j0,k0);
//	CheckIndex(i0,j1,k0);
//	CheckIndex(i1,j0,k0);
//	CheckIndex(i1,j1,k0);
//	CheckIndex(i0,j0,k1);
//	CheckIndex(i0,j1,k1);
//	CheckIndex(i1,j0,k1);
//	CheckIndex(i1,j1,k1);

#ifdef USE_MESH
    AssignEightVelocities(i0, j0, k0, i1, j1, k1, vvel, f);
	vel_index = 2;
	HandleEightRaysForVelocityInterpolation(voxel,
					i0, j0, k0, i1, j1, k1,
					vel_index, pi, f);

	vi =  r0*(s0*(t0*f[0]+t1*f[1])+
		      s1*(t0*f[2]+t1*f[3])) +
		  r1*(s0*(t0*f[4]+t1*f[5])+
			  s1*(t0*f[6]+t1*f[7]));
#else
    vi = r0*(s0*(t0*vvel[INDEX(i0,j0,k0)]+t1*vvel[INDEX(i0,j1,k0)])+
		     s1*(t0*vvel[INDEX(i1,j0,k0)]+t1*vvel[INDEX(i1,j1,k0)])) +
		 r1*(s0*(t0*vvel[INDEX(i0,j0,k1)]+t1*vvel[INDEX(i0,j1,k1)])+
		     s1*(t0*vvel[INDEX(i1,j0,k1)]+t1*vvel[INDEX(i1,j1,k1)]));
#endif

    k0=floor(zp)-1; if(k0 < 0) k0 = 0;
    k1 = k0+1;
    r1 = zp - 1 - k0;
    r0 = 1 - r1;
    i0 = floor(xp); x0 = i0;
    if(xp-i0 < 0.5f) i0--; if(i0 < 0) i0=0;
    i1=i0+1;
    xp0 = xp - x0;
    s1 = xp0 < 0.5f ? 0.5f+xp0 : xp0-0.5f;
    s0 = 1-s1;
    j0 = floor(yp); y0 = j0;
    if(yp-j0 < 0.5f) j0--; if(j0 < 0) j0=0;
    j1=j0+1;
    yp0 = yp-y0;
    t1 = yp0 < 0.5f ? 0.5f+yp0 : yp0-0.5f;
    t0 = 1-t1;

//	CheckIndex(i0,j0,k0);
//	CheckIndex(i0,j1,k0);
//	CheckIndex(i1,j0,k0);
//	CheckIndex(i1,j1,k0);
//	CheckIndex(i0,j0,k1);
//	CheckIndex(i0,j1,k1);
//	CheckIndex(i1,j0,k1);
//	CheckIndex(i1,j1,k1);

#ifdef USE_MESH
    AssignEightVelocities(i0, j0, k0, i1, j1, k1, wvel, f);
	vel_index = 3;
	HandleEightRaysForVelocityInterpolation(voxel,
				i0, j0, k0, i1, j1, k1,
				vel_index, pi, f);

	wi =  r0*(s0*(t0*f[0]+t1*f[1])+
			  s1*(t0*f[2]+t1*f[3])) +
		  r1*(s0*(t0*f[4]+t1*f[5])+
			  s1*(t0*f[6]+t1*f[7]));

#else
    wi = r0*(s0*(t0*wvel[INDEX(i0,j0,k0)]+t1*wvel[INDEX(i0,j1,k0)])+
	     	 s1*(t0*wvel[INDEX(i1,j0,k0)]+t1*wvel[INDEX(i1,j1,k0)])) +
	     r1*(s0*(t0*wvel[INDEX(i0,j0,k1)]+t1*wvel[INDEX(i0,j1,k1)])+
			 s1*(t0*wvel[INDEX(i1,j0,k1)]+t1*wvel[INDEX(i1,j1,k1)]));
#endif
//    if(i0 == I && j0 == J && k0 == K){
//		printf("particle at (%f, %f, %f) with w = %f \n ",
//			   	xp, yp, zp, wi);
//		printf("t0 = %f, t1 = %f, s0 = %f, s1 = %f, r0 = %f, r1 = %f \n",
//				t0, t1, s0, s1, r0, r1);
//		printf("w(i0,j0,k0) = %f,w(i0,j1,k0) = %f \n",
//				wvel[INDEX(i0,j0,k0)], wvel[INDEX(i0,j1,k0)]);
//		printf("w(i1,j0,k0) = %f,w(i1,j1,k0) = %f \n",
//				wvel[INDEX(i1,j0,k0)], wvel[INDEX(i1,j1,k0)]);
//		printf("w(i0,j0,k1) = %f,w(i0,j1,k1) = %f \n",
//				wvel[INDEX(i0,j0,k1)], wvel[INDEX(i0,j1,k1)]);
//		printf("w(i1,j0,k1) = %f,w(i1,j1,k1) = %f \n",
//				wvel[INDEX(i1,j0,k1)], wvel[INDEX(i1,j1,k1)]);
//
//	}
}

float VectorField3D::
	UpdateTentativeValue(const KnownPoint &p, u_int index,
						const char *mask, char MASK,
						 const Voxel &voxel,
						 float *phi){
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


		if(mask[indadj] & MASK){
//			printf("found in band %d \n", indadj);
			t[i] = phi[indadj];
			if(t[i] < 0.f)
				printf("i = %d, ii = %d, jj = %d, kk = %d, t[%d] = %f\n",
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

	if(p.ii == I && p.jj == J && p.kk == K)
		printf("p.ii = %d, p.jj = %d, p.kk = %d, "
				"t0=%f, t1=%f, t2=%f, t3=%f, t4=%f, t5=%f \n",
			p.ii, p.jj, p.kk, t[0], t[1], t[2], t[3], t[4], t[5]);

	if(axis[0] && axis[1] && axis[2]){ // 3 non-zero terms
		if(t[0] <= t[1]){     // quardrants: 0,1,2,3
			if(t[2] <= t[3]) { // quardrants: 0,2
				if(t[4] <= t[5]){ // quardrants: 0
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[2], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//			 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 2
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[2], t[5], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//								 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
			}
			else{               // quardrants: 1,3
				if(t[4] <= t[5]){ // quardrants: 1
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[3], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 3
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[3], t[5], voxel);
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
					float solution = ProcessThreeTerms(0, 1, 2, t[1], t[2], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 6
					float solution =  ProcessThreeTerms(0, 1, 2, t[1], t[2], t[5], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
			}
			else{               // quardrants: 5,7
				if(t[4] <= t[5]){ // quardrants: 5
					float solution = ProcessThreeTerms(0, 1, 2, t[1], t[3], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 7
					float solution =  ProcessThreeTerms(0, 1, 2, t[1], t[3], t[5], voxel);
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
				return ProcessTwoTerms(0, 1, t[0], t[2], voxel);
			else
				return ProcessTwoTerms(0, 1, t[0], t[3], voxel);
		}
		else{
			if(t[2] <= t[3])
				return ProcessTwoTerms(0, 1, t[1], t[2], voxel);
			else
				return ProcessTwoTerms(0, 1, t[1], t[3], voxel);
		}
	}
	else if(axis[0] && !axis[1] && axis[2]){ // 2 non-zero terms
		if( t[0] <= t[1]){
			if(t[4] <= t[5])
				return ProcessTwoTerms(0, 2, t[0], t[4], voxel);
			else
				return ProcessTwoTerms(0, 2, t[0], t[5], voxel);
		}
		else{
			if(t[4] <= t[5])
				return ProcessTwoTerms(0, 2, t[1], t[4], voxel);
			else
				return ProcessTwoTerms(0, 2, t[1], t[5], voxel);
		}
	}
	else if(!axis[0] && axis[1] && axis[2]){ // 2 non-zero terms
		if( t[2] <= t[3]){
			if(t[4] <= t[5])
				return ProcessTwoTerms(1, 2, t[2], t[4], voxel);
			else
				return ProcessTwoTerms(1, 2, t[2], t[5], voxel);
		}
		else{
			if(t[4] <= t[5])
				return ProcessTwoTerms(1, 2, t[3], t[4], voxel);
			else
				return ProcessTwoTerms(1, 2, t[3], t[5], voxel);
		}
	}
	else if(axis[0] && !axis[1] && !axis[2]){ // 1 non-zero term
		if(t[0] <= t[1])
			return ProcessOneTerm(0, t[0], voxel);
		else
			return ProcessOneTerm(0, t[1], voxel);
	}
	else if(!axis[0] && axis[1] && !axis[2]){ // 1 non-zero term
		if(t[2] <= t[3])
			return ProcessOneTerm(1, t[2], voxel);
		else
			return ProcessOneTerm(1, t[3], voxel);
	}
	else if(!axis[0] && !axis[1] && axis[2]){ // 1 non-zero term
		if(t[4] <= t[5])
			return ProcessOneTerm(2, t[4], voxel);
		else
			return ProcessOneTerm(2, t[5], voxel);
	}
	else
		return INFINITY;

}

float VectorField3D::
	UpdateTentativeValue(const KnownPoint &p, u_int index,
						 map<u_int, KnownPoint> &band,
						 const Voxel &voxel,
						 float *phi){
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
		 	indadj = INDEX(p.ii-1, p.jj, p.kk);
		}
		else if(i==1) {	 // 0,1 ==> x-axis
			if(p.ii==DimX-1){
				t[1] = INFINITY;
				continue;
			}
		 	indadj = INDEX(p.ii+1, p.jj, p.kk);
		}
		else if(i==4){ // 4,5 ==> z-axis
			if(p.kk==0){
				t[4] = INFINITY;
				continue;
			}
		    indadj = INDEX(p.ii, p.jj, p.kk-1);
//		   	printf("i=%d, p.kk = %d, p.kk+-1 = %d \n", i, p.kk, p.kk+((i%2==0) ? (-1) : 1) );
		}
		else if(i==5){ // 4,5 ==> z-axis
			if(p.kk==DimZ-1){
				t[5] = INFINITY;
				continue;
			}
		    indadj = INDEX(p.ii, p.jj, p.kk+1);
//		   	printf("i=%d, p.kk = %d, p.kk+-1 = %d \n", i, p.kk, p.kk+((i%2==0) ? (-1) : 1) );
		}
		else if(i==2){         // 2,3 ==> y-axis
			if(p.jj==0){
				t[2] = INFINITY;
				continue;
			}
			indadj = INDEX(p.ii, p.jj-1, p.kk);
		}
		else{         // 2,3 ==> y-axis
			if(p.jj==DimY-1){
				t[3] = INFINITY;
				continue;
			}
		    indadj = INDEX(p.ii, p.jj+1, p.kk);
		}
//		printf("trying to find in band i=%d ii=%d, jj=%d, kk=%d, %d \n",
//					i,	p.ii, p.jj, p.kk, indadj);
		map<u_int, KnownPoint>::iterator found = band.find(indadj);
		if(found != band.end()){
//			printf("found in band %d \n", indadj);
			t[i] = phi[indadj];
			if(t[i] < 0.f)
				printf("i = %d, ii=%d, jj=%d, kk=%d, t[%d]=%f, bandsize=%d\n",
						i, p.ii, p.jj, p.kk, i, t[i], band.size());
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

	if(p.ii == I && p.jj == J && p.kk == K)
		printf("p.ii = %d, p.jj = %d, p.kk = %d, "
				"t0=%f, t1=%f, t2=%f, t3=%f, t4=%f, t5=%f \n",
			p.ii, p.jj, p.kk, t[0], t[1], t[2], t[3], t[4], t[5]);

	if(axis[0] && axis[1] && axis[2]){ // 3 non-zero terms
		if(t[0] <= t[1]){     // quardrants: 0,1,2,3
			if(t[2] <= t[3]) { // quardrants: 0,2
				if(t[4] <= t[5]){ // quardrants: 0
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[2], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//			 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 2
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[2], t[5], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//								 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
			}
			else{               // quardrants: 1,3
				if(t[4] <= t[5]){ // quardrants: 1
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[3], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 3
					float solution = ProcessThreeTerms(0, 1, 2, t[0], t[3], t[5], voxel);
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
					float solution = ProcessThreeTerms(0, 1, 2, t[1], t[2], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 6
					float solution =  ProcessThreeTerms(0, 1, 2, t[1], t[2], t[5], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
			}
			else{               // quardrants: 5,7
				if(t[4] <= t[5]){ // quardrants: 5
					float solution = ProcessThreeTerms(0, 1, 2, t[1], t[3], t[4], voxel);
//					if(solution < 0.f)
//						printf("Wrong! phi < 0 occurs in Three Terms t0 = %f, t1 = %f, t2 = %f, t3 = %f, t4 = %f, t5 = %f \n",
//		 					t[0], t[1], t[2], t[3], t[4], t[5]);
					return solution;
				}
				else{ 			  // quardrants: 7
					float solution =  ProcessThreeTerms(0, 1, 2, t[1], t[3], t[5], voxel);
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
				return ProcessTwoTerms(0, 1, t[0], t[2], voxel);
			else
				return ProcessTwoTerms(0, 1, t[0], t[3], voxel);
		}
		else{
			if(t[2] <= t[3])
				return ProcessTwoTerms(0, 1, t[1], t[2], voxel);
			else
				return ProcessTwoTerms(0, 1, t[1], t[3], voxel);
		}
	}
	else if(axis[0] && !axis[1] && axis[2]){ // 2 non-zero terms
		if( t[0] <= t[1]){
			if(t[4] <= t[5])
				return ProcessTwoTerms(0, 2, t[0], t[4], voxel);
			else
				return ProcessTwoTerms(0, 2, t[0], t[5], voxel);
		}
		else{
			if(t[4] <= t[5])
				return ProcessTwoTerms(0, 2, t[1], t[4], voxel);
			else
				return ProcessTwoTerms(0, 2, t[1], t[5], voxel);
		}
	}
	else if(!axis[0] && axis[1] && axis[2]){ // 2 non-zero terms
		if( t[2] <= t[3]){
			if(t[4] <= t[5])
				return ProcessTwoTerms(1, 2, t[2], t[4], voxel);
			else
				return ProcessTwoTerms(1, 2, t[2], t[5], voxel);
		}
		else{
			if(t[4] <= t[5])
				return ProcessTwoTerms(1, 2, t[3], t[4], voxel);
			else
				return ProcessTwoTerms(1, 2, t[3], t[5], voxel);
		}
	}
	else if(axis[0] && !axis[1] && !axis[2]){ // 1 non-zero term
		if(t[0] <= t[1])
			return ProcessOneTerm(0, t[0], voxel);
		else
			return ProcessOneTerm(0, t[1], voxel);
	}
	else if(!axis[0] && axis[1] && !axis[2]){ // 1 non-zero term
		if(t[2] <= t[3])
			return ProcessOneTerm(1, t[2], voxel);
		else
			return ProcessOneTerm(1, t[3], voxel);
	}
	else if(!axis[0] && !axis[1] && axis[2]){ // 1 non-zero term
		if(t[4] <= t[5])
			return ProcessOneTerm(2, t[4], voxel);
		else
			return ProcessOneTerm(2, t[5], voxel);
	}
	else
		return INFINITY;

}

float  VectorField3D::ProcessOneTerm(int axis, float f, const Voxel &voxel){
	float delta = voxel.VoxelDelta();
	if(f + delta < 0.f)
		printf("Wrong! phi < 0 occurs in One Term, phi = %f \n", f + delta);
	return f + delta;
}

float  VectorField3D::ProcessTwoTerms(int axis1, int axis2,
									  float phi1, float phi2,
									  const Voxel &voxel){
	float delta = voxel.VoxelDelta();
	float delta_inv = 1.f / delta;
	float delta_squ = delta * delta;

	float phimax = max(phi1, phi2);
	float phi = phi1 <= phi2 ? phi1 : phi2;
	int axis = phi1 <= phi2 ? axis1 : axis2;
	float P1 = (phimax - phi1)*delta_inv;
	float P2 = (phimax - phi2)*delta_inv;
	float P  = P1*P1 + P2*P2;
	if(P > 1)
	 	return ProcessOneTerm(axis, phi, voxel);
	else{
		float A, B, C, t0, t1;
	 	A = 2.f;
	 	B = -2*(phi1 + phi2);
	 	C = phi1*phi1 + phi2*phi2 - delta_squ;
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

float  VectorField3D::ProcessThreeTerms(int axis1, int axis2, int axis3,
									    float phi1, float phi2, float phi3,
									    const Voxel &voxel){
	float delta = voxel.VoxelDelta();
	float delta_inv = 1.f / delta;
	float delta_squ = delta * delta;

	float phimax = max(phi1, phi2);
	phimax = max(phimax, phi3);
	float P1 = (phimax - phi1)*delta_inv;
	float P2 = (phimax - phi2)*delta_inv;
	float P3 = (phimax - phi3)*delta_inv;
	float P  = P1*P1 + P2*P2 + P3*P3;
	if(P > 1){
		if(fabsf(phi1-phimax) < E_EPSIL)
	 		return ProcessTwoTerms(axis2, axis3, phi2, phi3, voxel);
	 	else if(fabsf(phi2-phimax) < E_EPSIL)
	 		return ProcessTwoTerms(axis1, axis3, phi1, phi3, voxel);
	 	else
	 		return ProcessTwoTerms(axis1, axis2, phi1, phi2, voxel);
	}
	else{
		float A, B, C, t0, t1;
	 	A = 3.f;
	 	B = -2 * (phi1 + phi2 + phi3);
	 	C = phi1*phi1 + phi2*phi2 + phi3*phi3 - delta_squ;
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
