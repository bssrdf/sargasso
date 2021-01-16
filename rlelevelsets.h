#ifndef RLELEVELSETS_H_
#define RLELEVELSETS_H_

#include <algorithm>
#include <stdio.h>

//#include "CIsoSurface.h"

#ifdef CUDA
#include "reinit_cuda.h"
#endif

#define INDEX(i,j,k) (k)*DimX*DimY+(j)*DimX+(i)

#define INDEX_RUN(j,k) (k)*DimY+(j)

#define FOR_EACH_CELL for(int k=0;k<DimZ;k++) { \
							for(int j=0;j<DimY;j++) { \
								for(int i=0;i<DimX;i++){
#define END_FOR }}}

#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}

#define E_EPSIL  1.e-6f

struct RLERunSegment{

	vector<int> runcode;
	vector<int> runstart;
	vector<float> phi;

};

#define NEGATIVE_RUN  -9999999
#define POSITIVE_RUN   9999999


static inline float HJ_WENO_Coefficients(float q1, float q2, float q3, float q4, float q5){
	float one_sixth = 1.f / 6;
	float one_third = 1.f / 3;
	float thirteen_twelve = 13.f / 12.f;
//	float epsil = 1.e-6*max(q1,max(q2,max(q3,max(q4,q5))))+1.e-99;
	float s1 = thirteen_twelve*(q1-2*q2+q3)*(q1-2*q2+q3) + 0.25f*(q1-4*q2+3*q3)*(q1-4*q2+3*q3);
	float s2 = thirteen_twelve*(q2-2*q3+q4)*(q2-2*q3+q4) + 0.25f*(q2-q4)*(q2-q4);
	float s3 = thirteen_twelve*(q3-2*q4+q5)*(q3-2*q4+q5) + 0.25f*(3*q3-4*q4+q5)*(3*q3-4*q4+q5);
	float a1 = 0.1f/((E_EPSIL+s1)*(E_EPSIL+s1));
	float a2 = 0.6f/((E_EPSIL+s2)*(E_EPSIL+s2));
	float a3 = 0.3f/((E_EPSIL+s3)*(E_EPSIL+s3));
	float at = 1.f / (a1 + a2 + a3);
	float b1 = a1 * at;
	float b2 = a2 * at;
	float b3 = a3 * at;
	return b1*(q1 * one_third - 7.f * one_sixth * q2 + 11.f * one_sixth * q3) +
		   b2*(-q2 * one_sixth + 5.f * one_sixth * q3 + q4 * one_third) +
		   b3*(q3 * one_third + 5.f * one_sixth * q4 - q5 * one_sixth);
};

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
 		return fabsf(phi_p) > fabsf(phi_m)? phi_p : phi_m;
 	else{
 		printf("Error! Shouldn't reach here \n");
 		exit(1);
 	}
};

static float WENOCoeffs(int i, int j, int k, int DimX, int DimY, int DimZ,
						float *x0, float inv_delta, char axis, char s){
	float q1, q2, q3, q4, q5;
	switch(axis){
	case 1:
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
		break;
	case 2:
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
		break;
	case 3:
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
		break;
	}
	return HJ_WENO_Coefficients(q1, q2, q3, q4, q5);
};

class RLELevelSet{

public:
	RLELevelSet(int nx, int ny, int nz,
			float _x0, float _y0, float _z0,
			float _x1, float _y1, float _z1,
			float _h, float _thresh, float *value)
		: DimX(nx), DimY(ny), DimZ(nz),
		x0(_x0), y0(_y0), z0(_z0),
		x1(_x1), y1(_y1), z1(_z1),
		h(_h), thresh(_thresh), nthresh(-_thresh){

	runs = new RLERunSegment[nz * ny];
	for(int k = 0; k < DimZ; ++k)
		for(int j = 0; j < DimY; ++j){
			vector<float> vals;
			vector<int> index;
			RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
			for(int i = 0; i < DimX; ++i){
				if(fabs(value[INDEX(i,j,k)])+E_EPSIL > thresh){  // undefined band
					if( i == 0 ){
						if(value[INDEX(i,j,k)] > 0.f)
							run->runcode.push_back(POSITIVE_RUN);
						else
							run->runcode.push_back(NEGATIVE_RUN);
					}
					else if( run->runcode[run->runcode.size()-1] != POSITIVE_RUN &&
							 run->runcode[run->runcode.size()-1] != NEGATIVE_RUN ){
						if(run->runcode.size() > 0)
							run->runstart.push_back(i);
						if(value[INDEX(i,j,k)] > 0.f)
							run->runcode.push_back(POSITIVE_RUN);
						else
							run->runcode.push_back(NEGATIVE_RUN);
					}
				}
				else{
					run->phi.push_back(value[INDEX(i,j,k)]);
					if( i == 0 ){
						run->runcode.push_back(run->phi.size()-1);
//						run->runstart.push_back(i);
					}
					else if( run->runcode[run->runcode.size()-1] == POSITIVE_RUN ||
						     run->runcode[run->runcode.size()-1] == NEGATIVE_RUN ){
						if(run->runcode.size() > 0)
							run->runstart.push_back(i);
						run->runcode.push_back(run->phi.size()-1);
					}
				}

			}
		}
	}

	RLELevelSet(int nx, int ny, int nz,
				float _x0, float _y0, float _z0,
				float _x1, float _y1, float _z1,
				float _h, float _thresh)
			: DimX(nx), DimY(ny), DimZ(nz),
			x0(_x0), y0(_y0), z0(_z0),
			x1(_x1), y1(_y1), z1(_z1),
			h(_h), thresh(_thresh), nthresh(-_thresh){

		runs = new RLERunSegment[nz * ny];

	}

	RLELevelSet(const char *filename){
		   FILE *fp=fopen(filename, "rb");
		   if(!fp){
		      printf("Couldn't open gridimplicit file \"%s\" for reading\n", filename);
		      return;
		   }
		   // read in dimensions

		   fread(&DimX, sizeof(int), 1, fp);
		   fread(&DimY, sizeof(int), 1, fp);
		   fread(&DimZ, sizeof(int), 1, fp);

		   fread(&x0, sizeof(float), 1, fp);
		   fread(&y0, sizeof(float), 1, fp);
		   fread(&z0, sizeof(float), 1, fp);

		   fread(&x1, sizeof(float), 1, fp);
		   fread(&y1, sizeof(float), 1, fp);
		   fread(&z1, sizeof(float), 1, fp);

		   fread(&h, sizeof(float), 1, fp);

		   fread(&thresh, sizeof(float), 1, fp);

		   nthresh = -thresh;

		   runs = new RLERunSegment[DimY * DimZ];

		   // read in samples

		   //for(int i=0; i<dim[0]*dim[1]*dim[2]; ++i){
		   for(int k=0; k<DimZ; ++k)
	   		 for(int j=0; j<DimY; ++j){
	   			RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
	   			unsigned int siz;
	   			fread(&siz, sizeof(unsigned int), 1, fp);
	//   			printf("A. at (%d, %d), siz = %u \n", j, k, siz);
	   			if(siz > 0){
	   				for(int n=0; n < siz; ++n){
	   					int code;
	   					fread(&code, sizeof(int), 1, fp);
	   					run->runcode.push_back(code);
	   				}
	   			}
	   			fread(&siz, sizeof(unsigned int), 1, fp);
	//   			printf("B. at (%d, %d), siz = %u \n", j, k, siz);
	   			if(siz > 0){
	   				for(int n=0; n < siz; ++n){
	   					int code;
	   					fread(&code, sizeof(int), 1, fp);
	   					run->runstart.push_back(code);
	   				}
	   			}
	   			fread(&siz, sizeof(unsigned int), 1, fp);
	//   			printf("C. at (%d, %d), siz = %u \n", j, k, siz);
	   			if(siz > 0){
	   				for(int n=0; n < siz; ++n){
	   					float code;
	   					fread(&code, sizeof(float), 1, fp);
	   					run->phi.push_back(code);
	   				}
	   			}
		   }
		   fclose(fp);
		}

	~RLELevelSet(){
		delete [] runs;
	}

	void FillRunData(int i, int j, int k, float data){
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		if(fabs(data)+E_EPSIL > thresh){  // undefined band
			if( i == 0 ){
				if(data > 0.f)
					run->runcode.push_back(POSITIVE_RUN);
				else
					run->runcode.push_back(NEGATIVE_RUN);
			}
			else if( run->runcode[run->runcode.size()-1] != POSITIVE_RUN &&
					 run->runcode[run->runcode.size()-1] != NEGATIVE_RUN ){
				if(run->runcode.size() > 0)
					run->runstart.push_back(i);
				if(data > 0.f)
					run->runcode.push_back(POSITIVE_RUN);
				else
					run->runcode.push_back(NEGATIVE_RUN);
			}
		}
		else{
			run->phi.push_back(data);
			if( i == 0 ){
				run->runcode.push_back(run->phi.size()-1);
//						run->runstart.push_back(i);
			}
			else if( run->runcode[run->runcode.size()-1] == POSITIVE_RUN ||
					 run->runcode[run->runcode.size()-1] == NEGATIVE_RUN ){
				if(run->runcode.size() > 0)
					run->runstart.push_back(i);
				run->runcode.push_back(run->phi.size()-1);
			}
		}

	}

	void SetRunData(int i, int j, int k, float val) {
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		vector<int>::iterator pos = run->runstart.end();
		pos = upper_bound(run->runstart.begin(), run->runstart.end(), i);
//		printf(" i = %d, pos = %f \n", i, (int)pos);
		if(pos == run->runstart.begin()){ // i belongs to the first run
			int code = run->runcode[0];
			if(code == POSITIVE_RUN)
				;
			else if(code == NEGATIVE_RUN)
				;
			else{
				int index = code + i;
//				printf("pos = %d, index = %d \n", *pos, index);
				run->phi[index] = val;
			}
		}
		else if(pos == run->runstart.end()){ // i belongs to the last run
			int code = run->runcode[run->runcode.size()-1];
			if(code == POSITIVE_RUN)
				;
			else if(code == NEGATIVE_RUN)
				;
			else{
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *pos, index);
				run->phi[index] = val;
			}
		}
		else{
			int dist = distance(run->runstart.begin(), pos);
			int code = run->runcode[dist];
//			printf("dist = %d, code = %d \n", dist, code);
			if(code == POSITIVE_RUN)
				;
			else if(code == NEGATIVE_RUN)
				;
			else{
//				printf("pos = %d \n", *pos);
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *(pos-1), index);
				run->phi[index] = val;
			}
		}
	}


	void Dimensions(int *dims) const{
		dims[0] = DimX;
		dims[1] = DimY;
		dims[2] = DimZ;
	}

	Point GetP0() const{
		return Point(x0, y0, z0);
	}

	Point GetP1() const{
		return Point(x1, y1, z1);
	}

	float GetGridCellSize() const{
		return h;
	}

	void PrintRunStart(int j, int k) const{
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		for(int n=0; n < run->runstart.size(); ++n)
			printf("runstart at [%d] = %d \n", n, run->runstart[n]);
	}

	void PrintRunCode(int j, int k) const{
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		for(int n=0; n < run->runcode.size(); ++n)
			printf("runcode at [%d] = %d \n", n, run->runcode[n]);
	}

	void PrintPhi(int j, int k) const{
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		for(int n=0; n < run->phi.size(); ++n)
			printf("value at [%d] = %f \n", n, run->phi[n]);
	}

	void OutputBinaryGridData(const char *filename) const{
		FILE *fp = fopen(filename, "wb");
		if(!fp){
	      printf("Couldn't open file [ %s ] to write \n", filename);
	      exit(-1);
	   	}
	//   	fprintf(fp,"%d %d %d \n", DimX, DimY, DimZ);
	   	fwrite(&DimX, sizeof(int), 1, fp);
	   	fwrite(&DimY, sizeof(int), 1, fp);
	   	fwrite(&DimZ, sizeof(int), 1, fp);

	//   	fprintf(fp,"%f %f %f \n",p0.x, p0.y, p0.z);
	   	fwrite(&x0, sizeof(float), 1, fp);
	   	fwrite(&y0, sizeof(float), 1, fp);
	   	fwrite(&z0, sizeof(float), 1, fp);

	//    fprintf(fp,"%f %f %f \n",p1.x, p1.y, p1.z);
	    fwrite(&x1, sizeof(float), 1, fp);
	    fwrite(&y1, sizeof(float), 1, fp);
	    fwrite(&z1, sizeof(float), 1, fp);

	    fwrite(&h, sizeof(float), 1, fp);

	    fwrite(&thresh, sizeof(float), 1, fp);


	    for(int k=0;k<DimZ;k++)
			for(int j=0;j<DimY;j++){
					RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
					unsigned int siz = run->runcode.size();
					fwrite(&siz, sizeof(unsigned int), 1, fp);
					for(int n=0; n < siz; ++n)
						fwrite(&(run->runcode[n]), sizeof(int), 1, fp);
					siz = run->runstart.size();
					fwrite(&siz, sizeof(unsigned int), 1, fp);
					for(int n=0; n < siz; ++n)
						fwrite(&(run->runstart[n]), sizeof(int), 1, fp);
					siz = run->phi.size();
					fwrite(&siz, sizeof(unsigned int), 1, fp);
					for(int n=0; n < siz; ++n)
						fwrite(&(run->phi[n]), sizeof(float), 1, fp);

			}
	//				fprintf(fp,"%f\n", phi[k*(DimX*DimY)+j*DimX+i]);
		fclose(fp);

	}

	/*float Query(int i, int j, int k) const{
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		if(run->runcode.size() == 1){
			if(run->runcode[0] == POSITIVE_RUN)
				return thresh;
			else if(run->runcode[0] == NEGATIVE_RUN)
				return -thresh;
			else{
				printf("Error 1\n");
				exit(1);
			}
		}
		vector<int>::iterator pos = run->runstart.end();
		pos = upper_bound(run->runstart.begin(), run->runstart.end(), i);
//		printf(" i = %d, pos = %f \n", i, (int)pos);
		if(pos == run->runstart.begin()){ // i belongs to the first run
			int code = run->runcode[0];
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return -thresh;
			else{
				int index = code + i;
//				printf("pos = %d, index = %d \n", *pos, index);
				return run->phi[index];
			}
		}
		else if(pos == run->runstart.end()){ // i belongs to the last run
			int code = run->runcode[run->runcode.size()-1];
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return -thresh;
			else{
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *pos, index);
				return run->phi[index];
			}
		}
		else{
			int dist = distance(run->runstart.begin(), pos);
			int code = run->runcode[dist];
//			printf("dist = %d, code = %d \n", dist, code);
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return -thresh;
			else{
//				printf("pos = %d \n", *pos);
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *(pos-1), index);
				return run->phi[index];
			}
		}
	}*/

	float &QueryPoint(int i, int j, int k){
		return Query(i,j,k);
	}

	float& Query(int i, int j, int k) {
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		if(run->runcode.size() == 1){
			if(run->runcode[0] == POSITIVE_RUN)
				return thresh;
			else if(run->runcode[0] == NEGATIVE_RUN)
				return nthresh;
			else{
				printf("Error 1\n");
				exit(1);
			}
		}
		vector<int>::iterator pos = run->runstart.end();
		pos = upper_bound(run->runstart.begin(), run->runstart.end(), i);
//		printf(" i = %d, pos = %f \n", i, (int)pos);
		if(pos == run->runstart.begin()){ // i belongs to the first run
			int code = run->runcode[0];
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return nthresh;
			else{
				int index = code + i;
//				printf("pos = %d, index = %d \n", *pos, index);
				return run->phi[index];
			}
		}
		else if(pos == run->runstart.end()){ // i belongs to the last run
			int code = run->runcode[run->runcode.size()-1];
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return nthresh;
			else{
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *pos, index);
				return run->phi[index];
			}
		}
		else{
			int dist = distance(run->runstart.begin(), pos);
			int code = run->runcode[dist];
//			printf("dist = %d, code = %d \n", dist, code);
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return nthresh;
			else{
//				printf("pos = %d \n", *pos);
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *(pos-1), index);
				return run->phi[index];
			}
		}
	}

	const float& Query(int i, int j, int k) const {
		RLERunSegment *run = &(runs[INDEX_RUN(j,k)]);
		if(run->runcode.size() == 1){
			if(run->runcode[0] == POSITIVE_RUN)
				return thresh;
			else if(run->runcode[0] == NEGATIVE_RUN)
				return nthresh;
			else{
				printf("Error 1\n");
				exit(1);
			}
		}
		vector<int>::iterator pos = run->runstart.end();
		pos = upper_bound(run->runstart.begin(), run->runstart.end(), i);
//		printf(" i = %d, pos = %f \n", i, (int)pos);
		if(pos == run->runstart.begin()){ // i belongs to the first run
			int code = run->runcode[0];
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return nthresh;
			else{
				int index = code + i;
//				printf("pos = %d, index = %d \n", *pos, index);
				return run->phi[index];
			}
		}
		else if(pos == run->runstart.end()){ // i belongs to the last run
			int code = run->runcode[run->runcode.size()-1];
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return nthresh;
			else{
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *pos, index);
				return run->phi[index];
			}
		}
		else{
			int dist = distance(run->runstart.begin(), pos);
			int code = run->runcode[dist];
//			printf("dist = %d, code = %d \n", dist, code);
			if(code == POSITIVE_RUN)
				return thresh;
			else if(code == NEGATIVE_RUN)
				return nthresh;
			else{
//				printf("pos = %d \n", *pos);
				int index = code + i - *(pos-1);
//				printf("pos = %d, index = %d \n", *(pos-1), index);
				return run->phi[index];
			}
		}
	}

	void Reinitialize(){

		float *phi = new float[DimX*DimY*DimZ];
		char *needed = new char[DimX*DimY*DimZ];
		memset(needed, 0, sizeof(char)*DimX*DimY*DimZ);
		FOR_EACH_CELL
			phi[INDEX(i,j,k)] = Query(i,j,k);
			if(fabsf(phi[INDEX(i,j,k)])+E_EPSIL < thresh)
				needed[INDEX(i,j,k)] = 1;
		END_FOR
		float dtau = h / 8;
		int iterations = int(round(thresh/dtau));
#ifdef CUDA
		int dims[3];
		dims[0] = DimX; dims[1] = DimY; dims[2] = DimZ;
		ReinitializeCUDANeededRLE(phi,  needed, iterations, 20, 20, 20,
		 						  h, dtau, dims);
#else
		u_int totalCells = DimX*DimY*DimZ;
		float delta2 = h * h;
		float one_third = 1.f / 3.f;
		float two_third = 2.f / 3.f;
		printf(" Reinitializing with dtau = %f for %d iterations \n", dtau, iterations);
		float *phiPos = new float[totalCells];
		float *phiNeg = new float[totalCells];
		float *phi0   = new float[totalCells];
		float *S      = new float[totalCells];
		memcpy(phiPos, phi, sizeof(float)*totalCells);
		memcpy(phiNeg, phi, sizeof(float)*totalCells);
		memcpy(phi0,   phi, sizeof(float)*totalCells);
		memcpy(S,      phi, sizeof(float)*totalCells);
		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);
			if(needed[pos])
				S[pos] = phi[pos] / sqrtf(phi[pos] * phi[pos] + delta2);
		END_FOR

		for(int n = 0; n < iterations; ++n){
			SWAP(phi, phi0);
			EulerStep(dtau, needed, S, phiPos, phi0);
			EulerStep(dtau, needed, S, phiNeg, phiPos);
			FOR_EACH_CELL
				u_int pos = INDEX(i,j,k);
				if(needed[pos])
					phiPos[pos] = 0.75f*phi0[pos] + 0.25f * phiNeg[pos];
			END_FOR
			EulerStep(dtau, needed, S, phiNeg, phiPos);
			FOR_EACH_CELL
				u_int pos = INDEX(i,j,k);
				if(needed[pos])
					phi[pos] = one_third * phi0[pos] + two_third * phiNeg[pos];
			END_FOR
		}

		delete [] S;
		delete [] phiPos;
		delete [] phiNeg;
		delete [] phi0;
#endif

		FOR_EACH_CELL
			SetRunData(i,j,k,phi[INDEX(i,j,k)]);
		END_FOR

		delete [] phi;
		delete [] needed;

	}

	void EulerStep(float dtau, char *needed, const float *S, float *value, float *value0){


		float s_xp, s_yp, s_zp;
		float s_xm, s_ym, s_zm;
		float phi_x, phi_y, phi_z;
		float phi_xm, phi_ym, phi_zm;
		float phi_xp, phi_yp, phi_zp;

		float delta = h;
		float delta2 = delta * delta;
		float inv_delta = 1.f / delta;


		float *x0 = value0;

		FOR_EACH_CELL
			u_int pos = INDEX(i,j,k);

			if(needed[pos]){

				float s = S[pos];

				phi_xm = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1, -1);
				phi_xp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 1,  1);
				phi_ym = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 2, -1);
				phi_yp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 2,  1);
				phi_zm = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 3, -1);
				phi_zp = WENOCoeffs(i, j, k, DimX, DimY, DimZ, x0, inv_delta, 3,  1);

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
				float cfl = dtau / delta * (fabsf(s*phi_x)+fabsf(s*phi_y)+fabsf(s*phi_z))/sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z);
				if(cfl > 1.f){
	//			if(old_phi*value[pos]< 0.f)
	//				if(i == I && j == J && k == K){
					printf(" (reinit) phi changes sign at i=%d, j=%d, k=%d, old = %f, new = %f cfl = %f\n",
							i, j, k, value0[pos], value[pos], cfl);
					printf("S = %f, phi_x = %f, phi_y = %f, phi_z = %f, sqrt = %f \n ",
							s, phi_x, phi_y, phi_z, sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z));
					printf(" reinit cfl > 1, exiting \n ");
					exit(1);
				}
				value[pos] = value0[pos] - dtau * s * ( sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) - 1.f );

//			if(old_phi*value[pos]< 0.f)
//				if(i == 150 && j == 86 && k == 67){
//					printf(" (reinit) phi changes sign at i=%d, j=%d, k=%d, old = %f, new = %f \n",
//							i, j, k, value0[pos], value[pos]);
//					printf("S = %f, phi_x = %f, phi_y = %f, phi_z = %f \n ", s, phi_x, phi_y, phi_z);
//					printf("S = %f, phi_x = %f, phi_y = %f, phi_z = %f, sqrt = %f \n ",
//							s, phi_x, phi_y, phi_z, sqrtf(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z));
//				}
//			if(isnan(value[pos])){
//				printf(" (reinit) phi changes sign at i=%d, j=%d, k=%d, old = %f, new = %f s = %f\n",
//								i, j, k, value0[pos], value[pos], s);
//				exit(1);
//			}
		}
	END_FOR

	}

//	void ConstructMesh(){
//		mc = new CIsoSurface<float>();
////		float&(*ptrFunc)(int, int, int) = &QueryPoint;
//		float *phi = new float[DimX*DimY*DimZ];
//		FOR_EACH_CELL
//			phi[INDEX(i,j,k)] = Query(i,j,k);
//		END_FOR
//		mc->GenerateSurface(phi, 0.f, DimX-1, DimY-1, DimZ-1, h, h, h);
//		delete [] phi;
//	}
//	void DestroyMesh(){
//		mc->DeleteSurface();
//		delete mc;
//	}
//	void OutputGridDataRib(char *fname){
//		mc->RIBDump(fname);
//	}


private:

	int DimX;
	int DimY;
	int DimZ;
	float x0, y0, z0;
	float x1, y1, z1;
	float h; // uniform grid size in all 3 dimensions
	float thresh;
	float nthresh;
	RLERunSegment *runs;
//	CIsoSurface<float> *mc;

};

#endif /*RLELEVELSETS_H_*/
