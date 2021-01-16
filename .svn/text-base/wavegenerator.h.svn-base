/*
 * wavegenerator.h
 *
 *  Created on: Nov 16, 2010
 *      Author: Administrator
 */

#ifndef WAVEGENERATOR_H_
#define WAVEGENERATOR_H_

class WaveGenerator{
public:
	WaveGenerator(float s, float trans, const char *phifile, const char *velfile)
	: scale(s), translate_vel(trans*s) {

	   FILE *fp=fopen(phifile, "rb");
	   if(!fp){
		  printf("Couldn't open wave initial level set file \"%s\" for reading\n", phifile);
		  return;
	   }

	   fread(&DimX, sizeof(int), 1, fp);
	   fread(&DimZ, sizeof(int), 1, fp);

	   phi = new float[DimX*DimZ];
	   for(int k=0; k<DimZ; ++k){
		   for(int i=0; i<DimX; ++i){
				if(1!=fread(&phi[k*DimX+i], sizeof(float), 1, fp)){
				  printf("Problem reading wave initial level set file \"%s\" (at value %d, %d)\n",
						  phifile, i, k);
				  delete [] phi;
				  phi=0;
				}
		   }
	   }
	   fclose(fp);

	   fp=fopen(velfile, "rb");
	   if(!fp){
		  printf("Couldn't open wave initial level set file \"%s\" for reading\n", phifile);
		  return;
	   }

	   int DimXvel;
	   int DimZvel;
	   fread(&DimXvel, sizeof(int), 1, fp);
	   fread(&DimZvel, sizeof(int), 1, fp);
	   if(DimXvel == DimX && DimZvel == DimZ){
		   u = new float[DimX*DimZ];
		   w = new float[DimX*DimZ];
		   for(int k=0; k<DimZ; ++k){
			   for(int i=0; i<DimX; ++i){
					if(1!=fread(&u[k*DimX+i], sizeof(float), 1, fp)){
					  printf("Problem reading wave initial vel file \"%s\" (at value %d, %d)\n",
							  velfile, i, k);
					  delete [] u;
					  u=0;
					}
					if(1!=fread(&w[k*DimX+i], sizeof(float), 1, fp)){
					  printf("Problem reading wave initial vel file \"%s\" (at value %d, %d)\n",
							  velfile, i, k);
					  delete [] w;
					  w=0;
					}
			   }
		   }
	   }
	   else{
		   printf(" dimensions for wave init vel file differs from that of wave level set \n");
		   printf(" DimX for level set is %d for vel is %d \n", DimX, DimXvel);
		   printf(" DimZ for level set is %d for vel is %d \n", DimZ, DimZvel);
	   }
	   fclose(fp);
	}

	~WaveGenerator(){
		DeleteMemory();
	}
	void SetInitialWaveLevelset(const Voxel *voxel, int nx, int ny, int nz,
					int offx, int offy, int offz, float *f) const{

//		for(int k=offz; k<nz-offz; ++k){
//			for(int j=offy; j<ny-offy; ++j){
//				for(int i=offx; i<nx-offx; ++i){
//						f[k*nx*ny+j*nx+i] = scale*phi[(k-offz)*DimX+i-offx];
//				}
//			}
//		}
		for(int k=0; k<nz; ++k){
			for(int j=0; j<ny; ++j){
				for(int i=0; i<nx; ++i){
						f[k*nx*ny+j*nx+i] = scale*phi[(k)*DimX+i];
				}
			}
		}
	}

	void SetInitialWaveVelocity(const Voxel *voxel, int nx, int ny, int nz,
					int offx, int offy, int offz, float *ui, float *wi) const{

//		for(int k=offz; k<nz-offz; ++k){
//			for(int j=offy; j<ny-offy; ++j){
//				for(int i=offx; i<nx-offx; ++i){
//					if(!voxel->InSolid(i,j,k) && phi[(k-offz)*DimX+i-offx] <= 0.f){
//						ui[k*nx*ny+j*nx+i] = scale*u[(k-offz)*DimX+i-offx];
//						wi[k*nx*ny+j*nx+i] = scale*w[(k-offz)*DimX+i-offx];
//					}
//				}
//			}
//		}
		for(int k=0; k<nz; ++k){
			for(int j=0; j<ny; ++j){
				for(int i=0; i<nx; ++i){
					if(!voxel->InSolid(i,j,k) && phi[(k)*DimX+i] <= 0.f){
						ui[k*nx*ny+j*nx+i] = scale*u[(k)*DimX+i];
						wi[k*nx*ny+j*nx+i] = scale*w[(k)*DimX+i];
					}
				}
			}
		}

	}

	void SetBoundaryVelocity(const Voxel *voxel,
					int nx, int ny, int nz,
					int offx, int offy, int offz,
					float *u, float *v, float *w) const{

		int i = offx-1;
		for(int k=0; k<nz; ++k){
			for(int j=0; j<ny; ++j){
				if( voxel->InSolid(i,j,k)  &&
				    voxel->InLiquid(i+1,j,k) &&
				   !voxel->InSurface(i+1,j,k) ){
					u[k*nx*ny+j*nx+i] = translate_vel;
				}
			}
		}

		i = nx-offx;
		for(int k=0; k<nz; ++k){
			for(int j=0; j<ny; ++j){
				if( voxel->InSolid(i,j,k)    &&
				    voxel->InLiquid(i-1,j,k) &&
				   !voxel->InSurface(i-1,j,k) ){
					u[k*nx*ny+j*nx+i-1] = translate_vel;
				}
			}
		}

	}

	void DeleteMemory(){
		if(phi)
			delete [] phi;
		if(u)
			delete [] u;
		if(w)
			delete [] w;

	}

private:
	int DimX;
	int DimZ;
	float *phi;
	float *u, *w;
	float scale;
	float translate_vel;
};

#endif /* WAVEGENERATOR_H_ */
