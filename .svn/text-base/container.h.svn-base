/*
 * container.h
 *
 *  Created on: Sep 25, 2008
 *      Author: bzhao
 */

#ifndef CONTAINER_H_
#define CONTAINER_H_

#include "voxel.h"
#include "transform.h"
#include "physicsworld.h"

class Container{
public:
	Container(int X, int Y, int Z, float f):
	DimX(X), DimY(Y), DimZ(Z),
	friction(f), mHasSignedDistanceMap(false){
		phi  = new float[X*Y*Z];
		phiu = new float[X*Y*Z];
		phiv = new float[X*Y*Z];
		phiw = new float[X*Y*Z];
		memset(phi,  0, X*Y*Z*sizeof(float));
		memset(phiu, 0, X*Y*Z*sizeof(float));
		memset(phiv, 0, X*Y*Z*sizeof(float));
		memset(phiw, 0, X*Y*Z*sizeof(float));

	}
	virtual ~Container(){
		delete [] phi;
		delete [] phiu;
		delete [] phiv;
		delete [] phiw;
	}

	virtual bool IsSignedDistance() const = 0;

	bool InContainer(int i, int j, int k) const{
		if(phi[INDEX(i,j,k)] >= 0.f)
			return true;
		else
			return false;
	}

	virtual bool IsBorderCell(const Voxel &voxel, int i, int j, int k) const = 0;

	virtual void SetBorderCellPhi(const Voxel &voxel, int i, int j, int k, int n, float *f) const = 0;
	virtual void SetNormal(const Voxel &voxel, Vector &N, int i, int j, int k) const = 0;

	void MergeWith(const float *phi0, float *phi1){
		FOR_EACH_CELL
			u_int index = INDEX(i,j,k);
			if(phi0[index] >= 0.f && phi1[index] >= 0.f)
				phi1[index] = max(phi0[index], phi1[index]);
			else if(phi0[index] >= 0.f && phi1[index] < 0.f)
				phi1[index] = phi0[index];
			else if(phi0[index] < 0.f && phi1[index] < 0.f)
				phi1[index] = fabs(phi0[index]) < fabs(phi1[index]) ?
								phi0[index] : phi1[index];
		END_FOR
	}

	void BuilSignedDistanceMap(const Voxel *voxel, const PhysicsWorld *pw){
		float *phitmp  = new float[DimX*DimY*DimZ];
		memset(phitmp,  0, DimX*DimY*DimZ*sizeof(float));
		pw->BuildSignedDistanceMap(voxel, 0, phitmp);
		MergeWith(phitmp, phi);
		memset(phitmp,  0, DimX*DimY*DimZ*sizeof(float));
		pw->BuildSignedDistanceMap(voxel, 1, phitmp);
		MergeWith(phitmp, phiu);
		memset(phitmp,  0, DimX*DimY*DimZ*sizeof(float));
		pw->BuildSignedDistanceMap(voxel, 2, phitmp);
		MergeWith(phitmp, phiv);
		memset(phitmp,  0, DimX*DimY*DimZ*sizeof(float));
		pw->BuildSignedDistanceMap(voxel, 3, phitmp);
		MergeWith(phitmp, phiw);

		mHasSignedDistanceMap = false;
		delete [] phitmp;
	}

	void Friction(const Voxel &voxel, int vel_index, int i, int j, int k, float *f) const{
		if(vel_index == 1){
			if(phiu[INDEX(i,j,k)] >= 0.f)
				*f = friction;
		}
		else if(vel_index == 2){
			if(phiv[INDEX(i,j,k)] >= 0.f)
				*f = friction;
		}
		else if(vel_index == 3){
			if(phiw[INDEX(i,j,k)] >= 0.f)
				*f = friction;
		}
	}
	float *GetScalarPhi(int axis){
		if(axis == 0)
			return phi;
		else if(axis == 1)
			return phiu;
		else if(axis == 2)
			return phiv;
		else if(axis == 3)
			return phiw;
		else{
			printf("(Container::GetScalarPhi) Wrong! axis  = %d,  must be 0-3 \n", axis);
			exit(1);
		}
	}

	void OutputBinaryGridData(char *filename){
		FILE *fp = fopen(filename, "wb");
		if(!fp){
		  printf("Couldn't open file [ %s ] to write \n", filename);
		  exit(1);
		}
		fwrite(&DimX, sizeof(int), 1, fp);
		fwrite(&DimY, sizeof(int), 1, fp);
		fwrite(&DimZ, sizeof(int), 1, fp);

		for(int k=0;k<DimZ;k++)
			for(int j=0;j<DimY;j++)
				for(int i=0;i<DimX;i++)
					fwrite(&phi[INDEX(i,j,k)], sizeof(float), 1, fp);
		for(int k=0;k<DimZ;k++)
			for(int j=0;j<DimY;j++)
				for(int i=0;i<DimX;i++)
					fwrite(&phiu[INDEX(i,j,k)], sizeof(float), 1, fp);
		for(int k=0;k<DimZ;k++)
			for(int j=0;j<DimY;j++)
				for(int i=0;i<DimX;i++)
					fwrite(&phiv[INDEX(i,j,k)], sizeof(float), 1, fp);
		for(int k=0;k<DimZ;k++)
			for(int j=0;j<DimY;j++)
				for(int i=0;i<DimX;i++)
					fwrite(&phiw[INDEX(i,j,k)], sizeof(float), 1, fp);

		fclose(fp);

	}

	void LoadPhi(const char *file){

		int dim[3];

		FILE *fp=fopen(file, "rb");
	   if(!fp){
		  printf("Couldn't open file \"%s\" for reading\n", file);
		  exit(1);
	   }

	   // read in dimensions
	   fread(&dim[0], sizeof(int), 1, fp);
	   fread(&dim[1], sizeof(int), 1, fp);
	   fread(&dim[2], sizeof(int), 1, fp);
	   if(dim[0]<=1 || dim[1]<=1 || dim[2]<=1){
		  printf("Bad dimensions (%d, %d, %d) in gridimplicit file \"%s\"\n", dim[0], dim[1], dim[2], file);
		  dim[0]=dim[1]=dim[2]=0;
		  exit(1);
	   }
	   if(dim[0]!=DimX){
		   printf("dim[%d] != DimX, dim[%d] = %d, DimX = %d \n", 0, 0, dim[0], DimX);
		   exit(1);
	   }
		if(dim[1]!=DimY){
		   printf("dim[%d] != DimY, dim[%d] = %d, DimY = %d \n", 1, 1, dim[1], DimY);
		   exit(1);
	   }
		if(dim[2]!=DimZ){
		   printf("dim[%d] != DimZ, dim[%d] = %d, DimZ = %d \n", 2, 2, dim[2], DimZ);
		   exit(1);
	   }

	   printf("Read in DimX = %d, DimY = %d, DimZ = %d \n", DimX, DimY, DimZ);

	  // read in samples
	  //for(int i=0; i<dim[0]*dim[1]*dim[2]; ++i){
	  for(int k=0; k<dim[2]; ++k)
		for(int j=0; j<dim[1]; ++j)
			for(int i=0; i<dim[0]; ++i){
			   if(1!=fread(&phi[INDEX(i,j,k)], sizeof(float), 1, fp)){
					 printf("Problem reading gridimplicit file \"%s\" at value (%d, %d, %d)\n",
							 file, i, j, k);
					 exit(1);
			   }
			   //assert(phi[i]==phi[i]);
	  }
	  for(int k=0; k<dim[2]; ++k)
		for(int j=0; j<dim[1]; ++j)
			for(int i=0; i<dim[0]; ++i){
			   if(1!=fread(&phiu[INDEX(i,j,k)], sizeof(float), 1, fp)){
					 printf("Problem reading gridimplicit file \"%s\" at value (%d, %d, %d)\n",
							 file, i, j, k);
					 exit(1);
			   }
			//   assert(phiu[i]==phiu[i]);
	  }
	  for(int k=0; k<dim[2]; ++k)
		for(int j=0; j<dim[1]; ++j)
			for(int i=0; i<dim[0]; ++i){
			   if(1!=fread(&phiv[INDEX(i,j,k)], sizeof(float), 1, fp)){
					 printf("Problem reading gridimplicit file \"%s\" at value (%d, %d, %d)\n",
							 file, i, j, k);
					 exit(1);
			   }
			  // assert(phiw[i]==phiw[i]);
	  }
	  for(int k=0; k<dim[2]; ++k)
		for(int j=0; j<dim[1]; ++j)
			for(int i=0; i<dim[0]; ++i){
			   if(1!=fread(&phiw[INDEX(i,j,k)], sizeof(float), 1, fp)){
					 printf("Problem reading gridimplicit file \"%s\" at value (%d, %d, %d)\n",
							 file, i, j, k);
					 exit(1);
			   }
		//	   assert(phi[i]==phi[i]);
	  }
	  fclose(fp);
}

protected:
	float *phi;
	float *phiu, *phiv, *phiw;
	int DimX, DimY, DimZ;
	float friction;
	bool mHasSignedDistanceMap;
};


#endif /* CONTAINER_H_ */
