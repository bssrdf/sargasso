/*
 * renderlevelset.cpp
 *
 *  Created on: Nov 11, 2010
 *      Author: bzhao
 */

#include "renderlevelset.h"

void RenderLevelSet::
	OutputBinaryGridData(const Voxel &voxel, const char *filename) const{
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
