#include <stdio.h>
#include "geometry.h"
#include "voxel.h"
#include "particle.h"

#define INDEX(i,j,k) (k)*dim[0]*dim[1]+(j)*dim[0]+(i)
/* some macros to keep things cleaner */
#define SQ(val) ((val) * (val))
#define CUBE(val) ((val) * (val) * (val))
#define DIST(x1, y1, z1, x2, y2, z2) (SQ(x2 - x1) + SQ(y2 - y1) + SQ(z2-z1))
/* coefficients of influence equation */
#define A -.444444
#define B 1.888889
#define C -2.444444

static int dim[3];   
static float delta;
static float *phi;
static float thresh; 

void set_buf(const Voxel &voxel, float x, float y, float z, float strength)
{
   int i, j, k;
   float dist, effect;
   for(int k=0; k<dim[2]; ++k)
  		 for(int j=0; j<dim[1]; ++j)
   			 for(int i=0; i<dim[0]; ++i){	   	   	   			 
		   	   	  Point center = voxel.VoxelCenterPosition(i,j,k); 
/* only need to update things if distance from point(i,j,k) to metaball is
   within its sphere of influence */
		if((dist = DIST(center.x, center.y, center.z, x, y, z)) <= strength)	{
/* compute influence and update field */
		   effect = A * CUBE(dist)/CUBE(strength) + 
			B * SQ(dist)/SQ(strength) + C * dist/strength + 1.0;
	      	   phi[INDEX(i,j,k)] += effect;
		   }
   	}
}
void output(const Voxel &voxel, char *filename) {
	FILE *fp = fopen(filename, "wb");
	if(!fp){
      printf("Couldn't open file [ %s ] to write \n", filename);
      exit(-1);
   	}
//   	fprintf(fp,"%d %d %d \n", DimX, DimY, DimZ);
   	fwrite(&dim[0], sizeof(int), 1, fp);
   	fwrite(&dim[1], sizeof(int), 1, fp);
   	fwrite(&dim[2], sizeof(int), 1, fp);
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

    for(int k=0;k<dim[2];k++)
		for(int j=0;j<dim[1];j++)
			for(int i=0;i<dim[0];i++)
				fwrite(&phi[INDEX(i,j,k)], sizeof(float), 1, fp);
//				fprintf(fp,"%f\n", phi[k*(DimX*DimY)+j*DimX+i]);
	fclose(fp);

}

int main(int argc, char ** argv){

	char *filename = argv[1];
	thresh = atoi(argv[2]);	
	
#ifdef WIN32 // NOBOOK
	printf("WIN32 is defined!\n");
#else // NOBOOK
	printf("WIN32 is not defined!\n");
#endif // NOBOOK
	
	FILE *fp=fopen(filename, "rb");
	   if(!fp){
	      printf("Couldn't open gridimplicit file \"%s\" for reading\n", filename);
	      return;
	   }
	   
	   fread(&dim[0], sizeof(int), 1, fp);
	   fread(&dim[1], sizeof(int), 1, fp);
	   fread(&dim[2], sizeof(int), 1, fp);
	   if(dim[0]<=1 || dim[1]<=1 || dim[2]<=1){
	      printf("Bad dimensions (%d, %d, %d) in gridimplicit file \"%s\"\n", dim[0], dim[1], dim[2], filename);
	      dim[0]=dim[1]=dim[2]=0;
	      return;
	   }	   
	   printf("dimx = %d, dimy = %d, dimz = %d \n", dim[0], dim[1], dim[2]);
	   Point p0;
	   
	   fread(&p0.x, sizeof(float), 1, fp);
	   fread(&p0.y, sizeof(float), 1, fp);
	   fread(&p0.z, sizeof(float), 1, fp);
	   Point p1;
	   
	   fread(&p1.x, sizeof(float), 1, fp);
	   fread(&p1.y, sizeof(float), 1, fp);
	   fread(&p1.z, sizeof(float), 1, fp);
	  
	   fread(&delta, sizeof(float), 1, fp);
	   printf("p0 at (%f,%f,%f) p1 at (%f,%f,%f) delta = %f\n",
			   p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, delta);
	   Voxel voxel(p0, p1, delta);	  
	   // read in samples
	   phi=new float[dim[0]*dim[1]*dim[2]];
	   //for(int i=0; i<dim[0]*dim[1]*dim[2]; ++i){
	   for(int k=0; k<dim[2]; ++k)
	   		 for(int j=0; j<dim[1]; ++j)
	   			 for(int i=0; i<dim[0]; ++i){	      
	   		if(1!=fread(&phi[INDEX(i,j,k)], sizeof(float), 1, fp)){
			  printf("Problem reading gridimplicit file \"%s\" (at value %d)\n", filename, i);
			  dim[0]=dim[1]=dim[2]=0;
			  delete[] phi;
			  phi=0;
			  return;
		  }
	      assert(phi[i]==phi[i]);
	   }
	   fclose(fp);  
	   for(int k=0; k<dim[2]; ++k)
	    		 for(int j=0; j<dim[1]; ++j)
	     			 for(int i=0; i<dim[0]; ++i){	
	     				 phi[INDEX(i,j,k)] = 0.f;
	     			 }
	   set_buf(voxel, 10.5f, 10.5, 10.5, SQ(thresh));
	   
	   set_buf(voxel, 10.5f, 10.5, 15.5, SQ(thresh));
	   for(int k=0; k<dim[2]; ++k)
		 for(int j=0; j<dim[1]; ++j)
 			 for(int i=0; i<dim[0]; ++i){	
 				 phi[INDEX(i,j,k)] = -1 * phi[INDEX(i,j,k)];
 				 if(phi[INDEX(i,j,k)] == 0.f)
 					 phi[INDEX(i,j,k)] = 0.1f;
 			 }
//	   for(int k=0; k<dim[2]; ++k)
//   		 for(int j=0; j<dim[1]; ++j)
//    			 for(int i=0; i<dim[0]; ++i){	
////    				 if(phi[INDEX(i,j,k)] != 99.0)
//    				 if(fabs(phi[INDEX(i,j,k)]- 0.f) > 1.e-3)
//    					 printf(" phi = %f at (%d, %d, %d)\n", phi[INDEX(i,j,k)],i,j,k );
//    				 
//    			 }
	   char particle[] = "metaball.f0";
	   output(voxel, particle);
	   
	   delete [] phi;
	
	
}