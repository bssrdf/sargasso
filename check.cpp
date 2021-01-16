#include <stdio.h>
#include "geometry.h"
#include "scalarfield3d.h"
#include "physicsworld.h"

#include "rlelevelsets.h"

#define INDEX(i,j,k) (k)*dim[0]*dim[1]+(j)*dim[0]+(i)

int main(int argc, char ** argv){

	char *filename = argv[1];
	char *rlename  = argv[4];
	int J = atoi(argv[2]);
	int K = atoi(argv[3]);
	int dim[3];   
	float delta;
	float *phi;
	int I = 41;
	
	char *objectfile = "container";
	
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

	   // read in dimensions
	   /*if(3!=fscanf(fp, "%d %d %d", dim, dim+1, dim+2)){
	      Error("Problem reading gridimplicit file \"%s\" (dimensions)\n", filename);
	      dim[0]=dim[1]=dim[2]=0;
	      return;
	   }*/
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
	   fread(&delta, sizeof(float), 1, fp);
	   printf("p0 at (%f,%f,%f) p1 at (%f,%f,%f) delta = %f\n",
			   p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, delta);
	   
	   // read in samples
	   phi=new float[dim[0]*dim[1]*dim[2]];
	   //for(int i=0; i<dim[0]*dim[1]*dim[2]; ++i){
	   for(int k=0; k<dim[2]; ++k)
	   		 for(int j=0; j<dim[1]; ++j)
	   			 for(int i=0; i<dim[0]; ++i){
	      /*if(1!=fscanf(fp, "%f", phi+i)){
	         Error("Problem reading gridimplicit file \"%s\" (at value %d)\n", filename, i);
	         dim[0]=dim[1]=dim[2]=0;
	         delete[] phi;
	         phi=0;
	         return;
	      }*/
		  //if(1!=fread(phi+i, sizeof(float), 1, fp)){
	   		if(1!=fread(&phi[INDEX(i,j,k)], sizeof(float), 1, fp)){
			  printf("Problem reading gridimplicit file \"%s\" (at value %d)\n", filename, i);
			  dim[0]=dim[1]=dim[2]=0;
			  delete[] phi;
			  phi=0;
			  return;
		  }
	      assert(phi[i]==phi[i]);
	   }
	   PhysicsWorld *physWorld = new PhysicsWorld(p0, p1);
	   physWorld->BuildPhysicsWorld("scene.desc");

	   ScalarField3D object(physWorld, dim[0], dim[1], dim[2], I, J, K);
	   object.ReadBinaryGridData(objectfile);
	   float *phiobj = object.getScalar();	   
	   
//	   int k = K;
//	   for(int j=0; j<dim[1]; ++j){
//		   for(int i=0; i<dim[0]; ++i){		   
////			   if(phiobj[INDEX(i,j,k)] < 0.f && phi[INDEX(i,j,k)] > 0.f){
//			if(i == I && j == J && k == K){
//				   printf(" at (%d, %d, %d), phi = %f, phi[k-1] = %f \n", i, j, k, 
//					 phi[INDEX(i,j,k)],  phi[INDEX(i,j,k-1)]);
//			   printf(" at (%d, %d, %d), phi[k+1] = %f, phi[i-1] = %f \n", i, j, k, 
//			   			phi[INDEX(i,j,k+1)],  phi[INDEX(i-1,j,k)]);
//			   printf(" at (%d, %d, %d), phi[j+1] = %f, phi[i+1] = %f \n", i, j, k, 
//			   			   phi[INDEX(i,j+1,k)],  phi[INDEX(i+1,j,k)]);
//			   printf(" at (%d, %d, %d), phi[j-1] = %f  \n", i, j, k, 
//			   			   phi[INDEX(i,j-1,k)] );
//			   }
//		   }
//	   }
	   int k = K;
	   int j = J;
	      
//	   RLELevelSet rle(dim[0], dim[1], dim[2], 
//			   p0.x, p0.y, p0.z,
//			   p1.x, p1.y, p1.z, 
//			   delta, 5*delta, phi);
	   
	   	   
	   RLELevelSet rle(rlename);  
	   	   
	   printf(" ************* \n");
	   rle.PrintRunStart(j,k);
	   printf(" ************* \n");
	   rle.PrintRunCode(j,k);
	   printf(" ************* \n");
	   rle.PrintPhi(j,k);	   	   
	   for(int i=0; i<dim[0]; ++i)
	   		printf(" at (%d, %d, %d), phi = %f rle_phi = %f\n", i, j, k, 
	     				 phi[INDEX(i,j,k)], rle.Query(i,j,k));
//	   		printf(" at (%d, %d, %d), phi = %f \n", i, j, k, 
//	   	     				 phi[INDEX(i,j,k)]);
	   
	   char rhead[20]= "";
	   char rtail[] = ".rle";
	   strcat(rhead, filename);
	   strcat(rhead, rtail);
	   rle.OutputBinaryGridData(rhead);
   
//	   int k = K;
//	   int i = 114;
//	   int j = 5;
//   	   for(int k=0; k<dim[2]; ++k){
//   		   if(phi[INDEX(i,j,k)] < 0.f)
//   			   printf(" at (%d, %d, %d), phi = %f \n", i, j, k, phi[INDEX(i,j,k)]);
//   	   }   
//	   float phimax = -999.f;
//	   float phimin = 999.f;
//	   for(int k=0; k<dim[2]; ++k){
//		   for(int j=0; j<dim[1]; ++j){
//	   		   for(int i=0; i<dim[0]; ++i){		   
//	   			   if(phimax < phi[INDEX(i,j,k)])
////	   				   printf(" at (%d, %d, %d), phi = %f, phi[k-1] = %f \n", i, j, k, 
////	   					 phi[INDEX(i,j,k)],  phi[INDEX(i,j,k-1)]);
//	   				   phimax = phi[INDEX(i,j,k)];
//	   			   if(phimin > phi[INDEX(i,j,k)])
//	   				   phimin = phi[INDEX(i,j,k)];
//	   				
//	   		   }
//	   	   }
//	   }
//	   printf("max phi = %f, min phi = %f \n", phimax, phimin);
	   // finished with the file   
	   
	   fclose(fp);
	   
	   delete [] phi;
	   delete physWorld;

	
}
