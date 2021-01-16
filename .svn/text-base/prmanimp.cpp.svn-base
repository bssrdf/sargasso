#include <stdio.h>
#include <math.h>
#include <ImplicitField.h>
//#include "geometry.h"
//#include "voxel.h"

#define INDEX(i,j,k) (k)*DimX*DimY+(j)*DimX+(i)

	class PrmanGridImplicit: public ImplicitField{
	public:
		PrmanGridImplicit(char *);
		virtual ~PrmanGridImplicit();
		virtual RtFloat Eval(const RtPoint p);
		virtual void GradientEval(RtPoint grad, const RtPoint p);
		//virtual void Range(RtInterval r, const RtPoint corners[8],
		//	RtVolumeHandle h);
		float &val(int i, int j, int k){
			return phi[i+DimX*(j+DimY*k)]; 
		}

		const float &val(int i, int j, int k) const
		  { return phi[i+DimX*(j+DimY*k)]; }
	private:
	//	Voxel *voxel;
		float *phi;
		float delta;
		int DimX, DimY, DimZ;
	};
	PrmanGridImplicit::PrmanGridImplicit(char *filename){
		int dim[3];

	 //printf("Couldn't open gridimplicit file \"%s\" for reading\n", filename);
	   FILE *fp=fopen(filename, "rb");
	   if(!fp){
	      printf("Couldn't open gridimplicit file \"%s\" for reading\n", filename);
	      return;
	   }
	   else{
		  printf("Successfully opened gridimplicit file \"%s\" for reading\n", filename);
	   }
	   // read in dimensions
	   
	   fread(&dim[0], sizeof(int), 1, fp);
	   fread(&dim[1], sizeof(int), 1, fp);
	   fread(&dim[2], sizeof(int), 1, fp);
	   if(dim[0]<=1 || dim[1]<=1 || dim[2]<=1){
	      printf("Bad dimensions (%d, %d, %d) in gridimplicit file \"%s\"\n", dim[0], dim[1], dim[2], filename);
	      dim[0]=dim[1]=dim[2]=0;
	      return;
	   }
	   printf("dimx = %d, dimy = %d, dimz = %d \n", dim[0], dim[1], dim[2]);
	   DimX = dim[0];
	   DimY = dim[1];
	   DimZ = dim[2];
	  // Point p0;
	   
	   fread(&bbox[0], sizeof(float), 1, fp);
	   fread(&bbox[2], sizeof(float), 1, fp);
	   fread(&bbox[4], sizeof(float), 1, fp);
	   //Point p1;
	   
	   fread(&bbox[1], sizeof(float), 1, fp);
	   fread(&bbox[3], sizeof(float), 1, fp);
	   fread(&bbox[5], sizeof(float), 1, fp);
	   
	   fread(&delta, sizeof(float), 1, fp);
       printf("p0 at (%f,%f,%f) p1 at (%f,%f,%f) delta = %f\n",
			  bbox[0], bbox[2], bbox[4], bbox[1], bbox[3], bbox[5], delta);

	   // read in samples
	   phi=new float[dim[0]*dim[1]*dim[2]];
	   //for(int i=0; i<dim[0]*dim[1]*dim[2]; ++i){
	   for(int k=0; k<dim[2]; ++k)
	   		 for(int j=0; j<dim[1]; ++j)
	   			 for(int i=0; i<dim[0]; ++i){
	      
		  //if(1!=fread(phi+i, sizeof(float), 1, fp)){
	   		if(1!=fread(&phi[INDEX(i,j,k)], sizeof(float), 1, fp)){
			  printf("Problem reading gridimplicit file \"%s\" (at value %d)\n", filename, i);
			  dim[0]=dim[1]=dim[2]=0;
			  delete[] phi;
			  phi=0;
			  return;
		  }
	      //assert(phi[i]==phi[i]);
	   }

	   for(int k=0; k<dim[2]; ++k)
	   		 for(int j=0; j<dim[1]; ++j)
	   			 for(int i=0; i<dim[0]; ++i){
			phi[INDEX(i,j,k)] *= -1;
	   }
	   for(int k=0; k<dim[2]; ++k)
	   		 for(int j=0; j<dim[1]; ++j)
	   			 for(int i=0; i<dim[0]; ++i){
			phi[INDEX(i,j,k)] = phi[INDEX(i,j,k)] / 4 + 0.5f;
	   }


		//voxel = new Voxel(p0, p1, delta);
		//bbox[0]=p0.x;
		//bbox[1]=p0.y;
		//bbox[2]=p0.z;
		//bbox[3]=p1.x;
		//bbox[4]=p1.y;
		//bbox[5]=p1.z;
		
		fclose(fp);
		int I = 66;
		int J = 22;
		int K = 9;
		printf(" at (%d, %d, %d) phi = %f \n", I, J, K, phi[INDEX(I,J,K)]);
		/*bbox[0]=-1.;
		bbox[1]=1.;
		bbox[2]=-1.;
		bbox[3]=1.;
		bbox[4]=-1.;
		bbox[5]=1.;*/
	}

	static float geoff(float r2){
		if(r2>=1.f) return 0.f;
		return ((3.f-r2)*r2-3.f)*r2+1.f;
	}
	/*
	 * d geoff(r2)
	 * -----------
	 *    d r2
	 */
	static float dgeoff(float r2){
		if(r2>=1.f) return 0.f;
		return (6.f-3.f*r2)*r2-3.f;
	}
	static float lerp(float phi1, float phi2, float u){
		return (1.f-u) * phi1 + u * phi2;
	}

	static inline float bilerp(float phi1, float phi2, float phi3, float phi4, float u, float v){
		return lerp(lerp(phi1, phi2, u), lerp(phi3, phi4, u), v);
	}

	static inline float trilerp(float phi1, float phi2, float phi3, float phi4,
					  float phi5, float phi6, float phi7, float phi8,
					  float u, float v, float w){
		return bilerp(lerp(phi1, phi2, u), lerp(phi3, phi4, u),
				  lerp(phi5, phi6, u), lerp(phi7, phi8, u),
				  w, v);
	}

	static	void VoxelCornerPosition(int i, int j, int k, float delta, int index,
				const RtBound bbox, RtFloat *pos) {

	if(index == 1){
		pos[0] = bbox[0] +  delta * i;
		pos[1] = bbox[2] +  delta * j;
		pos[2] = bbox[4] +  delta * k;
	}

	else{
		fprintf(stderr, "Error! parameter [index] must be in range [1...8]\n ");
		//exit(1);
	}
	return;
}

	static float TriInterp(float delta, int  DimX, int DimY, int DimZ,
			                 const RtPoint p, const RtBound bbox, float *phi) {

		float x = p[0]/delta - 0.5f;
	    float y = p[1]/delta - 0.5f;
	    float z = p[2]/delta - 0.5f;
	    if(x < 0) x = 0;
	    if(y < 0) y = 0;
	    if(z < 0) z = 0;
		int xi = floor(x);
		int yi = floor(y);
		int zi = floor(z);
		//Point point = voxel.VoxelCenterPosition(xi,yi,zi);
		RtPoint point;
		VoxelCornerPosition(xi,yi,zi,delta,1,bbox,point);
		float u = x - point[0]/delta;
	    float v = y - point[1]/delta;
	    float w = z - point[2]/delta;

	//    printf("\npoint at (%f,%f,%f) \n", p.x, p.y, p.z);
	//    printf(" u = %f, v = %f, w = %f\n", u, v, w);
	//	printf("xi = %d, yi = %d, zi = %d \n\n", xi, yi, zi);
		return trilerp(phi[zi*(DimX*DimY)+yi*DimX+xi], phi[zi*(DimX*DimY)+yi*DimX+xi+1],
					   phi[(zi+1)*(DimX*DimY)+yi*DimX+xi], phi[(zi+1)*(DimX*DimY)+yi*DimX+xi+1],
					   phi[zi*(DimX*DimY)+(yi+1)*DimX+xi], phi[zi*(DimX*DimY)+(yi+1)*DimX+xi+1],
					   phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi], phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi+1],
					   u, v, w);
	}

	float PrmanGridImplicit::Eval(const RtPoint p){
		//Point pos(p[0], p[1], p[2]);
		if(p[0] <= bbox[0] || p[0] >= bbox[1] || p[1] <= bbox[2] || p[1] >= bbox[3] || p[2] <= bbox[4] || p[2] >= bbox[5] )
			return 0.f;
		else
		    return TriInterp(delta, DimX, DimY, DimZ, p, bbox, phi);
		/*RtPoint sq;
		float r2;
	
		sq[0]=p[0]*p[0];
		sq[1]=p[1]*p[1];
		sq[2]=p[2]*p[2];
		if(sq[0]>sq[1]) r2=sq[0]>sq[2]?sq[0]:sq[2];
		else r2=sq[1]>sq[2]?sq[1]:sq[2];
		return geoff(r2);*/
	}

	void PrmanGridImplicit::GradientEval(RtPoint grad, const RtPoint p){
		/*RtPoint sq;
	
		grad[0]=0.;
		grad[1]=0.;
		grad[2]=0.;
		sq[0]=p[0]*p[0];
		sq[1]=p[1]*p[1];
		sq[2]=p[2]*p[2];
		if(sq[0]>sq[1]){
			if(sq[0]>sq[2]) grad[0]=2.*p[0]*dgeoff(sq[0]);
			else grad[2]=2.*p[2]*dgeoff(sq[2]);
		}
		else if(sq[1]>sq[2])
			grad[1]=2.*p[1]*dgeoff(sq[1]);
		else
			grad[2]=2.*p[2]*dgeoff(sq[2]);*/


	   if(p[0] <= bbox[0] || p[0] >= bbox[1] || p[1] <= bbox[2] || p[1] >= bbox[3] || p[2] <= bbox[4] || p[2] >= bbox[5] ){

			grad[0] = 0.f;
			grad[1] = 0.f;
			grad[2] = 0.f;
		}
		else{
		int i=int(p[0]/delta);
		int j=int(p[1]/delta);
		int k=int(p[2]/delta);
	   // we rely on the user providing a band of empty cells around surface
	   if(i<1) i=1; else if(i>DimX-3) i=DimX-3;
	   if(j<1) j=1; else if(j>DimY-3) j=DimY-3;
	   if(k<1) k=1; else if(k>DimZ-3) k=DimZ-3;
	   // figure out coefficients for trilinear interpolation
	   float xm=p[0]/delta-i, xr=1.f-xm;
	   float ym=p[1]/delta-j, yr=1.f-ym;
	   float zm=p[2]/delta-k, zr=1.f-zm;
	   float c0=xr*yr*zr, c1=xr*yr*zm, c2=xr*ym*zr, c3=xr*ym*zm,
	         c4=xm*yr*zr, c5=xm*yr*zm, c6=xm*ym*zr, c7=xm*ym*zm;
	   grad[0]=c0*(val(i+1,j,  k  ) - val(i-1,j,  k  ))
	                   +c1*(val(i+1,j,  k+1) - val(i-1,j,  k+1))
	                   +c2*(val(i+1,j+1,k  ) - val(i-1,j+1,k  ))
	                   +c3*(val(i+1,j+1,k+1) - val(i-1,j+1,k+1))
	                   +c4*(val(i+2,j,  k  ) - val(i,  j,  k  ))
	                   +c5*(val(i+2,j,  k+1) - val(i,  j,  k+1))
	                   +c6*(val(i+2,j+1,k  ) - val(i,  j+1,k  ))
	                   +c7*(val(i+2,j+1,k+1) - val(i,  j+1,k+1));
	   grad[1]=c0*(val(i,  j+1,k  ) - val(i,  j-1,k  ))
	                   +c1*(val(i,  j+1,k+1) - val(i,  j-1,k+1))
	                   +c2*(val(i,  j+2,k  ) - val(i,  j,  k  ))
	                   +c3*(val(i,  j+2,k+1) - val(i,  j,  k+1))
	                   +c4*(val(i+1,j+1,k  ) - val(i+1,j-1,k  ))
	                   +c5*(val(i+1,j+1,k+1) - val(i+1,j-1,k+1))
	                   +c6*(val(i+1,j+2,k  ) - val(i+1,j,  k  ))
	                   +c7*(val(i+1,j+2,k+1) - val(i+1,j,  k+1));
	   grad[2]=c0*(val(i,  j,  k+1) - val(i,  j,  k-1))
	                   +c1*(val(i,  j,  k+2) - val(i,  j,  k  ))
	                   +c2*(val(i,  j+1,k+1) - val(i,  j+1,k-1))
	                   +c3*(val(i,  j+1,k+2) - val(i,  j+1,k  ))
	                   +c4*(val(i+1,j,  k+1) - val(i+1,j,  k-1))
	                   +c5*(val(i+1,j,  k+2) - val(i+1,j,  k  ))
	                   +c6*(val(i+1,j+1,k+1) - val(i+1,j+1,k-1))
	                   +c7*(val(i+1,j+1,k+2) - val(i+1,j+1,k  ));
		}
	}
	void isq(RtInterval sq, RtInterval x){
		if(x[0]>=0){
			sq[0]=x[0]*x[0];
			sq[1]=x[1]*x[1];
		}
		else if(x[1]<=0){
			sq[0]=x[1]*x[1];
			sq[1]=x[0]*x[0];
		}
		else{
			sq[0]=0;
			sq[1]=-x[0]>x[1]?x[0]*x[0]:x[1]+x[1];
		}
	}
	void imax(RtInterval max, RtInterval a, RtInterval b){
		max[0]=b[0]>a[0]?b[0]:a[0];
		max[1]=b[1]>a[1]?b[1]:a[1];
	}
	void igeoff(RtInterval g, RtInterval r2){
		g[0]=geoff(r2[1]);
		g[1]=geoff(r2[0]);
	}
	/*void PrmanGridImplicit::Range(RtInterval val, const RtPoint corners[8],
			RtVolumeHandle h){
		RtInterval r, x, y, z, xsq, ysq, zsq, maxxy, maxxyz;
		int i;

		x[0]=x[1]=corners[0][0];
		y[0]=y[1]=corners[0][0];
		z[0]=z[1]=corners[0][0];
		for(i=0;i!=8;i++){
			if(corners[i][0]<x[0]) x[0]=corners[i][0];
			if(corners[i][0]>x[1]) x[1]=corners[i][0];
			if(corners[i][1]<y[0]) y[0]=corners[i][1];
			if(corners[i][1]>y[1]) y[1]=corners[i][1];
			if(corners[i][2]<z[0]) z[0]=corners[i][2];
			if(corners[i][2]>z[1]) z[1]=corners[i][2];
		}
		isq(xsq, x);
		isq(ysq, y);
		isq(zsq, z);
		imax(maxxy, xsq, ysq);
		imax(maxxyz, maxxy, zsq);
		igeoff(val, maxxyz);
	}*/
	PrmanGridImplicit::~PrmanGridImplicit(){
		//delete voxel;
		delete [] phi;
	}
	FIELDCREATE{
		return new PrmanGridImplicit(string[0]);
	}

