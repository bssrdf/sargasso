#include "particle.h"
#include "kdtree.h"
#include "rlelevelsets.h"
#include "CIsoSurface.h"
//#include "reinit_cuda.h"

#define POS(i,j,k) (k)*dim[0]*dim[1]+(j)*dim[0]+(i)
#define POS1(i,j,k) (k)*nx*ny+(j)*nx+(i)
#define INDEX(i,j,k) (k)*DimX*DimY+(j)*DimX+(i)

static int dim[3];
static float delta;
static float *phi = NULL;
static Point p0, p1;

struct FloatPair{
	FloatPair(){
		nFoundParticles = 0;
		smallestR = INFINITY;
		value = INFINITY;
	}
	mutable vector<WaterParticle> foundParticles;
	mutable u_int nFoundParticles;
	mutable float smallestR, value;
	void operator()(const WaterParticle &p,
					 float distSquare,
					 float &maxDistSquared) const{
		foundParticles.push_back(p);
		nFoundParticles++;
		float dist = sqrt(distSquare);
		if( dist <= smallestR){
			smallestR = dist;
			value = smallestR;
		}
		return;
	}
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

void ReadPhi(char *filename, float *&phi){
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

		   /*if(3!=fscanf(fp, "%f %f %f", &p0.x, &p0.y, &p0.z)){
			     Error("Problem reading gridimplicit file \"%s\" (dimensions)\n", filename);
			     p0.x=p0.y=p0.z=0;
			     return;
		   }*/
		   fread(&p0.x, sizeof(float), 1, fp);
		   fread(&p0.y, sizeof(float), 1, fp);
		   fread(&p0.z, sizeof(float), 1, fp);
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
		   		if(1!=fread(&phi[POS(i,j,k)], sizeof(float), 1, fp)){
				  printf("Problem reading gridimplicit file \"%s\" (at value %d)\n", filename, i);
				  dim[0]=dim[1]=dim[2]=0;
				  delete[] phi;
				  phi=0;
				  return;
			  }
		      assert(phi[i]==phi[i]);
		   }
		   fclose(fp);
}

void OutputMergedData(char *filename,
		int nx, int ny, int nz, float h, float *phi){

	FILE *fp = fopen(filename, "wb");
	if(!fp){
      printf("Couldn't open file [ %s ] to write \n", filename);
      exit(-1);
   	}
//   	fprintf(fp,"%d %d %d \n", DimX, DimY, DimZ);
   	fwrite(&nx, sizeof(int), 1, fp);
   	fwrite(&ny, sizeof(int), 1, fp);
   	fwrite(&nz, sizeof(int), 1, fp);
//   	fprintf(fp,"%f %f %f \n",p0.x, p0.y, p0.z);
   	fwrite(&(p0.x), sizeof(float), 1, fp);
   	fwrite(&(p0.y), sizeof(float), 1, fp);
   	fwrite(&(p0.z), sizeof(float), 1, fp);
//    fprintf(fp,"%f %f %f \n",p1.x, p1.y, p1.z);
    fwrite(&(p1.x), sizeof(float), 1, fp);
    fwrite(&(p1.y), sizeof(float), 1, fp);
    fwrite(&(p1.z), sizeof(float), 1, fp);
     fwrite(&h, sizeof(float), 1, fp);
//    fprintf(fp,"%f \n",voxel.VoxelDelta());

    for(int k=0;k<nz;k++)
		for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++)
				fwrite(&phi[POS1(i,j,k)], sizeof(float), 1, fp);
//				fprintf(fp,"%f\n", phi[k*(DimX*DimY)+j*DimX+i]);
	fclose(fp);

}

static float Kernel(const Point &p0, const Point &p1, float R){
	float dist = (p1 - p0).Length() / R;
	float x = 1.f - dist * dist;
	float w = max(0.f, x*x*x);
	return w;
}

static float Kernel(const Point &p0, const Point &p1, const Vector &vel, float max_vel, float R){
	Vector veln = Normalize(vel);
	float dpara = AbsDot((p1-p0),  veln);
	Vector v2 = (p1 - p0) - dpara * veln;
	float mag_vel = vel.Length();
	float dperp = v2.Length();
	float t = mag_vel/(2*max_vel);
	float dist = ((1.f - t) * dpara + (1.f + t) * dperp ) / R;
	float x = 1.f - dist * dist;
	float w = max(0.f, x*x*x);
	return w;
}

static float Kernel(const Point &p0, const Point &p1, float rad, const Vector &vel, float max_vel){
	Vector veln = Normalize(vel);
	float dpara = AbsDot((p1-p0),  veln);
	Vector v2 = (p1 - p0) - dpara * veln;
	float mag_vel = vel.Length();
	float dperp = v2.Length();
	float t = mag_vel/(2*max_vel);
	float dist = ((1.f - t) * dpara + (1.f + t) * dperp ) / rad;
	float w = 1.f / (1.f + dist * dist);
	return w;
}

static float TriInterp(float delta, int  DimX, int DimY, int DimZ,
			           const Point &p, const Point &P0, const float *phi) {

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
	Point point = Point(P0.x + delta*xi, P0.y + delta*yi, P0.z + delta*zi);
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

//    printf("\npoint at (%f,%f,%f) \n", p.x, p.y, p.z);
//    printf(" u = %f, v = %f, w = %f\n", u, v, w);
//	printf("xi = %d, yi = %d, zi = %d \n\n", xi, yi, zi);
//	return trilerp(phi[zi*(DimX*DimY)+yi*DimX+xi], phi[zi*(DimX*DimY)+yi*DimX+xi+1],
//				   phi[(zi+1)*(DimX*DimY)+yi*DimX+xi], phi[(zi+1)*(DimX*DimY)+yi*DimX+xi+1],
//				   phi[zi*(DimX*DimY)+(yi+1)*DimX+xi], phi[zi*(DimX*DimY)+(yi+1)*DimX+xi+1],
//				   phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi], phi[(zi+1)*(DimX*DimY)+(yi+1)*DimX+xi+1],
//				   u, v, w);
	return trilerp(phi[INDEX(xi, yi, zi)], phi[INDEX(xj, yi, zi)],
				   phi[INDEX(xi, yi, zj)], phi[INDEX(xj, yi, zj)],
				   phi[INDEX(xi, yj, zi)], phi[INDEX(xj, yj, zi)],
				   phi[INDEX(xi, yj, zj)], phi[INDEX(xj, yj, zj)],
				   r, s, t);

}


int main(int argc, char ** argv){
  int frame0 = atoi(argv[1]);
  int frame1 = atoi(argv[2]);
  int refine = atoi(argv[3]);
  float kerRad = atof(argv[4]);
  float parRad = atof(argv[5]);
  int num;



  int nx, ny, nz;

  for(int m=frame0; m<=frame1; ++m){
	  list<WaterParticle> waterParticles;
	  char fr[6];
	  sprintf(fr, ".f%d", m);
	  char filename1[20] = "phi";
	  strcat(filename1, fr);
//	  char filename1[20]= "cylinder_container";
	  ReadPhi(filename1, phi);

	  char filename[20] = "particles";
	  strcat(filename, fr);
	  FILE *fp=fopen(filename, "rb");
	   if(!fp){
			 printf("Couldn't open gridimplicit file \"%s\" for reading\n", filename);
			   return -1;
	   }
	   printf("******************************** \n");
	   printf("Processing [%s] \n", filename);
	   printf("Processing [%s] \n", filename1);
	   printf("******************************** \n");
	   fread(&p0.x, sizeof(float), 1, fp);
	   fread(&p0.y, sizeof(float), 1, fp);
	   fread(&p0.z, sizeof(float), 1, fp);
	   fread(&p1.x, sizeof(float), 1, fp);
	   fread(&p1.y, sizeof(float), 1, fp);
	   fread(&p1.z, sizeof(float), 1, fp);
	   fread(&delta, sizeof(float), 1, fp);
//	   printf(" delta = %f \n ", delta);
	   fread(&num, sizeof(int), 1, fp);
	  printf("there are %d particles \n", num);
	  float x, y, z, u, v, w, r;
	  for (int i=0; i<num; i++){
		   fread(&x, sizeof(float), 1, fp);
		   fread(&y, sizeof(float), 1, fp);
		   fread(&z, sizeof(float), 1, fp);
		   fread(&u, sizeof(float), 1, fp);
		   fread(&v, sizeof(float), 1, fp);
		   fread(&w, sizeof(float), 1, fp);
		   fread(&r, sizeof(float), 1, fp);
//	       WaterParticle wps(Point(x,y,z), r, -1, 0, 1, 1, 1);
		   WaterParticle  wps(Point(x,y,z), r, -1, delta, u, v, w);
//         printf("particle at (%f, %f, %f) with r = %f \n", x,y,z,r);
		   waterParticles.push_back(wps);
	  }
	  fclose(fp);

	  float h = delta / refine;
	  float minRad = 0.1f * delta;

	  printf(" h = %f, delta = %f \n ", h, delta);
	  printf("p0 = (%f, %f, %f) p1 = (%f, %f, %f) \n",
			  p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
	  nx = int(round((p1-p0).x / h));
	  ny = int(round((p1-p0).y / h));
	  nz = int(round((p1-p0).z / h));
	  float fnx = round((p1-p0).x / h);
	  float fny = round((p1-p0).y / h);
	  float fnz = round((p1-p0).z / h);
	  printf("nx = %d, ny = %d, nz = %d \n", nx, ny, nz);
	  printf("fnx = %f, fny = %f, fnz = %f \n", fnx, fny, fnz);

	  float threshold = (8.f / 9) * (8.f / 9) * (8.f / 9);

	  RLELevelSet rle(nx, ny, nz,
			         p0.x, p0.y, p0.z,
			         p1.x, p1.y, p1.z,
			         h, 3*h);

	  list<WaterParticle>::iterator iter_particle;
	  	vector<WaterParticle> wp;
	  	float max_vel = 0.f;
	  	float min_vel = INFINITY;
	  	for(iter_particle = waterParticles.begin();
	  		iter_particle != waterParticles.end();
	  		){
	  			WaterParticle &wps = *iter_particle;
	  			Vector vel = wps.Velocity();
	  			float mag_vel = vel.LengthSquared();
	  			if(mag_vel == 0.f)
	  				iter_particle = waterParticles.erase(iter_particle);
	  			else{
	  				if(max_vel < mag_vel)
	  					max_vel = mag_vel;
	  				if(min_vel > mag_vel)
	  					min_vel = mag_vel;
	  				wp.push_back(wps);
	  				++iter_particle;
	  			}
	  	}
//	  	int T = 0;
//	  	for(iter_particle = waterParticles.begin();
//	  		iter_particle != waterParticles.end();
//	  		){
//	  			WaterParticle &wps = *iter_particle;
//				wp.push_back(wps);
//				++iter_particle;
//				++T;
//				if(T > 3)
//					break;
//	  	}
	  	max_vel = sqrt(max_vel);
	  	min_vel = sqrt(min_vel);
	  	printf("ratio max_vel/min_vel = %f \n", max_vel/min_vel);
//	  	printf("got here (1.2)\n");
//	  	for(int n=0; n< wp.size(); ++n){
//	  		WaterParticle &wps = wp[n];
//	  		float r = wps.Radius();
//	  		printf("wps at (%f, %f, %f) with r = %f \n", wps.p.x, wps.p.y, wps.p.z, r);
//	  	}
	  //	KdTree<WaterParticle, FloatPair> dist(waterParticles);
	  	if(wp.size() == 0){
	  		for(int k=0;k<nz;k++){
	  			for(int j=0;j<ny;j++){
	  				for(int i=0;i<nx;i++){
	  					rle.FillRunData(i, j, k, 4*h);
	  				}
	  			}
	  		}
	  	}
	  	else{
//	  	 	printf("got here (1.2.5) wp.size = %u\n", wp.size());
			KdTree<WaterParticle, FloatPair> dist(wp);
//		  	printf("got here (1.3)\n");
			FloatPair fp;
			float R = kerRad * h;
			float radiusSquared = R * R;
//		  	printf("got here (2)\n");

//			printf("h = %f\n", h);

			 for(int k=0;k<nz;k++){
					for(int j=0;j<ny;j++){
						for(int i=0;i<nx;i++){
		  //	if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
				Point p(p0.x + h*(i+0.5), p0.y + h*(j+0.5),	p0.z + h*(k+0.5));
				dist.Lookup(p, fp, radiusSquared);
				float totalWeight = 0.f;
				float *weight = NULL;
				u_int N = 0;
//		  	    printf("got here (3) %d, %d, %d\n", i,j,k);

				if(fp.nFoundParticles > 0){
					N = fp.Size();
					weight = new float[N];
					for(int n=0; n < N; ++n){
						WaterParticle &wps = fp.foundParticles[n];
						Point pos = wps.Position();
//						Vector vel = wps.Velocity();
						weight[n] = Kernel(pos, p, R);
//						weight[n] = Kernel(pos, p, vel, max_vel, R);
						totalWeight += weight[n];
					}
				}
//		  	    printf("got here (4) %d, %d, %d\n", i,j,k);
//#if 0
				float cdata;
				if(totalWeight != 0.f){
					for(int n=0; n < N; ++n)
						weight[n] /= totalWeight;
					Point posAve(0.f, 0.f, 0.f);
					float rAve = 0.f;
					for(int n=0; n < N; ++n){
						WaterParticle &wps = fp.foundParticles[n];
						Point pos = wps.Position();
						float rad = wps.Radius();
						float pari = rad / minRad;
						rad = pari * parRad * h; // increase particle radii to be at least
														 // 1/2 the particle spacing
						posAve += weight[n]*pos;
						rAve   += weight[n]*rad;
					}
					cdata = (p - posAve).Length() - rAve;
				}
//#endif
#if 0
				float cdata;
				if(fp.nFoundParticles > 0){
					N = fp.Size();
					weight = new float[N];
					cdata = 0.f;
					for(int n=0; n < N; ++n){
						WaterParticle &wps = fp.foundParticles[n];
						Point pos = wps.Position();
						float rad = wps.Radius();
						Vector vel = wps.Velocity();
						weight[n] = Kernel(pos, p, rad, vel, max_vel);
						cdata += weight[n];
					}
					cdata -= 2.f;
				}
#endif
				else{
					cdata = 4 * h;
				}
//#endif
//				float cdata = threshold - totalWeight;

				// blend with main levelset
//				float fdata = TriInterp(delta, dim[0], dim[1], dim[2], p, p0, phi);
//				cdata = min(cdata, fdata);

//				cdata = phi[POS(i,j,k)];

				rle.FillRunData(i, j, k, cdata);
//		  	    printf("got here (5) %d, %d, %d\n", i,j,k);
				if(weight != NULL)
					delete [] weight;
				fp.Clear();
		  //	    printf("got here (6) %d, %d, %d\n", i,j,k);
		  //		printf(" at (%d, %d, %d) phi = %f\n", i,j,k,phi[INDEX(i,j,k)]);
		  //	}
		  //	else
		  //		phi[INDEX(i,j,k)] = threshold;
			}
				}
			 }
	  	}

	  	rle.Reinitialize();

	  	RLELevelSet rlenew(nx, ny, nz,
					 p0.x, p0.y, p0.z,
					 p1.x, p1.y, p1.z,
					 h, 3*h);

//	  	int I = 121, J = 165, K = 103;
	  	for(int k=0;k<nz;k++){
			for(int j=0;j<ny;j++){
				for(int i=0;i<nx;i++){
					Point p(p0.x + h*(i+0.5), p0.y + h*(j+0.5),	p0.z + h*(k+0.5));
					float fdata = TriInterp(delta, dim[0], dim[1], dim[2], p, p0, phi);
					float cdata = rle.Query(i,j,k);
					cdata = min(cdata, fdata);
					rlenew.FillRunData(i, j, k, cdata);
//					rlenew.FillRunData(i, j, k, fdata);
					float updata = rle.Query(i,j,k);
//					if(fdata == 0.f){
//					if(i == I && j == J && k == K){
//						printf("at (%d, %d, %d), fdata = %f, cdata = %f updata = %f\n",
//								i,j,k, fdata, cdata, updata);
//					}
				}
			}
	  	}

	  //printf(" retrieve data = %f \n", rlenew.Query(I,J,K));

	  char mphiname[20] = "spray";
	  strcat(mphiname, fr);
	//  OutputMergedData(voxel, mphiname,phi);
//	  rle.Reinitialize();
	  rlenew.OutputBinaryGridData(mphiname);

//	  rle.ConstructMesh();
//	  rle.OutputGridDataRib(mphiname);
//	  rle.DestroyMesh();
	  char frrib[11];
	  sprintf(frrib, ".f%d.obj", m);
//	  sprintf(frrib, ".f%d.rib", m);
	  char ribname[20] = "spray";
	  strcat(ribname, frrib);
	  CIsoSurface<float> *mc = new CIsoSurface<float>();
	  mc->GenerateSurface(&rlenew, 0.f, nx-1, ny-1, nz-1, h, h, h);
//	  mc->GenerateSurface(&rle, 0.f, nx-1, ny-1, nz-1, h, h, h);
//	  mc->RIBDump(ribname);
	  mc->OBJDump(ribname);
	  mc->DeleteSurface();
	  delete mc;

	//  object.OutputBinaryGridData(voxel, objectfile);
	  if(phi!=NULL)
		  delete [] phi;
	 }

   return 0;

}


/*int main(int argc, char ** argv){



  int frame0 = atoi(argv[1]);
  int frame1 = atoi(argv[2]);
  int num;



  printf("Here \n");

//  object.EvaluatePhi(voxel, WALL_THICKNESS+1, WALL_THICKNESS+1, WALL_THICKNESS+1,
//  				dim[0]-(WALL_THICKNESS+1)-1,  dim[1]-(WALL_THICKNESS+1)-1, dim[2]-(WALL_THICKNESS+1)-1);
//  printf("Object Levelset Initialization and Fast Marching Begins! \n");
//  object.InitializeObject(voxel, band, negband);
//  object.FastMarching(band, negband, voxel, true);
  printf("Object Levelset Initialization and  Fast Marching Finished! \n");
  int nx, ny, nz;

  for(int n=frame0; n<=frame1; ++n){
	  char fr[6];
	  char phiname[20] = "phi";
	  sprintf(fr, ".f%d", n);
	  strcat(phiname, fr);
	  printf("reading file name %s\n",phiname);
	  ReadPhi(phiname, phi);

	  char filename[20] = "particles";
	  strcat(filename, fr);
	  float h = delta / 4;
	  Voxel voxel(p0, p1, h);
	  voxel.GetDimensions(nx, ny, nz);
	  BlendParticles water(nx, ny, nz, 2*h);
	  FILE *fp=fopen(filename, "rb");
	   if(!fp){
	         printf("Couldn't open gridimplicit file \"%s\" for reading\n", filename);
	           return;
	   }
	   fread(&p0.x, sizeof(float), 1, fp);
	   fread(&p0.y, sizeof(float), 1, fp);
	   fread(&p0.z, sizeof(float), 1, fp);
	   fread(&p1.x, sizeof(float), 1, fp);
	   fread(&p1.y, sizeof(float), 1, fp);
	   fread(&p1.z, sizeof(float), 1, fp);
	   fwrite(&delta, sizeof(float), 1, fp);
	   fread(&num, sizeof(int), 1, fp);
	  printf("there are %d particles \n", num);
	  float x, y, z, u, v, w, r;
	  for (int i=0; i<num; i++){
	   fread(&x, sizeof(float), 1, fp);
	   fread(&y, sizeof(float), 1, fp);
	   fread(&z, sizeof(float), 1, fp);
	   fread(&u, sizeof(float), 1, fp);
   	   fread(&v, sizeof(float), 1, fp);
   	   fread(&w, sizeof(float), 1, fp);
	   fread(&r, sizeof(float), 1, fp);
//	   WaterParticle wps(Point(x,y,z), r, -1, 0, 1, 1, 1);
	   WaterParticle  wps(Point(x,y,z), r, -1, delta, u, v, w);
	//   printf("particle at (%f, %f, %f) with r = %f \n", x,y,z,r);
	   water.AddWaterParticles(wps);
	  }
	  fclose(fp);

	  water.EvaluatePhi(voxel);

	//  char fspray[18]= "spray";
	//  strcat(fspray, fr);
	//  water.OutputBinaryGridData(voxel, fspray);
	//  printf("before Phi merge \n");
	//  water.MergeWith(voxel, phi);
	//  printf("Phi has merged \n");
	  RLELevelSet(nx, ny, nz, h, 2*h, water.phi);


	  char mphiname[20] = "spray";
	  strcat(mphiname, fr);
	//  OutputMergedData(voxel, mphiname,phi);
	  OutputMergedData(voxel, mphiname,
	  		nx, ny, nz, h, water.phi);

	//  object.OutputBinaryGridData(voxel, objectfile);

	  if(phi!=NULL)
	  	  delete [] phi;
	 }

}*/
