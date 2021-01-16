/*
 * statictarget.cpp
 *
 *  Created on: Mar 6, 2009
 *      Author: bzhao
 */

#include "statictarget.h"
#include "pcg_solver.h"

void StaticTargetShape::ReadinPhi(const char *file){

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

	  fclose(fp);
}

void StaticTargetShape::VelocityFeedbackForce(const Voxel &voxel, const float *u, const float *v, const float *w,
				  float *u0, float *v0, float *w0) const{
	FOR_EACH_CELL
		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) ){
			u0[INDEX(i,j,k)] += -beta * u[INDEX(i,j,k)];
			v0[INDEX(i,j,k)] += -beta * v[INDEX(i,j,k)];
			w0[INDEX(i,j,k)] += -beta * w[INDEX(i,j,k)];
		}
		else if((voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
			u0[INDEX(i,j,k)] += -beta * u[INDEX(i,j,k)];
		else if((voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
			v0[INDEX(i,j,k)] += -beta * v[INDEX(i,j,k)];
		else if((voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
			w0[INDEX(i,j,k)] += -beta * w[INDEX(i,j,k)];
	END_FOR
}

void StaticTargetShape::GeometricPotentialForce(const Voxel &voxel, float *u0, float *v0, float *w0) const{
	  float h = voxel.VoxelDelta();
	  float *pot = new float[DimX*DimY*DimZ];
	  memset(pot, 0, DimX*DimY*DimZ*sizeof(float));

	  FOR_EACH_CELL
		  u_int index = INDEX(i,j,k);
		  if(gamma == 2)
			  pot[index] = C * (phi[index] > 0 ? 1 : -1) * fabsf(phi[index]) * fabsf(phi[index]);
		  if(gamma == 1)
			  pot[index] = C * (phi[index] > 0 ? 1 : -1) * fabsf(phi[index]);
		  if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) ){
				u0[index] += (pot[INDEX(i+1,j,k)] - pot[index]) / h;
				v0[index] += (pot[INDEX(i,j+1,k)] - pot[index]) / h;
				w0[index] += (pot[INDEX(i,j,k+1)] - pot[index]) / h;
		  }
		  else if((voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
			    u0[index] += (pot[INDEX(i+1,j,k)] - pot[index]) / h;
		  else if((voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
			    v0[index] += (pot[INDEX(i,j+1,k)] - pot[index]) / h;
		  else if((voxel.InSolid(i,j,k) || voxel.InAir(i,j,k) || voxel.InSurface(i,j,k)) && voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
				w0[index] += (pot[INDEX(i,j,k+1)] - pot[index]) / h;
	  END_FOR

	  delete [] pot;
}

void StaticTargetShape::ShapeFeedbackForce(const Voxel &voxel, const float *phic, float *u0, float *v0, float *w0) const{

	float *fx = new float[DimX*DimY*DimZ];
	memset(fx, 0, DimX*DimY*DimZ*sizeof(float));
	float *fy = new float[DimX*DimY*DimZ];
	memset(fy, 0, DimX*DimY*DimZ*sizeof(float));
	float *fz = new float[DimX*DimY*DimZ];
	memset(fz, 0, DimX*DimY*DimZ*sizeof(float));
	float *H = new float[DimX*DimY*DimZ];
	memset(H, 0, DimX*DimY*DimZ*sizeof(float));
	float h = voxel.VoxelDelta();



	  float multiplier = 0.f;
	  u_int M = 0;
	  FOR_EACH_CELL
		  u_int index = INDEX(i,j,k);
		  float phiu = 0.5f * (phi[index] + phi[INDEX(i+1,j,k)]);
		  float phiv = 0.5f * (phi[index] + phi[INDEX(i,j+1,k)]);
		  float phiw = 0.5f * (phi[index] + phi[INDEX(i,j,k+1)]);
		  if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) &&
		    (voxel.InSurface(i+1,j,k) || voxel.InAir(i+1,j,k) || voxel.InSolid(i+1,j,k))){
			  if(phiu > E_EPSIL) // outside target shape
					fx[index] = -alpha * phi[index] * (phi[INDEX(i+1,j,k)] - phi[index]) / h;
			  else
					fx[index] = -alpha * phi[index] * (phic[INDEX(i+1,j,k)] - phic[index]) / h;
			  multiplier += fx[index];
			  ++M;
		  }
		  if((voxel.InSurface(i,j,k) || voxel.InAir(i,j,k) || voxel.InSolid(i,j,k)) &&
			  voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)){
			  if(phiu > E_EPSIL) // outside target shape
					fx[index] = -alpha * phi[index] * (phi[INDEX(i+1,j,k)] - phi[index]) / h;
			  else
					fx[index] = -alpha * phi[index] * (phic[INDEX(i+1,j,k)] - phic[index]) / h;
			  multiplier += fx[index];
			  ++M;
		  }
		  if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) &&
			(voxel.InSurface(i,j+1,k) || voxel.InAir(i,j+1,k) || voxel.InSolid(i,j+1,k))){
			  if(phiv > E_EPSIL) // outside target shape
					fy[index] = -alpha * phi[index] * (phi[INDEX(i,j+1,k)] - phi[index]) / h;
			  else
					fy[index] = -alpha * phi[index] * (phic[INDEX(i,j+1,k)] - phic[index]) / h;
			  multiplier += fy[index];
			  ++M;
		  }
		  if((voxel.InSurface(i,j,k) || voxel.InAir(i,j,k) || voxel.InSolid(i,j,k)) &&
			  voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)){
			  if(phiv > E_EPSIL) // outside target shape
					fy[index] = -alpha * phi[index] * (phi[INDEX(i,j+1,k)] - phi[index]) / h;
			  else
					fy[index] = -alpha * phi[index] * (phic[INDEX(i,j+1,k)] - phic[index]) / h;
			  multiplier += fy[index];
			  ++M;
		  }
		  if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) &&
			(voxel.InSurface(i,j,k+1) || voxel.InAir(i,j,k+1) || voxel.InSolid(i,j,k+1))){
			  if(phiw > E_EPSIL) // outside target shape
					fz[index] = -alpha * phi[index] * (phi[INDEX(i,j,k+1)] - phi[index]) / h;
			  else
					fz[index] = -alpha * phi[index] * (phic[INDEX(i,j,k+1)] - phic[index]) / h;
			  multiplier += fz[index];
			  ++M;
		  }
		  if((voxel.InSurface(i,j,k) || voxel.InAir(i,j,k) || voxel.InSolid(i,j,k)) &&
			  voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)){
			  if(phiw > E_EPSIL) // outside target shape
					fz[index] = -alpha * phi[index] * (phi[INDEX(i,j,k+1)] - phi[index]) / h;
			  else
					fz[index] = -alpha * phi[index] * (phic[INDEX(i,j,k+1)] - phic[index]) / h;
			  multiplier += fz[index];
			  ++M;
		  }
	  END_FOR

	  FOR_EACH_CELL
		  u_int index = INDEX(i,j,k);
		  fx[index] -= multiplier / M;
		  fy[index] -= multiplier / M;
		  fz[index] -= multiplier / M;
	  END_FOR

	  map<u_int, float> B;
	  map<u_int, u_int> A;
	  map<u_int, u_int>::iterator found;
	  u_int N = 0;

	  	// compute the r.h.s. of the pressure equation
	  FOR_EACH_CELL
		  float rhs = 0.f;
		   if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
			   if(!(voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)))
				   rhs += -fx[INDEX(i,j,k)];
			   if(!(voxel.InLiquid(i-1,j,k) && !voxel.InSurface(i-1,j,k)))
			   	   rhs += -fx[INDEX(i-1,j,k)];
			   if(!(voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)))
			   	   rhs += -fy[INDEX(i,j,k)];
			   if(!(voxel.InLiquid(i,j-1,k) && !voxel.InSurface(i,j-1,k)))
			   	   rhs += -fy[INDEX(i,j-1,k)];
			   if(!(voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)))
			   	   rhs += -fz[INDEX(i,j,k)];
			   if(!(voxel.InLiquid(i,j,k-1) && !voxel.InSurface(i,j,k-1)))
				   rhs += -fz[INDEX(i,j,k)];
			   rhs *= h;
			   B.insert( make_pair(INDEX(i,j,k), rhs) );
			   A.insert( make_pair(INDEX(i,j,k), N) );
			   ++N;
		   }
	  END_FOR
	  float *p = new float[N];
	  memset(p, 0, N*sizeof(float));
	  LinearSolver(voxel, p, A, B, N);
	  found = A.begin();
	  for(; found != A.end(); ++found){
	  		H[found->first] = p[found->second];
	  }
	  delete [] p;

	  FOR_EACH_CELL
	  		if(voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k)){
	  			if(!(voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k)))
	  				u0[INDEX(i,j,k)] += fx[INDEX(i,j,k)];
	  			else
	  				u0[INDEX(i,j,k)] += (H[INDEX(i+1,j,k)] - H[INDEX(i,j,k)]) / h;
	  			if(!(voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k)))
					v0[INDEX(i,j,k)] += fy[INDEX(i,j,k)];
				else
					v0[INDEX(i,j,k)] += (H[INDEX(i,j+1,k)] - H[INDEX(i,j,k)]) / h;
	  			if(!(voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1)))
					w0[INDEX(i,j,k)] += fz[INDEX(i,j,k)];
				else
					w0[INDEX(i,j,k)] += (H[INDEX(i,j,k+1)] - H[INDEX(i,j,k)]) / h;
	  		}
	  		else{
	  			if(voxel.InLiquid(i+1,j,k) && !voxel.InSurface(i+1,j,k))
					u0[INDEX(i,j,k)] += fx[INDEX(i,j,k)];
				if(voxel.InLiquid(i,j+1,k) && !voxel.InSurface(i,j+1,k))
					v0[INDEX(i,j,k)] += fy[INDEX(i,j,k)];
				if(voxel.InLiquid(i,j,k+1) && !voxel.InSurface(i,j,k+1))
					w0[INDEX(i,j,k)] += fz[INDEX(i,j,k)];
	  		}
	  	END_FOR

	delete [] H;
	delete [] fx;
	delete [] fy;
	delete [] fz;
}

void StaticTargetShape::LinearSolver( const Voxel &voxel,
					 float * x, map<u_int, u_int> &a,
					 map<u_int, float> &x0, u_int n ) const{

	float delta = voxel.VoxelDelta();
	float delta2 = delta * delta;
	float delta3 = delta * delta2;
	map<u_int, u_int>::iterator found;
	map<u_int, float>::iterator found1;
	PCGSolver<double> solver;
	solver.set_solver_parameters(1e-13, 100, 0.97, 0.25);
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
		 if( voxel.InLiquid(i,j,k) && !voxel.InSurface(i,j,k) ){
			 found = a.find(index);
			 diag = found->second;
			 found1 = x0.find(index);
			 rhs.push_back((double)found1->second);
			 TentativeGridPoint tp(0.f,i,j,k);
			 vector<TentativeGridPoint> neighbors;
			 tp.Neigbor(neighbors, DimX, DimY, DimZ);
			 int solidNeighbors=0;
			 float total_theta = 0.f;
			 float coeff[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
			 for(u_int m=0; m<neighbors.size(); m++){
					 if(voxel.InSolid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)){
						 ++solidNeighbors;
					 }
					 else if(voxel.InAir(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ||
							 voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) ){
						 ++solidNeighbors;
					 }

					 else if(voxel.InLiquid(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk) &&
							!voxel.InSurface(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk)
							){
						 found = a.find(INDEX(neighbors[m].ii,neighbors[m].jj, neighbors[m].kk));
//						 A.set_element(diag,found->second,-delta3);
						 A.set_element(diag,found->second,-1);
					 }
			 }
			 if(solidNeighbors == 6){
				printf("at (%d, %d, %d) Newmann b.c. for all directions \n", i,j,k);
				exit(1);
			 }
//			 printf("at (%d, %d, %d) total theta = %f \n", i,j,k,total_theta);

//             A.set_element(diag, diag, (double)diag_total);
			 A.set_element(diag,diag,(float)(neighbors.size()-solidNeighbors));
//			 printf("non solid neighbors= %d at(%d, %d,%d) \n",
//					 (neighbors.size()-solidNeighbors), i,j,k);
		 }


	END_FOR
	//A.write_matlab(cout,"A");
//	printf("\nmatrix A is ready \n\n ");

	solver.solve(A, rhs, p, residual, iterations);
	for(u_int i=0; i<n; ++i)
			x[i] = (float)p[i];
	printf("\nPCG iterations = %d, and residual = %16.13f \n\n ",
			iterations, residual);


}

