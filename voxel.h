#ifndef VOXEL_H_
#define VOXEL_H_

#include "geometry.h"
#include "particle.h"
#include "celltype.h"
//#include "scalarfield3d.h"
//class ScalarField3D;

#ifdef SPMDMPI
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"
extern int num_ghosts;
extern int imt, jmt, kmt;
extern int coords[2];
extern int nproc[2];
extern int pid;
extern int iphys_b, iphys_e;
extern int jphys_b, jphys_e;
extern bool comm_both;
extern char comm_dir;
extern int global_imt, global_jmt, global_kmt;

#define INDEX(i,j,k) (k)*imt*jmt+(j)*imt+(i)
#define INDEX_GLOBAL(i,j,k) (k)*global_imt*global_jmt+(j)*global_imt+(i)

#define FOR_EACH_CELL for(int k=0;k<DimZ;k++) { \
						 for(int j=jphys_b;j<=jphys_e;j++) { \
							for(int i=iphys_b;i<=iphys_e;i++){
#define END_FOR }}}

#else
#define INDEX(i,j,k) (k)*DimX*DimY+(j)*DimX+(i)

#define FOR_EACH_CELL for(int k=0;k<DimZ;k++) { \
							for(int j=0;j<DimY;j++) { \
								for(int i=0;i<DimX;i++){
#define END_FOR } \
				} \
	            }
#endif

#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}

#define WALL_THICKNESS 4

#define E_EPSIL  1.e-6f

class Container;

class Voxel{
public:
	Voxel(const Point &p0, const Point &p1, float delta )
#ifdef SPMDMPI
	{
		if(comm_both){
			P0.x = p0.x + coords[0]*(p1.x - p0.x)/nproc[0];
			P0.y = p0.y + coords[1]*(p1.y - p0.y)/nproc[1];
			P1.x = p0.x + (coords[0]+1)*(p1.x - p0.x)/nproc[0];
			P1.y = p0.y + (coords[1]+1)*(p1.y - p0.y)/nproc[1];
		}
		else{
			if(comm_dir == 1){
				P0.x = p0.x + coords[0]*(p1.x - p0.x)/nproc[0];
				P0.y = p0.y;
				P1.x = p0.x + (coords[0]+1)*(p1.x - p0.x)/nproc[0];
				P1.y = p1.y;
			}
			else{
				P0.x = p0.x;
				P0.y = p0.y + coords[1]*(p1.y - p0.y)/nproc[1];
				P1.x = p1.x;
				P1.y = p0.y + (coords[1]+1)*(p1.y - p0.y)/nproc[1];
			}
		}
		P0.z = p0.z;
		P1.z = p1.z;
		h = delta;
		DimX = imt;
		DimY = jmt;
		DimZ = kmt;
		if(!pid){
			printf(" the domain on this processor [%d] is p0(%f, %f, %f) to p1(%f, %f, %f) \n",
				pid, P0.x, P0.y, P0.z, P1.x, P1.y, P1.z);
			fflush(stdout);
		}
		cellType = new CellType(DimX*DimY*DimZ);
		cellFastMarchType = new CellType(DimX*DimY*DimZ);
	}
#else
	{
		P0 = p0;
		P1 = p1;
		h = delta;
		DimX = int(round((p1-p0).x / h));
		DimY = int(round((p1-p0).y / h));
		DimZ = int(round((p1-p0).z / h));
		printf("float dimx = %f, dimy = %f dimz= %f \n",
				(p1-p0).x / h , (p1-p0).y / h, (p1-p0).z / h);
		cellType = new CellType(DimX*DimY*DimZ);
		cellFastMarchType = new CellType(DimX*DimY*DimZ);

	}
#endif
	~Voxel(){
		delete cellType;
		delete cellFastMarchType;
	}
	void PlaceParticles(int ii, int jj, int kk,
						 Particle &ps, int face) const;
	void PlacePoints(int ii, int jj, int kk,
							 Point &ps, int face) const;
	Point VoxelCenterPosition(int i, int j, int k) const;
	Point SubVoxelCenterPosition(int i, int j, int k, int level, int index) const;
	Point VoxelCornerPosition(int i, int j, int k, int index) const;
	float VoxelDelta() const { return h;}
	float VoxelMaxExtent() const{
		return Distance(P0, P1);
	}
	void Initialize(int x0, int y0, int z0,
				  	int x1, int y1, int z1){
		FOR_EACH_CELL

			if( i <= x1 && i >= x0 &&
				j <= y1 && j >= y0 &&
				k <= z1 && k >= z0 )
//				k == z0  )
//				UpdateCellType(i, j, k, SET, LIQUID);
				UpdateCellType(i, j, k, SET, EMPTY);
			else if( i <= WALL_THICKNESS          ||
					 j <= WALL_THICKNESS          ||
					 k <= WALL_THICKNESS          ||
//					 k < z0         ||
					 i >= DimX-(WALL_THICKNESS+1) ||
					 j >= DimY-(WALL_THICKNESS+1) ||
					 k >= DimZ-(WALL_THICKNESS+1) )
				UpdateCellType(i, j, k, SET, SOLID);
			else
				UpdateCellType(i, j, k, SET, EMPTY);

		END_FOR
	}

	void Initialize(const Container *c);

	bool InLiquid(int i,int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
			printf("Error in (InLiquid())! index out of range i=%d, j=%d, k=%d \n", i, j, k);
			exit(1);
		}
		return cellType->InLiquid((u_int)INDEX(i,j,k));
	}
	bool InLiquid(const TentativeGridPoint &tp) const{
		Error(tp,"InLquid()");
		return cellType->InLiquid((u_int)INDEX(tp.ii,tp.jj,tp.kk));
	}
	bool InSurface(int i, int j, int k) const{
		Error(i,j,k,"InSurface()");
		return cellType->InSurface((u_int)INDEX(i,j,k));
	}
	bool InSurface(const TentativeGridPoint &tp) const{
		Error(tp,"InSurface()");
		return cellType->InSurface((u_int)INDEX(tp.ii,tp.jj,tp.kk));
	}
	bool InSolid(int i, int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
			printf("Error in (InSolid())! index out of range i=%d, j=%d, k=%d \n", i, j, k);
			printf("Error in (InSolid())! index out of range DimX=%d, DimY=%d, DimZ=%d \n", DimX, DimY, DimZ);
//			return true;
			exit(1);
		}
		return cellType->InSolid((u_int)INDEX(i,j,k));
	}
	bool InSolid(const TentativeGridPoint &tp) const{
		if( tp.ii < 0 || tp.ii >= DimX || tp.jj < 0 || tp.jj >= DimY  || tp.kk < 0 || tp.kk >= DimZ){
			printf("Error in (InSolid())! index out of range i=%d, j=%d, k=%d \n", tp.ii, tp.jj, tp.kk);
			printf("Error in (InSolid())! index out of range i=%d, j=%d, k=%d \n", DimX, DimY, DimZ);
//			return true;
			exit(1);
		}
		return cellType->InSolid((u_int)INDEX(tp.ii,tp.jj,tp.kk));
	}
	bool InMovingSolid(int i, int j, int k) const{
		Error(i,j,k,"InMovingSolid()");
		return cellType->InMovingSolid((u_int)INDEX(i,j,k));
	}
	bool InMovingSolid(const TentativeGridPoint &tp) const{
		Error(tp,"InMovingSolid()");
		return cellType->InMovingSolid((u_int)INDEX(tp.ii,tp.jj,tp.kk));
	}
	bool InAir(int i, int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
			printf("Error in (InAir())! index out of range i=%d, j=%d, k=%d \n", i, j, k);
			exit(1);
		}
		return cellType->InAir((u_int)INDEX(i,j,k));
	}
	bool InAir(const TentativeGridPoint &tp) const{
		Error(tp,"InAir()");
		return cellType->InAir((u_int)INDEX(tp.ii,tp.jj,tp.kk));
	}
	bool InSource(int i, int j, int k) const{
		Error(i,j,k,"InSource()");
		return cellType->InSource((u_int)INDEX(i,j,k));
	}
	bool InSource(const TentativeGridPoint &tp) const{
		Error(tp,"InSource()");
		return cellType->InSource((u_int)INDEX(tp.ii,tp.jj,tp.kk));
	}
	bool IsDone(int i, int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
			printf("Error in (InAir())! index out of range i=%d, j=%d, k=%d \n", i, j, k);
			exit(1);
		}
		return cellFastMarchType->IsDone((u_int)INDEX(i,j,k));
	}
	bool IsDone(u_int index) const{
		return cellFastMarchType->IsDone(index);
	}
	bool IsDone(const TentativeGridPoint &p) const{
		if( p.ii < 0 || p.ii >= DimX || p.jj < 0 || p.jj >= DimY  || p.kk < 0 || p.kk >= DimZ){
			printf("Error in (IsDone())! index out of range i=%d, j=%d, k=%d \n", p.ii, p.jj, p.kk);
			exit(1);
		}
		return cellFastMarchType->IsDone((u_int)INDEX(p.ii, p.jj, p.kk));
	}
	bool IsClose(int i, int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
			printf("Error in (InAir())! index out of range i=%d, j=%d, k=%d \n", i, j, k);
			exit(1);
		}
		return cellFastMarchType->IsClose((u_int)INDEX(i,j,k));
	}
	bool IsClose(const TentativeGridPoint &p) const{
		if( p.ii < 0 || p.ii >= DimX || p.jj < 0 || p.jj >= DimY  || p.kk < 0 || p.kk >= DimZ){
			printf("Error in (IsDone())! index out of range i=%d, j=%d, k=%d \n", p.ii, p.jj, p.kk);
			exit(1);
		}
		return cellFastMarchType->IsClose((u_int)INDEX(p.ii, p.jj, p.kk));
	}
	bool IsFar(int i, int j, int k) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
			printf("Error in (InAir())! index out of range i=%d, j=%d, k=%d \n", i, j, k);
			exit(1);
		}
		return cellFastMarchType->IsFar((u_int)INDEX(i,j,k));
	}
	bool IsFar(const TentativeGridPoint &p) const{
		if( p.ii < 0 || p.ii >= DimX || p.jj < 0 || p.jj >= DimY  || p.kk < 0 || p.kk >= DimZ){
			printf("Error in (IsDone())! index out of range i=%d, j=%d, k=%d \n", p.ii, p.jj, p.kk);
			exit(1);
		}
		return cellFastMarchType->IsFar((u_int)INDEX(p.ii, p.jj, p.kk));
	}
	void UpdateCellType(int i, int j, int k, bool set, char t){
		if(set)  // set a tag
			cellType->UpdateTag((u_int)INDEX(i,j,k),t);
		else     // clear a tag
			cellType->ClearTag((u_int)INDEX(i,j,k),t);
	}
	void UpdateCellFastMarchType(int i, int j, int k, bool set, char t){
		if(set)  // set a tag
			cellFastMarchType->UpdateTag((u_int)INDEX(i,j,k),t);
		else     // clear a tag
			cellFastMarchType->ClearTag((u_int)INDEX(i,j,k),t);
	}

	void Error(int i, int j, int k, const char *str) const{
		if( i < 0 || i >= DimX || j < 0 || j >= DimY  || k < 0 || k >= DimZ){
#ifdef SPMDMPI
			printf("pid = %d, Error in (%s)! index out of range i=%d, j=%d, k=%d \n", pid, str, i, j, k);
			fflush(stdout);
			MPI_Finalize();
#else
			printf("Error in (%s)! index out of range i=%d, j=%d, k=%d \n", i, j, k, str);
#endif
			exit(1);

		}
	}
	void Error(const TentativeGridPoint &tp, const char *str) const {
		if( tp.ii < 0 || tp.ii >= DimX || tp.jj < 0 || tp.jj >= DimY  || tp.kk < 0 || tp.kk >= DimZ){
#ifdef SPMDMPI
			printf("pid = %d, Error in (%s)! index out of range i=%d, j=%d, k=%d \n", pid, str, tp.ii, tp.jj, tp.kk);
			fflush(stdout);
			MPI_Finalize();
#else
			printf("Error in (%s)! index out of range i=%d, j=%d, k=%d \n", tp.ii, tp.jj, tp.kk, str);
#endif
			exit(1);
		}
	}
	void GetDimensions(int &X, int &Y, int &Z) const{
		X = DimX;
		Y = DimY;
		Z = DimZ;
		return;
	}

	TentativeGridPoint ContainsPoint( const Point &p ) const;
	TentativeGridPoint ContainsPointInVelCell( const Point &p, unsigned char vel_index ) const;

	float FacePostion(int axis, bool far, int i, int j, int k) const {
		if( axis == 1 ){
			if(far)
				return P0.x+(i+1)*h;
			else
				return P0.x+(i)*h;
		}
		else if( axis == 2 ){
			if(far)
				return P0.y+(j+1)*h;
			else
				return P0.y+(j)*h;
		}
		else if( axis == 3 ){
			if(far)
				return P0.z+(k+1)*h;
			else
				return P0.z+(k)*h;
		}
		else{
			fprintf(stderr, "Error! parameter [axis] must be in range [1...3]\n ");
			exit(1);
		}
	}

	Point SubFacePosition(int axis, int n, int i, int j, int k) const{
		if(axis == 0){
			if(n == 0)
				return Point(P0.x+(i)*h, P0.y+(j+0.25f)*h, P0.z+(k+0.25f)*h);
			else if(n == 1)
				return Point(P0.x+(i)*h, P0.y+(j+0.75f)*h, P0.z+(k+0.25f)*h);
			else if(n == 2)
				return Point(P0.x+(i)*h, P0.y+(j+0.75f)*h, P0.z+(k+0.75f)*h);
			else if(n == 3)
				return Point(P0.x+(i)*h, P0.y+(j+0.25f)*h, P0.z+(k+0.75f)*h);
		}
		else if(axis == 1){
			if(n == 0)
				return Point(P0.x+(i+1)*h, P0.y+(j+0.25f)*h, P0.z+(k+0.25f)*h);
			else if(n == 1)
				return Point(P0.x+(i+1)*h, P0.y+(j+0.75f)*h, P0.z+(k+0.25f)*h);
			else if(n == 2)
				return Point(P0.x+(i+1)*h, P0.y+(j+0.75f)*h, P0.z+(k+0.75f)*h);
			else if(n == 3)
				return Point(P0.x+(i+1)*h, P0.y+(j+0.25f)*h, P0.z+(k+0.75f)*h);

		}
		else if(axis == 2){
			if(n == 0)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j)*h, P0.z+(k+0.25f)*h);
			else if(n == 1)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j)*h, P0.z+(k+0.25f)*h);
			else if(n == 2)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j)*h, P0.z+(k+0.75f)*h);
			else if(n == 3)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j)*h, P0.z+(k+0.75f)*h);
		}
		else if(axis == 3){
			if(n == 0)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j+1)*h, P0.z+(k+0.25f)*h);
			else if(n == 1)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j+1)*h, P0.z+(k+0.25f)*h);
			else if(n == 2)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j+1)*h, P0.z+(k+0.75f)*h);
			else if(n == 3)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j+1)*h, P0.z+(k+0.75f)*h);
		}
		else if(axis == 4){
			if(n == 0)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j+0.25)*h, P0.z+(k)*h);
			else if(n == 1)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j+0.25)*h, P0.z+(k)*h);
			else if(n == 2)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j+0.75)*h, P0.z+(k)*h);
			else if(n == 3)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j+0.75)*h, P0.z+(k)*h);
		}
		else if(axis == 5){
			if(n == 0)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j+0.25)*h, P0.z+(k+1)*h);
			else if(n == 1)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j+0.25)*h, P0.z+(k+1)*h);
			else if(n == 2)
				return Point(P0.x+(i+0.75f)*h, P0.y+(j+0.75)*h, P0.z+(k+1)*h);
			else if(n == 3)
				return Point(P0.x+(i+0.25f)*h, P0.y+(j+0.75)*h, P0.z+(k+1)*h);
		}
		else{
				fprintf(stderr, "Error! parameter [axis] must be in range [0...5]\n ");
				exit(1);
			}
	}

	Point VelPosition(int axis, int i, int j, int k) const{
		if(axis == 1)
			return Point(P0.x+(i+1)*h, P0.y+(j+0.5f)*h, P0.z+(k+0.5f)*h);
		else if(axis == 2)
			return Point(P0.x+(i+0.5f)*h, P0.y+(j+1)*h, P0.z+(k+0.5f)*h);
		else if(axis == 3)
			return Point(P0.x+(i+0.5f)*h, P0.y+(j+0.5f)*h, P0.z+(k+1)*h);
		else{
				fprintf(stderr, "Error! parameter [axis] must be in range [1...3]\n ");
				exit(1);
			}
	}

	Point VelClosePosition(int m, int vel_index, int i, int j, int k) const{
		if(m == 0){
			if(vel_index == 1)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+0.5f)*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 2)
				return Point(P0.x+i*h, P0.y+(j+1.f)*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 3)
				return Point(P0.x+i*h, P0.y+(j+0.5f)*h, P0.z+(k+1.f)*h);
		}
		else if(m == 1){
			if(vel_index == 1)
				return Point(P0.x+(i+1.5f)*h, P0.y+(j+0.5f)*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 2)
				return Point(P0.x+(i+1)*h, P0.y+(j+1)*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 3)
				return Point(P0.x+(i+1)*h, P0.y+(j+0.5f)*h, P0.z+(k+1.f)*h);
		}
		else if(m == 2){
			if(vel_index == 1)
				return Point(P0.x+(i+1)*h, P0.y+j*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 2)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+0.5f)*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 3)
				return Point(P0.x+(i+0.5f)*h, P0.y+j*h, P0.z+(k+1.f)*h);
		}
		else if(m == 3){
			if(vel_index == 1)
				return Point(P0.x+(i+1)*h, P0.y+(j+1)*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 2)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+1.5f)*h, P0.z+(k+0.5f)*h);
			else if(vel_index == 3)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+1)*h, P0.z+(k+1)*h);
		}
		else if(m == 4){
			if(vel_index == 1)
				return Point(P0.x+(i+1)*h, P0.y+(j+0.5f)*h, P0.z+k*h);
			else if(vel_index == 2)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+1)*h, P0.z+k*h);
			else if(vel_index == 3)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+0.5f)*h, P0.z+(k+0.5f)*h);
		}
		else if(m == 5){
			if(vel_index == 1)
				return Point(P0.x+(i+1)*h, P0.y+(j+0.5f)*h, P0.z+(k+1)*h);
			else if(vel_index == 2)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+1)*h, P0.z+(k+1)*h);
			else if(vel_index == 3)
				return Point(P0.x+(i+0.5f)*h, P0.y+(j+0.5f)*h, P0.z+(k+1.5f)*h);
		}
		else{
			fprintf(stderr, "Error! parameter [m] must be in range [0...5]\n ");
			exit(1);
		}
	}

	bool CloseToSolid(const TentativeGridPoint &tp) const{
		vector<TentativeGridPoint> neighbors;
		tp.Neigbor(neighbors, DimX, DimY, DimZ);
		for(u_int m=0; m<neighbors.size(); m++){
			 if(InSolid(neighbors[m]))
				 return true;
		}
		return false;
	}
	
	bool CloseToSolidTwoRings(int xi, int yi, int zi) const{		
		for(int i=xi-1; i<=xi+1; ++i){
			for(int j=yi-1; j<=yi+1; ++j){
				for(int k=zi-1; k<=zi+1; ++k){
					if(InSolid(i, j, k))
						return true;
				}
			}
		}
		return false;
	}
	
	bool CloseToNonSolid(const TentativeGridPoint &tp, int &n) const{
		n = 0;
		vector<TentativeGridPoint> neighbors;
		tp.Neigbor(neighbors, DimX, DimY, DimZ);
		for(u_int m=0; m<neighbors.size(); m++){
			 if(!InSolid(neighbors[m]))
				 ++n;
		}
		if(n > 0)
			return true;
		else
			return false;
	}
	bool CloseToLiquid(const TentativeGridPoint &tp) const{
		int n = 0;
		vector<TentativeGridPoint> neighbors;
		tp.Neigbor(neighbors, DimX, DimY, DimZ);
		for(u_int m=0; m<neighbors.size(); m++){
				 if(InLiquid(neighbors[m]) && !InSurface(neighbors[m]))
						 ++n;
		}
		if(n > 0)
				return true;
		else
				return false;
	}

	bool CloseToAir(const TentativeGridPoint &tp) const{
		int n = 0;
		vector<TentativeGridPoint> neighbors;
		tp.Neigbor(neighbors, DimX, DimY, DimZ);
		for(u_int m=0; m<neighbors.size(); m++){
			 if(InAir(neighbors[m]) || InSurface(neighbors[m]))
				++n;
		}
		if(n > 0)
			return true;
		else
			return false;
	}

	TentativeGridPoint ClosestObjectCell(const TentativeGridPoint &tp) const{
		vector<TentativeGridPoint> neighbors;
		tp.Neigbor(neighbors, DimX, DimY, DimZ);
		for(u_int m=0; m<neighbors.size(); m++){
			if(InSolid(neighbors[m]))
				return neighbors[m];
		}
		printf("Error! not supposed to get here!");
		exit(1);
	}

	Point GetP0() const{
#ifdef SPMDMPI
		if(comm_both)
			return Point(P0.x-num_ghosts*h, P0.y-num_ghosts*h, P0.z);
		else{
			if(comm_dir == 1)
				return Point(P0.x-num_ghosts*h, P0.y, P0.z);
			else
				return Point(P0.x, P0.y-num_ghosts*h, P0.z);
		}
#else
		return P0;
#endif
	}

	Point GetP1() const{
#ifdef SPMDMPI
	if(comm_both)
		return Point(P1.x+num_ghosts*h, P1.y+num_ghosts*h, P1.z);
	else{
		if(comm_dir == 1)
			return Point(P1.x+num_ghosts*h, P1.y, P1.z);
		else
			return Point(P1.x, P1.y+num_ghosts*h, P1.z);
	}
#else
		return P1;
#endif
	}

#ifdef SPMDMPI
	void UpdateGhosts(){
		cellType->UpdateGhosts();
	}

	bool PosOutofEastBound(const Point &pos) const{
		if(pos.x > P1.x && (pos.y <= P1.y && pos.y >= P0.y))
			return true;
		else
			return false;
	}
	bool PosOutofWestBound(const Point &pos) const{
		if(pos.x < P0.x && (pos.y <= P1.y && pos.y >= P0.y))
			return true;
		else
			return false;
	}
	bool PosOutofNorthBound(const Point &pos) const{
		if(pos.y > P1.y && (pos.x <= P1.x && pos.x >= P0.x))
			return true;
		else
			return false;
	}
	bool PosOutofSouthBound(const Point &pos) const{
		if(pos.y < P0.y && (pos.x <= P1.x && pos.x >= P0.x))
			return true;
		else
			return false;
	}
	bool PosOutofNorthEastBound(const Point &pos) const{
		if(pos.x > P1.x && pos.y > P1.y)
			return true;
		else
			return false;
	}
	bool PosOutofNorthWestBound(const Point &pos) const{
		if(pos.x < P0.x && pos.y > P1.y)
			return true;
		else
			return false;
	}
	bool PosOutofSouthEastBound(const Point &pos) const{
		if(pos.x > P1.x && pos.y < P0.y)
			return true;
		else
			return false;
	}
	bool PosOutofSouthWestBound(const Point &pos) const{
		if(pos.x < P0.x && pos.y < P0.y)
			return true;
		else
			return false;
	}
#endif

private:
	Point P0, P1;
	//float DeltaX,DeltaY,DeltaZ;
	float h; // uniform grid size in all 3 dimensions
	int DimX;
	int DimY;
	int DimZ;
//	bool *SolidCell; // constant unless there is a moving object
//	bool *SurfcCell;
//	bool *EmptyCell;
	CellType *cellType;
	CellType *cellFastMarchType;
};

#endif /*VOXEL_H_*/
