#include <stdlib.h>
#include "mtrand.h"
#include "voxel.h"
#include "container.h"
// randomly place np particles in the voxel with index ii, jj, kk
// the position and radius of each particle are modified
void Voxel::PlaceParticles(int ii, int jj, int kk,
						 Particle &ps, int face) const{
	Point p0, p1;
#ifdef SPMDMPI
	if(comm_both){
		p0.x = P0.x + h*(ii-num_ghosts);
		p0.y = P0.y + h*(jj-num_ghosts);
	}
	else{
		if(comm_dir == 1){
			p0.x = P0.x + h*(ii-num_ghosts);
			p0.y = P0.y + h*(jj);
		}
		else{
			p0.x = P0.x + h*(ii);
			p0.y = P0.y + h*(jj-num_ghosts);
		}
	}
#else
	p0.x = P0.x + h*ii;
	p0.y = P0.y + h*jj;
#endif
	p0.z = P0.z + h*kk;
	MTRand_open mt;
	double rn1 = mt();
	double rn2 = mt();
	double rn3 = mt();
	double rn4 = mt();
//	printf(" %f, %f, %f,  %f \n", rn1, rn2, rn3, rn4);
	if(face == 1){
		p1.x = p0.x;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + (float)rn3 * h;
	}
	else if( face == 2 ){
		p1.x = p0.x + h;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + (float)rn3 * h;
	}
	else if( face == 3 ){
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y;
		p1.z = p0.z + (float)rn3 * h;
	}
	else if( face == 4 ){
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y + h;
		p1.z = p0.z + (float)rn3 * h;
	}
	else if( face == 5 ){
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z;
	}
	else if( face == 6 ){
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + h;
	}
	else{
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + (float)rn3 * h;
	}
	ps.SetParticle(p1, (0.1 + 0.4*(float)rn4) * h);
	return;
}

void Voxel::PlacePoints(int ii, int jj, int kk,
						 Point &ps, int face) const{
	Point p0, &p1 = ps;
	p0.x = P0.x + h*ii;
	p0.y = P0.y + h*jj;
	p0.z = P0.z + h*kk;
	MTRand mt;
	//for(int i=0; i<np; i++){
	double rn1 = mt();
	double rn2 = mt();
	double rn3 = mt();
	if(face == 1){
		//p1.x = p0.x + (float)rn1 * h;
		//p1.x = p0.x + 0.25f*h;
		p1.x = p0.x + 0.5f * h + (float)rn1 * h;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + (float)rn3 * h;
	}
	else if( face == 2 ){
		//p1.x = p0.x + (float)rn1 * h;
		//p1.x = p0.x + h - 0.25f*h;
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y +  0.5f * h + (float)rn2 * h;
		p1.z = p0.z + (float)rn3 * h;
	}
	else if( face == 3 ){
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + 0.5f * h + (float)rn3 * h;
	}
	else if( face == 4 ){
		p1.x = p0.x + (float)rn1 * h;
		//p1.y = p0.y + h - 0.25f*h;
		p1.y = p0.y + h;
		//p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + (float)rn3 * h;
	}
	else if( face == 5 ){
		p1.x = p0.x + (float)rn1 * h;
		//p1.y = p0.y + h;
		p1.y = p0.y + (float)rn2 * h;
		//p1.z = p0.z + 0.25f*h;
		p1.z = p0.z;
	}
	else if( face == 6 ){
		p1.x = p0.x + (float)rn1 * h;
		//p1.y = p0.y + h;
		p1.y = p0.y + (float)rn2 * h;
		//p1.z = p0.z + h- 0.25f*h;
		p1.z = p0.z + h;
	}
	else{
		p1.x = p0.x + (float)rn1 * h;
		p1.y = p0.y + (float)rn2 * h;
		p1.z = p0.z + (float)rn3 * h;
	}

	return;
}

Point Voxel::VoxelCenterPosition(int i, int j, int k) const{
#ifdef SPMDMPI
	if(comm_both)
		return Point(P0.x + h*(i-num_ghosts+0.5f),
				     P0.y + h*(j-num_ghosts+0.5f),
			  		 P0.z + h*(k+0.5));
	else{
		if(comm_dir == 1)
			return Point(P0.x + h*(i-num_ghosts+0.5f),
						 P0.y + h*(j+0.5f),
						 P0.z + h*(k+0.5));
		else
			return Point(P0.x + h*(i+0.5f),
					     P0.y + h*(j-num_ghosts+0.5f),
						 P0.z + h*(k+0.5));
	}
#else
	return Point(P0.x + h*(i+0.5),
		    	 P0.y + h*(j+0.5),
	  			 P0.z + h*(k+0.5));
#endif
}

Point Voxel::SubVoxelCenterPosition(int i, int j, int k, int level, int index) const{

	if(level == 1){
		switch(index){
		case 1:
			return Point(P0.x + h*(i+0.25f),
						 P0.y + h*(j+0.25f),
						 P0.z + h*(k+0.25f));
		case 2:
			return Point(P0.x + h*(i+0.75f),
						 P0.y + h*(j+0.25f),
						 P0.z + h*(k+0.25f));
		case 3:
			return Point(P0.x + h*(i+0.25f),
						 P0.y + h*(j+0.75f),
						 P0.z + h*(k+0.25f));
		case 4:
			return Point(P0.x + h*(i+0.75f),
						 P0.y + h*(j+0.75f),
						 P0.z + h*(k+0.25f));
		case 5:
			return Point(P0.x + h*(i+0.25f),
						 P0.y + h*(j+0.25f),
						 P0.z + h*(k+0.75f));
		case 6:
			return Point(P0.x + h*(i+0.75f),
						 P0.y + h*(j+0.25f),
						 P0.z + h*(k+0.75f));
		case 7:
			return Point(P0.x + h*(i+0.25f),
						 P0.y + h*(j+0.75f),
						 P0.z + h*(k+0.75f));
		case 8:
			return Point(P0.x + h*(i+0.75f),
						 P0.y + h*(j+0.75f),
						 P0.z + h*(k+0.75f));
		}
	}
	else{
		fprintf(stderr, "Error! now only level 1 is supported\n ");
		exit(1);
	}

}


Point Voxel::VoxelCornerPosition(int i, int j, int k, int index) const{
#ifdef SPMDMPI
	if(index == 1){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts),
						 P0.y + h*(j-num_ghosts),
					     P0.z + h*k);
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts),
							 P0.y + h*j,
							 P0.z + h*k);
			else
				return Point(P0.x + h*(i),
							 P0.y + h*(j-num_ghosts),
							 P0.z + h*k);
		}
	}
	else if( index  == 2 ){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts+1),
						 P0.y + h*(j-num_ghosts),
						 P0.z + h*(k));
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts+1),
							 P0.y + h*(j),
							 P0.z + h*(k));
			else
				return Point(P0.x + h*(i+1),
							 P0.y + h*(j-num_ghosts),
							 P0.z + h*(k));
		}
	}
	else if( index  == 3 ){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts+1),
						 P0.y + h*(j-num_ghosts+1),
						 P0.z + h*(k));
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts+1),
							 P0.y + h*(j+1),
							 P0.z + h*(k));
			else
				return Point(P0.x + h*(i+1),
							 P0.y + h*(j-num_ghosts+1),
							 P0.z + h*(k));
		}
	}
	else if( index  == 4 ){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts),
						 P0.y + h*(j-num_ghosts+1),
						 P0.z + h*(k));
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts),
							 P0.y + h*(j+1),
							 P0.z + h*(k));
			else
				return Point(P0.x + h*(i),
							 P0.y + h*(j-num_ghosts+1),
							 P0.z + h*(k));
		}
	}
	else if( index == 5 ){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts),
						 P0.y + h*(j-num_ghosts),
						 P0.z + h*(k+1));
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts),
							 P0.y + h*(j),
							 P0.z + h*(k+1));
			else
				return Point(P0.x + h*(i),
							 P0.y + h*(j-num_ghosts),
							 P0.z + h*(k+1));
		}
	}
	else if( index == 6 ){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts+1),
						 P0.y + h*(j-num_ghosts),
						 P0.z + h*(k+1));
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts+1),
							 P0.y + h*(j),
							 P0.z + h*(k+1));
			else
				return Point(P0.x + h*(i+1),
							 P0.y + h*(j-num_ghosts),
							 P0.z + h*(k+1));
		}
	}
	else if( index == 7 ){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts+1),
						 P0.y + h*(j-num_ghosts+1),
						 P0.z + h*(k+1));
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts+1),
							 P0.y + h*(j+1),
							 P0.z + h*(k+1));
			else
				return Point(P0.x + h*(i+1),
							 P0.y + h*(j-num_ghosts+1),
							 P0.z + h*(k+1));
		}
	}
	else if( index == 8 ){
		if(comm_both)
			return Point(P0.x + h*(i-num_ghosts),
						 P0.y + h*(j-num_ghosts+1),
						 P0.z + h*(k+1));
		else{
			if(comm_dir == 1)
				return Point(P0.x + h*(i-num_ghosts),
							 P0.y + h*(j+1),
							 P0.z + h*(k+1));
			else
				return Point(P0.x + h*(i),
							 P0.y + h*(j-num_ghosts+1),
							 P0.z + h*(k+1));
		}
	}
	else{
		if(!pid){
			printf("Error! parameter [axis] must be in range [1...3]\n ");
			fflush(stdout);
		}
		MPI_Finalize();
		exit(1);
	}
#else
	if(index == 1)
		return Point(P0.x + h*i,
					 P0.y + h*j,
					 P0.z + h*k);
	else if( index  == 2 )
		return Point(P0.x + h*(i+1),
					P0.y + h*(j),
					P0.z + h*(k));
	else if( index  == 3 )
			return Point(P0.x + h*(i+1),
						P0.y + h*(j+1),
						P0.z + h*(k));
	else if( index  == 4 )
			return Point(P0.x + h*(i),
						P0.y + h*(j+1),
						P0.z + h*(k));
	else if( index  == 5 )
			return Point(P0.x + h*(i),
						P0.y + h*(j),
						P0.z + h*(k+1));
	else if( index  == 6 )
			return Point(P0.x + h*(i+1),
						P0.y + h*(j),
						P0.z + h*(k+1));
	else if( index  == 7 )
			return Point(P0.x + h*(i+1),
						P0.y + h*(j+1),
						P0.z + h*(k+1));
	else if( index  == 8 )
			return Point(P0.x + h*(i),
						P0.y + h*(j+1),
						P0.z + h*(k+1));
	else{
		fprintf(stderr, "Error! parameter [index] must be in range [1...8]\n ");
		exit(1);
	}
#endif

}


TentativeGridPoint Voxel::ContainsPoint( const Point &p ) const{
	int ii, jj;
#ifdef SPMDMPI
	if(comm_both){
		ii = floor((p.x - P0.x) / h) + num_ghosts;
		jj = floor((p.y - P0.y) / h) + num_ghosts;
	}
	else{
		if(comm_dir == 1){
			ii = floor((p.x - P0.x) / h) + num_ghosts;
			jj = floor((p.y - P0.y) / h);
		}
		else{
			ii = floor((p.x - P0.x) / h);
			jj = floor((p.y - P0.y) / h) + num_ghosts;
		}
	}
#else
	ii = floor(p.x / h);
	jj = floor(p.y / h);
#endif
	int kk = floor(p.z / h);
	return TentativeGridPoint(0.f, ii, jj, kk);

}

void Voxel::Initialize(const Container *c){
	u_int solid = 0, nonsolid = 0;
	FOR_EACH_CELL
			if( c->InContainer(i,j,k) ){
				UpdateCellType(i, j, k, SET, SOLID);
				++solid;
			}
			else{
				UpdateCellType(i, j, k, SET, EMPTY);
				++nonsolid;
			}
	END_FOR
#ifdef SPMDMPI
	UpdateGhosts();
	printf("pid = %d, solid = %u, nonsolid = %u \n", pid, solid, nonsolid);
	fflush(stdout);
#else
	printf("solid = %u, nonsolid = %u \n", solid, nonsolid);
#endif
}

TentativeGridPoint Voxel::ContainsPointInVelCell( const Point &p,
		unsigned char vel_index ) const{
	int ii, jj, kk;
#ifdef SPMDMPI
	if(vel_index == 1){
		if(comm_both){
			ii = floor((p.x-P0.x) / h - 0.5f) + num_ghosts;
			jj = floor((p.y-P0.y) / h) + num_ghosts;
		}
		else{
			if(comm_dir == 1){
				ii = floor((p.x-P0.x) / h - 0.5f) + num_ghosts;
				jj = floor((p.y-P0.y) / h);
			}
			else{
				ii = floor((p.x-P0.x) / h - 0.5f);
				jj = floor((p.y-P0.y) / h) + num_ghosts;
			}
		}
		kk = floor((p.z-P0.z) / h);
	}
	else if(vel_index == 2){
		if(comm_both){
			ii = floor((p.x-P0.x) / h) + num_ghosts;
			jj = floor((p.y-P0.y) / h - 0.5f) + num_ghosts;
		}
		else{
			if(comm_dir == 1){
				ii = floor((p.x-P0.x) / h) + num_ghosts;
				jj = floor((p.y-P0.y) / h - 0.5f);
			}
			else{
				ii = floor((p.x-P0.x) / h);
				jj = floor((p.y-P0.y) / h - 0.5f) + num_ghosts;
			}
		}
		kk = floor((p.z-P0.z) / h);
	}
	else if(vel_index == 3){
		if(comm_both){
			ii = floor((p.x-P0.x) / h) + num_ghosts;
			jj = floor((p.y-P0.y) / h) + num_ghosts;
		}
		else{
			if(comm_dir == 1){
				ii = floor((p.x-P0.x) / h) + num_ghosts;
				jj = floor((p.y-P0.y) / h);
			}
			else{
				ii = floor((p.x-P0.x) / h);
				jj = floor((p.y-P0.y) / h) + num_ghosts;
			}
		}
		kk = floor((p.z-P0.z) / h - 0.5f);
	}
#else
	if(vel_index == 1){
		ii = floor(p.x / h  - 0.5f);
		jj = floor(p.y / h);
		kk = floor(p.z / h);
	}
	else if(vel_index == 2){
		ii = floor(p.x / h);
		jj = floor(p.y / h - 0.5f);
		kk = floor(p.z / h);
	}
	else if(vel_index == 3){
		ii = floor(p.x / h);
		jj = floor(p.y / h);
		kk = floor(p.z / h - 0.5f);
	}
#endif
#ifdef SPMDMPI
	else{
		printf("Error! parameter [vel_index] must be in range [1...3]\n ");
		fflush(stdout);
		MPI_Finalize();
		exit(1);
	}
#else
	else{
		fprintf(stderr, "Error! parameter [vel_index] must be in range [1...3]\n ");
		exit(1);
	}
#endif
	return TentativeGridPoint(0.f, ii, jj, kk);
}
