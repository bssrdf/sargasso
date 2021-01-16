#ifndef MPIUTILITIES_H_
#define MPIUTILITIES_H_

#ifdef SPMDMPI
#ifdef WIN32
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#endif

#include "mpi.h"

#define INDEX_ALL(i,j,k) (k)*imt*jmt+(j)*imt+(i)
#define INDEX(i,j,k) (k)*imt*jmt+(j)*imt+(i)

#define FOR_EACH_CELL for(int k=0;k<DimZ;k++) { \
						 for(int j=jphys_b;j<=jphys_e;j++) { \
							for(int i=iphys_b;i<=iphys_e;i++){
#define END_FOR }}}

const int ndims = 2;    // dimensions of the virtual process grid

int num_ghosts = 5;  // # of gohst cells in each direction

const int reorder = 0;

const int
      mpitag_nshift = 1,        // MPI tags for various
      mpitag_sshift = 2,         // communication patterns
      mpitag_eshift = 3,
      mpitag_wshift = 4,
      mpitag_io     = 5,
      mpitag_gthr   = 6,
      mpitag_neshift = 7,
      mpitag_swshift = 8,
      mpitag_nwshift = 9,
      mpitag_seshift =10;

int nproc_s;  // total # of pes

int nproc[ndims];  // pes in each dimension

int periodic[ndims];

int cellx, celly;
int imt, jmt, kmt;
int global_imt, global_jmt, global_kmt;
int iphys_b, iphys_e;

int jphys_b, jphys_e;



int nbr_east, nbr_west;

int nbr_south, nbr_north;

int nbr_ne, nbr_se, nbr_nw, nbr_sw;

int coords[ndims];  // coordinates of the current pe within the virtual grid

int coords_ne[ndims], coords_se[ndims], coords_nw[ndims], coords_sw[ndims];

int pid;   // rank of the current pe

MPI_Comm this_cart;

bool comm_both = false;
char comm_dir = 0;  // 1: x-dir
                    // 2: y-dir

MPI_Group world, pres_solver_group;
MPI_Comm pres_solver_comm;
int non_pres_solver_procs;
int pres_solver_procs;

static float *f_sb_w2e, *f_sb_e2w;
static float *f_rb_w2e, *f_rb_e2w;
static float *f_sb_n2s, *f_sb_s2n;
static float *f_rb_n2s, *f_rb_s2n;
static float *f_sb_sw2ne, *f_sb_ne2sw;
static float *f_rb_sw2ne, *f_rb_ne2sw;
static float *f_sb_nw2se, *f_sb_se2nw;
static float *f_rb_nw2se, *f_rb_se2nw;

static char *c_sb_w2e, *c_sb_e2w;
static char *c_rb_w2e, *c_rb_e2w;
static char *c_sb_n2s, *c_sb_s2n;
static char *c_rb_n2s, *c_rb_s2n;
static char *c_sb_sw2ne, *c_sb_ne2sw;
static char *c_rb_sw2ne, *c_rb_ne2sw;
static char *c_sb_nw2se, *c_sb_se2nw;
static char *c_rb_nw2se, *c_rb_se2nw;

static unsigned int *ui_sb_w2e, *ui_sb_e2w;
static unsigned int *ui_rb_w2e, *ui_rb_e2w;
static unsigned int *ui_sb_n2s, *ui_sb_s2n;
static unsigned int *ui_rb_n2s, *ui_rb_s2n;
static unsigned int *ui_sb_sw2ne, *ui_sb_ne2sw;
static unsigned int *ui_rb_sw2ne, *ui_rb_ne2sw;
static unsigned int *ui_sb_nw2se, *ui_sb_se2nw;
static unsigned int *ui_rb_nw2se, *ui_rb_se2nw;
static int l_ew, l_sn, l_corner;

void init_mpi(int c, char *v[], int nx, int ny, int dimx, int dimy, int dimz){
	int pes;
	nproc[0] = nx;
	nproc[1] = ny;
	nproc_s  = nx * ny;

	global_imt = dimx;
	global_jmt = dimy;
	global_kmt = dimz;

	kmt = dimz;

	MPI_Init (&c, &v);
	MPI_Comm_rank (MPI_COMM_WORLD, &pid);
	MPI_Comm_size (MPI_COMM_WORLD, &pes);
	if(pes != nproc_s){
		if(!pid){
			printf(" pes in x-dir [%d] * pes in y-dir [%d] != total pes [%d] \n",
		          nx, ny, pes);
			fflush(stdout);
		}
		MPI_Finalize();
		exit(1);
	}
	if(dimx % nproc[0] != 0){
		if(!pid){
			printf(" dimx [%d] must be muiltiples of pes in x-dir [%d] \n",
		          dimx, nproc[0]);
			fflush(stdout);
		}
		MPI_Finalize();
		exit(1);
	}
	if(dimy % nproc[1] != 0){
		if(!pid){
			printf(" dimy [%d] must be muiltiples of pes in y-dir [%d] \n",
		          dimy, nproc[1]);
			fflush(stdout);
		}
		MPI_Finalize();
		exit(1);
	}

	cellx = dimx / nproc[0];
	celly = dimy / nproc[1];

	periodic[0] = periodic[1] = 1;

	// create Cartitian virtual process grid
	MPI_Dims_create(pes, ndims, nproc);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, nproc, periodic, reorder, &this_cart);

	MPI_Cart_coords(this_cart, pid, ndims, coords);
	MPI_Cart_shift(this_cart, 0, 1, &nbr_west, &nbr_east);
	MPI_Cart_shift(this_cart, 1, 1, &nbr_south, &nbr_north);

    if(nx > 1 && ny > 1){

    	comm_both = true;

		coords_ne[0] = coords[0] + 1;
	    coords_ne[1] = coords[1] + 1;
	    if(coords_ne[0] > nproc[0]-1)
	       coords_ne[0] = 0;
	    if(coords_ne[1] > nproc[1]-1)
	       coords_ne[1] = 0;
	    MPI_Cart_rank(this_cart, coords_ne, &nbr_ne);

	    coords_se[0] = coords[0] + 1;
	    coords_se[1] = coords[1] - 1;
	    if(coords_se[0] > nproc[0]-1)
	       coords_se[0] = 0;
	    if(coords_se[1] < 0)
	       coords_se[1] = nproc[1]-1;
	    MPI_Cart_rank(this_cart, coords_se, &nbr_se);

	    coords_nw[0] = coords[0] - 1;
	    coords_nw[1] = coords[1] + 1;
	    if(coords_nw[0] < 0)
	       coords_nw[0] = nproc[0]-1;
	    if(coords_nw[1] > nproc[1]-1)
	       coords_nw[1] = 0;
	    MPI_Cart_rank(this_cart, coords_nw, &nbr_nw);

	    coords_sw[0] = coords[0] - 1;
	    coords_sw[1] = coords[1] - 1;
	    if(coords_sw[0] < 0)
	       coords_sw[0] = nproc[0]-1;
	    if(coords_sw[1] < 0)
	       coords_sw[1] = nproc[1]-1;
	    MPI_Cart_rank(this_cart, coords_sw, &nbr_sw);

    }
    else{
    	if(nx > 1)
    		comm_dir = 1;
    	else
    		comm_dir = 2;
    }
    if(comm_both){
    	imt = cellx + 2 * num_ghosts;
    	jmt = celly + 2 * num_ghosts;
    	iphys_b = num_ghosts;
    	iphys_e = imt - num_ghosts - 1;
    	jphys_b = num_ghosts;
    	jphys_e = jmt - num_ghosts - 1;
    	l_ew = num_ghosts*(jmt-2*num_ghosts)*kmt;
		f_sb_w2e = new float[l_ew];
		f_sb_e2w = new float[l_ew];
		f_rb_w2e = new float[l_ew];
		f_rb_e2w = new float[l_ew];
		c_sb_w2e = new  char[l_ew];
		c_sb_e2w = new  char[l_ew];
		c_rb_w2e = new  char[l_ew];
		c_rb_e2w = new  char[l_ew];
		ui_sb_w2e = new  unsigned int[l_ew];
		ui_sb_e2w = new  unsigned int[l_ew];
		ui_rb_w2e = new  unsigned int[l_ew];
		ui_rb_e2w = new  unsigned int[l_ew];
		l_sn = num_ghosts*(imt-2*num_ghosts)*kmt;
		f_sb_n2s = new float[l_sn];
		f_sb_s2n = new float[l_sn];
		f_rb_n2s = new float[l_sn];
		f_rb_s2n = new float[l_sn];
		c_sb_n2s = new  char[l_sn];
		c_sb_s2n = new  char[l_sn];
		c_rb_n2s = new  char[l_sn];
		c_rb_s2n = new  char[l_sn];
		ui_sb_n2s = new unsigned int[l_sn];
		ui_sb_s2n = new unsigned int[l_sn];
		ui_rb_n2s = new unsigned int[l_sn];
		ui_rb_s2n = new unsigned int[l_sn];
		l_corner = num_ghosts*num_ghosts*kmt;
		f_sb_ne2sw = new float[l_corner];
		f_sb_sw2ne = new float[l_corner];
		f_rb_sw2ne = new float[l_corner];
		f_rb_ne2sw = new float[l_corner];
		f_sb_nw2se = new float[l_corner];
		f_sb_se2nw = new float[l_corner];
		f_rb_nw2se = new float[l_corner];
		f_rb_se2nw = new float[l_corner];
		c_sb_ne2sw = new  char[l_corner];
		c_sb_sw2ne = new  char[l_corner];
		c_rb_sw2ne = new  char[l_corner];
		c_rb_ne2sw = new  char[l_corner];
		c_sb_nw2se = new  char[l_corner];
		c_sb_se2nw = new  char[l_corner];
		c_rb_nw2se = new  char[l_corner];
		c_rb_se2nw = new  char[l_corner];
		ui_sb_ne2sw = new  unsigned int[l_corner];
		ui_sb_sw2ne = new  unsigned int[l_corner];
		ui_rb_sw2ne = new  unsigned int[l_corner];
		ui_rb_ne2sw = new  unsigned int[l_corner];
		ui_sb_nw2se = new  unsigned int[l_corner];
		ui_sb_se2nw = new  unsigned int[l_corner];
		ui_rb_nw2se = new  unsigned int[l_corner];
		ui_rb_se2nw = new  unsigned int[l_corner];
    }
    else{
    	if(comm_dir == 1){
    		imt = cellx + 2 * num_ghosts;
    		jmt = dimy;
    		iphys_b = num_ghosts;
    		iphys_e = imt - num_ghosts - 1;
    		jphys_b = 0;
    		jphys_e = jmt - 1;
    		l_ew = num_ghosts*jmt*kmt;
			f_sb_w2e = new float[l_ew];
			f_sb_e2w = new float[l_ew];
			f_rb_w2e = new float[l_ew];
			f_rb_e2w = new float[l_ew];
			c_sb_w2e = new  char[l_ew];
			c_sb_e2w = new  char[l_ew];
			c_rb_w2e = new  char[l_ew];
			c_rb_e2w = new  char[l_ew];
			ui_sb_w2e = new  unsigned int[l_ew];
			ui_sb_e2w = new  unsigned int[l_ew];
			ui_rb_w2e = new  unsigned int[l_ew];
			ui_rb_e2w = new  unsigned int[l_ew];
    	}
    	else{
    		imt = dimx;
    		jmt = celly + 2 * num_ghosts;
    		iphys_b = 0;
	    	iphys_e = imt - 1;
	    	jphys_b = num_ghosts;
	    	jphys_e = jmt - num_ghosts - 1;
	    	l_sn = num_ghosts*imt*kmt;
			f_sb_n2s = new float[l_sn];
			f_sb_s2n = new float[l_sn];
			f_rb_n2s = new float[l_sn];
			f_rb_s2n = new float[l_sn];
			c_sb_n2s = new  char[l_sn];
			c_sb_s2n = new  char[l_sn];
			c_rb_n2s = new  char[l_sn];
			c_rb_s2n = new  char[l_sn];
			ui_sb_n2s = new  unsigned int[l_sn];
			ui_sb_s2n = new  unsigned int[l_sn];
			ui_rb_n2s = new  unsigned int[l_sn];
			ui_rb_s2n = new  unsigned int[l_sn];
    	}
    }

}

void boundary3d_scalar(float *array){

	MPI_Request request[2];
	MPI_Status  status[2];

	if(comm_both){

		// send east-west boundary info

		float *sb_w2e = f_sb_w2e, *sb_e2w = f_sb_e2w;
		float *rb_w2e = f_rb_w2e, *rb_e2w = f_rb_e2w;
		int l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int i = 0; i < num_ghosts; ++i){
				int i1 = iphys_e-num_ghosts+i+1;
				int i2 = iphys_b+i;
				for(int j = jphys_b; j <= jphys_e; ++j){
				    sb_w2e[l] = array[INDEX_ALL(i1,j,k)];
				    sb_e2w[l] = array[INDEX_ALL(i2,j,k)];
				    ++l;
			}
		}
		MPI_Isend(sb_w2e, l, MPI_FLOAT, nbr_east, mpitag_wshift, this_cart, &request[0]);
		MPI_Isend(sb_e2w, l, MPI_FLOAT, nbr_west, mpitag_eshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_w2e, l_ew, MPI_FLOAT, nbr_west, mpitag_wshift, this_cart, &status[0]);
		MPI_Recv(rb_e2w, l_ew, MPI_FLOAT, nbr_east, mpitag_eshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int i = 0; i < num_ghosts; ++i){
				int i1 = i;
				int i2 = iphys_e+i+1;
				for(int j = jphys_b; j <= jphys_e; ++j){
				    array[INDEX_ALL(i1,j,k)] = rb_w2e[l];
				    array[INDEX_ALL(i2,j,k)] = rb_e2w[l];
				    ++l;
			}
		}

		// receive north-south boundary info

		float *sb_n2s = f_sb_n2s, *sb_s2n = f_sb_s2n;
		float *rb_n2s = f_rb_n2s, *rb_s2n = f_rb_s2n;
		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int j = 0; j < num_ghosts; ++j){
				int j1 = jphys_e-num_ghosts+j+1;
				int j2 = jphys_b+j;
				for(int i = iphys_b; i <= iphys_e; ++i){
				    sb_s2n[l] = array[INDEX_ALL(i,j1,k)];
				    sb_n2s[l] = array[INDEX_ALL(i,j2,k)];
				    ++l;
			}
		}
		MPI_Isend(sb_n2s, l, MPI_FLOAT, nbr_south, mpitag_nshift, this_cart, &request[0]);
		MPI_Isend(sb_s2n, l, MPI_FLOAT, nbr_north, mpitag_sshift, this_cart, &request[1]);

		MPI_Recv(rb_n2s, l_sn, MPI_FLOAT, nbr_north, mpitag_nshift, this_cart, &status[0]);
		MPI_Recv(rb_s2n, l_sn, MPI_FLOAT, nbr_south, mpitag_sshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int j = 0; j < num_ghosts; ++j){
				int j1 = j;
				int j2 = jphys_e+j+1;
				for(int i = iphys_b; i <= iphys_e; ++i){
				    array[INDEX_ALL(i,j1,k)] = rb_s2n[l];
				    array[INDEX_ALL(i,j2,k)] = rb_n2s[l];
				    ++l;
			}
		}

		// send-recv ne-sw boundary information
		float *sb_ne2sw = f_sb_ne2sw, *sb_sw2ne = f_sb_sw2ne;
		float *rb_ne2sw = f_rb_ne2sw, *rb_sw2ne = f_rb_sw2ne;
		int l1 = 0, l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = iphys_b; i < iphys_b + num_ghosts; ++i)
				for(int j = jphys_b; j < jphys_b + num_ghosts; ++j){
				    sb_ne2sw[l1] = array[INDEX_ALL(i,j,k)];
				    ++l1;
				}
			for(int i = iphys_e-num_ghosts+1; i <= iphys_e; ++i)
				for(int j = jphys_e-num_ghosts+1; j <= jphys_e; ++j){
				    sb_sw2ne[l2] = array[INDEX_ALL(i,j,k)];
				    ++l2;
				}
		}
		MPI_Isend(sb_ne2sw, l1, MPI_FLOAT, nbr_sw, mpitag_neshift, this_cart, &request[0]);
		MPI_Isend(sb_sw2ne, l2, MPI_FLOAT, nbr_ne, mpitag_swshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_ne2sw, l_corner, MPI_FLOAT, nbr_ne, mpitag_neshift, this_cart, &status[0]);
		MPI_Recv(rb_sw2ne, l_corner, MPI_FLOAT, nbr_sw, mpitag_swshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l1 = 0; l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = 0; i < num_ghosts; ++i)
				for(int j = 0; j < num_ghosts; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_sw2ne[l1];
				    ++l1;
				}
			for(int i = iphys_e+1; i < imt; ++i)
				for(int j = jphys_e+1; j < jmt; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_ne2sw[l2];
				    ++l2;
				}
		}

		// send-recv nw-se boundary information
		float *sb_nw2se = f_sb_nw2se, *sb_se2nw = f_sb_se2nw;
		float *rb_nw2se = f_sb_nw2se, *rb_se2nw = f_rb_se2nw;

		l1 = 0, l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = iphys_e-num_ghosts+1; i <= iphys_e; ++i)
				for(int j = jphys_b; j < jphys_b + num_ghosts; ++j){
				    sb_nw2se[l1] = array[INDEX_ALL(i,j,k)];
				    ++l1;
				}
			for(int i = iphys_b; i < iphys_b + num_ghosts; ++i)
				for(int j = jphys_e-num_ghosts+1; j <= jphys_e; ++j){
				    sb_se2nw[l2] = array[INDEX_ALL(i,j,k)];
				    ++l2;
				}
		}
		MPI_Isend(sb_nw2se, l1, MPI_FLOAT, nbr_se, mpitag_nwshift, this_cart, &request[0]);
		MPI_Isend(sb_se2nw, l2, MPI_FLOAT, nbr_nw, mpitag_seshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_nw2se, l_corner, MPI_FLOAT, nbr_nw, mpitag_nwshift, this_cart, &status[0]);
		MPI_Recv(rb_se2nw, l_corner, MPI_FLOAT, nbr_se, mpitag_seshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l1 = 0; l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = 0; i < num_ghosts; ++i)
				for(int j = jphys_e+1; j < jmt; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_nw2se[l1];
				    ++l1;
				}
			for(int i = iphys_e+1; i < imt; ++i)
				for(int j = 0; j < num_ghosts; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_se2nw[l2];
				    ++l2;
				}
		}

	}
	else{
		if(comm_dir == 1){

			// send east-west boundary info

			float *sb_w2e = f_sb_w2e, *sb_e2w = f_sb_e2w;
			float *rb_w2e = f_rb_w2e, *rb_e2w = f_rb_e2w;
			int l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int i = 0; i < num_ghosts; ++i){
					int i1 = iphys_e-num_ghosts+i+1;
					int i2 = iphys_b+i;
					for(int j = jphys_b; j <= jphys_e; ++j){
					    sb_w2e[l] = array[INDEX_ALL(i1,j,k)];
					    sb_e2w[l] = array[INDEX_ALL(i2,j,k)];
					    ++l;
				}
			}
//			if(!pid){
//				printf(" l_ew = %d, l = %d \n", l_ew, l);
//				fflush(stdout);
//			}

			MPI_Isend(sb_w2e, l, MPI_FLOAT, nbr_east, mpitag_wshift, this_cart, &request[0]);
			MPI_Isend(sb_e2w, l, MPI_FLOAT, nbr_west, mpitag_eshift, this_cart, &request[1]);

			MPI_Recv(rb_w2e, l_ew, MPI_FLOAT, nbr_west, mpitag_wshift, this_cart, &status[0]);
			MPI_Recv(rb_e2w, l_ew, MPI_FLOAT, nbr_east, mpitag_eshift, this_cart, &status[1]);

			MPI_Waitall(2, request, status);

			l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int i = 0; i < num_ghosts; ++i){
					int i1 = i;
					int i2 = iphys_e+i+1;
					for(int j = jphys_b; j <= jphys_e; ++j){
					    array[INDEX_ALL(i1,j,k)] = rb_w2e[l];
					    array[INDEX_ALL(i2,j,k)] = rb_e2w[l];
					    ++l;
				}
			}
//			if(!pid){
//				printf(" l_ew = %d, l = %d, iphys_b = %d, iphys_e = %d ",
//						l_ew, l, iphys_b, iphys_e);
//				fflush(stdout);
//			}

		}
		else{
			// receive north-south boundary info

			float *sb_n2s = f_sb_n2s, *sb_s2n = f_sb_s2n;
			float *rb_n2s = f_rb_n2s, *rb_s2n = f_rb_s2n;

			int l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int j = 0; j < num_ghosts; ++j){
					int j1 = jphys_e-num_ghosts+j+1;
					int j2 = jphys_b+j;
					for(int i = iphys_b; i <= iphys_e; ++i){
					    sb_s2n[l] = array[INDEX_ALL(i,j1,k)];
					    sb_n2s[l] = array[INDEX_ALL(i,j2,k)];
					    ++l;
				}
			}
			MPI_Isend(sb_n2s, l, MPI_FLOAT, nbr_south, mpitag_nshift, this_cart, &request[0]);
			MPI_Isend(sb_s2n, l, MPI_FLOAT, nbr_north, mpitag_sshift, this_cart, &request[1]);

			MPI_Recv(rb_n2s, l_sn, MPI_FLOAT, nbr_north, mpitag_nshift, this_cart, &status[0]);
			MPI_Recv(rb_s2n, l_sn, MPI_FLOAT, nbr_south, mpitag_sshift, this_cart, &status[1]);

			MPI_Waitall(2, request, status);

			l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int j = 0; j < num_ghosts; ++j){
					int j1 = j;
					int j2 = jphys_e+j+1;
					for(int i = iphys_b; i <= iphys_e; ++i){
					    array[INDEX_ALL(i,j1,k)] = rb_s2n[l];
					    array[INDEX_ALL(i,j2,k)] = rb_n2s[l];
//					    if(!pid && k == 0){
//					    	printf(" wrong i = %d, j1= %d, j2 = %d\n", i, j1, j2);
//					    	fflush(stdout);
//					    }
//					    if(!pid && i == 21 && j2 == 10 && k == 0){
//					    	printf(" wrong 2 \n");
//					    	fflush(stdout);
//					    }
					    ++l;
				}
			}
		}
	}
}

void boundary3d_scalarc(char *array){

	MPI_Request request[2];
	MPI_Status  status[2];

	if(comm_both){

		// send east-west boundary info

		char *sb_w2e = c_sb_w2e, *sb_e2w = c_sb_e2w;
		char *rb_w2e = c_rb_w2e, *rb_e2w = c_rb_e2w;
		int l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int i = 0; i < num_ghosts; ++i){
				int i1 = iphys_e-num_ghosts+i+1;
				int i2 = iphys_b+i;
				for(int j = jphys_b; j <= jphys_e; ++j){
				    sb_w2e[l] = array[INDEX_ALL(i1,j,k)];
				    sb_e2w[l] = array[INDEX_ALL(i2,j,k)];
				    ++l;
			}
		}
		MPI_Isend(sb_w2e, l, MPI_CHAR, nbr_east, mpitag_wshift, this_cart, &request[0]);
		MPI_Isend(sb_e2w, l, MPI_CHAR, nbr_west, mpitag_eshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_w2e, l_ew, MPI_CHAR, nbr_west, mpitag_wshift, this_cart, &status[0]);
		MPI_Recv(rb_e2w, l_ew, MPI_CHAR, nbr_east, mpitag_eshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int i = 0; i < num_ghosts; ++i){
				int i1 = i;
				int i2 = iphys_e+i+1;
				for(int j = jphys_b; j <= jphys_e; ++j){
				    array[INDEX_ALL(i1,j,k)] = rb_w2e[l];
				    array[INDEX_ALL(i2,j,k)] = rb_e2w[l];
				    ++l;
			}
		}

		// receive north-south boundary info

		char *sb_n2s = c_sb_n2s, *sb_s2n = c_sb_s2n;
		char *rb_n2s = c_rb_n2s, *rb_s2n = c_rb_s2n;
		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int j = 0; j < num_ghosts; ++j){
				int j1 = jphys_e-num_ghosts+j+1;
				int j2 = jphys_b+j;
				for(int i = iphys_b; i <= iphys_e; ++i){
				    sb_s2n[l] = array[INDEX_ALL(i,j1,k)];
				    sb_n2s[l] = array[INDEX_ALL(i,j2,k)];
				    ++l;
			}
		}
		MPI_Isend(sb_n2s, l, MPI_CHAR, nbr_south, mpitag_nshift, this_cart, &request[0]);
		MPI_Isend(sb_s2n, l, MPI_CHAR, nbr_north, mpitag_sshift, this_cart, &request[1]);

		MPI_Recv(rb_n2s, l_sn, MPI_CHAR, nbr_north, mpitag_nshift, this_cart, &status[0]);
		MPI_Recv(rb_s2n, l_sn, MPI_CHAR, nbr_south, mpitag_sshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int j = 0; j < num_ghosts; ++j){
				int j1 = j;
				int j2 = jphys_e+j+1;
				for(int i = iphys_b; i <= iphys_e; ++i){
				    array[INDEX_ALL(i,j1,k)] = rb_s2n[l];
				    array[INDEX_ALL(i,j2,k)] = rb_n2s[l];
				    ++l;
			}
		}

		// send-recv ne-sw boundary information
		char *sb_ne2sw = c_sb_ne2sw, *sb_sw2ne = c_sb_sw2ne;
		char *rb_ne2sw = c_rb_ne2sw, *rb_sw2ne = c_rb_sw2ne;

		int l1 = 0, l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = iphys_b; i < iphys_b + num_ghosts; ++i)
				for(int j = jphys_b; j < jphys_b + num_ghosts; ++j){
				    sb_ne2sw[l1] = array[INDEX_ALL(i,j,k)];
				    ++l1;
				}
			for(int i = iphys_e-num_ghosts+1; i <= iphys_e; ++i)
				for(int j = jphys_e-num_ghosts+1; j <= jphys_e; ++j){
				    sb_sw2ne[l2] = array[INDEX_ALL(i,j,k)];
				    ++l2;
				}
		}
		MPI_Isend(sb_ne2sw, l1, MPI_CHAR, nbr_sw, mpitag_neshift, this_cart, &request[0]);
		MPI_Isend(sb_sw2ne, l2, MPI_CHAR, nbr_ne, mpitag_swshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_ne2sw, l_corner, MPI_CHAR, nbr_ne, mpitag_neshift, this_cart, &status[0]);
		MPI_Recv(rb_sw2ne, l_corner, MPI_CHAR, nbr_sw, mpitag_swshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l1 = 0; l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = 0; i < num_ghosts; ++i)
				for(int j = 0; j < num_ghosts; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_sw2ne[l1];
				    ++l1;
				}
			for(int i = iphys_e+1; i < imt; ++i)
				for(int j = jphys_e+1; j < jmt; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_ne2sw[l2];
				    ++l2;
				}
		}

		// send-recv nw-se boundary information
		char *sb_nw2se = c_sb_nw2se, *sb_se2nw = c_sb_se2nw;
		char *rb_nw2se = c_rb_nw2se, *rb_se2nw = c_rb_se2nw;
		l1 = 0, l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = iphys_e-num_ghosts+1; i <= iphys_e; ++i)
				for(int j = jphys_b; j < jphys_b + num_ghosts; ++j){
				    sb_nw2se[l1] = array[INDEX_ALL(i,j,k)];
				    ++l1;
				}
			for(int i = iphys_b; i < iphys_b + num_ghosts; ++i)
				for(int j = jphys_e-num_ghosts+1; j <= jphys_e; ++j){
				    sb_se2nw[l2] = array[INDEX_ALL(i,j,k)];
				    ++l2;
				}
		}
		MPI_Isend(sb_nw2se, l1, MPI_CHAR, nbr_se, mpitag_nwshift, this_cart, &request[0]);
		MPI_Isend(sb_se2nw, l2, MPI_CHAR, nbr_nw, mpitag_seshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_nw2se, l_corner, MPI_CHAR, nbr_nw, mpitag_nwshift, this_cart, &status[0]);
		MPI_Recv(rb_se2nw, l_corner, MPI_CHAR, nbr_se, mpitag_seshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l1 = 0; l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = 0; i < num_ghosts; ++i)
				for(int j = jphys_e+1; j < jmt; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_nw2se[l1];
				    ++l1;
				}
			for(int i = iphys_e+1; i < imt; ++i)
				for(int j = 0; j < num_ghosts; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_se2nw[l2];
				    ++l2;
				}
		}
	}
	else{
		if(comm_dir == 1){

			// send east-west boundary info

			char *sb_w2e = c_sb_w2e, *sb_e2w = c_sb_e2w;
			char *rb_w2e = c_rb_w2e, *rb_e2w = c_rb_e2w;
			int l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int i = 0; i < num_ghosts; ++i){
					int i1 = iphys_e-num_ghosts+i+1;
					int i2 = iphys_b+i;
					for(int j = jphys_b; j <= jphys_e; ++j){
					    sb_w2e[l] = array[INDEX_ALL(i1,j,k)];
					    sb_e2w[l] = array[INDEX_ALL(i2,j,k)];
					    ++l;
				}
			}
			MPI_Isend(sb_w2e, l, MPI_CHAR, nbr_east, mpitag_wshift, this_cart, &request[0]);
			MPI_Isend(sb_e2w, l, MPI_CHAR, nbr_west, mpitag_eshift, this_cart, &request[1]);

			MPI_Recv(rb_w2e, l_ew, MPI_CHAR, nbr_west, mpitag_wshift, this_cart, &status[0]);
			MPI_Recv(rb_e2w, l_ew, MPI_CHAR, nbr_east, mpitag_eshift, this_cart, &status[1]);

			MPI_Waitall(2, request, status);

			l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int i = 0; i < num_ghosts; ++i){
					int i1 = i;
					int i2 = iphys_e+i+1;
					for(int j = jphys_b; j <= jphys_e; ++j){
					    array[INDEX_ALL(i1,j,k)] = rb_w2e[l];
					    array[INDEX_ALL(i2,j,k)] = rb_e2w[l];
					    ++l;
				}
			}
		}
		else{
			// receive north-south boundary info

			char *sb_n2s = c_sb_n2s, *sb_s2n = c_sb_s2n;
			char *rb_n2s = c_rb_n2s, *rb_s2n = c_rb_s2n;
			int l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int j = 0; j < num_ghosts; ++j){
					int j1 = jphys_e-num_ghosts+j+1;
					int j2 = jphys_b+j;
					for(int i = iphys_b; i <= iphys_e; ++i){
					    sb_s2n[l] = array[INDEX_ALL(i,j1,k)];
					    sb_n2s[l] = array[INDEX_ALL(i,j2,k)];
					    ++l;
				}
			}
			MPI_Isend(sb_n2s, l, MPI_CHAR, nbr_south, mpitag_nshift, this_cart, &request[0]);
			MPI_Isend(sb_s2n, l, MPI_CHAR, nbr_north, mpitag_sshift, this_cart, &request[1]);

			MPI_Recv(rb_n2s, l_sn, MPI_CHAR, nbr_north, mpitag_nshift, this_cart, &status[0]);
			MPI_Recv(rb_s2n, l_sn, MPI_CHAR, nbr_south, mpitag_sshift, this_cart, &status[1]);

			MPI_Waitall(2, request, status);

			l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int j = 0; j < num_ghosts; ++j){
					int j1 = j;
					int j2 = jphys_e+j+1;
					for(int i = iphys_b; i <= iphys_e; ++i){
					    array[INDEX_ALL(i,j1,k)] = rb_s2n[l];
					    array[INDEX_ALL(i,j2,k)] = rb_n2s[l];
					    ++l;
				}
			}
		}
	}
}

void boundary3d_scalarui(unsigned int *array){

	MPI_Request request[2];
	MPI_Status  status[2];

	if(comm_both){

		// send east-west boundary info

		unsigned int *sb_w2e = ui_sb_w2e, *sb_e2w = ui_sb_e2w;
		unsigned int *rb_w2e = ui_rb_w2e, *rb_e2w = ui_rb_e2w;
		int l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int i = 0; i < num_ghosts; ++i){
				int i1 = iphys_e-num_ghosts+i+1;
				int i2 = iphys_b+i;
				for(int j = jphys_b; j <= jphys_e; ++j){
				    sb_w2e[l] = array[INDEX_ALL(i1,j,k)];
				    sb_e2w[l] = array[INDEX_ALL(i2,j,k)];
				    ++l;
			}
		}
		MPI_Isend(sb_w2e, l, MPI_UNSIGNED, nbr_east, mpitag_wshift, this_cart, &request[0]);
		MPI_Isend(sb_e2w, l, MPI_UNSIGNED, nbr_west, mpitag_eshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_w2e, l_ew, MPI_UNSIGNED, nbr_west, mpitag_wshift, this_cart, &status[0]);
		MPI_Recv(rb_e2w, l_ew, MPI_UNSIGNED, nbr_east, mpitag_eshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int i = 0; i < num_ghosts; ++i){
				int i1 = i;
				int i2 = iphys_e+i+1;
				for(int j = jphys_b; j <= jphys_e; ++j){
				    array[INDEX_ALL(i1,j,k)] = rb_w2e[l];
				    array[INDEX_ALL(i2,j,k)] = rb_e2w[l];
				    ++l;
			}
		}

		// receive north-south boundary info

		unsigned int *sb_n2s = ui_sb_n2s, *sb_s2n = ui_sb_s2n;
		unsigned int *rb_n2s = ui_rb_n2s, *rb_s2n = ui_rb_s2n;
		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int j = 0; j < num_ghosts; ++j){
				int j1 = jphys_e-num_ghosts+j+1;
				int j2 = jphys_b+j;
				for(int i = iphys_b; i <= iphys_e; ++i){
				    sb_s2n[l] = array[INDEX_ALL(i,j1,k)];
				    sb_n2s[l] = array[INDEX_ALL(i,j2,k)];
				    ++l;
			}
		}
		MPI_Isend(sb_n2s, l, MPI_UNSIGNED, nbr_south, mpitag_nshift, this_cart, &request[0]);
		MPI_Isend(sb_s2n, l, MPI_UNSIGNED, nbr_north, mpitag_sshift, this_cart, &request[1]);

		MPI_Recv(rb_n2s, l_sn, MPI_UNSIGNED, nbr_north, mpitag_nshift, this_cart, &status[0]);
		MPI_Recv(rb_s2n, l_sn, MPI_UNSIGNED, nbr_south, mpitag_sshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l = 0;
		for(int k = 0; k < kmt; ++k)
			for(int j = 0; j < num_ghosts; ++j){
				int j1 = j;
				int j2 = jphys_e+j+1;
				for(int i = iphys_b; i <= iphys_e; ++i){
				    array[INDEX_ALL(i,j1,k)] = rb_s2n[l];
				    array[INDEX_ALL(i,j2,k)] = rb_n2s[l];
				    ++l;
			}
		}

		// send-recv ne-sw boundary information
		unsigned int *sb_ne2sw = ui_sb_ne2sw, *sb_sw2ne = ui_sb_sw2ne;
		unsigned int *rb_ne2sw = ui_rb_ne2sw, *rb_sw2ne = ui_rb_sw2ne;
		int l1 = 0, l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = iphys_b; i < iphys_b + num_ghosts; ++i)
				for(int j = jphys_b; j < jphys_b + num_ghosts; ++j){
				    sb_ne2sw[l1] = array[INDEX_ALL(i,j,k)];
				    ++l1;
				}
			for(int i = iphys_e-num_ghosts+1; i <= iphys_e; ++i)
				for(int j = jphys_e-num_ghosts+1; j <= jphys_e; ++j){
				    sb_sw2ne[l2] = array[INDEX_ALL(i,j,k)];
				    ++l2;
				}
		}
		MPI_Isend(sb_ne2sw, l1, MPI_UNSIGNED, nbr_sw, mpitag_neshift, this_cart, &request[0]);
		MPI_Isend(sb_sw2ne, l2, MPI_UNSIGNED, nbr_ne, mpitag_swshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_ne2sw, l_corner, MPI_UNSIGNED, nbr_ne, mpitag_neshift, this_cart, &status[0]);
		MPI_Recv(rb_sw2ne, l_corner, MPI_UNSIGNED, nbr_sw, mpitag_swshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l1 = 0; l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = 0; i < num_ghosts; ++i)
				for(int j = 0; j < num_ghosts; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_sw2ne[l1];
				    ++l1;
				}
			for(int i = iphys_e+1; i < imt; ++i)
				for(int j = jphys_e+1; j < jmt; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_ne2sw[l2];
				    ++l2;
				}
		}

		// send-recv nw-se boundary information
		unsigned int *sb_nw2se = ui_sb_nw2se, *sb_se2nw = ui_sb_se2nw;
		unsigned int *rb_nw2se = ui_rb_nw2se, *rb_se2nw = ui_rb_se2nw;
		l1 = 0, l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = iphys_e-num_ghosts+1; i <= iphys_e; ++i)
				for(int j = jphys_b; j < jphys_b + num_ghosts; ++j){
				    sb_nw2se[l1] = array[INDEX_ALL(i,j,k)];
				    ++l1;
				}
			for(int i = iphys_b; i < iphys_b + num_ghosts; ++i)
				for(int j = jphys_e-num_ghosts+1; j <= jphys_e; ++j){
				    sb_se2nw[l2] = array[INDEX_ALL(i,j,k)];
				    ++l2;
				}
		}
		MPI_Isend(sb_nw2se, l1, MPI_UNSIGNED, nbr_se, mpitag_nwshift, this_cart, &request[0]);
		MPI_Isend(sb_se2nw, l2, MPI_UNSIGNED, nbr_nw, mpitag_seshift, this_cart, &request[1]);

		// recv east-west boundary info

		MPI_Recv(rb_nw2se, l_corner, MPI_UNSIGNED, nbr_nw, mpitag_nwshift, this_cart, &status[0]);
		MPI_Recv(rb_se2nw, l_corner, MPI_UNSIGNED, nbr_se, mpitag_seshift, this_cart, &status[1]);

		MPI_Waitall(2, request, status);

		l1 = 0; l2 = 0;
		for(int k = 0; k < kmt; ++k){
			for(int i = 0; i < num_ghosts; ++i)
				for(int j = jphys_e+1; j < jmt; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_nw2se[l1];
				    ++l1;
				}
			for(int i = iphys_e+1; i < imt; ++i)
				for(int j = 0; j < num_ghosts; ++j){
				    array[INDEX_ALL(i,j,k)] = rb_se2nw[l2];
				    ++l2;
				}
		}

	}
	else{
		if(comm_dir == 1){

			// send east-west boundary info

			unsigned int *sb_w2e = ui_sb_w2e, *sb_e2w = ui_sb_e2w;
			unsigned int *rb_w2e = ui_rb_w2e, *rb_e2w = ui_rb_e2w;
			int l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int i = 0; i < num_ghosts; ++i){
					int i1 = iphys_e-num_ghosts+i+1;
					int i2 = iphys_b+i;
					for(int j = jphys_b; j <= jphys_e; ++j){
					    sb_w2e[l] = array[INDEX_ALL(i1,j,k)];
					    sb_e2w[l] = array[INDEX_ALL(i2,j,k)];
					    ++l;
				}
			}
			MPI_Isend(sb_w2e, l, MPI_UNSIGNED, nbr_east, mpitag_wshift, this_cart, &request[0]);
			MPI_Isend(sb_e2w, l, MPI_UNSIGNED, nbr_west, mpitag_eshift, this_cart, &request[1]);

			MPI_Recv(rb_w2e, l_ew, MPI_UNSIGNED, nbr_west, mpitag_wshift, this_cart, &status[0]);
			MPI_Recv(rb_e2w, l_ew, MPI_UNSIGNED, nbr_east, mpitag_eshift, this_cart, &status[1]);

			MPI_Waitall(2, request, status);

			l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int i = 0; i < num_ghosts; ++i){
					int i1 = i;
					int i2 = iphys_e+i+1;
					for(int j = jphys_b; j <= jphys_e; ++j){
					    array[INDEX_ALL(i1,j,k)] = rb_w2e[l];
					    array[INDEX_ALL(i2,j,k)] = rb_e2w[l];
					    ++l;
				}
			}
		}
		else{
			// receive north-south boundary info

			unsigned int *sb_n2s = ui_sb_n2s, *sb_s2n = ui_sb_s2n;
			unsigned int *rb_n2s = ui_rb_n2s, *rb_s2n = ui_rb_s2n;
			int l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int j = 0; j < num_ghosts; ++j){
					int j1 = jphys_e-num_ghosts+j+1;
					int j2 = jphys_b+j;
					for(int i = iphys_b; i <= iphys_e; ++i){
					    sb_s2n[l] = array[INDEX_ALL(i,j1,k)];
					    sb_n2s[l] = array[INDEX_ALL(i,j2,k)];
					    ++l;
				}
			}
			MPI_Isend(sb_n2s, l, MPI_UNSIGNED, nbr_south, mpitag_nshift, this_cart, &request[0]);
			MPI_Isend(sb_s2n, l, MPI_UNSIGNED, nbr_north, mpitag_sshift, this_cart, &request[1]);

			MPI_Recv(rb_n2s, l_sn, MPI_UNSIGNED, nbr_north, mpitag_nshift, this_cart, &status[0]);
			MPI_Recv(rb_s2n, l_sn, MPI_UNSIGNED, nbr_south, mpitag_sshift, this_cart, &status[1]);

			MPI_Waitall(2, request, status);

			l = 0;
			for(int k = 0; k < kmt; ++k)
				for(int j = 0; j < num_ghosts; ++j){
					int j1 = j;
					int j2 = jphys_e+j+1;
					for(int i = iphys_b; i <= iphys_e; ++i){
					    array[INDEX_ALL(i,j1,k)] = rb_s2n[l];
					    array[INDEX_ALL(i,j2,k)] = rb_n2s[l];
					    ++l;
				}
			}
		}
	}
}

void LocalToGlobalIndex(int il, int jl, int kl, int &ig, int &jg, int &kg){
	kg = kl;
	if(comm_both){
		ig = coords[0] * cellx + il - num_ghosts;
		jg = coords[1] * celly + jl - num_ghosts;
	}
	else{
		if(comm_dir == 1){
			ig = coords[0] * cellx + il - num_ghosts;
			jg = jl;
		}
		else{
			ig = il;
			jg = coords[1] * celly + jl - num_ghosts;
		}
	}
}

void GlobalToLocalIndex(int ig, int jg, int kg, int &il, int &jl, int &kl){
	kl = kg;
	il = ig + num_ghosts - coords[0] * cellx;
	jl = jg + num_ghosts - coords[1] * celly;
}

bool CellIsGhost(int i, int j, int k){
	if(comm_both){
		if(i < iphys_b || i > iphys_e || j < jphys_b || j > jphys_e)
			return true;
		else
			return false;
	}
	else{
		if(comm_dir == 1){
			if(i < iphys_b || i > iphys_e)
				return true;
			else
				return false;
		}
		else{
			if(j < jphys_b || j > jphys_e)
				return true;
			else
				return false;
		}
	}
}

void form_pres_solver_group(unsigned int N){
	   MPI_Comm_group(MPI_COMM_WORLD, &world);
	   int gsize = 0;
	   MPI_Group_size(world, &gsize);
//	   if(!pid){
//		   printf("There are %d processes in global group \n", gsize);
//		   fflush(stdout);
//	   }
	   int rank;
	   int *masks = (int *)malloc(nproc_s*sizeof(int));
	   MPI_Group_rank(world, &rank);
	   int mask;
	   if( N > 0 )
		  mask = 1;
	   else
		  mask = 0;
	   MPI_Allgather(&mask, 1, MPI_INT, masks, 1, MPI_INT, MPI_COMM_WORLD);
	   non_pres_solver_procs = 0;
	   for(int i=0; i < nproc_s; ++i){
		   if(masks[i] == 0)
			   non_pres_solver_procs++;
	   }
	   if(!pid){
		   printf("non_pres_solver_procs = %d gsize = %d\n", non_pres_solver_procs, gsize);
		   fflush(stdout);
	   }
	   if(non_pres_solver_procs == gsize)
		   return;
//	   if( N == 0 ) return;
	   if(non_pres_solver_procs > 0){
	     int *exl = (int *)malloc(non_pres_solver_procs*sizeof(int));
	     int l = 0;
	     for(int i=0; i < nproc_s; ++i){
	    	 if(masks[i] == 0){
	    		 exl[l] = i;
	    		 ++l;
	    	 }
	     }
//	     if(!pid){
//	    	 for(int i=0; i < exl_procs;++i)
//	    		 printf(" exl = %d \n", exl[i]);
//	    	 fflush(stdout);
//	     }
	     MPI_Group_excl(world, non_pres_solver_procs, exl, &pres_solver_group);
	     free(exl);
	   }
//	   printf("B. pid = %d non_pres_solver_procs = %d gsize = %d\n", pid, non_pres_solver_procs, gsize);
//	   fflush(stdout);
	 // MPI_Group_size(pres_solver_group, &gsize);
	  if(non_pres_solver_procs > 0){
		  MPI_Comm_create(MPI_COMM_WORLD, pres_solver_group, &pres_solver_comm);
	  }
	  else{
		  pres_solver_group = world;
		  MPI_Comm_create(MPI_COMM_WORLD, pres_solver_group, &pres_solver_comm);
	  }
	  pres_solver_procs = gsize - non_pres_solver_procs;
//	  printf("C. pid = %d non_pres_solver_procs = %d gsize = %d\n", pid, non_pres_solver_procs, gsize);
//	  fflush(stdout);
	  //MPI_Comm_size(pres_solver_comm, &pres_solver_procs);

	  free(masks);
}

void finalize_mpi(){
	if(comm_both){
		delete [] f_sb_w2e;
		delete [] f_sb_e2w;
		delete [] f_rb_w2e;
		delete [] f_rb_e2w;
		delete [] f_sb_n2s;
		delete [] f_sb_s2n;
		delete [] f_rb_n2s;
		delete [] f_rb_s2n;
		delete [] f_sb_ne2sw;
		delete [] f_sb_sw2ne;
		delete [] f_rb_ne2sw;
		delete [] f_rb_sw2ne;
		delete [] f_sb_nw2se;
		delete [] f_sb_se2nw;
		delete [] f_rb_nw2se;
		delete [] f_rb_se2nw;
		delete [] c_sb_w2e;
		delete [] c_sb_e2w;
		delete [] c_rb_w2e;
		delete [] c_rb_e2w;
		delete [] c_sb_n2s;
		delete [] c_sb_s2n;
		delete [] c_rb_n2s;
		delete [] c_rb_s2n;
		delete [] c_sb_ne2sw;
		delete [] c_sb_sw2ne;
		delete [] c_rb_ne2sw;
		delete [] c_rb_sw2ne;
		delete [] c_sb_nw2se;
		delete [] c_sb_se2nw;
		delete [] c_rb_nw2se;
		delete [] c_rb_se2nw;
		delete [] ui_sb_w2e;
		delete [] ui_sb_e2w;
		delete [] ui_rb_w2e;
		delete [] ui_rb_e2w;
		delete [] ui_sb_n2s;
		delete [] ui_sb_s2n;
		delete [] ui_rb_n2s;
		delete [] ui_rb_s2n;
		delete [] ui_sb_ne2sw;
		delete [] ui_sb_sw2ne;
		delete [] ui_rb_ne2sw;
		delete [] ui_rb_sw2ne;
		delete [] ui_sb_nw2se;
		delete [] ui_sb_se2nw;
		delete [] ui_rb_nw2se;
		delete [] ui_rb_se2nw;
	}
	else{
		if(comm_dir == 1){
			delete [] f_sb_w2e;
			delete [] f_sb_e2w;
			delete [] f_rb_w2e;
			delete [] f_rb_e2w;
			delete [] c_sb_w2e;
			delete [] c_sb_e2w;
			delete [] c_rb_w2e;
			delete [] c_rb_e2w;
			delete [] ui_sb_w2e;
			delete [] ui_sb_e2w;
			delete [] ui_rb_w2e;
			delete [] ui_rb_e2w;
		}
		else{
			delete [] f_sb_n2s;
			delete [] f_sb_s2n;
			delete [] f_rb_n2s;
			delete [] f_rb_s2n;
			delete [] c_sb_n2s;
			delete [] c_sb_s2n;
			delete [] c_rb_n2s;
			delete [] c_rb_s2n;
			delete [] ui_sb_n2s;
			delete [] ui_sb_s2n;
			delete [] ui_rb_n2s;
			delete [] ui_rb_s2n;
		}
	}
}

#endif

#endif /*MPIUTILITIES_H_*/
