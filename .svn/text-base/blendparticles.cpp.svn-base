#include <map>
#include <list>
using std::map;
using std::list;

#include "blendparticles.h"

#define EPS 1.e-7f

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

static float Kernel(const Point &p0, const Point &p1, float R){
	float dist = (p1 - p0).Length() / R;
	float x = 1.f - dist * dist;
	float w = max(0.f, x*x*x);
	return w;
}

/*class DtCell {
public:
  int level;          // depth in the octree of this cell
  Point corners[8];   // coordinates of the lower and upper corners of the cell
  float phi[8];    // indices into the octree's pts array of the vertices of this cell
  DtCell *parent;         // index of the parent of this cell
  DtCell *child[8];       // indices of the children of this cell
  float h;
  DtCell(float _h)
  : h(_h) {
	  parent = NULL;
	  for(int n=0; n<8; ++n)
		  child[n] = NULL;
  }
  ~DtCell(){
	for(int n=0; n < 8; ++n)
		if(child[n]!=NULL)
			delete child[n];
  	}

  void Refine(const Point &input_lc, const Point &input_uc, int level, float h, const DtCell *par,
		  list<WaterParticle> &particles, const KdTree<WaterParticle, FloatPair> &dist);

   bool contains(const Point &p); // does x lie in this cell?

};

bool DtCell::contains(const Point &x){
	if (x[0] < corners[0].x-EPS || x[0] > corners[7].x+EPS) return false;
    if (x[1] < corners[0].y-EPS || x[1] > corners[7].y+EPS) return false;
	if (x[2] < corners[0].z-EPS || x[2] > corners[7].z+EPS) return false;
	return true;
}

static void SetCorner(const Point &lc, float _h, DtCell &cell){
	cell.corners[0] = Point(lc.x, lc.y, lc.z);
	cell.corners[1] = Point(lc.x+_h, lc.y, lc.z);
	cell.corners[2] = Point(lc.x, lc.y+_h, lc.z);
	cell.corners[3] = Point(lc.x+_h, lc.y+_h, lc.z);
	cell.corners[4] = Point(lc.x, lc.y, lc.z+_h);
	cell.corners[5] = Point(lc.x+_h, lc.y, lc.z+_h);
	cell.corners[6] = Point(lc.x, lc.y+_h, lc.z+_h);
	cell.corners[7] = Point(lc.x+_h, lc.y+_h, lc.z+_h);
}

static void SetPhi(const KdTree<WaterParticle, FloatPair> &dist, float h, DtCell &cell){
	float radiusSquared = h * h;
	FloatPair fp;
	for(int n=0; n<8; ++n){
		Point p = cell.corners[n];
		dist.Lookup(p, fp, radiusSquared);
		float totalWeight = 0.f;
		float *weight = NULL;
		u_int N = 0;
		if(fp.nFoundParticles > 0){
			N = fp.Size();
			weight = new float[N];
			for(int m=0; m < N; ++m){
				WaterParticle &wps = fp.foundParticles[m];
				Point pos = wps.Position();
				float R = wps.Radius();
				weight[m] = Kernel(pos, p, R);
				totalWeight += weight[m];
			}
		}
		if(totalWeight != 0.f){
			for(int m=0; m < N; ++m)
				weight[m] /= totalWeight;
			Point posAve(0.f, 0.f, 0.f);
			float rAve = 0.f;
			for(int m=0; m < N; ++m){
				WaterParticle &wps = fp.foundParticles[m];
				Point pos = wps.Position();
//				float rad = wps.Radius()+ 0.5f * R;
				float rad = wps.Radius(); // increase particle radii to be at least
													 // 1/2 the particle spacing
				posAve += weight[m]*pos;
				rAve += weight[m]*rad;
			}
			cell.phi[n] = (p - posAve).Length() - rAve;
		}
		else
			cell.phi[n] = h;
	}
}

void DtCell::Refine(const Point &lc, const Point &uc, int l, float _h, const DtCell *par, list<WaterParticle> &particles,
			const KdTree<WaterParticle, FloatPair> &dist){
	for(int n=0; n<8; ++n){
		child[n] = new DtCell(_h);
		child[n].parent = par;
		child[n].level = l;
	}
	int n = 0;
	SetCorner(Point(lc.x, lc.y, lc.z), _h, child[n]);
	n = 1;
	SetCorner(Point(lc.x+_h, lc.y, lc.z), _h, child[n]);
	n = 2;
	SetCorner(Point(lc.x, lc.y+_h, lc.z), _h, child[n]);
	n = 3;
	SetCorner(Point(lc.x+_h, lc.y+_h, lc.z), _h, child[n]);
	n = 4;
	SetCorner(Point(lc.x, lc.y, lc.z+_h), _h, child[n]);
	n = 5;
	SetCorner(Point(lc.x+_h, lc.y, lc.z+_h), _h, child[n]);
	n = 6;
	SetCorner(Point(lc.x, lc.y+_h, lc.z+_h), _h, child[n]);
	n = 7;
	SetCorner(Point(lc.x+_h, lc.y+_h, lc.z+_h), _h, child[n]);

	for(int n=0; n<8; ++n){
		SetPhi(dist, _h, child[n]);
	}
	for(int n=0; n<8; ++n){
		list<WaterParticle>::iterator iter_particle;
		bool hasParticle = false;
		floar rmin = INFINITY;
		for(iter_particle = particles.begin();
			iter_particle != particles.end();
			++iter_particle){
			WaterParticle &wps = *iter_particle;
			Point pos = wps.Position();
			float r = wps.Radius();
			if(child[n].contains(pos)){
				hasParticle = true;
				if(rmin > r)
					rmin = r;
			}
		}
		if(hasParticle && _h > rmin)
			child[n].Refine(child[n].corners[0], child[n].corners[7], l+1, _h/2,
					&(child[n]), particles, dist);
	}

}

class DtTree {
public:
  Point lc, uc; // upper and lower coordinates of the tree
	               // (the domain is (lc,lc,lc)x(uc,uc,uc)

  int max_level; // how deep can the tree go?

  DtTree(int nx, int ny, int nz, float _h)
  : DimX(nx), DimY(ny), DimZ(nz), h(_h){
  		cells = new DtCell[nx*ny*nz](h);
  };

  ~DtTree() {
	  delete [] cells;
  };

	// functions to build the tree.  One takes a phi function, the other
	// a list of triangles & vertices.
	void buildTree(const Point &input_lc, const Point &input_uc, list<WaterParticle> &particles, const KdTree<WaterParticle, FloatPair> &dist);

private:
	int DimX, DimY, DimZ;
	float h;
	DtCell *cells;

};

void DtTree::buildTree(const Point &input_lc, const Point &input_uc, list<WaterParticle> &particles,
						const KdTree<WaterParticle, FloatPair> &dist){
	lc = input_lc;
	uc = input_uc;
	FloatPair fp;
	float radiusSquared = h * h;
	FOR_EACH_CELL
		cells[INDEX(i,j,k)].level = 0;
		cells[INDEX(i,j,k)].corners[0] = Point(lc.x+h*i, lc.y+h*j, lc.z+h*k);
		cells[INDEX(i,j,k)].corners[1] = Point(lc.x+h*(i+1), lc.y+h*j, lc.z+h*k);
		cells[INDEX(i,j,k)].corners[2] = Point(lc.x+h*i, lc.y+h*(j+1), lc.z+h*k);
		cells[INDEX(i,j,k)].corners[3] = Point(lc.x+h*(i+1), lc.y+h*(j+1), lc.z+h*k);
		cells[INDEX(i,j,k)].corners[4] = Point(lc.x+h*i, lc.y+h*j, lc.z+h*(k+1));
		cells[INDEX(i,j,k)].corners[5] = Point(lc.x+h*(i+1), lc.y+h*j, lc.z+h*(k+1));
		cells[INDEX(i,j,k)].corners[6] = Point(lc.x+h*i, lc.y+h*(j+1), lc.z+h*(k+1));
		cells[INDEX(i,j,k)].corners[7] = Point(lc.x+h*(i+1), lc.y+h*(j+1), lc.z+h*(k+1));
		for(int n=0; n<8; ++n){
			Point p = cells[INDEX(i,j,k)].corners[n];
			dist.Lookup(p, fp, radiusSquared);
			float totalWeight = 0.f;
			float *weight = NULL;
			u_int N = 0;
			if(fp.nFoundParticles > 0){
				N = fp.Size();
				weight = new float[N];
				for(int m=0; m < N; ++m){
					WaterParticle &wps = fp.foundParticles[m];
					Point pos = wps.Position();
					float R = wps.Radius();
					weight[m] = Kernel(pos, p, R);
					totalWeight += weight[m];
				}
			}
			if(totalWeight != 0.f){
				for(int m=0; m < N; ++m)
					weight[m] /= totalWeight;
				Point posAve(0.f, 0.f, 0.f);
				float rAve = 0.f;
				for(int m=0; m < N; ++m){
					WaterParticle &wps = fp.foundParticles[m];
					Point pos = wps.Position();
	//				float rad = wps.Radius()+ 0.5f * R;
					float rad = wps.Radius(); // increase particle radii to be at least
														 // 1/2 the particle spacing
					posAve += weight[m]*pos;
					rAve += weight[m]*rad;
				}
				cells[INDEX(i,j,k)].phi[n] = (p - posAve).Length() - rAve;
			}
			else
				cells[INDEX(i,j,k)].phi[n] = h;
		}
		list<WaterParticle>::iterator iter_particle;
		bool hasParticle = false;
		for(iter_particle = particles.begin();
			iter_particle != particles.end();
			++iter_particle){
			WaterParticle &wps = *iter_particle;
			Point pos = wps.Position();
			if(cells[INDEX(i,j,k)].contains(pos)){
				hasParticle = true;
				break;
			}
		}
		if(hasParticle)
			cells[INDEX(i,j,k)].Refine(cells[INDEX(i,j,k)].corners[0], cells[INDEX(i,j,k)].corners[7], level+1, h/2,
					&(cells[INDEX(i,j,k)]), particles, dist);
	END_FOR

}*/





//void FloatPair::operator()(const Particle &p,
//							float distSquare,
//							float &maxDistSquared) const {
//	foundParticles.push_back(p);
//	nFoundParticles++;
//	float dist = sqrt(distSquare);
//	if( dist <= smallestR){
//		smallestR = dist;
//		value = p.EvaluatePhiP(smallestR);
//	}
//	return;
//}

void BlendParticles::AddWaterParticles(const WaterParticle &p){
	waterParticles.push_back(p);
}

//void BlendParticles::EvaluatePhi(const Voxel &voxel){
//
//	KdTree<Particle, FloatPair> dist(waterParticles);
//	FloatPair fp;
//	float R = voxel.VoxelDelta();
//	float radiusSquared = R * R;
//	FOR_EACH_CELL
//		Point p = voxel.VoxelCenterPosition(i,j,k);
//		dist.Lookup(p, fp, radiusSquared);
//		if(fp.nFoundParticles > 0){
//			float *weight = new float[fp.Size()];
//			for(int i=0; i < fp.Size(); ++i){
//				Particle &particle = fp.foundParticles[i];
//				Point pos = particle.Position();
//				weight[i] = Kernel(pos, p, R);
//			}
//			float totalWeight = 0.f;
//			for(int i=0; i < fp.Size(); ++i)
//				totalWeight += weight[i];
//			if(totalWeight != 0.f){
//				for(int i=0; i < fp.Size(); ++i)
//					weight[i] /= totalWeight;
//				Point posAve(0.f, 0.f, 0.f);
//				float rAve = 0.f;
//				for(int i=0; i < fp.Size(); ++i){
//					Particle &particle = fp.foundParticles[i];
//					Point pos = particle.Position();
//					float rad = particle.Radius();
//					posAve += weight[i]*pos;
//					rAve += weight[i]*rad;
//				}
//				phi[INDEX(i,j,k)] = (p - posAve).Length() - rAve;
//			}
//			else{  // what if all weights = 0?
//
//			}
//			delete [] weight;
//		}
//		else{ // what if there are no particles nearby?
//
//		}
//	END_FOR
//}

void BlendParticles::EvaluatePhi(const Voxel &voxel){

//	printf("got here (1)\n");
//	list<WaterParticle>::iterator iter_particle;
	WaterParticleList::iterator iter_particle;
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
	max_vel = sqrt(max_vel);
	min_vel = sqrt(min_vel);
	printf("ratio max_vel/min_vel = %f \n", max_vel/min_vel);
//	printf("got here (1.2)\n");
//	KdTree<WaterParticle, FloatPair> dist(waterParticles);
	if(wp.size() == 0){
		FOR_EACH_CELL
			phi[INDEX(i,j,k)] = threshold;
		END_FOR
	}
	else{
	KdTree<WaterParticle, FloatPair> dist(wp);
//	printf("got here (1.3)\n");
	FloatPair fp;
	float R = voxel.VoxelDelta();
	float radiusSquared = R * R;
//	printf("got here (2)\n");
	FOR_EACH_CELL
//	if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
		Point p = voxel.VoxelCenterPosition(i,j,k);
	    dist.Lookup(p, fp, radiusSquared);
	    float totalWeight = 0.f;
	    float *weight = NULL;
	    u_int N = 0;
//	    printf("got here (3) %d, %d, %d\n", i,j,k);
	    if(fp.nFoundParticles > 0){
	    	N = fp.Size();
	    	weight = new float[N];
			for(int n=0; n < N; ++n){
				WaterParticle &wps = fp.foundParticles[n];
				Point pos = wps.Position();
				weight[n] = Kernel(pos, p, R);
				totalWeight += weight[n];
			}
		}
//	    printf("got here (4) %d, %d, %d\n", i,j,k);
	    if(totalWeight != 0.f){
			for(int n=0; n < N; ++n)
				weight[n] /= totalWeight;
			Point posAve(0.f, 0.f, 0.f);
			float rAve = 0.f;
			for(int n=0; n < N; ++n){
				WaterParticle &wps = fp.foundParticles[n];
				Point pos = wps.Position();
				float rad = wps.Radius();
//				float rad = wps.Radius()+ 0.25f * R; // increase particle radii to be at least
				                                     // 1/2 the particle spacing
				posAve += weight[n]*pos;
				rAve += weight[n]*rad;
			}
			phi[INDEX(i,j,k)] = (p - posAve).Length() - rAve;
	    }
		else
			phi[INDEX(i,j,k)] = threshold;
//	    printf("got here (5) %d, %d, %d\n", i,j,k);
	    if(weight != NULL)
	    	delete [] weight;
	    fp.Clear();
//	    printf("got here (6) %d, %d, %d\n", i,j,k);
//		printf(" at (%d, %d, %d) phi = %f\n", i,j,k,phi[INDEX(i,j,k)]);
//	}
//	else
//		phi[INDEX(i,j,k)] = threshold;
	END_FOR
	}

//	ReInitialize(10, voxel, phi);

	u_int NM = 0;
	FOR_EACH_CELL
		if(phi[INDEX(i,j,k)] < 0.f)
			++NM;
	END_FOR
	printf("\nThere are %u points in water spray \n\n", NM);
	return;

}

/*void BlendParticles::EvaluatePhi(const Voxel &voxel){

//	printf("got here (1)\n");
	list<WaterParticle>::iterator iter_particle;
	vector<WaterParticle> wp;
	float max_vel = 0.f;
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
				wp.push_back(wps);
				++iter_particle;
			}
	}
	max_vel = sqrt(max_vel);

//	KdTree<WaterParticle, FloatPair> dist(waterParticles);
	KdTree<WaterParticle, FloatPair> dist(wp);
	FloatPair fp;
	float R = voxel.VoxelDelta();
	float radiusSquared = R * R;
//	printf("got here (2)\n");
	FOR_EACH_CELL
		Point p = voxel.VoxelCenterPosition(i,j,k);
	    float sum = 0.f;
	    dist.Lookup(p, fp, radiusSquared);
	    if(fp.nFoundParticles > 0){
//		for(iter_particle = waterParticles.begin();
//			iter_particle != waterParticles.end();
//			){
			for(int n=0; n < fp.Size(); ++n){
				WaterParticle &wps = fp.foundParticles[n];
	//				WaterParticle &wps = *iter_particle;
				Vector vel = wps.Velocity();
				Point pos = wps.Position();
				float r = wps.Radius();
				float mag_vel = vel.Length();
				float dpara, dperp;
	//			printf(" u = %f, v = %f, w = %f \n",vel.x, vel.y, vel.z);
				dpara = AbsDot((p-pos), Normalize(vel));
				Vector v2 = (p - pos) - dpara * Normalize(vel);
				dperp = v2.Length();
				float t = mag_vel/(2*max_vel);
				float d = (1.f - t) * dpara + (1.f + t) * dperp;
				float dd = d / r;
				sum += 1.f/(1.f+ dd * dd);
			}
//			++iter_particle;
		}
	    fp.Clear();
		phi[INDEX(i,j,k)] = threshold + sum;
//		printf(" at (%d, %d, %d) phi = %f\n", i,j,k,phi[INDEX(i,j,k)]);
	END_FOR
	FOR_EACH_CELL
		phi[INDEX(i,j,k)] = -1 * phi[INDEX(i,j,k)];
	END_FOR
//	printf("got here (3)\n");
	return;
}*/


void BlendParticles::OutputBinaryGridData(const Voxel &voxel, char *filename) const{

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

void BlendParticles::OutputWaterParticles(const Voxel &voxel, char *name) {
//	list<WaterParticle>::iterator iter_particle;
	WaterParticleList::iterator iter_particle;
	FILE *filename = fopen(name, "w");
	if(!filename){
      printf("Couldn't open file [ %s ] to write \n", name);
      exit(-1);
   	}
	for(iter_particle = waterParticles.begin();
		iter_particle != waterParticles.end();
		){
		WaterParticle &wps = *iter_particle;
		Vector vel = wps.Velocity();
//		float mag_vel = vel.LengthSquared();
		float r = wps.Radius();
		Point pos = wps.Position();
//		if(mag_vel == 0.f)
//			iter_particle = waterParticles.erase(iter_particle);
//		else{
//			vel = Normalize(vel);
			fprintf(filename, "AttributeBegin \n");
			fprintf(filename, "Translate %f  %f  %f \n", pos.x, pos.y, pos.z);
//			fprintf(filename, "Scale %f  %f  %f \n", vel.x, vel.y, vel.z);
			fprintf(filename, "Shape \"sphere\" \"float radius\"  %f \n", r);
			fprintf(filename, "AttributeEnd \n");
			++iter_particle;
//		}
	}
	fclose(filename);
}

void BlendParticles::OutputWaterParticlesBinary(const Voxel &voxel, char *name) {

	float delta = voxel.VoxelDelta();

//	list<WaterParticle>::iterator iter_particle;
	WaterParticleList::iterator iter_particle;

	FILE *fp = fopen(name, "wb");
	if(!fp){
      printf("Couldn't open file [ %s ] to write \n", name);
      exit(-1);
   	}
	u_int num = Size();
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

    fwrite(&delta, sizeof(float), 1, fp);
	fwrite(&num, sizeof(u_int), 1, fp);
	for(iter_particle = waterParticles.begin();
		iter_particle != waterParticles.end();
		){
		WaterParticle &wps = *iter_particle;
		float r = wps.Radius();
		Point pos = wps.Position();
		Vector vel = wps.Velocity();
		fwrite(&(pos.x), sizeof(float), 1, fp);
		fwrite(&(pos.y), sizeof(float), 1, fp);
		fwrite(&(pos.z), sizeof(float), 1, fp);
		fwrite(&(vel.x), sizeof(float), 1, fp);
		fwrite(&(vel.y), sizeof(float), 1, fp);
		fwrite(&(vel.z), sizeof(float), 1, fp);
		fwrite(&r, sizeof(float), 1, fp);
		++iter_particle;
	}
	fclose(fp);
}

void BlendParticles::MergeWith(const Voxel &voxel, float *phi0) const{
	FOR_EACH_CELL
		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k)){
//			printf("spray phi < 0 at (%d, %d, %d)\n", i,j,k);
			phi0[INDEX(i,j,k)] = phi0[INDEX(i,j,k)] < phi[INDEX(i,j,k)] ?
								 phi0[INDEX(i,j,k)] : phi[INDEX(i,j,k)];
		}
	END_FOR
}

void BlendParticles::ReInitialize(int steps, const Voxel &voxel, float *d){

	float *value = new float[DimX*DimY*DimZ];
	FOR_EACH_CELL
		value[INDEX(i,j,k)] = d[INDEX(i,j,k)];
	END_FOR

	float delta = voxel.VoxelDelta();
	float dx = voxel.VoxelDelta();
	float dy = voxel.VoxelDelta();
	float dz = voxel.VoxelDelta();
	float dt = .1f;

	int ii, jj, kk;
	u_int N=0;

	float old_phi;

	FOR_EACH_CELL

//         printf("at (%d, %d, %d) phi = %f \n", i,j,k,value[INDEX(i,j,k)]);
//		if(!voxel.InSolid(i,j,k) && !voxel.InSource(i,j,k) && fabs(value[INDEX(i,j,k)]) < 5*delta){
//			TentativeGridPoint tp(0.f, i,j,k);
//		 if(!voxel.CloseToSolid(tp)){
			N++;
			if(i==DimX-1) ii = i - 1;
			else ii = i;
			int posxp1 = k*DimX*DimY+j*DimX+ii+1;
			if(i==0) ii = i + 1;
			else ii = i;
			int posxm1 = k*DimX*DimY+j*DimX+ii-1;

			if(j==DimY-1) jj = j - 1;
			else jj = j;
			int posyp1 = k*DimX*DimY+(jj+1)*DimX+i;
			if(j==0) jj = j + 1;
			else jj = j;
			int posym1 = k*DimX*DimY+(jj-1)*DimX+i;

			if(k==DimZ-1) kk = k - 1;
			else kk = k;
			int poszp1 = (kk+1)*DimX*DimY+j*DimX+i;
			if(k==0) kk = k + 1;
			else kk = k;
			int poszm1 = (kk-1)*DimX*DimY+j*DimX+i;

			float t = sqrt(value[INDEX(i,j,k)]*value[INDEX(i,j,k)] + delta*delta);
			float s = value[INDEX(i,j,k)] / t;

			dt = 3.f / 4 * delta / fabs(s); // 3/4 * (CFL allowed)

			old_phi = value[INDEX(i,j,k)];

			for(int m=0; m<steps; m++){

				float phixp = (value[posxp1]-value[INDEX(i,j,k)]) / dx;
				float phixm = (value[INDEX(i,j,k)]-value[posxm1]) / dx;
				float phiyp = (value[posyp1]-value[INDEX(i,j,k)]) / dy;
				float phiym = (value[INDEX(i,j,k)]-value[posym1]) / dy;
				float phizp = (value[poszp1]-value[INDEX(i,j,k)]) / dz;
				float phizm = (value[INDEX(i,j,k)]-value[poszm1]) / dz;
				float phix2, phiy2, phiz2, H;

				if(s > 0) {
					phix2 = max( max(phixm,0.f)*max(phixm,0.f),
								 min(phixp,0.f)*min(phixp,0.f) );
					phiy2 = max( max(phiym,0.f)*max(phiym,0.f),
								 min(phiyp,0.f)*min(phiyp,0.f) );
					phiz2 = max( max(phizm,0.f)*max(phizm,0.f),
								 min(phizp,0.f)*min(phizp,0.f) );
					H = s*sqrt(phix2+phiy2+phiz2);
				}
				else if (s < 0){
					phix2 = max( min(phixm,0.f)*min(phixm,0.f),
								 max(phixp,0.f)*max(phixp,0.f) );
					phiy2 = max( min(phiym,0.f)*min(phiym,0.f),
								 max(phiyp,0.f)*max(phiyp,0.f) );
					phiz2 = max( min(phizm,0.f)*min(phizm,0.f),
								 max(phizp,0.f)*max(phizp,0.f) );
					H = s*sqrt(phix2+phiy2+phiz2);
				}
				else
					H = 0.f;
//				if(i == I & j == J && k == K){
//					printf("i=%d, j=%d, k=%d, m = %d, dt = %f, phi = %f, H = %f, S = %f, dt*(-H+s) = %f, grad(phi) = %f \n",
//								i,j, k, m, dt, value[pos], H, s, dt*(-H+s), sqrt(phix2+phiy2+phiz2));
//					printf("i=%d, j=%d, k=%d, m = %d, phixp = %f, phixm = %f, phiyp = %f, phiym = %f, phizp = %f, phizm = %f \n",
//										i,j, k, m, phixp, phixm, phiyp, phiym, phizp, phizm );
//					printf("i=%d, j=%d, k=%d, m = %d, phixp = %f, phixm = %f, phiyp = %f, phiym = %f, phizp = %f, phizm = %f \n",
//							i,j, k, m, value[posxp1],value[posxm1],value[posyp1],value[posym1],value[poszp1],value[poszm1]);
//					//print_grad_phi(I,J,K);
//				}
				d[INDEX(i,j,k)] += dt*(-H+s);
			}
//			if(old_phi*value[pos]< 0.f)
//				printf(" (reinit) phi changes sign at i=%d, j=%d, k=%d, old = %f, new = %f \n",
//						i, j, k, old_phi, value[pos]);
//		 }
//		}
	END_FOR
	delete [] value;
	printf(" There are %d points in total got reinitialized \n", N);


//	for(int i=0; i<DimX*DimY*DimZ; ++i){
//			phi[i] = phi_tmp[i];
//	}
//	delete [] phi_tmp;

}
