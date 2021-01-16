/*
 * octree.h
 *
 *  Created on: Sep 2, 2008
 *      Author: bzhao
 */

#ifndef OCTREE_H_
#define OCTREE_H_

#define MAX_LEVELS 8

struct OcNode{
	OcNode(float _h, unsigned char l)
	: level(l), h(_h){
		for(int n=0; n < 8; ++n)
			child[n] = NULL;
	}
	~OcNode(){
		for(int n=0; n < 8; ++n)
			if(child[n]!=NULL)
				delete child[n];
	}
	bool IsLeaf() const{
		bool isleaf = true;
		for(int n=0; n < 8; ++n)
			if(child[n]!=NULL){
				isleaf = false;
				break;
			}
		return isleaf;
	}
	void Refine(unsigned char l){
		if(level < l){
			for(int n=0; n < 8; ++n){
				child[n] = new OcNode(h/2, ++level);
				child[n].Refine(l);
			}
		}
	}
	void Coarsen(){
		if(!IsLeaf()){
			float sum = 0.f;
			for(int n=0; n < 8; ++n){
				phi[n] = child[n].phi[n];
				sum += child[n].p;
//				phi[0] = child[0].phi[0];
//				phi[1] = child[1].phi[1];
//				phi[2] = child[2].phi[2];
//				phi[3] = child[3].phi[3];
			}
			p = sum / 8;
			u = (child[1].u + child[2].u + child[5].u + child[6].u)/4;
			v = (child[2].v + child[3].v + child[6].v + child[7].v)/4;
			w = (child[4].w + child[5].w + child[6].w + child[7].w)/4;
			for(int n=0; n < 8; ++n){
				delete [] child;
				child = NULL;
			}
		}
	}

	float phi[8];
	float p;
	float u, v, w;
	unsigned char level;
	float h;
	OcNode *child[8];
};

class OcTree{
public:
	OcTree(int nx, int ny, int nz, float _h):
	DimX(nx), DimY(ny), DimZ(nz), h(_h){
		cells = new OcNode[nx*ny*nz](h, 0);
	}
	~OcTree(){
		delete [] cells;
	}
private:
	int DimX, DimY, DimZ;
	float h;
	OcNode *cells;
};


#endif /* OCTREE_H_ */
