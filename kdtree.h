
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef KDTREE_H
#define KDTREE_H
// kdtree.h*
//#include "pbrt.h"
#include <vector>
using namespace std;
#include "geometry.h"
// KdTree Declarations
struct KdNode {
	void init(float p, u_int a) {
		splitPos = p;
		splitAxis = a;
		rightChild = ~(0U);
		hasLeftChild = 0;
	}
	void initLeaf() {
		splitAxis = 3;
		rightChild = ~(0U);
		hasLeftChild = 0;
	}
	// KdNode Data
	float splitPos;
	u_int splitAxis:2;
	u_int hasLeftChild:1;
	u_int rightChild:29;
};
template <class NodeData, class LookupProc> class KdTree {
public:
	// KdTree Public Methods
	KdTree(const vector<NodeData> &data);
//	KdTree(list<NodeData> &data);
	~KdTree() {
//		printf("free memory \n");
		FreeAligned(nodes);
//		printf("freealigned memory finished\n");
		delete[] nodeData;
//		printf("free memory finished\n");
	}
	void recursiveBuild(u_int nodeNum, int start, int end,
		vector<const NodeData *> &buildNodes);
	void Lookup(const Point &p, const LookupProc &process,
			float &maxDistSquared) const;
private:
	// KdTree Private Methods
	void privateLookup(u_int nodeNum, const Point &p,
		const LookupProc &process, float &maxDistSquared) const;
	// KdTree Private Data
	KdNode *nodes;
	NodeData *nodeData;
	u_int nNodes, nextFreeNode;	
};

template<class NodeData> struct CompareNode {
	CompareNode(int a) { axis = a; }
	int axis;
	bool operator()(const NodeData *d1,
			const NodeData *d2) const {
		return d1->p[axis] == d2->p[axis] ? (d1 < d2) :
			d1->p[axis] < d2->p[axis];
	}
};

// KdTree Method Definitions
template <class NodeData, class LookupProc>
KdTree<NodeData,
       LookupProc>::KdTree(const vector<NodeData> &d) {
	nNodes = d.size();
	nextFreeNode = 1;	
	nodes = (KdNode *)AllocAligned(nNodes * sizeof(KdNode));
	nodeData = new NodeData[nNodes];
	vector<const NodeData *> buildNodes;
//	printf("ok kdtree (1) \n");
	for (u_int i = 0; i < nNodes; ++i)
		buildNodes.push_back(&d[i]);
//	printf("ok kdtree (2) \n");	
	// Begin the KdTree building process
	recursiveBuild(0, 0, nNodes, buildNodes);
}
//template <class NodeData, class LookupProc>
//KdTree<NodeData,
//       LookupProc>::KdTree(list<NodeData> &d) {
//	nNodes = d.size();
//	nextFreeNode = 1;
//	nodes = (KdNode *)AllocAligned(nNodes * sizeof(KdNode));
////	if(nodes == NULL)
////		printf("Kdtree memory allocation error \n");
//	nodeData = new NodeData[nNodes];
//	vector<const NodeData *> buildNodes;
//	list<NodeData>::iterator iter_particle;
//	for(iter_particle = d.begin();
//		iter_particle != d.end();
//		++iter_particle){
//	//		Point pos = particles[m].Position();
//	//		float radius = particles[m].Radius();
//		NodeData &ps = *iter_particle;
//		buildNodes.push_back(&ps);
//	}
//	// Begin the KdTree building process
//	recursiveBuild(0, 0, nNodes, buildNodes);
//}
template <class NodeData, class LookupProc> void
KdTree<NodeData, LookupProc>::recursiveBuild(u_int nodeNum,
		int start, int end,
		vector<const NodeData *> &buildNodes) {
	// Create leaf node of kd-tree if we've reached the bottom	
	if (start + 1 == end) {
		nodes[nodeNum].initLeaf();
		nodeData[nodeNum] = *buildNodes[start];
		return;
	}
	// Choose split direction and partition data
	// Compute bounds of data from _start_ to _end_
	BBox bound;
	for (int i = start; i < end; ++i)
		bound = Union(bound, buildNodes[i]->p);
//	printf("Bbox = (%f, %f, %f) - (%f, %f, %f) \n",
//			bound.pMin.x, bound.pMin.y, bound.pMin.z,
//			bound.pMax.x, bound.pMax.y, bound.pMax.z);
	int splitAxis = bound.MaximumExtent();
	int splitPos = (start+end)/2;
//	printf("splitAxis = %d, splitPos = %d \n", splitAxis, splitPos);
//	printf("start = %d, end = %d \n", start, end);
	if(end == nNodes)
		std::nth_element(&buildNodes[start], &buildNodes[splitPos],
				&buildNodes[end-1], CompareNode<NodeData>(splitAxis));
	else
		std::nth_element(&buildNodes[start], &buildNodes[splitPos],
				&buildNodes[end], CompareNode<NodeData>(splitAxis));
//	printf("kdtree Z. \n");
	// Allocate kd-tree node and continue recursively
	nodes[nodeNum].init(buildNodes[splitPos]->p[splitAxis],
		splitAxis);
//	printf("kdtree A. \n");
	nodeData[nodeNum] = *buildNodes[splitPos];
//	printf("kdtree B. \n");
	if (start < splitPos) {
		nodes[nodeNum].hasLeftChild = 1;
		u_int childNum = nextFreeNode++;
		recursiveBuild(childNum, start, splitPos, buildNodes);
//		printf("kdtree C. \n");
	}
	if (splitPos+1 < end) {
		nodes[nodeNum].rightChild = nextFreeNode++;
		recursiveBuild(nodes[nodeNum].rightChild, splitPos+1,
		               end, buildNodes);
//		printf("kdtree D. \n");
	}
}
template <class NodeData, class LookupProc> void
KdTree<NodeData, LookupProc>::Lookup(const Point &p,
		const LookupProc &proc,
		float &maxDistSquared) const {
	privateLookup(0, p, proc, maxDistSquared);
}
template <class NodeData, class LookupProc> void
KdTree<NodeData, LookupProc>::privateLookup(u_int nodeNum,
		const Point &p,	const LookupProc &process,
		float &maxDistSquared) const {
	KdNode *node = &nodes[nodeNum];
	// Process kd-tree node's children
	int axis = node->splitAxis;
	if (axis != 3) {
		float dist2 = (p[axis] - node->splitPos) *
			(p[axis] - node->splitPos);
		if (p[axis] <= node->splitPos) {
			if (node->hasLeftChild)
				privateLookup(nodeNum+1, p,
				              process, maxDistSquared);
			if (dist2 < maxDistSquared &&
			    node->rightChild < nNodes)
				privateLookup(node->rightChild,
				              p,
							  process,
							  maxDistSquared);
		}
		else {
			if (node->rightChild < nNodes)
				privateLookup(node->rightChild,
				              p,
							  process,
							  maxDistSquared);
			if (dist2 < maxDistSquared && node->hasLeftChild)
				privateLookup(nodeNum+1,
				              p,
							  process,
							  maxDistSquared);
		}
	}
	// Hand kd-tree node to processing function
	float dist2 = DistanceSquared(nodeData[nodeNum].p, p);
	if (dist2 < maxDistSquared)
		process(nodeData[nodeNum], dist2, maxDistSquared);
}

#endif // PBRT_KDTREE_H
