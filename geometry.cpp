
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// geometry.cpp*
#include "geometry.h"
// BBox Method Definitions
COREDLL BBox Union(const BBox &b, const Point &p) {
	BBox ret = b;
	ret.pMin.x = min(b.pMin.x, p.x);
	ret.pMin.y = min(b.pMin.y, p.y);
	ret.pMin.z = min(b.pMin.z, p.z);
	ret.pMax.x = max(b.pMax.x, p.x);
	ret.pMax.y = max(b.pMax.y, p.y);
	ret.pMax.z = max(b.pMax.z, p.z);
	return ret;
}
COREDLL BBox Union(const BBox &b, const BBox &b2) {
	BBox ret;
	ret.pMin.x = min(b.pMin.x, b2.pMin.x);
	ret.pMin.y = min(b.pMin.y, b2.pMin.y);
	ret.pMin.z = min(b.pMin.z, b2.pMin.z);
	ret.pMax.x = max(b.pMax.x, b2.pMax.x);
	ret.pMax.y = max(b.pMax.y, b2.pMax.y);
	ret.pMax.z = max(b.pMax.z, b2.pMax.z);
	return ret;
}
void BBox::BoundingSphere(Point *c, float *rad) const {
	*c = .5f * pMin + .5f * pMax;
	*rad = Inside(*c) ? Distance(*c, pMax) : 0.f;
}
