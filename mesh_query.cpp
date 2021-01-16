#include <cfloat>
#include "mesh_query.h"
#include "bounding_box_tree.h"
#include "predicates.h"

//==============================================================================
// helper functions for ray casts along positive z axis

// returns true if the ray cast from p along positive z axis hits box
static bool
box_zcast(const Vec3d& p,
          const BoundingBox& box)
{
    return p[0]>=box.xmin[0] && p[0]<=box.xmax[0]
        && p[1]>=box.xmin[1] && p[1]<=box.xmax[1]
        && p[2]<=box.xmax[2];
}

// helper function for tri_zcast below...
// Here we are given a robust 2d orientation of qrs in xy projection, which
// must be positive.
static bool
tri_zcast_inner(const Vec3d& p,
                double qrs,
                const Vec3d& q,
                const Vec3d& r,
                const Vec3d& s)
{
    assert(qrs>0);

    // first check if point is above or below triangle in z
    double pqrs=orient3d(p.v, q.v, r.v, s.v);
    if(pqrs>=0) return false; // point is on or above triangle - no intersection

    // then check if point lies outside triangle in 2D xy projection
    double pqr=orient2d(p.v, q.v, r.v);
    if(pqr<0) return false;
    double prs=orient2d(p.v, r.v, s.v);
    if(prs<0) return false;
    double psq=orient2d(p.v, s.v, q.v);
    if(psq<0) return false;

    // note: the following tests are somewhat redundant, but it's a pretty
    // tiny optimization to eliminate the redundancy compared to the loss in
    // clarity.

    // check if point is strictly inside the triangle in xy
    if(pqr>0 && prs>0 && psq>0) return true;

    // check if point is strictly on edge qr
    if(pqr==0 && prs>0 && psq>0){
        if(q[1]<r[1]) return false;
        if(q[1]>r[1]) return true;
        if(q[0]<r[0]) return true;
        assert(q[0]>r[0]); // q!=r because triangle is not degenerate
        return false;
    }

    // check if point is strictly on edge rs
    if(prs==0 && pqr>0 && psq>0){
        if(r[1]<s[1]) return false;
        if(r[1]>s[1]) return true;
        if(r[0]<s[0]) return true;
        assert(r[0]>s[0]); // r!=s because triangle is not degenerate
        return false;
    }

    // check if point is strictly on edge sq
    if(psq==0 && pqr>0 && prs>0){
        if(s[1]<q[1]) return false;
        if(s[1]>q[1]) return true;
        if(s[0]<q[0]) return true;
        assert(s[0]>q[0]); // r!=s because triangle is not degenerate
        return false;
    }

    // check if point is on vertex q
    if(p[0]==q[0] && p[1]==q[1]){
        return q[1]>=r[1] && q[1]<s[1];
    }

    // check if point is on vertex r
    if(p[0]==r[0] && p[1]==r[1]){
        return r[1]>=s[1] && r[1]<q[1];
    }

    // check if point is on vertex s
    if(p[0]==s[0] && p[1]==s[1]){
        return s[1]>=q[1] && s[1]<r[1];
    }

    assert(false); // we should have covered all cases at this point
    return false; // just to quiet compiler warnings
}

// returns true if the ray cast from p along positive z axis hits triangle in
// exactly one spot, with edge cases handled appropriately
static bool
tri_zcast(const Vec3d& p,
          const Vec3d& q,
          const Vec3d& r,
          const Vec3d& s)
{
    // robustly find orientation of qrs in 2D xy projection
    double qrs=orient2d(q.v, r.v, s.v);
    if(qrs>0)
        return tri_zcast_inner(p, qrs, q, r, s);
    else if(qrs<0)
        return tri_zcast_inner(p, -qrs, q, s, r); // flip triangle to reorient
    else
        return false; // triangle is degenerate in 2D projection - ignore
}

//==============================================================================
// helper functions for segment intersecting triangles

// returns false only if the segment for sure can't intersect the box
static bool
segment_box_intersect(const Vec3d& p,
                      const Vec3d& q,
                      const BoundingBox& box)
{
    // these are conservative bounds on error factor from rounding
    const double lo=1-5*DBL_EPSILON, hi=1+5*DBL_EPSILON;
    double s=0, t=1; // bounds on parameter for intersection interval
    for(unsigned int i=0; i<3; ++i){
        if(p[i]<q[i]){
            double d=q[i]-p[i];
            double s0=lo*(box.xmin[i]-p[i])/d, t0=hi*(box.xmax[i]-p[i])/d;
            if(s0>s) s=s0;
            if(t0<t) t=t0;
        }else if(p[i]>q[i]){
            double d=q[i]-p[i];
            double s0=lo*(box.xmax[i]-p[i])/d, t0=hi*(box.xmin[i]-p[i])/d;
            if(s0>s) s=s0;
            if(t0<t) t=t0;
        }else{
            if(p[i]<box.xmin[i] || p[i]>box.xmax[i]) return false;
        }
        if(s>t) return false;
    }
    return true;
}

static bool
ray_box_intersect(const Vec3d& P,
				  const Vec3d& d,
				  const BoundingBox& box)
{
	double tNear = -DBL_MAX;
	double tFar  = DBL_MAX;

    for(unsigned int i=0; i<3; ++i){
    	if(fabs(d[i]) < DBL_EPSILON){
    		if(P[i] < box.xmin[i] || P[i] > box.xmax[i])
    			return false;
    		continue;
    	}
    	// The line is not parallel to the slab planes, so compute the
	    // parameters of intersection.  The line segment of intersection is
	    // P+t*d, where t0 <= t <= t1.
    	double t0 = (box.xmin[i] - P[i]) / d[i];
		double t1 = (box.xmax[i] - P[i]) / d[i];
		double tmp;
		if (t0 > t1){
		  tmp = t1;
		  t1 = t0;
		  t0 = tmp;
	     }
		 // Compare with current values.  The current slab may have increased
		 // tNear and/or decreased tFar.
		 if (t0 > tNear){
			  tNear = t0;
		 }
		 if (t1 < tFar){
			  tFar = t1;
		 }
		 // Check if the line misses the AABB entirely.
		 if (tNear > tFar){
			  return false;
		  }
	  }
      // The ray is (P+t*d, t >= 0).  We need to compute the intersection of
      // the interval tNear <= t <= tFar with the interval t >= 0.
      if (tFar < 0.0)
	  {
		  return false;
	  }
	  return true;

}

static bool
sphere_box_intersect(const Vec3d& p,
                     float d,
                     const BoundingBox& box)
{

	double distSquared = 0.0;

	for(unsigned int i=0; i<3; ++i){
	  if(p[i] < box.xmin[i])
		  distSquared += sqr(p[i]-box.xmin[i]);
	  else if(p[i] > box.xmax[i])
		  distSquared += sqr(p[i]-box.xmax[i]);
	}
	if(distSquared <= d * d)
		return true;
	else
		return false;
}


// determine if segment pq intersects triangle uvw, setting approximate
// barycentric coordinates if so.
static bool
segment_tri_intersect(const Vec3d& p,
                      const Vec3d& q,
                      const Vec3d& u,
                      const Vec3d& v,
                      const Vec3d& w,
                      double* s,
                      double* t,
                      double* a,
                      double* b,
                      double* c)
{
    // find where segment hits plane of triangle
    double puvw=orient3d(p.v, u.v, v.v, w.v),
           uvwq=orient3d(u.v, v.v, w.v, q.v);
    if((puvw<=0 && uvwq>=0) || (puvw>=0 && uvwq<=0))
        return false; // either no intersection, or a degenerate one
    if(puvw<0 || uvwq<0){
        double pqvw=orient3d(p.v, q.v, v.v, w.v);
        if(pqvw>0) return false;
        double puqw=orient3d(p.v, u.v, q.v, w.v);
        if(puqw>0) return false;
        double puvq=orient3d(p.v, u.v, v.v, q.v);
        if(puvq>0) return false;
        *s=uvwq/(puvw+uvwq);
        *t=puvw/(puvw+uvwq);
        *a=pqvw/(pqvw+puqw+puvq);
        *b=puqw/(pqvw+puqw+puvq);
        *c=puvq/(pqvw+puqw+puvq);
        return true;
    }else{ //(puvw>0 || uvwq>0)
        double pqvw=orient3d(p.v, q.v, v.v, w.v);
        if(pqvw<0) return false;
        double puqw=orient3d(p.v, u.v, q.v, w.v);
        if(puqw<0) return false;
        double puvq=orient3d(p.v, u.v, v.v, q.v);
        if(puvq<0) return false;
        *s=uvwq/(puvw+uvwq);
        *t=puvw/(puvw+uvwq);
        *a=pqvw/(pqvw+puqw+puvq);
        *b=puqw/(pqvw+puqw+puvq);
        *c=puvq/(pqvw+puqw+puvq);
        return true;
    }
}

static bool
ray_tri_intersect(const Vec3d& P,
				  const Vec3d& d,
				  const Vec3d& v0,
				  const Vec3d& v1,
				  const Vec3d& v2,
				  double* t,
				  double* a,
				  double* b,
				  double* c)
{
	Vec3d e1 = v1 - v0;
	Vec3d e2 = v2 - v1;
	Vec3d p = cross(d, e2);
    double tmp = dot(p, e1);
    if (tmp > -DBL_EPSILON && tmp < DBL_EPSILON) {
    	return false;
    }
    tmp = 1.0 / tmp;
    Vec3d s = P - v0;
    *a = tmp * dot(s, p);
    if (*a < 0.0 || *a > 1.0) {
    	return false;
    }
    Vec3d q = cross(s, e1);
    *b = tmp * dot(d, q);
    if (*b < 0.0 || *b > 1.0) {
    	return false;
    }
    if ( (*a) + (*b) > 1.0 )
         return false;
    *t = tmp * dot(e2, q);
    if( (*t) < 0.0)
    	return false;
    *c = 1.0 - (*a) - (*b);
    return true;
}

typedef double Real;

static double
point_tri_distance(const Vec3d& p,
				  const Vec3d& u,
				  const Vec3d& v,
				  const Vec3d& w,
				  double* barya,
				  double* baryb,
				  double* baryc)
{

	    Vec3d diff =  u - p;
	    Vec3d edge0 = v - u;
	    Vec3d edge1 = w - u;
	    Real a00 = mag2(edge0);
	    Real a01 = dot(edge0, edge1);
	    Real a11 = mag2(edge1);
	    Real b0 = dot(diff, edge0);
	    Real b1 = dot(diff, edge1);
	    Real c = mag2(diff);
	    Real det = fabs(a00*a11 - a01*a01);
	    Real s = a01*b1 - a11*b0;
	    Real t = a01*b0 - a00*b1;
	    Real sqrDistance;

	    if (s + t <= det)
	    {
	        if (s < (Real)0)
	        {
	            if (t < (Real)0)  // region 4
	            {
	                if (b0 < (Real)0)
	                {
	                    t = (Real)0;
	                    if (-b0 >= a00)
	                    {
	                        s = (Real)1;
	                        sqrDistance = a00 + ((Real)2)*b0 + c;
	                    }
	                    else
	                    {
	                        s = -b0/a00;
	                        sqrDistance = b0*s + c;
	                    }
	                }
	                else
	                {
	                    s = (Real)0;
	                    if (b1 >= (Real)0)
	                    {
	                        t = (Real)0;
	                        sqrDistance = c;
	                    }
	                    else if (-b1 >= a11)
	                    {
	                        t = (Real)1;
	                        sqrDistance = a11 + ((Real)2)*b1 + c;
	                    }
	                    else
	                    {
	                        t = -b1/a11;
	                        sqrDistance = b1*t + c;
	                    }
	                }
	            }
	            else  // region 3
	            {
	                s = (Real)0;
	                if (b1 >= (Real)0)
	                {
	                    t = (Real)0;
	                    sqrDistance = c;
	                }
	                else if (-b1 >= a11)
	                {
	                    t = (Real)1;
	                    sqrDistance = a11 + ((Real)2)*b1 + c;
	                }
	                else
	                {
	                    t = -b1/a11;
	                    sqrDistance = b1*t + c;
	                }
	            }
	        }
	        else if (t < (Real)0)  // region 5
	        {
	            t = (Real)0;
	            if (b0 >= (Real)0)
	            {
	                s = (Real)0;
	                sqrDistance = c;
	            }
	            else if (-b0 >= a00)
	            {
	                s = (Real)1;
	                sqrDistance = a00 + ((Real)2)*b0 + c;
	            }
	            else
	            {
	                s = -b0/a00;
	                sqrDistance = b0*s + c;
	            }
	        }
	        else  // region 0
	        {
	            // minimum at interior point
	            Real invDet = ((Real)1)/det;
	            s *= invDet;
	            t *= invDet;
	            sqrDistance = s*(a00*s + a01*t + ((Real)2)*b0) +
	                t*(a01*s + a11*t + ((Real)2)*b1) + c;
	        }
	    }
	    else
	    {
	        Real tmp0, tmp1, numer, denom;

	        if (s < (Real)0)  // region 2
	        {
	            tmp0 = a01 + b0;
	            tmp1 = a11 + b1;
	            if (tmp1 > tmp0)
	            {
	                numer = tmp1 - tmp0;
	                denom = a00 - ((Real)2)*a01 + a11;
	                if (numer >= denom)
	                {
	                    s = (Real)1;
	                    t = (Real)0;
	                    sqrDistance = a00 + ((Real)2)*b0 + c;
	                }
	                else
	                {
	                    s = numer/denom;
	                    t = (Real)1 - s;
	                    sqrDistance = s*(a00*s + a01*t + ((Real)2)*b0) +
	                        t*(a01*s + a11*t + ((Real)2)*b1) + c;
	                }
	            }
	            else
	            {
	                s = (Real)0;
	                if (tmp1 <= (Real)0)
	                {
	                    t = (Real)1;
	                    sqrDistance = a11 + ((Real)2)*b1 + c;
	                }
	                else if (b1 >= (Real)0)
	                {
	                    t = (Real)0;
	                    sqrDistance = c;
	                }
	                else
	                {
	                    t = -b1/a11;
	                    sqrDistance = b1*t + c;
	                }
	            }
	        }
	        else if (t < (Real)0)  // region 6
	        {
	            tmp0 = a01 + b1;
	            tmp1 = a00 + b0;
	            if (tmp1 > tmp0)
	            {
	                numer = tmp1 - tmp0;
	                denom = a00 - ((Real)2)*a01 + a11;
	                if (numer >= denom)
	                {
	                    t = (Real)1;
	                    s = (Real)0;
	                    sqrDistance = a11 + ((Real)2)*b1 + c;
	                }
	                else
	                {
	                    t = numer/denom;
	                    s = (Real)1 - t;
	                    sqrDistance = s*(a00*s + a01*t + ((Real)2)*b0) +
	                        t*(a01*s + a11*t + ((Real)2)*b1) + c;
	                }
	            }
	            else
	            {
	                t = (Real)0;
	                if (tmp1 <= (Real)0)
	                {
	                    s = (Real)1;
	                    sqrDistance = a00 + ((Real)2)*b0 + c;
	                }
	                else if (b0 >= (Real)0)
	                {
	                    s = (Real)0;
	                    sqrDistance = c;
	                }
	                else
	                {
	                    s = -b0/a00;
	                    sqrDistance = b0*s + c;
	                }
	            }
	        }
	        else  // region 1
	        {
	            numer = a11 + b1 - a01 - b0;
	            if (numer <= (Real)0)
	            {
	                s = (Real)0;
	                t = (Real)1;
	                sqrDistance = a11 + ((Real)2)*b1 + c;
	            }
	            else
	            {
	                denom = a00 - ((Real)2)*a01 + a11;
	                if (numer >= denom)
	                {
	                    s = (Real)1;
	                    t = (Real)0;
	                    sqrDistance = a00 + ((Real)2)*b0 + c;
	                }
	                else
	                {
	                    s = numer/denom;
	                    t = (Real)1 - s;
	                    sqrDistance = s*(a00*s + a01*t + ((Real)2)*b0) +
	                        t*(a01*s + a11*t + ((Real)2)*b1) + c;
	                }
	            }
	        }
	    }

	    // Account for numerical round-off error.
	    if (sqrDistance < (Real)0)
	    {
	        sqrDistance = (Real)0;
	    }

	    *barya = s;
	    *baryb = t;
	    *baryc = (Real)1 - s - t;
	    return sqrt(sqrDistance);

}

//==============================================================================
// the actual accelerated mesh class

struct MeshObject
{
    int n;
    const Vec3d *x;
    const Vec3d *N;
    int nt;
    const Vec3i *tri;
    BoundingBoxTree tree;

    MeshObject(int n_,
               const double *x_,
               int nt_,
               const int *tri_,
               const double *y_=NULL)
        : n(n_), x((const Vec3d*)x_), N((const Vec3d*)y_), nt(nt_), tri((const Vec3i*)tri_)
    {
        assert(x && tri && n>=0 && nt>=0);
        if(nt==0) return;
        std::vector<BoundingBox> box(nt);
        for(int t=0; t<nt; ++t){
            int i, j, k; assign(tri[t], i, j, k);
            box[t].build_from_points(x[i], x[j], x[k]);
        }
        tree.construct_from_leaf_boxes(nt, &box[0]);
    }

    bool
    inside(const double point[3]) const
    {
        Vec3d p(point);
        // quick check on root bounding box
        if(!tree.box.contains(p)) return false;
        // count intersections along a ray-cast, check parity for inside/outside
        int intersection_count=0;
        // we cast ray along positive z axis
        std::vector<const BoundingBoxTree*> stack;
        stack.push_back(&tree);
        while(!stack.empty()){
            const BoundingBoxTree *node=stack.back();
            stack.pop_back();
            // check any triangles in this node
            for(unsigned int i=0; i<node->index.size(); ++i){
                int t=node->index[i];
                if(tri_zcast(p, x[tri[t][0]], x[tri[t][1]], x[tri[t][2]]))
                    ++intersection_count;
            }
            // check any subtrees for this node
            for(unsigned int i=0; i<node->children.size(); ++i){
                if(box_zcast(p, node->children[i]->box))
                    stack.push_back(node->children[i]);
            }
        }
        return intersection_count%2;
    }

    bool
	distance_to_point(const Vec3d &p, double &rad, Vec3d &normal) const{
	    std::vector<const BoundingBoxTree*> stack;
//	    if(!sphere_box_intersect(p, rad, tree.box))
//	    	return false;
	    int triangle_index = 0;
	    bool hasCloser = false;
	    Vec3d barycentric;
		stack.push_back(&tree);
		while(!stack.empty()){
			const BoundingBoxTree *node=stack.back();
			stack.pop_back();
			if(!sphere_box_intersect(p, rad, node->box))
				continue; // no need to go further with this node
			// check any triangles in this node
			for(unsigned int i=0; i<node->index.size(); ++i){
				int u, v, w; assign(tri[node->index[i]], u, v, w);
				double a, b, c;
				double dist = point_tri_distance(p, x[u], x[v], x[w], &a, &b, &c);
				if(dist < rad){
					rad = dist;
					triangle_index = node->index[i];
					hasCloser = true;
					barycentric[0] = a;
					barycentric[1] = b;
					barycentric[2] = c;
				}
			}
			// check any subtrees for this node
			for(unsigned int i=0; i<node->children.size(); ++i)
				stack.push_back(node->children[i]);
		}
		if(hasCloser){
			Vec3i tri_id = tri[triangle_index];
			if(N){ // has per-vertex normal
					normal =  barycentric[0] * N[tri_id[0]] +
							  barycentric[1] * N[tri_id[1]] +
							  barycentric[2] * N[tri_id[2]];
			}
			else{
				int u, v, w; assign(tri_id, u, v, w);
				Vec3d vu = x[v] - x[u];
				Vec3d wu = x[w] - x[u];
				normal = cross(vu, wu);
				normalize(normal);
			}
		}
		return hasCloser;

	}

    bool
    intersects(const Vec3d &p,
               const Vec3d &q,
               int *triangle_index,
               double *s,
               double *t,
               double *a,
               double *b,
               double *c) const
    {
        std::vector<const BoundingBoxTree*> stack;
        stack.push_back(&tree);
        while(!stack.empty()){
            const BoundingBoxTree *node=stack.back();
            stack.pop_back();
            if(!segment_box_intersect(p, q, node->box))
                continue; // no need to go further with this node
            // check any triangles in this node
            for(unsigned int i=0; i<node->index.size(); ++i){
                int u, v, w; assign(tri[node->index[i]], u, v, w);
                if(segment_tri_intersect(p, q, x[u], x[v], x[w],
                                         s, t, a, b, c)){
                    *triangle_index=node->index[i];
                    return true;
                }
            }
            // check any subtrees for this node
            for(unsigned int i=0; i<node->children.size(); ++i)
                stack.push_back(node->children[i]);
        }
        return false;
    }

    bool
    ray_intersects(const Vec3d &p,
                   const Vec3d &q,
                   int *triangle_index,
                   double *t,
                   double *a,
                   double *b,
                   double *c,
                   bool &in) const
	{
		std::vector<const BoundingBoxTree*> stack;
		stack.push_back(&tree);
		*t = DBL_MAX;
		double tNear = DBL_MAX;
		int intersections = 0;
		while(!stack.empty()){
			const BoundingBoxTree *node=stack.back();
			stack.pop_back();
			if(!ray_box_intersect(p, q, node->box))
				continue; // no need to go further with this node
			// check any triangles in this node
			for(unsigned int i=0; i<node->index.size(); ++i){
				int u, v, w; assign(tri[node->index[i]], u, v, w);
				if(ray_tri_intersect(p, q, x[u], x[v], x[w],
										 &tNear, a, b, c)){
					*triangle_index=node->index[i];
					++intersections;
//					printf(" inters = %d, t = %f \n", intersections, *t);
					if( tNear < *t){
						*t = tNear;
						Vec3d normal;
						if(N){ // has per-vertex normal
								normal =  (*c) * N[u] +
										  (*a) * N[v] +
										  (*b) * N[w];
						}
						else{
							Vec3d vu = x[v] - x[u];
							Vec3d wu = x[w] - x[u];
							normal = cross(vu, wu);
							normalize(normal);
						}
						if(dot(q, normal) < 0.f)
							in = false;
						else
							in = true;
					}
				}
			}
			// check any subtrees for this node
			for(unsigned int i=0; i<node->children.size(); ++i)
				stack.push_back(node->children[i]);
		}
//		printf(" at the end of intersects, t = %f, ints = %d \n ", *t, intersections);
		in = intersections % 2;
		return intersections > 0;
	}
};

//==============================================================================
// the plain C wrapper API

MeshObject*
construct_mesh_object(int num_vertices,
                      const double *positions,
                      int num_triangles,
                      const int *triangles)
{
    exactinit(); // really only need to call this once, but safe to call always
    return new MeshObject(num_vertices, positions, num_triangles, triangles);
}

void
destroy_mesh_object(MeshObject *mesh)
{
    delete mesh;
}

bool
point_inside_mesh(const double point[3],
                  const MeshObject *mesh)
{
    return mesh->inside(point);
}

bool
point_to_mesh(const double point[3],
		      double &rad,
			  double normal[3],
              const MeshObject *mesh){
	Vec3d norm(0.0);
	Vec3d p(point);
	bool p2mesh = mesh->distance_to_point(p, rad, norm);
	normal[0] = norm[0];
	normal[1] = norm[1];
	normal[2] = norm[2];
	return p2mesh;
}

bool
segment_intersects_mesh(const double point0[3],
                        const double point1[3],
                        const MeshObject *mesh,
                        int *triangle_index,
                        double *s,
                        double *t,
                        double *a,
                        double *b,
                        double *c)
{
    return mesh->intersects(Vec3d(point0), Vec3d(point1),
                            triangle_index, s, t, a, b, c);
}

bool
ray_intersects_mesh(const double point0[3],
                    const double point1[3],
                    const MeshObject *mesh,
                    int *triangle_index,
					double *t,
					double *a,
					double *b,
					double *c,
					bool &in)
{
    return mesh->ray_intersects(Vec3d(point0), Vec3d(point1),
                            triangle_index, t, a, b, c, in);
}
