#ifndef GEOMETRY_H_
#define GEOMETRY_H_

// geometry.h*
//#include "pbrt.h"
#include <float.h>
#ifdef WIN32
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <assert.h>

#if !defined(__APPLE__)
#include <malloc.h> // for _alloca, memalign
#endif

#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
using std::ostream;
using std::min;
using std::max;
using std::vector;
using std::list;
using std::swap;


//#ifdef WIN32
//#ifdef CORE_SOURCE
//#define COREDLL __declspec(dllexport)
//#else
//#define COREDLL __declspec(dllimport)
//#endif
//#define DLLEXPORT __declspec(dllexport)
//#else
#define COREDLL
#define DLLEXPORT
//#endif
#ifdef __APPLE__
#define powf pow
#define sinf sin
#define cosf cos
#define tanf tan
#define asinf asin
#define acosf acos
#define atanf atan
#define atan2f atan2
#define logf log
#define log10f log10
#define expf exp
#define sqrtf sqrt
#if __GNUC__ == 3
extern "C" {
  int isinf(double);
  int isnan(double);
}
#endif // ONLY GCC 3
#endif // __APPLE__

#ifdef WIN32
#define cbrtf(a) pow(a,b)
#else
#define powf(a,b) pow(a,b)
#endif

//#define memalign(a,b) valloc(b)
#ifdef WIN32
#define memalign(a,b) _aligned_malloc(b, a)
#else
#define memalign(a,b) valloc(b)
#endif

COREDLL void *AllocAligned(size_t size);

COREDLL void FreeAligned(void *ptr);

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif


#ifndef M_PI
#define M_PI  3.14159265358979323846f
#endif

#ifndef EPSIL_C
#define EPSIL_C 1.e-3f
#endif
#define sign(a) ((a) >= 0.f ? 1 : -1)


typedef unsigned int u_int;

inline bool Quadratic(float A, float B, float C, float *t0,
		float *t1) {
	// Find quadratic discriminant
	float discrim = B * B - 4.f * A * C;
	if (discrim < 0.) return false;
	float rootDiscrim = sqrtf(discrim);
	// Compute quadratic _t_ values
	float q;
	if (B < 0) q = -.5f * (B - rootDiscrim);
	else       q = -.5f * (B + rootDiscrim);
	*t0 = q / A;
	*t1 = C / q;
	if (*t0 > *t1){
		 float tmp = *t0;
		 *t0 = *t1;
		 *t1 = tmp;
	}
	return true;
}

const float one_sixth = 1.f / 6;
const float one_third = 1.f / 3;
const float thirteen_twelve = 13.f / 12.f;
const float seven_sixth = 7.f * one_sixth;
const float five_sixth  = 5.f * one_sixth;
const float eleven_sixth  = 11.f * one_sixth;

#define POSINDONEBAND 0x01
#define NEGINDONEBAND 0x02
#define POSINTRIALBAND 0x04
#define NEGINTRIALBAND 0x08

class Vector;
class Point;
class Normal;

// Geometry Declarations
class COREDLL Vector {
public:
	// Vector Public Methods
	Vector(float _x=0, float _y=0, float _z=0)
		: x(_x), y(_y), z(_z) {
	}
	explicit Vector(const Point &p);
	Vector operator+(const Vector &v) const {
		return Vector(x + v.x, y + v.y, z + v.z);
	}

	Vector& operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Vector operator-(const Vector &v) const {
		return Vector(x - v.x, y - v.y, z - v.z);
	}

	Vector& operator-=(const Vector &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	bool operator==(const Vector &v) const {
		return x == v.x && y == v.y && z == v.z;
	}
	Vector operator*(float f) const {
		return Vector(f*x, f*y, f*z);
	}

	Vector &operator*=(float f) {
		x *= f; y *= f; z *= f;
		return *this;
	}
	Vector operator/(float f) const {
		assert(f!=0);
		float inv = 1.f / f;
		return Vector(x * inv, y * inv, z * inv);
	}

	Vector &operator/=(float f) {
		assert(f!=0);
		float inv = 1.f / f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	Vector operator-() const {
		return Vector(-x, -y, -z);
	}
	float operator[](int i) const {
		assert(i >= 0 && i <= 2);
		return (&x)[i];
	}

	float &operator[](int i) {
		assert(i >= 0 && i <= 2);
		return (&x)[i];
	}
	float LengthSquared() const { return x*x + y*y + z*z; }
	float Length() const { return sqrtf(LengthSquared()); }
	explicit Vector(const Normal &n);
	// Vector Public Data
	float x, y, z;
};
class COREDLL Point {
public:
	// Point Methods
	Point(float _x=0, float _y=0, float _z=0)
		: x(_x), y(_y), z(_z) {
	}
	Point operator+(const Vector &v) const {
		return Point(x + v.x, y + v.y, z + v.z);
	}

	Point &operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Vector operator-(const Point &p) const {
		return Vector(x - p.x, y - p.y, z - p.z);
	}

	Point operator-(const Vector &v) const {
		return Point(x - v.x, y - v.y, z - v.z);
	}

	Point &operator-=(const Vector &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	Point &operator+=(const Point &p) {
		x += p.x; y += p.y; z += p.z;
		return *this;
	}
	Point operator+(const Point &p) const {
		return Point(x + p.x, y + p.y, z + p.z);
	}
	Point operator* (float f) const {
		return Point(f*x, f*y, f*z);
	}
	Point &operator*=(float f) {
		x *= f; y *= f; z *= f;
		return *this;
	}
	Point operator/ (float f) const {
		float inv = 1.f/f;
		return Point(inv*x, inv*y, inv*z);
	}
	Point &operator/=(float f) {
		float inv = 1.f/f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	float operator[](int i) const { return (&x)[i]; }
	float &operator[](int i) { return (&x)[i]; }
	// Point Public Data
	float x,y,z;
};
class COREDLL Normal {
public:
	// Normal Methods
	Normal(float _x=0, float _y=0, float _z=0)
		: x(_x), y(_y), z(_z) {}
	Normal operator-() const {
		return Normal(-x, -y, -z);
	}
	Normal operator+ (const Normal &v) const {
		return Normal(x + v.x, y + v.y, z + v.z);
	}

	Normal& operator+=(const Normal &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Normal operator- (const Normal &v) const {
		return Normal(x - v.x, y - v.y, z - v.z);
	}

	Normal& operator-=(const Normal &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	Normal operator* (float f) const {
		return Normal(f*x, f*y, f*z);
	}

	Normal &operator*=(float f) {
		x *= f; y *= f; z *= f;
		return *this;
	}
	Normal operator/ (float f) const {
		float inv = 1.f/f;
		return Normal(x * inv, y * inv, z * inv);
	}

	Normal &operator/=(float f) {
		float inv = 1.f/f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	float LengthSquared() const { return x*x + y*y + z*z; }
	float Length() const        { return sqrtf(LengthSquared()); }

	explicit Normal(const Vector &v)
	  : x(v.x), y(v.y), z(v.z) {}
	float operator[](int i) const { return (&x)[i]; }
	float &operator[](int i) { return (&x)[i]; }
	// Normal Public Data
	float x,y,z;
};
class COREDLL BBox {
public:
	// BBox Public Methods
	BBox() {
		pMin = Point( INFINITY,  INFINITY,  INFINITY);
		pMax = Point( -INFINITY, -INFINITY, -INFINITY);
	}
	BBox(const Point &p) : pMin(p), pMax(p) { }
	BBox(const Point &p1, const Point &p2) {
		pMin = Point(min(p1.x, p2.x),
					 min(p1.y, p2.y),
					 min(p1.z, p2.z));
		pMax = Point(max(p1.x, p2.x),
					 max(p1.y, p2.y),
					 max(p1.z, p2.z));
	}
	friend inline ostream &
		operator<<(ostream &os, const BBox &b);
	friend COREDLL BBox Union(const BBox &b, const Point &p);
	friend COREDLL BBox Union(const BBox &b, const BBox &b2);
	bool Overlaps(const BBox &b) const {
		bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
		bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
		bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
		return (x && y && z);
	}
	bool Inside(const Point &pt) const {
		return (pt.x >= pMin.x && pt.x <= pMax.x &&
	            pt.y >= pMin.y && pt.y <= pMax.y &&
	            pt.z >= pMin.z && pt.z <= pMax.z);
	}
	void Expand(float delta) {
		pMin -= Vector(delta, delta, delta);
		pMax += Vector(delta, delta, delta);
	}
	float Volume() const {
		Vector d = pMax - pMin;
		return d.x * d.y * d.z;
	}
	int MaximumExtent() const {
		Vector diag = pMax - pMin;
		if (diag.x > diag.y && diag.x > diag.z)
			return 0;
		else if (diag.y > diag.z)
			return 1;
		else
			return 2;
	}
	void BoundingSphere(Point *c, float *rad) const;
	// BBox Public Data
	Point pMin, pMax;
};

#define SET true
#define CLEAR false

struct TentativeGridPoint;
struct KnownPoint{
	KnownPoint(int i=0, int j=0, int k=0)
	: ii(i), jj(j), kk(k) { }
//	KnownPoint(const TentativeGridPoint &r)
//	: ii(r.ii), jj(r.jj), kk(r.kk) { }
	int ii, jj, kk;

	bool operator==(const KnownPoint &p){
		if(ii == p.ii && jj == p.jj && kk == p.kk)
			return true;
		else
			return false;
	}
	bool operator!=(const KnownPoint &p){
		if(ii != p.ii || jj != p.jj || kk != p.kk)
			return true;
		else
			return false;
	}
	KnownPoint& operator=(const KnownPoint &rhs){
		if(this != &rhs){
			ii = rhs.ii;
			jj = rhs.jj;
			kk = rhs.kk;
		}
		return *this;
	}
};

struct TentativeGridPoint{
	TentativeGridPoint(float v, int i, int j, int k)
	: value(v), ii(i), jj(j), kk(k) { }
	TentativeGridPoint(float v, const KnownPoint &p)
	: value(v), ii(p.ii), jj(p.jj), kk(p.kk) { }
	TentativeGridPoint() { }
	KnownPoint Pos() const{
		return KnownPoint(ii,jj,kk);
	}
	bool LeftNeighbor(const TentativeGridPoint &n) const{
		if(ii == n.ii+1)
			return true;
		else
			return false;
	}
	bool RightNeighbor(const TentativeGridPoint &n) const{
		if(ii == n.ii-1)
			return true;
		else
			return false;
	}
	bool BackNeighbor(const TentativeGridPoint &n) const{
		if(jj == n.jj+1)
			return true;
		else
			return false;
	}
	bool FrontNeighbor(const TentativeGridPoint &n) const{
		if(jj == n.jj-1)
			return true;
		else
			return false;
	}
	bool BottomNeighbor(const TentativeGridPoint &n) const{
		if(kk == n.kk+1)
			return true;
		else
			return false;
	}
	bool TopNeighbor(const TentativeGridPoint &n) const{
		if(kk == n.kk-1)
			return true;
		else
			return false;
	}
	TentativeGridPoint& operator=(const TentativeGridPoint &rhs){
		if(this != &rhs){
			ii = rhs.ii;
			jj = rhs.jj;
			kk = rhs.kk;
		}
		return *this;
	}
	bool operator==(const TentativeGridPoint &v) const {
		return ii == v.ii && jj == v.jj && kk == v.kk;
	}
	void Neigbor(vector<TentativeGridPoint> &n, int X, int Y, int Z) const{
		if( ii != 0 )
			n.push_back(TentativeGridPoint(0,ii-1,jj,kk));
		if( ii+1 < X )
			n.push_back(TentativeGridPoint(0,ii+1,jj,kk));
		if( jj != 0 )
		 	n.push_back(TentativeGridPoint(0,ii,jj-1,kk));
		if( jj+1 < Y )
			n.push_back(TentativeGridPoint(0,ii,jj+1,kk));
		if( kk != 0 )
			n.push_back(TentativeGridPoint(0,ii,jj,kk-1));
		if( kk+1 < Z )
			n.push_back(TentativeGridPoint(0,ii,jj,kk+1));
		return;
	}
	void AllNeigbor(vector<TentativeGridPoint> &n, int X, int Y, int Z) const{
		n.push_back(TentativeGridPoint(0,ii-1,jj,kk));
		n.push_back(TentativeGridPoint(0,ii+1,jj,kk));
	 	n.push_back(TentativeGridPoint(0,ii,jj-1,kk));
		n.push_back(TentativeGridPoint(0,ii,jj+1,kk));
		n.push_back(TentativeGridPoint(0,ii,jj,kk-1));
		n.push_back(TentativeGridPoint(0,ii,jj,kk+1));
		return;
	}

	bool NeighborExists(int index, int X, int Y, int Z){

		switch(index){
			case 0:
				if( ii-1 < 0 )
					return false;
				else{
					return true;
				}
				break;
			case 1:
				if( ii+1 > X-1 )
					return false;
				else{
					return true;
				}
				break;
			case 2:
				if( jj-1 < 0 )
					return false;
				else{
					return true;
				}
				break;
			case 3:
				if( jj+1 > Y-1 )
					return false;
				else{
					return true;
				}
				break;
			case 4:
				if( kk-1 < 0 )
					return false;
				else
					return true;
				break;
			case 5:
				if( kk+1 > Z-1 )
					return false;
				else
					return true;
				break;
		}
	}
	float value;
	int ii, jj, kk;
};

struct COREDLL Matrix4x4 {
	// Matrix4x4 Public Methods
	Matrix4x4() {
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				if (i == j) m[i][j] = 1.;
				else m[i][j] = 0.;
	}
	Matrix4x4(float mat[4][4]);
	Matrix4x4(float t00, float t01, float t02, float t03,
	          float t10, float t11, float t12, float t13,
	          float t20, float t21, float t22, float t23,
	          float t30, float t31, float t32, float t33);
	Matrix4x4 Transpose() const;
	void Print(ostream &os) const {
		os << "[ ";
		for (int i = 0; i < 4; ++i) {
			os << "[ ";
			for (int j = 0; j < 4; ++j)  {
				os << m[i][j];
				if (j != 3) os << ", ";
			}
			os << " ] ";
		}
		os << " ] ";
	}
	static Matrix4x4
		Mul(const Matrix4x4 &m1,
	        const Matrix4x4 &m2) {
		float r[4][4];
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				r[i][j] = m1.m[i][0] * m2.m[0][j] +
				          m1.m[i][1] * m2.m[1][j] +
				          m1.m[i][2] * m2.m[2][j] +
				          m1.m[i][3] * m2.m[3][j];
		return Matrix4x4(r);
	}
	Matrix4x4 Inverse() const;
	float m[4][4];
};
// Geometry Inline Functions
inline Vector::Vector(const Point &p)
	: x(p.x), y(p.y), z(p.z) {
}
inline ostream &operator<<(ostream &os, const Vector &v) {
	os << v.x << ", " << v.y << ", " << v.z;
	return os;
}
inline Vector operator*(float f, const Vector &v) {
	return v*f;
}
inline float Dot(const Vector &v1, const Vector &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
inline float AbsDot(const Vector &v1, const Vector &v2) {
	return fabsf(Dot(v1, v2));
}
inline Vector Cross(const Vector &v1, const Vector &v2) {
	return Vector((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}
inline Vector Cross(const Vector &v1, const Normal &v2) {
	return Vector((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}
inline Vector Cross(const Normal &v1, const Vector &v2) {
	return Vector((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}
inline Vector Normalize(const Vector &v) {
	return v / v.Length();
}
inline void CoordinateSystem(const Vector &v1, Vector *v2, Vector *v3) {
	if (fabsf(v1.x) > fabsf(v1.y)) {
		float invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
		*v2 = Vector(-v1.z * invLen, 0.f, v1.x * invLen);
	}
	else {
		float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
		*v2 = Vector(0.f, v1.z * invLen, -v1.y * invLen);
	}
	*v3 = Cross(v1, *v2); //left handed
//	*v3 = Cross(*v2, v1); //right handed
}
inline float Distance(const Point &p1, const Point &p2) {
	return (p1 - p2).Length();
}
inline float DistanceSquared(const Point &p1, const Point &p2) {
	return (p1 - p2).LengthSquared();
}
inline float DistanceToLine(const Point &p1, const Point &p2, const Point &p0, float f){
	Vector d = p2 - p1;
	Vector d1 = p0 - p1;
	if(d.Length() == 0.f){
		printf(" (DistanceToLine): p0 = (%f, %f, %f), p1 = (%f, %f, %f), p2 = (%f, %f, %f)\n",
			p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
//		exit(1);
		return f;
	}
	else{
		d = Normalize(d);
		float t0 = Dot(d, d1);
		Point q = p1 + t0 * d;
		return Distance(p0, q);
	}
}
inline ostream &operator<<(ostream &os, const Point &v) {
	os << v.x << ", " << v.y << ", " << v.z;
	return os;
}
inline Point operator*(float f, const Point &p) {
	return p*f;
}
inline Normal operator*(float f, const Normal &n) {
	return Normal(f*n.x, f*n.y, f*n.z);
}
inline Normal Normalize(const Normal &n) {
	return n / n.Length();
}
inline Vector::Vector(const Normal &n)
  : x(n.x), y(n.y), z(n.z) { }
inline float Dot(const Normal &n1, const Vector &v2) {
	return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}
inline float Dot(const Vector &v1, const Normal &n2) {
	return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}
inline float Dot(const Normal &n1, const Normal &n2) {
	return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}
inline float AbsDot(const Normal &n1, const Vector &v2) {
	return fabsf(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}
inline float AbsDot(const Vector &v1, const Normal &n2) {
	return fabsf(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}
inline float AbsDot(const Normal &n1, const Normal &n2) {
	return fabsf(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}
inline ostream &operator<<(ostream &os, const Normal &v) {
	os << v.x << ", " << v.y << ", " << v.z;
	return os;
}
inline float lerp(float phi1, float phi2, float u){
	return (1.f-u) * phi1 + u * phi2;
}

inline float Radians(float deg) {
	return ((float)M_PI/180.f) * deg;
}

inline float FindMinimum(float *r, int m){
     float s = INFINITY;
     for(int i=0; i < m; ++i){
    	 if(s > r[i]) s = r[i];
     }
     return s;
}

inline float TwoPointDistance(float phi1, float phi2){
	return phi1 * phi2 / sqrt( phi1*phi1 + phi2*phi2 );
}


inline float bilerp(float phi1, float phi2, float phi3, float phi4, float u, float v){
	return lerp(lerp(phi1, phi2, u), lerp(phi3, phi4, u), v);
}

inline float trilerp(float phi1, float phi2, float phi3, float phi4,
					  float phi5, float phi6, float phi7, float phi8,
					  float u, float v, float w){
	return bilerp(lerp(phi1, phi2, u), lerp(phi3, phi4, u),
				  lerp(phi5, phi6, u), lerp(phi7, phi8, u),
				  w, v);
}

#endif /*GEOMETRY_H_*/
