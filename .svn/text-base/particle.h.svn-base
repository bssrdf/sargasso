#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "geometry.h"

class Particle{
public:
	Particle(const Point &_p=Point(0,0,0), float r=0, char s = -1, float cd = 0.f)
		: p(_p), rad(r), particle_sign(s), collisiondist(cd) {

	}
	virtual ~Particle(){}
	float EvaluatePhiP(const Point &p) const;
	float EvaluatePhiP(float d) const;
	inline void SetParticle(const Point &_p, float r){
		rad = r;
		p = _p;
	}
	inline void SetParticle(const Point &_p){
		p = _p;
	}
	void Print(int i) const;
	float DistanceTo(const Point &p) const;
	bool PointInside(const Point &pos) const{
		if(DistanceTo(pos) < rad)
			return true;
		else
			return false;
	}
	inline Point Position() const{
//		return pos;
		return p;
	}
	inline float Radius() const{
		return rad;
	}
	inline float CollisionDistance() const{
		return collisiondist;
	}
	void OutputPBRT(FILE *filename) const;
	void MoveUnderGravity(float dt);
	virtual float Influence(const Point &pos, float h) const{
		float r2 = h * h;
		float dist2 = DistanceSquared(p, pos);
		if(dist2 <= r2)
			return 15.f / (8 * M_PI * h * h * h) * (1.f - dist2 / r2);
		else
			return 0.f;
	}
	Point p;
protected:
	float rad;
private:
//	Point pos;
	char particle_sign;
	float collisiondist;
};

typedef list<Particle> ParticleList;
typedef vector<Particle> ParticleVector;

class WaterParticle : public Particle{

public:
	WaterParticle(const Point &_p=Point(0,0,0), float r=0, char s = -1, float cd = 0,
				float u=0.f, float v=0.f, float w = 0.f)
		: Particle(_p, r, s, cd), vel(u,v,w){
	}
	~WaterParticle(){}
	inline Vector Velocity() const{
		return vel;
	}
	inline float UVelocity() const{
		return vel.x;
	}
	inline float VVelocity() const{
		return vel.y;
	}
	inline float WVelocity() const{
		return vel.z;
	}
	inline void SetVelocity(float u, float v, float w){
		vel = Vector(u, v, w);
	}
	inline void SetVelocity(const Vector &v){
		vel = v;
	}

private:
	Vector vel;
};

typedef list<WaterParticle> WaterParticleList;
//typedef vector<WaterParticle> WaterParticleList;

class VortexParticle : public Particle{

public:
	VortexParticle(const Point &_p=Point(0,0,0), float r=0, char s = -1, float cd = 0,
				float m = 0, float u=0.f, float v=0.f, float w = 0.f)
		: Particle(_p, r, s, cd), magnitude(m), dir(u,v,w){
	}
	~VortexParticle(){}

	inline Vector Direction() const{
		return dir;
	}

	inline Vector Vorticity() const{
		return magnitude * dir;
	}

	inline float Magnitude() const{
			return magnitude;
		}

	inline void SetDirection(float u, float v, float w){
		dir = Vector(u, v, w);
	}

	inline void SetDirection(const Vector &v){
		dir = v;
	}

	inline void SetMagnitude(float m){
		magnitude = m;
	}

	inline float Influence(const Point &pos, float r){
		float r2 = r * r;
		float dist2 = DistanceSquared(p, pos);
		if(dist2 <= r2)
			return expf(-dist2/(2*r2))/(r*r2*powf(2*M_PI,1.5f));
		else
			return 0.f;
	}
private:
	float magnitude;
	Vector dir;

};
#endif /*PARTICLE_H_*/
