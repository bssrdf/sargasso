#include <iostream>

using std::cout;
using std::endl;


#include "particle.h"

#define g -9.8

float Particle::EvaluatePhiP(const Point &pos) const{
	return DistanceTo(pos) - rad;
}

float Particle::EvaluatePhiP(float d) const{
	return d - rad;
}

void Particle::Print(int i) const{
	cout << "Particle# = " << i << " (" << p << ") rad= " << rad << endl;
}

void Particle::OutputPBRT(FILE *filename) const{
	fprintf(filename, "AttributeBegin \n");
	fprintf(filename, "Translate %f  %f  %f \n", p.x, p.y, p.z);
	fprintf(filename, "Shape \"sphere\" \"float radius\"  %f \n", rad);
	fprintf(filename, "AttributeEnd \n");
}

void Particle::MoveUnderGravity(float dt) {
	p.z += dt*g;
	return;
}

float Particle::DistanceTo(const Point &pos) const{
	return Distance(p, pos);
}
