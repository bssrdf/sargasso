#ifndef SOURCEREGION_H_
#define SOURCEREGION_H_

#include "geometry.h" 

struct SourceRegion{
	
	SourceRegion(int _start_x0=0, int _start_y0=0, int _start_z0=0,
	int _end_x0=0, int _end_y0=0, int _end_z0=0, 
	int X = 0, int Y = 0, int Z = 0, float _speed = 0.f)
	: start_x0(_start_x0), start_y0(_start_y0), start_z0(_start_z0),
	  end_x0(_end_x0), end_y0(_end_y0), end_z0(_end_z0),
	  DimX(X), DimY(Y), DimZ(Z),
	  speed(_speed){ 
		
 			
	}
	
	void AddSource(int i, int j, int k){
		KnownPoint tmp(i,j,k);
		source.insert(make_pair(INDEX(i,j,k), tmp));
	}
	
	void AddSourceNeighbor(int i, int j, int k){
		KnownPoint tmp(i,j,k);
		sourceNeighbor.insert(make_pair(INDEX(i,j,k), tmp));
	}
	
	
	bool IsSourceNeighbor(int i, int j, int k) {		
		found = sourceNeighbor.find(INDEX(i,j,k));
		if(found != sourceNeighbor.end())
			return true;
		else
			return false;
	}
	
	bool IsSourceCell(int i, int j, int k) {		
		found = source.find(INDEX(i,j,k));
		if(found != source.end())
			return true;
		else
			return false;
	}
	
	void SourcePressure(float *p) const{
		int i, j, k;
		k =  end_z0;  
		for(j = start_y0; j <= end_y0; ++j)
			for(i = start_x0; i <= end_x0; ++i){
				p[INDEX(i,j,k)] = p[INDEX(i,j,k+1)];
			}
	}
	
	int start_x0, start_y0, start_z0;
	int end_x0, end_y0, end_z0;
	float speed;
	int DimX, DimY, DimZ;
	map<u_int, KnownPoint> source;
	map<u_int, KnownPoint> sourceNeighbor;
	map<u_int, KnownPoint>::iterator found;
};

#endif /*SOURCEREGION_H_*/