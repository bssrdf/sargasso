#ifndef CELLTYPE_H_
#define CELLTYPE_H_

#include "geometry.h"

#define LIQUID 0x01
#define SOLID  0x02
#define EMPTY  0x04
#define SURFACE 0x08
#define SOURCE 0x10
#define MOVSOLID  1<<5

#define DONE   0x01
#define CLOSE  0x02
#define FAR    0x04
#ifdef SPMDMPI
void boundary3d_scalarc(char *array);
#endif
class CellType{

public:
	CellType(u_int size){
		type = new char[size];
		memset(type, 0, size);
	}
	~CellType(){
		delete [] type;
	}
	bool InLiquid(u_int index) const{
		return type[index] & LIQUID;
	}
	bool InSurface(u_int index) const{
		return type[index] & SURFACE;
	}
	bool InSolid(u_int index) const{
		return type[index] & SOLID;
	}
	bool InMovingSolid(u_int index) const{
		return type[index] & MOVSOLID;
	}
	bool InAir(u_int index) const{
		return type[index] & EMPTY;
	}
	bool InSource(u_int index) const{
		return type[index] & SOURCE;
	}
	bool IsDone(u_int index) const{
		return type[index] & DONE;
	}
	bool IsClose(u_int index) const{
		return type[index] & CLOSE;
	}
	bool IsFar(u_int index) const{
		return type[index] & FAR;
	}
	void UpdateTag(u_int index, char t){
		type[index] |= t;
	}
#ifdef SPMDMPI
	void UpdateGhosts(){
		boundary3d_scalarc(type);
	}
#endif
	void ClearTag(u_int index, char t){
		type[index] ^= t;
	}
private:
	char *type; // 0x01 liquid
	            // 0x02 solid
	            // ox04 empty
	            // 0x08 surface cell
				// 0x10 source cell
};

#endif /*CELLTYPE_H_*/
