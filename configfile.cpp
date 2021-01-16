/*
 * configfile.cpp
 *
 *  Created on: Nov 5, 2010
 *      Author: bzhao
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "configfile.h"

#define CHARSTRSIZE 100
#define TRANSFORMMATRIX 16
#define TRANSFORMMATRIXROWS 4

void ConfigFile::Parse(const char *f){
	char token[CHARSTRSIZE], buf[CHARSTRSIZE], filename[CHARSTRSIZE];
	char tran[CHARSTRSIZE];
	float trans[TRANSFORMMATRIXROWS][TRANSFORMMATRIXROWS];
	FILE *scene = fopen(f,"r");
	if(scene){
		while(!feof(scene)){
			token[0] = NULL;
			fscanf(scene,"%s", token);
			if(!strcmp(token,"object")){
				++mNumObjects;
				fscanf(scene,"%s",buf);
				if(!strcmp(buf,"filename")){
					if(!strcmp(buf,"filename")){
						fscanf(scene,"%s",filename);
//						printf("filename is %s \n", filename);
						mObjectFileNames.push_back(string(filename));
					}
				}
				fscanf(scene,"%s",buf);
//				printf("buf is %s \n", buf);
				if(!strcmp(buf,"transform")){
					for(unsigned int i=0; i<TRANSFORMMATRIXROWS; ++i){
						for(unsigned int j=0; j<TRANSFORMMATRIXROWS; ++j){
							fscanf(scene,"%s",tran);
							trans[i][j] = (float)atof(tran);
//						printf(" %u component of transform is %s \n", n, tran[n]);
						}
					}
					mObjectTransforms.push_back(Transform(trans));
				}
			}

		}
		fclose(scene);
	}
	else{
		printf(" No scene file with name [%s] exists \n", f);
	}

}

bool ConfigFile::HasNextObject() {
	if(mNumObjects > 0){
		--mNumObjects;
		return true;
	}
	else
		return false;
}

Transform ConfigFile::ReadTransform(unsigned int n){
	return mObjectTransforms[n];
}

string ConfigFile::ReadObjectFileName(unsigned int n){
	return mObjectFileNames[n];
}

bool ConfigFile::IsAWObj(string objfile){
	char suffix[5]= ".obj" ;
	const char *of = objfile.c_str();
	if(NULL != strstr(of, suffix))
		return true;
	else
		return false;
}
