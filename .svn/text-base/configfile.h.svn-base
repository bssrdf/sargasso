/*
 * configfile.h
 *
 *  Created on: Nov 5, 2010
 *      Author: bzhao
 */

#ifndef CONFIGFILE_H_
#define CONFIGFILE_H_

#include <iostream>
#include <vector>
using namespace std;

#include "transform.h"

class ConfigFile{
public:
	ConfigFile(const char* f)
	: mNumObjects(0){
		Parse(f);
	}
	~ConfigFile(){
	}
	void Parse(const char *f);
	bool HasNextObject();
	Transform ReadTransform(unsigned int n);
	string ReadObjectFileName(unsigned int n);
	bool IsAWObj(string objfile);

private:
	int    mNumObjects;
	vector<string>    mObjectFileNames;
	vector<Transform> mObjectTransforms;


};

#endif /* CONFIGFILE_H_ */
