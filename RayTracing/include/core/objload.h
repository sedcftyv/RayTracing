#pragma once

#ifndef OBJLOAD_H
#define OBJLOAD_H

#include "pbrt.h"
#include<sstream>
#include<fstream>
#include "geometry.h"
using std::ifstream;
using std::string;
using std::istringstream;

inline bool objload(std::string path, std::vector<Point3f>&p,std::vector<int>&vi)
{
	ifstream fin(path);
	//ifstream fin("ratings1m.data");
	string in;
	while (getline(fin, in))
	{
		//std::cout << in << std::endl;
		istringstream line(in);
		char id;
		line >> id;
		switch (id)
		{
		case 'v':
			Float a, b, c;
			line >> a >> b >> c;
			p.push_back(Point3f(a, b, c));
			break;
		case 'f':
			int a1, b1, c1;
			line >> a1 >> b1 >> c1;
			vi.push_back(a1-1);
			vi.push_back(b1-1);
			vi.push_back(c1-1);
			break;
		}
	}
	return true;
}












#endif