#pragma once

#ifndef RAY_H
#define RAY_H

//#include "Vec3.h"
#include "core/geometry.h"

class ray {
public:
	ray() {}
	ray(const Vector3f& origin, const Vector3f& direction, double time = 0.0)
		: orig(origin), dir(direction), tm(time)
	{}

	Vector3f origin() const { return orig; }
	Vector3f direction() const { return dir; }
	double time() const { return tm; }

	Vector3f at(double t) const {
		return orig + t * dir;
	}

public:
	Vector3f orig;
	Vector3f dir;
	double tm;
};

#endif