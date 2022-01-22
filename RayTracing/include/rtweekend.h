#pragma once

#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <memory>
// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

inline double random_double() {
	// Returns a random real in [0,1).
	return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
	// Returns a random real in [min,max).
	return min + (max - min)*random_double();
}

inline double clamp(double x, double min, double max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

inline int random_int(int min, int max) {
	// Returns a random integer in [min,max].
	return static_cast<int>(random_double(min, max + 1));
}


// Common Headers

#include "ray.h"
#include "vec3.h"

inline vec3 random_cosine_direction() {
	auto r1 = random_double();
	auto r2 = random_double();
	auto z = sqrt(1 - r2);

	auto phi = 2 * pi*r1;
	auto x = cos(phi)*sqrt(r2);
	auto y = sin(phi)*sqrt(r2);

	return vec3(x, y, z);
}

#endif