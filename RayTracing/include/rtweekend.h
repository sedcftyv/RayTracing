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

//const double infinity = std::numeric_limits<double>::infinity();
//const double pi = 3.1415926535897932385;
//
//// Utility Functions
//
//inline double degrees_to_radians(double degrees) {
//	return degrees * pi / 180.0;
//}
//
//inline double random_double() {
//	// Returns a random real in [0,1).
//	return rand() / (RAND_MAX + 1.0);
//}
//
//inline double random_double(double min, double max) {
//	// Returns a random real in [min,max).
//	return min + (max - min)*random_double();
//}
//
//inline double clamp(double x, double min, double max) {
//	if (x < min) return min;
//	if (x > max) return max;
//	return x;
//}
//
//inline int random_int(int min, int max) {
//	// Returns a random integer in [min,max].
//	return static_cast<int>(random_double(min, max + 1));
//}
//
//
//// Common Headers
//
////#include "ray.h"
//#include "core/geometry.h"
//
//inline Vector3f random_cosine_direction() {
//	auto r1 = random_double();
//	auto r2 = random_double();
//	auto z = sqrt(1 - r2);
//
//	auto phi = 2 * pi*r1;
//	auto x = cos(phi)*sqrt(r2);
//	auto y = sin(phi)*sqrt(r2);
//
//	return Vector3f(x, y, z);
//}
//
//inline  Vector3f random() {
//	return Vector3f(random_double(), random_double(), random_double());
//}
//
//inline  Vector3f random(double min, double max) {
//	return Vector3f(random_double(min, max), random_double(min, max), random_double(min, max));
//}
//
//inline Vector3f random_in_unit_sphere() {
//	while (true) {
//		auto p = random(-1, 1);
//		if (p.LengthSquared() >= 1) continue;
//		return p;
//	}
//}
//inline Vector3f random_unit_vector() {
//	return Normalize(random_in_unit_sphere());
//}
//inline Vector3f random_to_sphere(double radius, double distance_squared) {
//	auto r1 = random_double();
//	auto r2 = random_double();
//	auto z = 1 + r2 * (sqrt(1 - radius * radius / distance_squared) - 1);
//
//	auto phi = 2 * pi*r1;
//	auto x = cos(phi)*sqrt(1 - z * z);
//	auto y = sin(phi)*sqrt(1 - z * z);
//
//	return Vector3f(x, y, z);
//}
//inline Vector3f random_in_unit_disk() {
//	while (true) {
//		auto p = Vector3f(random_double(-1, 1), random_double(-1, 1), 0);
//		if (p.LengthSquared() >= 1) continue;
//		return p;
//	}
//}
//
//inline Vector3f random_in_hemisphere(const Vector3f& normal) {
//	Vector3f in_unit_sphere = random_in_unit_sphere();
//	if (Dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
//		return in_unit_sphere;
//	else
//		return -in_unit_sphere;
//}
//
//inline Vector3f refract(const Vector3f& uv, const Vector3f& n, double etai_over_etat) {
//	auto cos_theta = fmin(Dot(-uv, n), 1.0);
//	Vector3f r_out_perp = etai_over_etat * (uv + cos_theta * n);
//	Vector3f r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.LengthSquared())) * n;
//	return r_out_perp + r_out_parallel;
//}
//
//inline Vector3f reflect(const Vector3f& v, const Vector3f& n) {
//	return v - 2 * Dot(v, n)*n;
//}

#endif