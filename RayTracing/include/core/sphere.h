#pragma once

#ifndef PBRT_SHAPES_SPHERE_H
#define PBRT_SHAPES_SPHERE_H

#include "pbrt.h"
#include "shape.h"

class Sphere : public Shape {
public:
	Sphere(const Transform *ObjectToWorld, const Transform *WorldToObject,bool reverseOrientation, Float radius)
		: Shape(ObjectToWorld, WorldToObject, reverseOrientation),radius(radius) {}
	Bounds3f ObjectBound() const;
	bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
		bool testAlphaTexture) const;
	bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
	Float Area() const;			
private:
	const Float radius;
};
















#endif