#pragma once

#ifndef PBRT_SHAPES_SPHERE_H
#define PBRT_SHAPES_SPHERE_H

#include "pbrt.h"
#include "shape.h"

class Sphere : public Shape {
public:
	// Sphere Public Methods
	Sphere(const Transform *ObjectToWorld, const Transform *WorldToObject,
		bool reverseOrientation, Float radius)
		: Shape(ObjectToWorld, WorldToObject, reverseOrientation),
		radius(radius) {}
	Bounds3f ObjectBound() const;
	bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
		bool testAlphaTexture) const;
	bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
	Float Area() const;
	Interaction Sample(const Point2f &u, Float *pdf) const;
	//Interaction Sample(const Interaction &ref, const Point2f &u,
	//	Float *pdf) const;
	//Float Pdf(const Interaction &ref, const Vector3f &wi) const;
	//Float SolidAngle(const Point3f &p, int nSamples) const;

private:
	// Sphere Private Data
	const Float radius;
	//const Float zMin, zMax;
	//const Float thetaMin, thetaMax, phiMax;
};
















#endif