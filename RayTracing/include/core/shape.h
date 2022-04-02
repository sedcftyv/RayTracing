#pragma once

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

#include "pbrt.h"
#include "geometry.h"
#include "transform.h"

class Shape {
public:
	// Shape Interface
	Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
		bool reverseOrientation);
	virtual ~Shape();
	virtual Bounds3f ObjectBound() const = 0;
	virtual Bounds3f WorldBound() const;
	virtual bool Intersect(const Ray &ray, Float *tHit,
		SurfaceInteraction *isect,
		bool testAlphaTexture = true) const = 0;
	virtual bool IntersectP(const Ray &ray,
		bool testAlphaTexture = true) const {
		return Intersect(ray, nullptr, nullptr, testAlphaTexture);
	}
	virtual Float Area() const = 0;
	virtual Interaction Sample(const Point2f &u, Float *pdf) const = 0;
	virtual Float Pdf(const Interaction &) const { return 1 / Area(); }
	virtual Interaction Sample(const Interaction &ref, const Point2f &u,
		Float *pdf) const;
	virtual Float Pdf(const Interaction &ref, const Vector3f &wi) const;
	//virtual Float SolidAngle(const Point3f &p, int nSamples = 512) const;

	// Shape Public Data
	const Transform *ObjectToWorld, *WorldToObject;
	const bool reverseOrientation;
	const bool transformSwapsHandedness;
};











#endif