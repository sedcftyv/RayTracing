#pragma once

#ifndef INTERACTION_H
#define INTERACTION_H

#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "material.h"

struct Interaction {
	Interaction() :time(0) { }
	Interaction(const Point3f &p, const Normal3f &n,const Vector3f &wo, Float time)
		: p(p), time(time), wo(wo), n(n) { }
	bool IsSurfaceInteraction() const {
		return n != Normal3f();
	}
	Ray SpawnRay(const Vector3f &d) const {
		Point3f  o = OffsetRayOrigin(p,  n, d);
		return Ray(o, d, Infinity, time);
	}
	Ray SpawnRayTo(const Point3f  &p2) const {
		Point3f  origin = OffsetRayOrigin(p, n, Vector3f(p2 - p));
		Vector3f d = Vector3f(p2 - origin);
		return Ray(origin, d, 1 - ShadowEpsilon, time );
	}
	Ray SpawnRayTo(const Interaction &it) const {
		Point3f  origin = OffsetRayOrigin(p,  n, it.p - p);
		Point3f  target = OffsetRayOrigin(it.p,  it.n, origin - it.p);
		Vector3f d = target - origin;
		return Ray(origin, d, 1 - ShadowEpsilon, time);
	}
	Interaction(const Point3f &p, const Vector3f &wo, Float time)
		: p(p), time(time), wo(wo){ }
	Interaction(const Point3f &p, Float time)
		: p(p), time(time){ }
	
	Point3f p;
	Float time;
		Vector3f wo;
	Normal3f n;
	};

class SurfaceInteraction : public Interaction {
public:
	SurfaceInteraction() { }
	SurfaceInteraction(const Point3f &p,const Point2f &uv,
		const Vector3f &wo, const Vector3f &dpdu, const Vector3f &dpdv,
		const Normal3f &dndu, const Normal3f &dndv,
		Float time, const Shape *sh, int faceIndex=0);
	~SurfaceInteraction();
	void SetShadingGeometry(const Vector3f &dpdu, const Vector3f &dpdv,
		const Normal3f &dndu, const Normal3f &dndv, bool orientationIsAuthoritative);
		void ComputeScatteringFunctions(
		const Ray &ray, bool allowMultipleLobes = false,
		TransportMode mode = TransportMode::Radiance);
		Spectrum Le(const Vector3f &w) const;

	Point2f uv;
	Vector3f dpdu, dpdv;
	Normal3f dndu, dndv;
	const Shape *shape = nullptr;
	struct {
		Normal3f n;
		Vector3f dpdu, dpdv;
		Normal3f dndu, dndv;
	} shading;
	const Primitive *primitive = nullptr;
	BSDF *bsdf = nullptr;
	mutable Vector3f dpdx, dpdy;
	mutable Float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0;
	int faceIndex = 0;
};


#endif