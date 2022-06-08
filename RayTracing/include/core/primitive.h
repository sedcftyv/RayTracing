#pragma once

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

#include "pbrt.h"
#include "shape.h"
#include "material.h"

class Primitive {
public:
	virtual ~Primitive();
	virtual Bounds3f WorldBound() const = 0;
	virtual bool Intersect(const Ray &r, SurfaceInteraction *) const = 0;
	virtual bool IntersectP(const Ray &r) const = 0;
	virtual const AreaLight *GetAreaLight() const = 0;
	virtual const Material *GetMaterial() const = 0;
	virtual void ComputeScatteringFunctions(SurfaceInteraction *isect, TransportMode mode, bool allowMultipleLobes) const = 0;
};

class GeometricPrimitive : public Primitive {
public:
	virtual Bounds3f WorldBound() const;
	virtual bool Intersect(const Ray &r, SurfaceInteraction *isect) const;
	virtual bool IntersectP(const Ray &r) const;
	GeometricPrimitive(const std::shared_ptr<Shape> &shape,
		const std::shared_ptr<Material> &material,
	const std::shared_ptr<AreaLight> &areaLight);
	const AreaLight *GetAreaLight() const;
	const Material *GetMaterial() const;
	void ComputeScatteringFunctions(SurfaceInteraction *isect, TransportMode mode, bool allowMultipleLobes) const;
			
private:
	std::shared_ptr<Shape> shape;
	std::shared_ptr<Material> material;
	std::shared_ptr<AreaLight> areaLight;
	};

class Aggregate : public Primitive {
public:
	const AreaLight *GetAreaLight() const;
	const Material *GetMaterial() const;
	void ComputeScatteringFunctions(SurfaceInteraction *isect, TransportMode mode, bool allowMultipleLobes) const {};
};

class SimpleAccel : public Aggregate {
public:
	SimpleAccel(std::vector<std::shared_ptr<Primitive>> p);
	~SimpleAccel();
	bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
	bool IntersectP(const Ray &ray) const;
	Bounds3f WorldBound() const;
private:
	std::vector<std::shared_ptr<Primitive>> primitives;
	Bounds3f bound;
};












#endif