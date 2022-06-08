#pragma once

#ifndef PBRT_CORE_SCENE_H
#define PBRT_CORE_SCENE_H

#include "pbrt.h"
#include "geometry.h"
#include "primitive.h"
#include "sampler.h"

class Scene {
public:
	Scene(std::shared_ptr<Primitive> aggregate, const std::vector<std::shared_ptr<Light>> &lights);
	const Bounds3f &WorldBound() const { return worldBound; }
	bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
	bool IntersectP(const Ray &ray) const;
	
	std::vector<std::shared_ptr<Light>> lights;
	std::vector<std::shared_ptr<Light>> infiniteLights;

private:
	std::shared_ptr<Primitive> aggregate;
	Bounds3f worldBound;
};






#endif