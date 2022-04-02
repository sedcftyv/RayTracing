#pragma once
#ifndef PBRT_CORE_LIGHTDISTRIB_H
#define PBRT_CORE_LIGHTDISTRIB_H

#include "pbrt.h"
#include "geometry.h"
#include "sampling.h"

class LightDistribution {
public:
	virtual ~LightDistribution();

	// Given a point |p| in space, this method returns a (hopefully
	// effective) sampling distribution for light sources at that point.
	virtual const Distribution1D *Lookup(const Point3f &p) const = 0;
};

std::unique_ptr<LightDistribution> CreateLightSampleDistribution(
	const std::string &name, const Scene &scene);

class UniformLightDistribution : public LightDistribution {
public:
	UniformLightDistribution(const Scene &scene);
	const Distribution1D *Lookup(const Point3f &p) const;

private:
	std::unique_ptr<Distribution1D> distrib;
};









#endif