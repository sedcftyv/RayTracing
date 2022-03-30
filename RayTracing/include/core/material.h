#pragma once
#ifndef PBRT_CORE_MATERIAL_H
#define PBRT_CORE_MATERIAL_H

#include "pbrt.h"

enum class TransportMode {Radiance, Importance};

class Material {
public:
	// Material Interface
	virtual void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,bool allowMultipleLobes) const = 0;
	virtual ~Material() {};
	//static void Bump(const std::shared_ptr<Texture<Float>> &d,SurfaceInteraction *si);
};

class MatteMaterial : public Material {
public:
	// MatteMaterial Public Methods
	MatteMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
		const std::shared_ptr<Texture<Float>> &sigma,
		const std::shared_ptr<Texture<Float>> &bumpMap)
		: Kd(Kd), sigma(sigma), bumpMap(bumpMap) {}
	void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,bool allowMultipleLobes) const;

private:
	// MatteMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> Kd;
	std::shared_ptr<Texture<Float>> sigma, bumpMap;
};







#endif