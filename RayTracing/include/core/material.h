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

class MirrorMaterial : public Material {
public:
	// MirrorMaterial Public Methods
	MirrorMaterial(const std::shared_ptr<Texture<Spectrum>> &r,
		const std::shared_ptr<Texture<Float>> &bump) {
		Kr = r;
		bumpMap = bump;
	}
	void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,
		bool allowMultipleLobes) const;

private:
	// MirrorMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> Kr;
	std::shared_ptr<Texture<Float>> bumpMap;
};

class GlassMaterial : public Material {
public:
	// GlassMaterial Public Methods
	GlassMaterial(const std::shared_ptr<Texture<Spectrum>> &Kr,
		const std::shared_ptr<Texture<Spectrum>> &Kt,
		//const std::shared_ptr<Texture<Float>> &uRoughness,
		//const std::shared_ptr<Texture<Float>> &vRoughness,
		const std::shared_ptr<Texture<Float>> &index,
		const std::shared_ptr<Texture<Float>> &bumpMap)
		: Kr(Kr),
		Kt(Kt),
		index(index),
		bumpMap(bumpMap){}
	void ComputeScatteringFunctions(SurfaceInteraction *si,
		TransportMode mode,
		bool allowMultipleLobes) const;

private:
	// GlassMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> Kr, Kt;
	//std::shared_ptr<Texture<Float>> uRoughness, vRoughness;
	std::shared_ptr<Texture<Float>> index;
	std::shared_ptr<Texture<Float>> bumpMap;
	//bool remapRoughness;
};



#endif