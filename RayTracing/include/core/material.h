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
		const std::shared_ptr<Texture<Float>> &sigma)
		//const std::shared_ptr<Texture<Float>> &bumpMap)
		: Kd(Kd), sigma(sigma) {}
	//bumpMap(bumpMap) 
	void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,bool allowMultipleLobes) const;

private:
	// MatteMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> Kd;
	std::shared_ptr<Texture<Float>> sigma;// bumpMap;
};

class MirrorMaterial : public Material {
public:
	// MirrorMaterial Public Methods
	MirrorMaterial(const std::shared_ptr<Texture<Spectrum>> &r)
		//const std::shared_ptr<Texture<Float>> &bump
		:Kr(r){}
	void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,
		bool allowMultipleLobes) const;

private:
	// MirrorMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> Kr;
	//std::shared_ptr<Texture<Float>> bumpMap;
};

class GlassMaterial : public Material {
public:
	// GlassMaterial Public Methods
	GlassMaterial(const std::shared_ptr<Texture<Spectrum>> &Kr,
		const std::shared_ptr<Texture<Spectrum>> &Kt,
		//const std::shared_ptr<Texture<Float>> &uRoughness,
		//const std::shared_ptr<Texture<Float>> &vRoughness,
		const std::shared_ptr<Texture<Float>> &index)
		//const std::shared_ptr<Texture<Float>> &bumpMap)
		: Kr(Kr),
		Kt(Kt),
		index(index) {}
		//bumpMap(bumpMap)
	void ComputeScatteringFunctions(SurfaceInteraction *si,
		TransportMode mode,
		bool allowMultipleLobes) const;

private:
	// GlassMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> Kr, Kt;
	//std::shared_ptr<Texture<Float>> uRoughness, vRoughness;
	std::shared_ptr<Texture<Float>> index;
	//std::shared_ptr<Texture<Float>> bumpMap;
	//bool remapRoughness;
};

class PlasticMaterial : public Material {
public:
	// PlasticMaterial Public Methods
	PlasticMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
		const std::shared_ptr<Texture<Spectrum>> &Ks,
		const std::shared_ptr<Texture<Float>> &roughness,
		//const std::shared_ptr<Texture<Float>> &bumpMap,
		bool remapRoughness)
		: Kd(Kd),
		Ks(Ks),
		roughness(roughness),
		//bumpMap(bumpMap),
		remapRoughness(remapRoughness) {}
	void ComputeScatteringFunctions(SurfaceInteraction *si,
		TransportMode mode,
		bool allowMultipleLobes) const;

private:
	// PlasticMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> Kd, Ks;
	std::shared_ptr<Texture<Float>> roughness;
	//std::shared_ptr<Texture<Float>> bumpMap;
	const bool remapRoughness;
};

class MetalMaterial : public Material {
public:
	// MetalMaterial Public Methods
	MetalMaterial(const std::shared_ptr<Texture<Spectrum>> &eta,
		const std::shared_ptr<Texture<Spectrum>> &k,
		const std::shared_ptr<Texture<Float>> &rough,
		const std::shared_ptr<Texture<Float>> &urough,
		const std::shared_ptr<Texture<Float>> &vrough,
		//const std::shared_ptr<Texture<Float>> &bump,
		bool remapRoughness):eta(eta),
		k(k),
		roughness(rough),
		uRoughness(urough),
		vRoughness(vrough),
		//bumpMap(bump),
		remapRoughness(remapRoughness) {}
	void ComputeScatteringFunctions(SurfaceInteraction *si, 
		TransportMode mode,
		bool allowMultipleLobes) const;

private:
	// MetalMaterial Private Data
	std::shared_ptr<Texture<Spectrum>> eta, k;
	std::shared_ptr<Texture<Float>> roughness, uRoughness, vRoughness;
	//std::shared_ptr<Texture<Float>> bumpMap;
	bool remapRoughness;
};

class MetalRoughnessMaterial : public Material {
public:

	MetalRoughnessMaterial(const std::shared_ptr<Texture<Spectrum>> &base_color,
		const std::shared_ptr<Texture<Spectrum>> &metalroughness) 
		:base_color(base_color),metalroughness(metalroughness){}
	void ComputeScatteringFunctions(SurfaceInteraction *si,
		TransportMode mode,
		bool allowMultipleLobes) const;

private:

	std::shared_ptr<Texture<Spectrum>> base_color,metalroughness;
};
#endif