#pragma once
#ifndef PBRT_CORE_MATERIAL_H
#define PBRT_CORE_MATERIAL_H

#include "pbrt.h"

enum class TransportMode {Radiance, Importance};

class Material {
public:
	virtual void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,bool allowMultipleLobes) const = 0;
	virtual ~Material() {};
};

class MatteMaterial : public Material {
public:
	MatteMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
	const std::shared_ptr<Texture<Float>> &sigma): Kd(Kd), sigma(sigma) {}
	void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,bool allowMultipleLobes) const;
private:
	std::shared_ptr<Texture<Spectrum>> Kd;
	std::shared_ptr<Texture<Float>> sigma;};

class MirrorMaterial : public Material {
public:
	MirrorMaterial(const std::shared_ptr<Texture<Spectrum>> &r):Kr(r){}
	void ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,bool allowMultipleLobes) const;
private:
	std::shared_ptr<Texture<Spectrum>> Kr;
};

class GlassMaterial : public Material {
public:
	GlassMaterial(const std::shared_ptr<Texture<Spectrum>> &Kr,
	const std::shared_ptr<Texture<Spectrum>> &Kt,const std::shared_ptr<Texture<Float>> &index)
		: Kr(Kr),Kt(Kt),index(index) {}
	void ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,
	bool allowMultipleLobes) const;
private:
	std::shared_ptr<Texture<Spectrum>> Kr, Kt;
	std::shared_ptr<Texture<Float>> index;
};

class PlasticMaterial : public Material {
public:
	PlasticMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,const std::shared_ptr<Texture<Spectrum>> &Ks,
		const std::shared_ptr<Texture<Float>> &roughness,bool remapRoughness)
		: Kd(Kd),Ks(Ks),roughness(roughness),remapRoughness(remapRoughness) {}
	void ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,
	bool allowMultipleLobes) const;

private:
	std::shared_ptr<Texture<Spectrum>> Kd, Ks;
	std::shared_ptr<Texture<Float>> roughness;
	const bool remapRoughness;
};

class MetalMaterial : public Material {
public:
	MetalMaterial(const std::shared_ptr<Texture<Spectrum>> &eta,
	const std::shared_ptr<Texture<Spectrum>> &k,
	const std::shared_ptr<Texture<Float>> &rough,
	const std::shared_ptr<Texture<Float>> &urough,
	const std::shared_ptr<Texture<Float>> &vrough,
			bool remapRoughness):eta(eta),
	k(k),
	roughness(rough),
	uRoughness(urough),
	vRoughness(vrough),
			remapRoughness(remapRoughness) {}
	void ComputeScatteringFunctions(SurfaceInteraction *si, TransportMode mode,
		bool allowMultipleLobes) const;

private:
	std::shared_ptr<Texture<Spectrum>> eta, k;
	std::shared_ptr<Texture<Float>> roughness, uRoughness, vRoughness;
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