#pragma once

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "scene.h"
//#include "light.h"
//#include "reflection.h"
#include "sampler.h"
//#include "material.h"
#include "stratified.h"
#include "geometry.h"
#include "camera.h"
#include "triangle.h"
#include "lightdistrib.h"

class Integrator {
public:
	// Integrator Interface
	virtual ~Integrator();
	virtual void Render(const Scene &scene) = 0;
};

Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene, Sampler &sampler,
	const std::vector<int> &nLightSamples,
	bool handleMedia = false);
Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene, Sampler &sampler,
	bool handleMedia = false,
	const Distribution1D *lightDistrib = nullptr);
Spectrum EstimateDirect(const Interaction &it, const Point2f &uShading,
	const Light &light, const Point2f &uLight,
	const Scene &scene, Sampler &sampler,bool handleMedia = false,
	bool specular = false);
//std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
//	const Scene &scene);

// SamplerIntegrator Declarations
class SamplerIntegrator : public Integrator {
public:
	// SamplerIntegrator Public Methods
	SamplerIntegrator(std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i &pixelBounds)
		: camera(camera), sampler(sampler), pixelBounds(pixelBounds) {}
	virtual void Preprocess(const Scene &scene, Sampler &sampler) {}
	void Render(const Scene &scene);
	virtual Spectrum Li(const Ray &ray, const Scene &scene,
		Sampler &sampler, int depth = 0) const = 0;
	Spectrum SpecularReflect(const Ray &ray,
		const SurfaceInteraction &isect,
		const Scene &scene, Sampler &sampler,int depth) const;
	Spectrum SpecularTransmit(const Ray &ray,
		const SurfaceInteraction &isect,
		const Scene &scene, Sampler &sampler, int depth) const;

protected:
	// SamplerIntegrator Protected Data
	std::shared_ptr<const Camera> camera;

private:
	// SamplerIntegrator Private Data
	std::shared_ptr<Sampler> sampler;
	const Bounds2i pixelBounds;
};

class WhittedIntegrator : public SamplerIntegrator {
public:
	// WhittedIntegrator Public Methods
	WhittedIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i &pixelBounds)
		: SamplerIntegrator(camera, sampler, pixelBounds), maxDepth(maxDepth) {}
	Spectrum Li(const Ray &ray, const Scene &scene,Sampler &sampler, int depth) const;

private:
	// WhittedIntegrator Private Data
	const int maxDepth;
};

class PathIntegrator : public SamplerIntegrator {
public:
	// PathIntegrator Public Methods
	PathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i &pixelBounds, Float rrThreshold = 1,
		const std::string &lightSampleStrategy = "spatial") :SamplerIntegrator(camera, sampler, pixelBounds),
		maxDepth(maxDepth),
		rrThreshold(rrThreshold),
		lightSampleStrategy(lightSampleStrategy) {}

	void Preprocess(const Scene &scene, Sampler &sampler);
	Spectrum Li(const Ray &ray, const Scene &scene,
		Sampler &sampler, int depth) const;

private:
	// PathIntegrator Private Data
	const int maxDepth;
	const Float rrThreshold;
	const std::string lightSampleStrategy;
	std::unique_ptr<LightDistribution> lightDistribution;
};





#endif