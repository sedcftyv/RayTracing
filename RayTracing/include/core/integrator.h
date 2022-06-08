#pragma once

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "scene.h"
#include "sampler.h"
#include "stratified.h"
#include "geometry.h"
#include "camera.h"
#include "triangle.h"
#include "lightdistrib.h"

class Integrator {
public:
	virtual ~Integrator();
	virtual void Render(const Scene &scene, int image_width, int image_height) = 0;
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
class SamplerIntegrator : public Integrator {
public:
	SamplerIntegrator(std::shared_ptr<const Camera> camera,std::shared_ptr<Sampler> sampler,
	const Bounds2i &pixelBounds)
	: camera(camera), sampler(sampler), pixelBounds(pixelBounds) {}
	virtual void Preprocess(const Scene &scene, Sampler &sampler) {}
	void Render(const Scene &scene, int image_width, int image_height);
	virtual Spectrum Li(const Ray &ray, const Scene &scene,
	Sampler &sampler, int depth = 0) const = 0;
	virtual Spectrum Li_re(const Ray &ray, const Scene &scene,
	Sampler &sampler, int depth = 0) const;
	Spectrum SpecularReflect(const Ray &ray,
	const SurfaceInteraction &isect,
	const Scene &scene, Sampler &sampler,int depth) const;
	Spectrum SpecularTransmit(const Ray &ray,
	const SurfaceInteraction &isect,
	const Scene &scene, Sampler &sampler, int depth) const;
	void render_pixel(int nowpos);
protected:
	std::shared_ptr<const Camera> camera;

private:
	std::shared_ptr<Sampler> sampler;
	const Bounds2i pixelBounds;
};

class WhittedIntegrator : public SamplerIntegrator {
public:
	WhittedIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
	std::shared_ptr<Sampler> sampler,
	const Bounds2i &pixelBounds): SamplerIntegrator(camera, sampler, pixelBounds), maxDepth(maxDepth) {}
	Spectrum Li(const Ray &ray, const Scene &scene,Sampler &sampler, int depth) const;
private:
		const int maxDepth;
};

class PathIntegrator : public SamplerIntegrator {
public:
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

	Spectrum Li_re(const Ray &ray, const Scene &scene,
	Sampler &sampler, int depth) const;

private:
	const int maxDepth;
	const Float rrThreshold;
	const std::string lightSampleStrategy;
	std::unique_ptr<LightDistribution> lightDistribution;
};





#endif