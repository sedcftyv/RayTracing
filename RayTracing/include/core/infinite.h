#pragma once

#ifndef PBRT_LIGHTS_INFINITE_H
#define PBRT_LIGHTS_INFINITE_H

#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"


class InfiniteAreaLight : public Light {
public:
	// InfiniteAreaLight Public Methods
	InfiniteAreaLight(const Transform &LightToWorld, const Spectrum &power,
		int nSamples, const std::string &texmap);
	void Preprocess(const Scene &scene) {
		scene.WorldBound().BoundingSphere(&worldCenter, &worldRadius);
	}
	Spectrum Power() const;
	Spectrum Le(const Ray &ray) const;
	Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi,
		Float *pdf, VisibilityTester *vis) const;
	Float Pdf_Li(const Interaction &, const Vector3f &) const;
	Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
		Ray *ray, Normal3f *nLight, Float *pdfPos,
		Float *pdfDir) const;
	void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
		Float *pdfDir) const;

private:
	// InfiniteAreaLight Private Data
	std::unique_ptr<MIPMap<RGBSpectrum>> Lmap;
	Point3f worldCenter;
	Float worldRadius;
	std::unique_ptr<Distribution2D> distribution;
};

//std::shared_ptr<InfiniteAreaLight> CreateInfiniteLight(
//	const Transform &light2world, const ParamSet &paramSet);










#endif