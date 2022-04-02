#pragma once
#ifndef PBRT_CORE_LIGHT_H
#define PBRT_CORE_LIGHT_H

#include "pbrt.h"
#include "interaction.h"
#include "spectrum.h"



enum class LightFlags : int {
	DeltaPosition = 1,
	DeltaDirection = 2,
	Area = 4,
	Infinite = 8
};

inline bool IsDeltaLight(int flags) {
	return flags & (int)LightFlags::DeltaPosition ||
		flags & (int)LightFlags::DeltaDirection;
}

class Light {
public:
	// Light Interface
	virtual ~Light();
	Light(int flags, const Transform &LightToWorld, int nSamples = 1);
	virtual Spectrum Sample_Li(const Interaction &ref, const Point2f &u,
		Vector3f *wi, Float *pdf, VisibilityTester *vis) const = 0;
	virtual Spectrum Power() const = 0;
	virtual void Preprocess(const Scene &scene) {}
	virtual Spectrum Le(const Ray &r) const {
		return Spectrum(0.f);
	}

	virtual Float Pdf_Li(const Interaction &ref, const Vector3f &wi) const = 0;
	virtual Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
		Ray *ray, Normal3f *nLight, Float *pdfPos,
		Float *pdfDir) const = 0;
	virtual void Pdf_Le(const Ray &ray, const Normal3f &nLight, Float *pdfPos,
		Float *pdfDir) const = 0;

	// Light Public Data
	const int flags;
	const int nSamples;
	//const MediumInterface mediumInterface;

protected:
	// Light Protected Data
	const Transform LightToWorld, WorldToLight;
};

class VisibilityTester {
public:
	VisibilityTester() {}
	// VisibilityTester Public Methods
	VisibilityTester(const Interaction &p0, const Interaction &p1)
		: p0(p0), p1(p1) {}
	const Interaction &P0() const { return p0; }
	const Interaction &P1() const { return p1; }
	bool Unoccluded(const Scene &scene) const;
	Spectrum Tr(const Scene &scene, Sampler &sampler) const {
		return Spectrum(0.f);
	}

private:
	Interaction p0, p1;
};

class AreaLight : public Light {
public:
	// AreaLight Interface
	AreaLight(const Transform &LightToWorld,int nSamples);
	virtual Spectrum L(const Interaction &intr, const Vector3f &w) const = 0;
};

class DiffuseAreaLight : public AreaLight {
public:
	// DiffuseAreaLight Public Methods
	DiffuseAreaLight(const Transform &LightToWorld, const Spectrum &Le,
		int nSamples, const std::shared_ptr<Shape> &shape,
		bool twoSided = false);
	Spectrum L(const Interaction &intr, const Vector3f &w) const {
		return (twoSided || Dot(intr.n, w) > 0) ? Lemit : Spectrum(0.f);
	}
	Spectrum Power() const;
	Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wo,
		Float *pdf, VisibilityTester *vis) const;
	Float Pdf_Li(const Interaction &, const Vector3f &) const;
	Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
		Ray *ray, Normal3f *nLight, Float *pdfPos,
		Float *pdfDir) const;
	void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
		Float *pdfDir) const;

protected:
	// DiffuseAreaLight Protected Data
	const Spectrum Lemit;
	std::shared_ptr<Shape> shape;
	// Added after book publication: by default, DiffuseAreaLights still
	// only emit in the hemimsphere around the surface normal.  However,
	// this behavior can now be overridden to give emission on both sides.
	const bool twoSided;
	const Float area;
};

class PointLight : public Light {
public:
	// PointLight Public Methods
	PointLight(const Transform &LightToWorld, const Spectrum &I)
		: Light((int)LightFlags::DeltaPosition, LightToWorld),
		pLight(LightToWorld(Point3f(0, 0, 0))),
		I(I) {}
	Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi,
		Float *pdf, VisibilityTester *vis) const;
	Spectrum Power() const;
	Float Pdf_Li(const Interaction &, const Vector3f &) const;
	Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
		Ray *ray, Normal3f *nLight, Float *pdfPos,
		Float *pdfDir) const;
	void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
		Float *pdfDir) const;

private:
	// PointLight Private Data
	const Point3f pLight;
	const Spectrum I;
};




#endif