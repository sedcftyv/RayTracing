#include "light.h"
#include "interaction.h"
#include "sampling.h"
#include "geometry.h"
#include "scene.h"
static long long numLights = 0;
static long long numAreaLights = 0;

Light::Light(int flags, const Transform &LightToWorld, int nSamples)
	: flags(flags),
	nSamples(std::max(1, nSamples)),
	LightToWorld(LightToWorld),
	WorldToLight(Inverse(LightToWorld)) {
	++numLights;
}

Light::~Light() {}

bool VisibilityTester::Unoccluded(const Scene &scene) const {
	return !scene.IntersectP(p0.SpawnRayTo(p1));
}

AreaLight::AreaLight(const Transform &LightToWorld, int nSamples)
	: Light((int)LightFlags::Area, LightToWorld, nSamples) {
	++numAreaLights;
}

DiffuseAreaLight::DiffuseAreaLight(const Transform &LightToWorld,
	const Spectrum &Lemit, int nSamples,
	const std::shared_ptr<Shape> &shape,
	bool twoSided)
	: AreaLight(LightToWorld, nSamples),
	Lemit(Lemit),
	shape(shape),
	twoSided(twoSided),
	area(shape->Area()) {
									}

Spectrum DiffuseAreaLight::Power() const {
	return (twoSided ? 2 : 1) * Lemit * area * Pi;
}

Spectrum DiffuseAreaLight::Sample_Li(const Interaction &ref, const Point2f &u,
	Vector3f *wi, Float *pdf,
	VisibilityTester *vis) const {
	Interaction pShape = shape->Sample(ref, u, pdf);
	if (*pdf == 0 || (pShape.p - ref.p).LengthSquared() == 0) {
	*pdf = 0;
	return 0.f;
	}
	*wi = Normalize(pShape.p - ref.p);
	*vis = VisibilityTester(ref, pShape);
	return L(pShape, -*wi);
}

Float DiffuseAreaLight::Pdf_Li(const Interaction &ref,
	const Vector3f &wi) const {
	return shape->Pdf(ref, wi);
}

Spectrum DiffuseAreaLight::Sample_Le(const Point2f &u1, const Point2f &u2,
	Float time, Ray *ray, Normal3f *nLight,
	Float *pdfPos, Float *pdfDir) const {
	Interaction pShape = shape->Sample(u1, pdfPos);
	*nLight = pShape.n;

	Vector3f w;
	if (twoSided) {
		Point2f u = u2;
		if (u[0] < .5) {
			u[0] = std::min(u[0] * 2, OneMinusEpsilon);
			w = CosineSampleHemisphere(u);
		}
		else {
			u[0] = std::min((u[0] - .5f) * 2, OneMinusEpsilon);
			w = CosineSampleHemisphere(u);
			w.z *= -1;
		}
		*pdfDir = 0.5f * CosineHemispherePdf(std::abs(w.z));
	}
	else {
		w = CosineSampleHemisphere(u2);
		*pdfDir = CosineHemispherePdf(w.z);
	}

	Vector3f v1, v2, n(pShape.n);
	CoordinateSystem(n, &v1, &v2);
	w = w.x * v1 + w.y * v2 + w.z * n;
	*ray = pShape.SpawnRay(w);
	return L(pShape, w);
}

void DiffuseAreaLight::Pdf_Le(const Ray &ray, const Normal3f &n, Float *pdfPos,
	Float *pdfDir) const {
	Interaction it(ray.o, n,  Vector3f(n), ray.time);
	*pdfPos = shape->Pdf(it);
	*pdfDir = twoSided ? (.5 * CosineHemispherePdf(AbsDot(n, ray.d))): CosineHemispherePdf(Dot(n, ray.d));
}

