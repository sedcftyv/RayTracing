#include "interaction.h"
#include "shape.h"
#include "reflection.h"
#include "primitive.h"
#include "light.h"

SurfaceInteraction::~SurfaceInteraction()
{
	if (bsdf)bsdf->~BSDF();
}
SurfaceInteraction::SurfaceInteraction(
	const Point3f &p,const Point2f &uv,
	const Vector3f &wo, const Vector3f &dpdu, const Vector3f &dpdv,
	const Normal3f &dndu, const Normal3f &dndv, Float time, const Shape *shape,
	int faceIndex)
	: Interaction(p, Normal3f(Normalize(Cross(dpdu, dpdv))), wo, time),
	uv(uv),
	dpdu(dpdu),
	dpdv(dpdv),
	dndu(dndu),
	dndv(dndv),
	shape(shape),
	faceIndex(faceIndex) {
	// Initialize shading geometry from true geometry
	shading.n = n;
	shading.dpdu = dpdu;
	shading.dpdv = dpdv;
	shading.dndu = dndu;
	shading.dndv = dndv;

	// Adjust normal based on orientation and handedness
	if (shape &&
		(shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
		n *= -1;
		shading.n *= -1;
	}
}

void SurfaceInteraction::ComputeScatteringFunctions(const Ray &ray,
	bool allowMultipleLobes,
	TransportMode mode) {
	//ComputeDifferentials(ray);
	primitive->ComputeScatteringFunctions(this, mode, allowMultipleLobes);
}

Spectrum SurfaceInteraction::Le(const Vector3f &w) const {
	const AreaLight *area = primitive->GetAreaLight();
	return area ? area->L(*this, w) : Spectrum(0.f);
}