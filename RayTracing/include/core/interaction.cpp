#include "interaction.h"
#include "shape.h"
#include "reflection.h"
#include "primitive.h"
#include "light.h"

SurfaceInteraction::~SurfaceInteraction(){
	if (bsdf)
	delete bsdf;
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
		shading.n = n;
		shading.dpdu = dpdu;
		shading.dpdv = dpdv;
		shading.dndu = dndu;
		shading.dndv = dndv;
		if (shape &&(shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
		n *= -1;
		shading.n *= -1;
	}
}

void SurfaceInteraction::SetShadingGeometry(const Vector3f &dpdus,const Vector3f &dpdvs,const Normal3f &dndus,const Normal3f &dndvs,
bool orientationIsAuthoritative) {
	shading.n = Normalize((Normal3f)Cross(dpdus, dpdvs));
	if (orientationIsAuthoritative)
		n = Faceforward(n, shading.n);
	else
		shading.n = Faceforward(shading.n, n);
		shading.dpdu = dpdus;
		shading.dpdv = dpdvs;
		shading.dndu = dndus;
		shading.dndv = dndvs;
}

void SurfaceInteraction::ComputeScatteringFunctions(const Ray &ray,bool allowMultipleLobes,TransportMode mode) {
	primitive->ComputeScatteringFunctions(this, mode, allowMultipleLobes);
}

Spectrum SurfaceInteraction::Le(const Vector3f &w) const {
	const AreaLight *area = primitive->GetAreaLight();
	return area ? area->L(*this, w) : Spectrum(0.f);
}