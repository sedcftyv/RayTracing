#include "sphere.h"
#include "sampling.h"
#include "interaction.h"

Bounds3f Sphere::ObjectBound() const {
	return Bounds3f(Point3f(-radius, -radius, -radius),
		Point3f(radius, radius, radius));
}

Float Sphere::Area() const { return 4 * Pi * radius * radius; }

bool Sphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
	bool testAlphaTexture) const {
	
	Point3f pHit;
	Ray ray = (*WorldToObject)(r);
	Vector3f oc = ray.o - Point3f(0.0f, 0.0f, 0.0f);
	// Compute quadratic sphere coefficients

	// Initialize _EFloat_ ray coordinate values
	Float a = Dot(ray.d,ray.d);
	Float b = 2.0*Dot(oc, ray.d);
	Float c = Dot(oc, oc) - radius * radius;
	Float discriminant = b * b - 4 * a*c;
	return discriminant > 0;
}

bool Sphere::IntersectP(const Ray &r, bool testAlphaTexture) const {

	Point3f pHit;
	Ray ray = (*WorldToObject)(r);
	Vector3f oc = ray.o - Point3f(0.0f, 0.0f, 0.0f);
	// Compute quadratic sphere coefficients

	// Initialize _EFloat_ ray coordinate values
	Float a = Dot(ray.d, ray.d);
	Float b = 2.0*Dot(oc, ray.d);
	Float c = Dot(oc, oc) - radius * radius;
	Float discriminant = b * b - 4 * a*c;
	return discriminant > 0;
}

Interaction Sphere::Sample(const Point2f &u, Float *pdf) const {
	Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
	Interaction it;
	it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
	if (reverseOrientation) it.n *= -1;
	// Reproject _pObj_ to sphere surface and compute _pObjError_
	pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
	Vector3f pObjError = gamma(5) * Abs((Vector3f)pObj);
	it.p = (*ObjectToWorld)(pObj);
	*pdf = 1 / Area();
	return it;
}

Interaction Sphere::Sample(const Interaction &ref, const Point2f &u,
	Float *pdf) const {
	Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));

	// Sample uniformly on sphere if Undefined control sequence \pt is inside it
	Point3f pOrigin =
		OffsetRayOrigin(ref.p,  ref.n, pCenter - ref.p);
	if (DistanceSquared(pOrigin, pCenter) <= radius * radius) {
		Interaction intr = Sample(u, pdf);
		Vector3f wi = intr.p - ref.p;
		if (wi.LengthSquared() == 0)
			*pdf = 0;
		else {
			// Convert from area measure returned by Sample() call above to
			// solid angle measure.
			wi = Normalize(wi);
			*pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
		}
		if (std::isinf(*pdf)) *pdf = 0.f;
		return intr;
	}

	// Sample sphere uniformly inside subtended cone

	// Compute coordinate system for sphere sampling
	Float dc = Distance(ref.p, pCenter);
	Float invDc = 1 / dc;
	Vector3f wc = (pCenter - ref.p) * invDc;
	Vector3f wcX, wcY;
	CoordinateSystem(wc, &wcX, &wcY);

	// Compute  and  values for sample in cone
	Float sinThetaMax = radius * invDc;
	Float sinThetaMax2 = sinThetaMax * sinThetaMax;
	Float invSinThetaMax = 1 / sinThetaMax;
	Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));

	Float cosTheta = (cosThetaMax - 1) * u[0] + 1;
	Float sinTheta2 = 1 - cosTheta * cosTheta;

	if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
		/* Fall back to a Taylor series expansion for small angles, where
		   the standard approach suffers from severe cancellation errors */
		sinTheta2 = sinThetaMax2 * u[0];
		cosTheta = std::sqrt(1 - sinTheta2);
	}

	// Compute angle ¦Á from center of sphere to sampled point on surface
	Float cosAlpha = sinTheta2 * invSinThetaMax +
		cosTheta * std::sqrt(std::max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
	Float sinAlpha = std::sqrt(std::max((Float)0.f, 1.f - cosAlpha * cosAlpha));
	Float phi = u[1] * 2 * Pi;

	// Compute surface normal and sampled point on sphere
	Vector3f nWorld =
		SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
	Point3f pWorld = pCenter + radius * Point3f(nWorld.x, nWorld.y, nWorld.z);

	// Return _Interaction_ for sampled point on sphere
	Interaction it;
	it.p = pWorld;
	//it.pError = gamma(5) * Abs((Vector3f)pWorld);
	it.n = Normal3f(nWorld);
	if (reverseOrientation) it.n *= -1;

	// Uniform cone PDF.
	*pdf = 1 / (2 * Pi * (1 - cosThetaMax));

	return it;
}

Float Sphere::Pdf(const Interaction &ref, const Vector3f &wi) const {
	Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
	// Return uniform PDF if point is inside sphere
	Point3f pOrigin =
		OffsetRayOrigin(ref.p, ref.n, pCenter - ref.p);
	if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
		return Shape::Pdf(ref, wi);

	// Compute general sphere PDF
	Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
	Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
	return UniformConePdf(cosThetaMax);
}