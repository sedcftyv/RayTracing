#include "sphere.h"
#include "sampling.h"
#include "interaction.h"

Bounds3f Sphere::ObjectBound() const {
	return Bounds3f(Point3f(-radius, -radius, -radius),
		Point3f(radius, radius, radius));
}

Float Sphere::Area() const { return 4 * Pi * radius * radius; }

inline bool Quadratic(Float A, Float B, Float C, Float *t0, Float *t1) {
	double discrim = (double)B * (double)B - 4. * (double)A * (double)C;
	if (discrim < 0.) return false;
	double rootDiscrim = std::sqrt(discrim);
	Float floatRootDiscrim=rootDiscrim;

	Float q;
	if ((float)B < 0)
		q = -.5 * (B - floatRootDiscrim);
	else
		q = -.5 * (B + floatRootDiscrim);
	*t0 = q / A;
	*t1 = C / q;
	if ((float)*t0 > (float)*t1) std::swap(*t0, *t1);
	return true;
}

bool Sphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
	bool testAlphaTexture) const {
	Ray ray = (*WorldToObject)(r);
		
	Vector3f oc = ray.o - Point3f(0.0f, 0.0f, 0.0f);
	Float a = ray.d.LengthSquared();
	Float half_b = Dot(oc, ray.d);
	Float c = oc.LengthSquared() - radius * radius;

	Float discriminant = half_b * half_b - a * c;
	if (discriminant < 0) return false;
		Float sqrtd = sqrt(discriminant);

		Float t0 = (-half_b - sqrtd) / a;
	Float t1 = (-half_b + sqrtd) / a;
	
	if (t0 > ray.tMax || t1 <= 0) return false;
	Float tShapeHit = t0;
	if (tShapeHit <= 0) {
		tShapeHit = t1;
		if (tShapeHit > ray.tMax) return false;
	}
	Point3f pHit = Point3f(ray((Float)tShapeHit));
		Float phi = std::atan2(pHit.y, pHit.x);
	if (phi < 0) phi += 2 * Pi;

	Float phiMax = 2 * Pi;
	Float u = phi / phiMax;
	Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
	Float v = theta / Pi;

	Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
	Float invZRadius = 1 / zRadius;
	Float cosPhi = pHit.x * invZRadius;
	Float sinPhi = pHit.y * invZRadius;
	Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
	Vector3f dpdv = Pi *Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));
	Normal3f n1 = Normal3f(Normalize(Vector3f(pHit)));
	*isect = (*ObjectToWorld)(SurfaceInteraction(pHit, Point2f(u, v),
		-ray.d, dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0),
		ray.time, this));
	isect->n = isect->shading.n = (*ObjectToWorld)(n1);
	*tHit = tShapeHit;
	return true;
}

bool Sphere::IntersectP(const Ray &r, bool testAlphaTexture) const {
	Ray ray = (*WorldToObject)(r);
	Vector3f oc = ray.o - Point3f(0.0f, 0.0f, 0.0f);
	Float a = ray.d.LengthSquared();
	Float half_b = Dot(oc, ray.d);
	Float c = oc.LengthSquared() - radius * radius;
	Float discriminant = half_b * half_b - a * c;
	if (discriminant < 0) return false;
	Float sqrtd = sqrt(discriminant);
	Float t0 = (-half_b - sqrtd) / a;
	Float t1 = (-half_b + sqrtd) / a;
	if (t0 > ray.tMax || t1 <= 0) return false;
	Float tShapeHit = t0;
	if (tShapeHit <= 0) {
		tShapeHit = t1;
		if (tShapeHit > ray.tMax) return false;
	}
	return true;
}


//}