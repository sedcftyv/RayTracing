#include "triangle.h"
#include "sampling.h"
static long long nTris = 0;
static long long nMeshes = 0;
static long long nHits = 0;
static long long nTests = 0;

TriangleMesh::TriangleMesh(
	const Transform &ObjectToWorld, int nTriangles, const int *vertexIndices,
	int nVertices, const Point3f *P, const Vector3f *S, const Normal3f *N,
	const Point2f *UV,const int *fIndices)
	: nTriangles(nTriangles),
	nVertices(nVertices),
	vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles){
	++nMeshes;
	nTris += nTriangles;
	triMeshBytes += sizeof(*this) + this->vertexIndices.size() * sizeof(int) +nVertices * (sizeof(*P) + (N ? sizeof(*N) : 0) +
		(S ? sizeof(*S) : 0) + (UV ? sizeof(*UV) : 0) +(fIndices ? sizeof(*fIndices) : 0));
	p.reset(new Point3f[nVertices]);
	for (int i = 0; i < nVertices; ++i) p[i] = ObjectToWorld(P[i]);
	if (UV) {
		uv.reset(new Point2f[nVertices]);
		memcpy(uv.get(), UV, nVertices * sizeof(Point2f));
	}
	if (N) {
		n.reset(new Normal3f[nVertices]);
		for (int i = 0; i < nVertices; ++i) n[i] = ObjectToWorld(N[i]);
	}
	if (S) {
		s.reset(new Vector3f[nVertices]);
		for (int i = 0; i < nVertices; ++i) s[i] = ObjectToWorld(S[i]);
	}

	if (fIndices)
		faceIndices = std::vector<int>(fIndices, fIndices + nTriangles);
}

Bounds3f Triangle::ObjectBound() const {
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),(*WorldToObject)(p2));
}

Bounds3f Triangle::WorldBound() const {
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	return Union(Bounds3f(p0, p1), p2);
}

bool Triangle::Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
	bool testAlphaTexture) const {		
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	Point3f p0t = p0 - Vector3f(ray.o);
	Point3f p1t = p1 - Vector3f(ray.o);
	Point3f p2t = p2 - Vector3f(ray.o);

	int kz = MaxDimension(Abs(ray.d));
	int kx = kz + 1;
	if (kx == 3) kx = 0;
	int ky = kx + 1;
	if (ky == 3) ky = 0;
	Vector3f d = Permute(ray.d, kx, ky, kz);
	p0t = Permute(p0t, kx, ky, kz);
	p1t = Permute(p1t, kx, ky, kz);
	p2t = Permute(p2t, kx, ky, kz);

	Float Sx = -d.x / d.z;
	Float Sy = -d.y / d.z;
	Float Sz = 1.f / d.z;
	p0t.x += Sx * p0t.z;
	p0t.y += Sy * p0t.z;
	p1t.x += Sx * p1t.z;
	p1t.y += Sy * p1t.z;
	p2t.x += Sx * p2t.z;
	p2t.y += Sy * p2t.z;

	Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
	Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
	Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

	if (sizeof(Float) == sizeof(float) &&
		(e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
		double p2txp1ty = (double)p2t.x * (double)p1t.y;
		double p2typ1tx = (double)p2t.y * (double)p1t.x;
		e0 = (float)(p2typ1tx - p2txp1ty);
		double p0txp2ty = (double)p0t.x * (double)p2t.y;
		double p0typ2tx = (double)p0t.y * (double)p2t.x;
		e1 = (float)(p0typ2tx - p0txp2ty);
		double p1txp0ty = (double)p1t.x * (double)p0t.y;
		double p1typ0tx = (double)p1t.y * (double)p0t.x;
		e2 = (float)(p1typ0tx - p1txp0ty);
	}

	if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
		return false;
	Float det = e0 + e1 + e2;
	if (det == 0) return false;

	p0t.z *= Sz;
	p1t.z *= Sz;
	p2t.z *= Sz;
	Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
	if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
		return false;
	else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
		return false;

	Float invDet = 1 / det;
	Float b0 = e0 * invDet;
	Float b1 = e1 * invDet;
	Float b2 = e2 * invDet;
	Float t = tScaled * invDet;

	if (t <= eps) return false;

	Vector3f dpdu, dpdv;
	Point2f uv[3];
	GetUVs(uv);

	Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
	Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
	Point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
	Point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;

	Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
	bool degenerateUV = std::abs(determinant) < 1e-8;
	if (!degenerateUV)
	{
		Float invdet = 1 / determinant;
		dpdu = (duv12[1] * dp02 - duv02[1] * dp12)*invdet;
		dpdv = (-duv12[0] * dp02 + duv02[0] * dp12)*invdet;
	}

	if (degenerateUV || Cross(dpdu, dpdv).LengthSquared() == 0) {
		Vector3f ng = Cross(p2 - p0, p1 - p0);
		if (ng.LengthSquared() == 0)
			return false;
		CoordinateSystem(Normalize(ng), &dpdu, &dpdv);
	}
	*isect = SurfaceInteraction(pHit,uvHit,-ray.d,dpdu,dpdv,Normal3f(0,0,0),Normal3f(0,0,0),ray.time,this,faceIndex);
	isect->n = isect->shading.n = Normal3f(Normalize(Cross(dp02, dp12)));																																								
	*tHit = t;
	return true;
}

bool Triangle::IntersectP(const Ray &ray, bool testAlphaTexture) const {
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];

	Point3f p0t = p0 - Vector3f(ray.o);
	Point3f p1t = p1 - Vector3f(ray.o);
	Point3f p2t = p2 - Vector3f(ray.o);

	int kz = MaxDimension(Abs(ray.d));
	int kx = kz + 1;
	if (kx == 3) kx = 0;
	int ky = kx + 1;
	if (ky == 3) ky = 0;
	Vector3f d = Permute(ray.d, kx, ky, kz);
	p0t = Permute(p0t, kx, ky, kz);
	p1t = Permute(p1t, kx, ky, kz);
	p2t = Permute(p2t, kx, ky, kz);

	Float Sx = -d.x / d.z;
	Float Sy = -d.y / d.z;
	Float Sz = 1.f / d.z;
	p0t.x += Sx * p0t.z;
	p0t.y += Sy * p0t.z;
	p1t.x += Sx * p1t.z;
	p1t.y += Sy * p1t.z;
	p2t.x += Sx * p2t.z;
	p2t.y += Sy * p2t.z;

	Float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
	Float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
	Float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

	if (sizeof(Float) == sizeof(float) &&
		(e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
		double p2txp1ty = (double)p2t.x * (double)p1t.y;
		double p2typ1tx = (double)p2t.y * (double)p1t.x;
		e0 = (float)(p2typ1tx - p2txp1ty);
		double p0txp2ty = (double)p0t.x * (double)p2t.y;
		double p0typ2tx = (double)p0t.y * (double)p2t.x;
		e1 = (float)(p0typ2tx - p0txp2ty);
		double p1txp0ty = (double)p1t.x * (double)p0t.y;
		double p1typ0tx = (double)p1t.y * (double)p0t.x;
		e2 = (float)(p1typ0tx - p1txp0ty);
	}

	if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
		return false;
	Float det = e0 + e1 + e2;
	if (det == 0) return false;

	p0t.z *= Sz;
	p1t.z *= Sz;
	p2t.z *= Sz;
	Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
	if (det < 0 && (tScaled >= 0 || tScaled < ray.tMax * det))
		return false;
	else if (det > 0 && (tScaled <= 0 || tScaled > ray.tMax * det))
		return false;

	Float invDet = 1 / det;
	Float b0 = e0 * invDet;
	Float b1 = e1 * invDet;
	Float b2 = e2 * invDet;
	Float t = tScaled * invDet;

	if (t <= eps) return false;

	return true;
}

Float Triangle::Area() const {
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	return 0.5 * Cross(p1 - p0, p2 - p0).Length();
}

Interaction Triangle::Sample(const Point2f &u, Float *pdf) const {
	Point2f b = UniformSampleTriangle(u);
	const Point3f &p0 = mesh->p[v[0]];
	const Point3f &p1 = mesh->p[v[1]];
	const Point3f &p2 = mesh->p[v[2]];
	Interaction it;
	it.p = b[0] * p0 + b[1] * p1 + (1 - b[0] - b[1]) * p2;
	it.n = Normalize(Normal3f(Cross(p1 - p0, p2 - p0)));
		if (mesh->n) {
		Normal3f ns(b[0] * mesh->n[v[0]] + b[1] * mesh->n[v[1]] +
			(1 - b[0] - b[1]) * mesh->n[v[2]]);
		it.n = Faceforward(it.n, ns);
	}
	else if (reverseOrientation ^ transformSwapsHandedness)
		it.n *= -1;
		Point3f pAbsSum =
		Abs(b[0] * p0) + Abs(b[1] * p1) + Abs((1 - b[0] - b[1]) * p2);
		*pdf = 1 / Area();
	return it;
}

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
	const Transform *ObjectToWorld, const Transform *WorldToObject,
	bool reverseOrientation, int nTriangles, const int *vertexIndices,
	int nVertices, const Point3f *p, const Vector3f *s, const Normal3f *n,
	const Point2f *uv,const int *faceIndices) 
{
	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
		*ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv,faceIndices);
	std::vector<std::shared_ptr<Shape>> tris;
	tris.reserve(nTriangles);
	for (int i = 0; i < nTriangles; ++i)
		tris.push_back(std::make_shared<Triangle>(ObjectToWorld, WorldToObject,
			reverseOrientation, mesh, i));
	return tris;
}