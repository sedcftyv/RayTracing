#pragma once

#ifndef PBRT_SHAPES_TRIANGLE_H
#define PBRT_SHAPES_TRIANGLE_H

#include "shape.h"
#include "interaction.h"

static long long triMeshBytes = 0;

struct TriangleMesh {
	// TriangleMesh Public Methods
	TriangleMesh(const Transform &ObjectToWorld, int nTriangles,
		const int *vertexIndices, int nVertices, const Point3f *P,
		const Vector3f *S, const Normal3f *N, const Point2f *uv,
		const int *faceIndices);

	// TriangleMesh Data
	const int nTriangles, nVertices;
	std::vector<int> vertexIndices;
	std::unique_ptr<Point3f[]> p;
	std::unique_ptr<Normal3f[]> n;
	std::unique_ptr<Vector3f[]> s;
	std::unique_ptr<Point2f[]> uv;
	std::vector<int> faceIndices;
};

class Triangle : public Shape {
public:
	// Triangle Public Methods
	Triangle(const Transform *ObjectToWorld, const Transform *WorldToObject,
		bool reverseOrientation, const std::shared_ptr<TriangleMesh> &mesh,
		int triNumber)
		: Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh) {
		v = &mesh->vertexIndices[3 * triNumber];
		triMeshBytes += sizeof(*this);
		faceIndex = mesh->faceIndices.size() ? mesh->faceIndices[triNumber] : 0;
	}
	Bounds3f ObjectBound() const;
	Bounds3f WorldBound() const;
	bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
		bool testAlphaTexture = true) const;
	bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const;
	Float Area() const;

	//using Shape::Sample;  // Bring in the other Sample() overload.
	//Interaction Sample(const Point2f &u, Float *pdf) const;

	// Returns the solid angle subtended by the triangle w.r.t. the given
	// reference point p.
	//Float SolidAngle(const Point3f &p, int nSamples = 0) const;

private:
	// Triangle Private Methods
	void GetUVs(Point2f uv[3]) const {
		if (mesh->uv) {
			uv[0] = mesh->uv[v[0]];
			uv[1] = mesh->uv[v[1]];
			uv[2] = mesh->uv[v[2]];
		}
		else {
			uv[0] = Point2f(0, 0);
			uv[1] = Point2f(1, 0);
			uv[2] = Point2f(1, 1);
		}
	}

	// Triangle Private Data
	std::shared_ptr<TriangleMesh> mesh;
	const int *v;
	int faceIndex;
};

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
	const Transform *o2w, const Transform *w2o, bool reverseOrientation,
	int nTriangles, const int *vertexIndices, int nVertices, const Point3f *p,
	const Vector3f *s, const Normal3f *n, const Point2f *uv,const int *faceIndices = nullptr);











#endif