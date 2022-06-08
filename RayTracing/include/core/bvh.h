#pragma once

#ifndef PBRT_ACCELERATORS_BVH_H
#define PBRT_ACCELERATORS_BVH_H

#include "pbrt.h"
#include "primitive.h"
#include "geometry.h"

struct BVHBuildNode;
struct BVHPrimitiveInfo;
struct LinearBVHNode;
struct BucketInfo;

class BVHAccel : public Aggregate {
public:
		enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };
		BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
		int maxPrimsInNode = 1,
		SplitMethod splitMethod = SplitMethod::Middle);
		Bounds3f WorldBound() const;
		~BVHAccel();
		bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
		bool IntersectP(const Ray &ray) const;

private:
		BVHBuildNode *recursiveBuild(
		std::vector<BVHPrimitiveInfo> &primitiveInfo,
		int start, int end, int *totalNodes,
		std::vector<std::shared_ptr<Primitive>> &orderedPrims);
		int flattenBVHTree(BVHBuildNode *node, int *offset);

		const int maxPrimsInNode;
		const SplitMethod splitMethod;
		std::vector<std::shared_ptr<Primitive>> primitives;
		LinearBVHNode *nodes = nullptr;
};













#endif