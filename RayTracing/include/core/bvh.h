#pragma once

#ifndef PBRT_ACCELERATORS_BVH_H
#define PBRT_ACCELERATORS_BVH_H

#include "pbrt.h"
#include "primitive.h"
#include "geometry.h"

struct BVHBuildNode;
// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;
//struct MortonPrimitive;
struct LinearBVHNode;
struct BucketInfo;

class BVHAccel : public Aggregate {
public:
	// BVHAccel Public Types
	enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

	// BVHAccel Public Methods
	BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
		int maxPrimsInNode = 1,
		SplitMethod splitMethod = SplitMethod::Middle);
	Bounds3f WorldBound() const;
	~BVHAccel();
	bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
	bool IntersectP(const Ray &ray) const;

private:
	// BVHAccel Private Methods
	BVHBuildNode *recursiveBuild(
		std::vector<BVHPrimitiveInfo> &primitiveInfo,
		int start, int end, int *totalNodes,
		std::vector<std::shared_ptr<Primitive>> &orderedPrims);
	//BVHBuildNode *HLBVHBuild(
	//	const std::vector<BVHPrimitiveInfo> &primitiveInfo,
	//	int *totalNodes,
	//	std::vector<std::shared_ptr<Primitive>> &orderedPrims) const;
	//BVHBuildNode *emitLBVH(
	//	BVHBuildNode *&buildNodes,
	//	const std::vector<BVHPrimitiveInfo> &primitiveInfo,
	//	MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
	//	std::vector<std::shared_ptr<Primitive>> &orderedPrims,
	//	std::atomic<int> *orderedPrimsOffset, int bitIndex) const;
	//BVHBuildNode *buildUpperSAH(MemoryArena &arena,
	//	std::vector<BVHBuildNode *> &treeletRoots,
	//	int start, int end, int *totalNodes) const;
	int flattenBVHTree(BVHBuildNode *node, int *offset);

	// BVHAccel Private Data
	const int maxPrimsInNode;
	const SplitMethod splitMethod;
	std::vector<std::shared_ptr<Primitive>> primitives;
	LinearBVHNode *nodes = nullptr;
};













#endif