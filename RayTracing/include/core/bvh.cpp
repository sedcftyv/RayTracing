#include "bvh.h"

static long long treeBytes = 0;
static long long totalPrimitives = 0;
static long long totalLeafNodes = 0;
static long long interiorNodes = 0;
static long long leafNodes = 0;

struct BVHPrimitiveInfo {
	BVHPrimitiveInfo() {}
	BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
		: primitiveNumber(primitiveNumber),
		bounds(bounds),
		centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}
	size_t primitiveNumber;
	Bounds3f bounds;
	Point3f centroid;
};

struct BVHBuildNode {
		void InitLeaf(int first, int n, const Bounds3f &b) {
		firstPrimOffset = first;
		nPrimitives = n;
		bounds = b;
		children[0] = children[1] = nullptr;
		++leafNodes;
		++totalLeafNodes;
		totalPrimitives += n;
	}
	void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
		children[0] = c0;
		children[1] = c1;
		bounds = Union(c0->bounds, c1->bounds);
		splitAxis = axis;
		nPrimitives = 0;
		++interiorNodes;
	}
	Bounds3f bounds;
	BVHBuildNode *children[2];
	int splitAxis, firstPrimOffset, nPrimitives;
};

struct LinearBVHNode {
	Bounds3f bounds;
	union {
		int primitivesOffset;   		int secondChildOffset;  	};
	uint16_t nPrimitives;  	uint8_t axis;          	uint8_t pad[1];        };

struct BucketInfo {
	int count = 0;
	Bounds3f bounds;
};

BVHAccel::~BVHAccel() { 
	if (!nodes)
		delete[]nodes;
}

Bounds3f BVHAccel::WorldBound() const {
	return nodes ? nodes[0].bounds : Bounds3f();
}

BVHAccel::BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
	int maxPrimsInNode, SplitMethod splitMethod)
	: maxPrimsInNode(std::min(255, maxPrimsInNode)),
	splitMethod(splitMethod),
	primitives(std::move(p)) {
	if (primitives.empty()) return;	
	std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
	for (size_t i = 0; i < primitives.size(); ++i)
		primitiveInfo[i] = { i, primitives[i]->WorldBound() };
	int totalNodes = 0;
	std::vector<std::shared_ptr<Primitive>> orderedPrims;
	orderedPrims.reserve(primitives.size());
	BVHBuildNode *root;
	root = recursiveBuild(primitiveInfo, 0, primitives.size(),&totalNodes, orderedPrims);
	primitives.swap(orderedPrims);
	primitiveInfo.resize(0);						
	treeBytes += totalNodes * sizeof(LinearBVHNode) + sizeof(*this) +
	primitives.size() * sizeof(primitives[0]);
	nodes = new LinearBVHNode[totalNodes];
	int offset = 0;
	flattenBVHTree(root, &offset);
}

BVHBuildNode *BVHAccel::recursiveBuild(
	std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
	int end, int *totalNodes,
	std::vector<std::shared_ptr<Primitive>> &orderedPrims) {
	BVHBuildNode *node = new BVHBuildNode;
	(*totalNodes)++;
		Bounds3f bounds;
	for (int i = start; i < end; ++i)
		bounds = Union(bounds, primitiveInfo[i].bounds);
	int nPrimitives = end - start;
	if (nPrimitives == 1) {
				int firstPrimOffset = orderedPrims.size();
		for (int i = start; i < end; ++i) {
			int primNum = primitiveInfo[i].primitiveNumber;
			orderedPrims.push_back(primitives[primNum]);
		}
		node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
		return node;
	}
	else {
				Bounds3f centroidBounds;
		for (int i = start; i < end; ++i)
			centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
		int dim = centroidBounds.MaximumExtent();
		int mid = (start + end) / 2;
		if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
						int firstPrimOffset = orderedPrims.size();
			for (int i = start; i < end; ++i) {
				int primNum = primitiveInfo[i].primitiveNumber;
				orderedPrims.push_back(primitives[primNum]);
			}
			node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
			return node;
		}
		else {
			switch (splitMethod) {
			case SplitMethod::Middle: {
								Float pmid =
					(centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
				BVHPrimitiveInfo *midPtr = std::partition(
					&primitiveInfo[start], &primitiveInfo[end - 1] + 1,
					[dim, pmid](const BVHPrimitiveInfo &pi) {
					return pi.centroid[dim] < pmid;
				});
				mid = midPtr - &primitiveInfo[0];
				if (mid != start && mid != end) break;
			}
			case SplitMethod::EqualCounts: {
				mid = (start + end) / 2;
				std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
					&primitiveInfo[end - 1] + 1,
					[dim](const BVHPrimitiveInfo &a,
						const BVHPrimitiveInfo &b) {
					return a.centroid[dim] < b.centroid[dim];
				});
				break;
			}
			case SplitMethod::SAH:
			default: {
					if (nPrimitives <= 2) {
					mid = (start + end) / 2;
					std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
						&primitiveInfo[end - 1] + 1,
						[dim](const BVHPrimitiveInfo &a,
							const BVHPrimitiveInfo &b) {
						return a.centroid[dim] <
							b.centroid[dim];
					});
				}
				else {
					PBRT_CONSTEXPR int nBuckets = 12;
					BucketInfo buckets[nBuckets];
					for (int i = start; i < end; ++i) {
					int b = nBuckets *
						centroidBounds.Offset(
							primitiveInfo[i].centroid)[dim];
					if (b == nBuckets) b = nBuckets - 1;
					buckets[b].count++;
					buckets[b].bounds =
						Union(buckets[b].bounds, primitiveInfo[i].bounds);
					}
					Float cost[nBuckets - 1];
					for (int i = 0; i < nBuckets - 1; ++i) {
						Bounds3f b0, b1;
						int count0 = 0, count1 = 0;
						for (int j = 0; j <= i; ++j) {
							b0 = Union(b0, buckets[j].bounds);
							count0 += buckets[j].count;
						}
						for (int j = i + 1; j < nBuckets; ++j) {
							b1 = Union(b1, buckets[j].bounds);
							count1 += buckets[j].count;
						}
						cost[i] = 1 +
							(count0 * b0.SurfaceArea() +
								count1 * b1.SurfaceArea()) /
							bounds.SurfaceArea();
					}

					Float minCost = cost[0];
					int minCostSplitBucket = 0;
					for (int i = 1; i < nBuckets - 1; ++i) {
						if (cost[i] < minCost) {
							minCost = cost[i];
							minCostSplitBucket = i;
						}
					}

					Float leafCost = nPrimitives;
					if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
						BVHPrimitiveInfo *pmid = std::partition(
							&primitiveInfo[start], &primitiveInfo[end - 1] + 1,
							[=](const BVHPrimitiveInfo &pi) {
							int b = nBuckets *
								centroidBounds.Offset(pi.centroid)[dim];
							if (b == nBuckets) b = nBuckets - 1;
							return b <= minCostSplitBucket;
						});
						mid = pmid - &primitiveInfo[0];
					}
					else {
						int firstPrimOffset = orderedPrims.size();
						for (int i = start; i < end; ++i) {
							int primNum = primitiveInfo[i].primitiveNumber;
							orderedPrims.push_back(primitives[primNum]);
						}
						node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
						return node;
					}
				}
				break;
			}
			}
			node->InitInterior(dim,
				recursiveBuild(primitiveInfo, start, mid,
					totalNodes, orderedPrims),
				recursiveBuild(primitiveInfo, mid, end,
					totalNodes, orderedPrims));
		}
	}
	return node;
}

int BVHAccel::flattenBVHTree(BVHBuildNode *node, int *offset) {
	LinearBVHNode *linearNode = &nodes[*offset];
	linearNode->bounds = node->bounds;
	int myOffset = (*offset)++;
	if (node->nPrimitives > 0) {
		linearNode->primitivesOffset = node->firstPrimOffset;
		linearNode->nPrimitives = node->nPrimitives;
	}
	else {
		linearNode->axis = node->splitAxis;
		linearNode->nPrimitives = 0;
		flattenBVHTree(node->children[0], offset);
		linearNode->secondChildOffset =
			flattenBVHTree(node->children[1], offset);
	}
	return myOffset;
}

bool BVHAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
	if (!nodes) return false;
		bool hit = false;
	Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
		int toVisitOffset = 0, currentNodeIndex = 0;
	int nodesToVisit[64];
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
				if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			if (node->nPrimitives > 0) {
					for (int i = 0; i < node->nPrimitives; ++i)
					if (primitives[node->primitivesOffset + i]->Intersect(
						ray, isect))
						hit = true;
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
			else {
					if (dirIsNeg[node->axis]) {
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				}
				else {
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex = currentNodeIndex + 1;
				}
			}
		}
		else {
			if (toVisitOffset == 0) break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
	return hit;
}

bool BVHAccel::IntersectP(const Ray &ray) const {
	if (!nodes) return false;
		Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	int nodesToVisit[64];
	int toVisitOffset = 0, currentNodeIndex = 0;
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
				if (node->nPrimitives > 0) {
				for (int i = 0; i < node->nPrimitives; ++i) {
					if (primitives[node->primitivesOffset + i]->IntersectP(
						ray)) {
						return true;
					}
				}
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
			else {
				if (dirIsNeg[node->axis]) {
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				}
				else {
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex = currentNodeIndex + 1;
				}
			}
		}
		else {
			if (toVisitOffset == 0) break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
	return false;
}
