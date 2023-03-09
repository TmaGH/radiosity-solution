#pragma once


#include "BvhNode.hpp"

#include <queue>
#include <vector>
#include <iostream>
#include <memory>
#include <stack>


namespace FW {

    struct NodeTraverse {
        const BvhNode &node;
        size_t childCount;
    };

	struct SplitPlane {
		int axis;
		int planeIndex;
	};

    constexpr int partitionsPerAxis = 24;

class Bvh {
public:

    Bvh();
    Bvh(std::istream& is);

    // move assignment for performance
    Bvh& operator=(Bvh&& other) {
        mode_ = other.mode_;
        std::swap(rootNode_, other.rootNode_);
        std::swap(indices_, other.indices_);
        std::swap(preOrderTraverse, other.preOrderTraverse);
        return *this;
    }

    BvhNode&			root() { return *rootNode_; }
    const BvhNode&		root() const { return *rootNode_; }

    void				save(std::ostream& os);

    uint32_t			getIndex(uint32_t index) const { return indices_[index]; };

    void                constructBvh(std::vector<RTTriangle>& triangles, SplitMode splitMode);

	std::vector<NodeTraverse>     preOrderTraverse;

private:

    size_t                          constructBvhNode(BvhNode* node, std::vector<RTTriangle>& triangles, SplitMode splitMode, int depth);

    size_t                          rebuildPreOrderTraverse(BvhNode* node);

    void                            spatialMedianSplit(BvhNode* node, std::vector<RTTriangle>& triangles);
    void                            objectMedianSplit(BvhNode *node, std::vector<RTTriangle>& triangles);
    void                            surfaceAreaHeuristicSplit(BvhNode *node, std::vector<RTTriangle>& triangles);

    bool                            hasOnlyOneGroup(size_t startPrim, size_t endPrim, std::vector<uint32_t>::iterator splitPrimIt);
    void                            buildNodeAABB(BvhNode* node, std::vector<RTTriangle>& triangles);
    void                            createChildNodes(BvhNode* node, size_t splitPrim, Vec3f &planeToSplitBy, std::vector<RTTriangle>& triangles, std::vector<std::shared_ptr<BvhNode>>& results, std::vector<std::pair<char, Vec3f>>& splitPlanes, char splitAxis);

    void                            calculateSahCostsForAxis(BvhNode *node, int axis, std::vector<RTTriangle>& triangles, std::vector<std::pair<float, SplitPlane>> &sahCosts);
    void                            extendBoundingBoxByTriangle(AABB& bb, RTTriangle& triangle);
    void                            extendBoundingBoxByBoundingBox(AABB& bb, AABB& bb2);
    int                             calculateTrianglePartition(RTTriangle& triangle, AABB rootBox, Vec3f rootBoxLength, int axis, std::vector<RTTriangle>& triangles);


    SplitMode						mode_;
    std::unique_ptr<BvhNode>		rootNode_;
    
	std::vector<uint32_t>			indices_; // triangle index list that will be sorted during BVH construction
};


}