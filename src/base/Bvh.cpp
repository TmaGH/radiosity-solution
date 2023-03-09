
#include "Bvh.hpp"
#include "filesaves.hpp"
#include <queue>
#include <algorithm>

namespace FW {


	Bvh::Bvh() { }


	// reconstruct from a file
	Bvh::Bvh(std::istream& is) {
		// Load file header.
		fileload(is, mode_);
		Statusbar nodeInfo("Loading nodes", 0);
		Loader loader(is, nodeInfo);

		// Load elements.
		{
			size_t size;
			fileload(is, size);

			for (size_t i = 0; i < size; ++i) {
				uint32_t idx;
				loader(idx);
				indices_.push_back(idx);
			}
		}

		// Load the rest.
		rootNode_.reset(new BvhNode(loader));

		rebuildPreOrderTraverse(rootNode_.get());
	}

	void Bvh::save(std::ostream& os) {
		// Save file header.
		filesave(os, mode_);
		Statusbar nodeInfo("Saving nodes", 0);
		Saver saver(os, nodeInfo);

		// Save elements.
		{
			filesave(os, (size_t)indices_.size());

			for (auto& i : indices_) {
				saver(i);
			}
		}

		// Save the rest.
		rootNode_->save(saver);
	}

	size_t Bvh::rebuildPreOrderTraverse(BvhNode* node) {
		if (!node) return 0;

		size_t nodeTraverseIndex = preOrderTraverse.size();
		size_t childCount = 0;
		preOrderTraverse.push_back({ *node, 0 });

		childCount += rebuildPreOrderTraverse(node->left.get());
		childCount += rebuildPreOrderTraverse(node->right.get());

		preOrderTraverse[nodeTraverseIndex].childCount = childCount;

		return childCount + 1;
	}

	bool Bvh::hasOnlyOneGroup(size_t startPrim, size_t endPrim, std::vector<uint32_t>::iterator splitPrimIt) {
		return std::distance(std::next(indices_.begin(), startPrim), splitPrimIt) == 0 || std::distance(std::next(indices_.begin(), startPrim), splitPrimIt) == endPrim;
	}

	void Bvh::buildNodeAABB(BvhNode* node, std::vector<RTTriangle>& triangles) {
		size_t startPrim = node->startPrim, endPrim = node->endPrim;
		Vec3f max, min;

		max = triangles[indices_[startPrim]].m_vertices[0].p;
		min = triangles[indices_[startPrim]].m_vertices[0].p;

		for (size_t i = startPrim; i < endPrim; i++) {
			auto triangle = triangles[indices_[i]];

			for (VertexPNTC& vertex : triangle.m_vertices) {
				auto vertexPos = vertex.p;
				max = Vec3f(FW::max(vertexPos.x, max.x), FW::max(vertexPos.y, max.y), FW::max(vertexPos.z, max.z));
				min = Vec3f(FW::min(vertexPos.x, min.x), FW::min(vertexPos.y, min.y), FW::min(vertexPos.z, min.z));
			}
		}

		node->bb.min = min;
		node->bb.max = max;
	}

	void Bvh::createChildNodes(BvhNode* node, size_t splitPrim, Vec3f& planeToSplitBy, std::vector<RTTriangle>& triangles, std::vector<std::shared_ptr<BvhNode>>& results, std::vector<std::pair<char, Vec3f>>& splitPlanes, char splitAxis) {
		size_t startPrim = node->startPrim, endPrim = node->endPrim;

		std::shared_ptr<BvhNode> L = std::make_shared<BvhNode>(startPrim, splitPrim);
		std::shared_ptr<BvhNode> R = std::make_shared<BvhNode>(splitPrim, endPrim);
		buildNodeAABB(L.get(), triangles);
		buildNodeAABB(R.get(), triangles);

		results.push_back(L);
		results.push_back(R);
		splitPlanes.push_back(std::make_pair(splitAxis, planeToSplitBy));
	}

	void Bvh::objectMedianSplit(BvhNode* node, std::vector<RTTriangle>& triangles) {
		// Not implemented yet
		return;
	}

	void Bvh::spatialMedianSplit(BvhNode* node, std::vector<RTTriangle>& triangles) {
		size_t startPrim = node->startPrim, endPrim = node->endPrim;
		Vec3f max, min, bbMax, bbMin;

		max = triangles[indices_[startPrim]].centroid();
		min = triangles[indices_[startPrim]].centroid();

		bbMax = triangles[indices_[startPrim]].m_vertices[0].p;
		bbMin = triangles[indices_[startPrim]].m_vertices[0].p;

		for (size_t i = startPrim; i < endPrim; i++) {
			auto triangle = triangles[indices_[i]];
			auto centroid = triangle.centroid();

			max = Vec3f(FW::max(centroid.x, max.x), FW::max(centroid.y, max.y), FW::max(centroid.z, max.z));
			min = Vec3f(FW::min(centroid.x, min.x), FW::min(centroid.y, min.y), FW::min(centroid.z, min.z));

			for (VertexPNTC& vertex : triangle.m_vertices) {
				auto vertexPos = vertex.p;
				bbMax = Vec3f(FW::max(vertexPos.x, bbMax.x), FW::max(vertexPos.y, bbMax.y), FW::max(vertexPos.z, bbMax.z));
				bbMin = Vec3f(FW::min(vertexPos.x, bbMin.x), FW::min(vertexPos.y, bbMin.y), FW::min(vertexPos.z, bbMin.z));
			}
		}

		node->bb.min = bbMin;
		node->bb.max = bbMax;

		if (endPrim - startPrim < 8) return;

		Vec3f center = (min + max) / 2;

		// std::cout << "Max: " << max << ", Min: " << min << ", Center: " << center << std::endl;

		std::vector<uint32_t>::iterator splitPrimIt;

		auto chooseSplit = [&triangles, &center](int index) {return triangles[index].centroid().x < center.x; };
		splitPrimIt = std::partition(std::next(indices_.begin(), startPrim), std::next(indices_.begin(), endPrim), chooseSplit);

		if (hasOnlyOneGroup(startPrim, endPrim, splitPrimIt)) {
			auto chooseSplit = [&triangles, &center](int index) { return triangles[index].centroid().y < center.y; };
			splitPrimIt = std::partition(std::next(indices_.begin(), startPrim), std::next(indices_.begin(), endPrim), chooseSplit);
		}

		if (hasOnlyOneGroup(startPrim, endPrim, splitPrimIt)) {
			auto chooseSplit = [&triangles, &center](int index) { return triangles[index].centroid().z < center.z; };
			splitPrimIt = std::partition(std::next(indices_.begin(), startPrim), std::next(indices_.begin(), endPrim), chooseSplit);
		}

		size_t splitPrim = std::distance(indices_.begin(), splitPrimIt);

		node->left = std::make_unique<BvhNode>(startPrim, splitPrim);
		node->right = std::make_unique<BvhNode>(splitPrim, endPrim);
	}

	void Bvh::extendBoundingBoxByTriangle(AABB& bb, RTTriangle& triangle) {

		if (bb.max == Vec3f(0, 0, 0) && bb.min == Vec3f(0, 0, 0)) {
			bb.max = triangle.m_vertices[0].p;
			bb.min = triangle.m_vertices[0].p;
		}

		for (VertexPNTC& vertex : triangle.m_vertices) {
			auto vertexPos = vertex.p;
			bb.max = Vec3f(FW::max(vertexPos.x, bb.max.x), FW::max(vertexPos.y, bb.max.y), FW::max(vertexPos.z, bb.max.z));
			bb.min = Vec3f(FW::min(vertexPos.x, bb.min.x), FW::min(vertexPos.y, bb.min.y), FW::min(vertexPos.z, bb.min.z));
		}
	}

	void Bvh::extendBoundingBoxByBoundingBox(AABB& bb, AABB& bb2) {
		if (bb.max == Vec3f(0, 0, 0) && bb.min == Vec3f(0, 0, 0)) {
			bb.max = bb2.max;
			bb.min = bb2.min;
		}

		bb.max = Vec3f(FW::max(bb2.max.x, bb.max.x), FW::max(bb2.max.y, bb.max.y), FW::max(bb2.max.z, bb.max.z));
		bb.min = Vec3f(FW::min(bb2.min.x, bb.min.x), FW::min(bb2.min.y, bb.min.y), FW::min(bb2.min.z, bb.min.z));
	}

	int Bvh::calculateTrianglePartition(RTTriangle& triangle, AABB rootBox, Vec3f rootBoxLength, int axis, std::vector<RTTriangle>& triangles) {
		float triangleCentroidAxis = triangle.centroid()[axis];
		float rootBoxMinAxis = rootBox.min[axis];
		float rootBoxLengthAxis = rootBoxLength[axis];

		if (rootBoxLengthAxis <= 0) rootBoxLengthAxis = 1;

		float triangleToBoxOffset = (triangleCentroidAxis - rootBoxMinAxis) / rootBoxLengthAxis;
		int partitionIndex = partitionsPerAxis * triangleToBoxOffset;
		if (partitionIndex == partitionsPerAxis) partitionIndex--;

		return partitionIndex;
	}

	void Bvh::calculateSahCostsForAxis(BvhNode* node, int axis, std::vector<RTTriangle>& triangles, std::vector<std::pair<float, SplitPlane>>& sahCosts) {
		size_t startPrim = node->startPrim, endPrim = node->endPrim;
		AABB rootBox = node->bb;

		Vec3f rootBoxLength = rootBox.max - rootBox.min;

		struct BoxPartitionInfo {
			int triangleCount = 0;
			AABB bb;
		};

		AABB bb;

		BoxPartitionInfo bbPartitions[partitionsPerAxis];

		for (int i = startPrim; i < endPrim; i++) {
			RTTriangle triangle = triangles[indices_[i]];
			int partitionIndex = calculateTrianglePartition(triangle, rootBox, rootBoxLength, axis, triangles);

			bbPartitions[partitionIndex].triangleCount++;

			extendBoundingBoxByTriangle(bbPartitions[partitionIndex].bb, triangle);
		}

		int L = 0, R = partitionsPerAxis - 1;
		BoxPartitionInfo leftPartitions[partitionsPerAxis - 1], rightPartitions[partitionsPerAxis - 1];

		int triangleCountLeft = 0, triangleCountRight = 0;
		AABB bbLeft, bbRight;

		while (L < partitionsPerAxis - 1 && R > 0) {
			triangleCountLeft += bbPartitions[L].triangleCount;
			triangleCountRight += bbPartitions[R].triangleCount;

			if (bbPartitions[L].triangleCount > 0) extendBoundingBoxByBoundingBox(bbLeft, bbPartitions[L].bb);
			if (bbPartitions[R].triangleCount > 0) extendBoundingBoxByBoundingBox(bbRight, bbPartitions[R].bb);

			leftPartitions[L] = { triangleCountLeft, bbLeft };
			rightPartitions[R - 1] = { triangleCountRight, bbRight };

			L++;
			R--;
		}

		for (int i = 0; i < partitionsPerAxis - 1; i++) {
			if (leftPartitions[i].triangleCount == 0 || rightPartitions[i].triangleCount == 0) continue;
			float sah = leftPartitions[i].triangleCount * leftPartitions[i].bb.area() + rightPartitions[i].triangleCount * rightPartitions[i].bb.area();
			sahCosts.push_back(std::make_pair(sah, SplitPlane{ axis, i }));
		}
	}

	void Bvh::surfaceAreaHeuristicSplit(BvhNode* node, std::vector<RTTriangle>& triangles) {
		size_t startPrim = node->startPrim, endPrim = node->endPrim;
		AABB rootBox = node->bb;
		Vec3f rootBoxLength = rootBox.max - rootBox.min;

		if (endPrim - startPrim < 10) return;

		size_t splitPrim;
		SplitPlane planeToSplitBy;

		std::vector<std::pair<float, SplitPlane>> sahCosts;

		calculateSahCostsForAxis(node, 0, triangles, sahCosts);
		calculateSahCostsForAxis(node, 1, triangles, sahCosts);
		calculateSahCostsForAxis(node, 2, triangles, sahCosts);

		int minPairIndex = 0;
		float minSah = INFINITY;

		if (sahCosts.size() == 0) return;

		for (int i = 0; i < sahCosts.size(); i++) {
			if (sahCosts[i].first < minSah) {
				minSah = sahCosts[i].first;
				minPairIndex = i;
			}
		}

		planeToSplitBy = sahCosts[minPairIndex].second;
		int axisToSplitBy = planeToSplitBy.axis;

		auto splitPrimPred = [&triangles, &rootBox, &rootBoxLength, &planeToSplitBy, &axisToSplitBy, this](int index) {
			int partitionIndex = calculateTrianglePartition(triangles[index], rootBox, rootBoxLength, axisToSplitBy, triangles);
			return partitionIndex <= planeToSplitBy.planeIndex;
		};

		splitPrim = std::partition(indices_.begin() + startPrim, indices_.begin() + endPrim, splitPrimPred) - indices_.begin();

		if (splitPrim == startPrim || splitPrim == endPrim) return;

		node->left = std::make_unique<BvhNode>(startPrim, splitPrim);
		buildNodeAABB(node->left.get(), triangles);
		node->right = std::make_unique<BvhNode>(splitPrim, endPrim);
		buildNodeAABB(node->right.get(), triangles);
	}

	size_t Bvh::constructBvhNode(BvhNode* node, std::vector<RTTriangle>& triangles, SplitMode splitMode, int depth) {
		if (!node) return 0;

		// std::cout << "Node: " << node->startPrim << " -> " << node->endPrim << ", depth: " << depth << ", size: " << node->endPrim - node->startPrim << std::endl;
		size_t nodeTraverseIndex = preOrderTraverse.size();
		size_t childCount = 0;
		preOrderTraverse.push_back({ *node, 0 });

		size_t startPrim = node->startPrim, endPrim = node->endPrim;
		size_t splitPrim = endPrim;

		switch (splitMode) {
		case SplitMode_SpatialMedian:
			spatialMedianSplit(node, triangles);
			break;
		case SplitMode_ObjectMedian:
			objectMedianSplit(node, triangles);
			break;
		case SplitMode_Sah:
			surfaceAreaHeuristicSplit(node, triangles);
			break;
		case SplitMode_Linear:
			break;
		case SplitMode_None:
			break;
		default:
			break;
		}

		childCount += constructBvhNode(node->left.get(), triangles, splitMode, depth + 1);
		childCount += constructBvhNode(node->right.get(), triangles, splitMode, depth + 1);

		preOrderTraverse[nodeTraverseIndex].childCount = childCount;

		return childCount + 1;
	}

	void Bvh::constructBvh(std::vector<RTTriangle>& triangles, SplitMode splitMode) {
		for (int i = 0; i < triangles.size(); i++) {
			indices_.push_back(i);
		}

		rootNode_ = std::make_unique<BvhNode>(0, indices_.size());

		buildNodeAABB(rootNode_.get(), triangles);

		size_t treeSize = constructBvhNode(rootNode_.get(), triangles, splitMode, 1);
	}

}
