#define _CRT_SECURE_NO_WARNINGS

#include "base/Defs.hpp"
#include "base/Math.hpp"
#include "RayTracer.hpp"
#include <stdio.h>
#include "rtIntersect.inl"
#include <fstream>

#include "rtlib.hpp"


// Helper function for hashing scene data for caching BVHs
extern "C" void MD5Buffer(void* buffer, size_t bufLen, unsigned int* pDigest);


namespace FW
{


	Vec2f getTexelCoords(Vec2f uv, const Vec2i size)
	{
		return Vec2f(fmod(1.0f + fmod(uv.x, 1.0f), 1.0f) * size.x, fmod(1.0f + fmod(uv.y, 1.0f), 1.0f) * size.y);
	}

	Mat3f formBasis(const Vec3f& n) {
		Mat3f R;

		Vec3f Q = Vec3f(n);
		float smallestAbs = min(abs(Q.x), abs(Q.y), abs(Q.z));
		if (smallestAbs == abs(Q.x)) Q.x = 1;
		else if (smallestAbs == abs(Q.y)) Q.y = 1;
		else if (smallestAbs == abs(Q.z)) Q.z = 1;

		Vec3f T = n.cross(Q).normalized();
		Vec3f B = n.cross(T);

		R.setCol(0, T);
		R.setCol(1, B);
		R.setCol(2, n);

		return R;
	}


	String RayTracer::computeMD5(const std::vector<Vec3f>& vertices)
	{
		unsigned char digest[16];
		MD5Buffer((void*)&vertices[0], sizeof(Vec3f) * vertices.size(), (unsigned int*)digest);

		// turn into string
		char ad[33];
		for (int i = 0; i < 16; ++i)
			::sprintf(ad + i * 2, "%02x", digest[i]);
		ad[32] = 0;

		return FW::String(ad);
	}


	// --------------------------------------------------------------------------


	RayTracer::RayTracer()
	{
	}

	RayTracer::~RayTracer()
	{
	}


	void RayTracer::loadHierarchy(const char* filename, std::vector<RTTriangle>& triangles)
	{
		std::ifstream ifs(filename, std::ios::binary);
		m_bvh = Bvh(ifs);

		m_triangles = &triangles;
	}

	void RayTracer::saveHierarchy(const char* filename, const std::vector<RTTriangle>& triangles) {
		(void)triangles; // Not used.

		std::ofstream ofs(filename, std::ios::binary);
		m_bvh.save(ofs);
	}

	void RayTracer::constructHierarchy(std::vector<RTTriangle>& triangles, SplitMode splitMode) {
		m_bvh.constructBvh(triangles, splitMode);
		m_triangles = &triangles;
	}

	RaycastResult RayTracer::raycast(const Vec3f& orig, const Vec3f& dir) const {
		++m_rayCount;
		RaycastResult result;
		Vec3f dirInverse = Vec3f(1.0f, 1.0f, 1.0f) / dir;

		const std::vector<NodeTraverse>& preOrderTraverse = m_bvh.preOrderTraverse;
		size_t i = 0;

		float tRayMax = 1;

		while (i < preOrderTraverse.size()) {
			NodeTraverse nodeTraverse = preOrderTraverse[i];
			const BvhNode& node = nodeTraverse.node;

			Vec3f tMin = (node.bb.min - orig) * dirInverse;
			Vec3f tMax = (node.bb.max - orig) * dirInverse;

			if (tMin.x > tMax.x) std::swap(tMin.x, tMax.x);
			if (tMin.y > tMax.y) std::swap(tMin.y, tMax.y);
			if (tMin.z > tMax.z) std::swap(tMin.z, tMax.z);

			float tStart = FW::max(tMin.x, tMin.y, tMin.z);
			float tEnd = FW::min(tMax.x, tMax.y, tMax.z);

			if (tStart > 0) {
				if (tStart > tRayMax || tEnd < tStart) {
					i++;
					i += nodeTraverse.childCount;
					continue;
				}
			}
			else {
				if (tEnd < 0) {
					i++;
					i += nodeTraverse.childCount;
					continue;
				}
			}

			if (!node.left && !node.right) {
				float closest_t = 1.0f, closest_u = 0.0f, closest_v = 0.0f;
				int closest_i = -1;

				for (int i = node.startPrim; i < node.endPrim; ++i)
				{
					float t, u, v;
					if ((*m_triangles)[m_bvh.getIndex(i)].intersect_woop(orig, dir, t, u, v))
					{
						if (t > 0.0f && t < closest_t)
						{
							closest_i = i;
							closest_t = t;
							closest_u = u;
							closest_v = v;
						}
					}
				}

				if (closest_i != -1) {
					if (!result.tri) {
						result = RaycastResult(&(*m_triangles)[m_bvh.getIndex(closest_i)], closest_t, closest_u, closest_v, orig + closest_t * dir, orig, dir);
					}
					else {
						RaycastResult newResult = RaycastResult(&(*m_triangles)[m_bvh.getIndex(closest_i)], closest_t, closest_u, closest_v, orig + closest_t * dir, orig, dir);
						if (newResult.t < result.t) result = newResult;
					}
					tRayMax = result.t;
				}
			}
			i++;
		}

		return result;
	}

} // namespace FW