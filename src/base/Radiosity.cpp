#include "Radiosity.hpp"
#include "AreaLight.hpp"
#include "RayTracer.hpp"
#include <cassert>


namespace FW {


	// --------------------------------------------------------------------------

	Radiosity::~Radiosity()
	{
		if (isRunning())
		{
			m_context.m_bForceExit = true;
			while (m_launcher.getNumTasks() > m_launcher.getNumFinished())
				Sleep(1);
			m_launcher.popAll();
		}
	}


	// --------------------------------------------------------------------------
	void Radiosity::vertexTaskFunc(MulticoreLauncher::Task& task)
	{
		RadiosityContext& ctx = *(RadiosityContext*)task.data;

		if (ctx.m_bForceExit)
			return;

		Random rnd;

		int v = task.idx;

		Vec3f n = ctx.m_scene->vertex(v).n.normalized();
		Vec3f o = ctx.m_scene->vertex(v).p + 0.01f * n;

		int samplingIndex = 0;

		// Simplified spherical harmonics constants from the provided slides
		// Comes from multiplying C4 * Y0 and 2 * C2 * Y1 based on the provided paper
		// Means below I should be calculating L coefficients multiplied C4/C2 constants to get irradiance coefficients 
		float AY0 = 0.25f;
		float AY1 = 0.50f;
		
		if (ctx.m_currentBounce == 0)
		{
			Vec3f E(0);
			Vec3f sphericalC(0);
			Vec3f sphericalX(0);
			Vec3f sphericalY(0);
			Vec3f sphericalZ(0);
			for (int r = 0; r < ctx.m_numDirectRays; ++r)
			{
				float pdf;
				Vec3f Pl;

				ctx.m_light->sample(pdf, Pl, 2, 3, samplingIndex, rnd);
				samplingIndex++;
				Vec3f d = Pl - o;

				RaycastResult shadowResult = ctx.m_rt->raycast(o, d);

				if (!shadowResult.tri) {
					float cosIncDir = n.dot(d.normalized());
					float cosYXDir = ctx.m_light->getNormal().normalized().dot(-d.normalized());

					cosIncDir = max(cosIncDir, 0.0f);
					cosYXDir = max(cosYXDir, 0.0f);


					Vec3f Ei = (ctx.m_light->getEmission() * cosYXDir * cosIncDir) / pow(d.length(), 2) * (1 / pdf);

					// I know the spherical harmonics here don't give the correct result.
					// It's already clear from the light bleeding into the vertices that should be in shadow because of the direct illumination calculations.
					// However, I spent enough time on it to get into diminishing returns with learning so I gave up.
					if (ctx.m_useSpherical) {
						Ei = (ctx.m_light->getEmission() * cosYXDir) / pow(d.length(), 2) * (1 / pdf);
					}

					E += Ei;

					if (ctx.m_useSpherical) {
						sphericalC += Ei;
						sphericalX += Ei * d.normalized().x;
						sphericalY += Ei * d.normalized().y;
						sphericalZ += Ei * d.normalized().z;
					}
				}
			}

			// According to the reference implementation, should be possible to calculate this at the same time as normal irradiance
			// but since I wasn't able to make it work correctly I didn't bother trying to implement it that way.
			if (ctx.m_useSpherical) {
				ctx.m_vecSphericalC[v] = sphericalC * (AY0 / ctx.m_numDirectRays);
				ctx.m_vecSphericalX[v] = sphericalX * (AY1 / ctx.m_numDirectRays);
				ctx.m_vecSphericalY[v] = sphericalY * (AY1 / ctx.m_numDirectRays);
				ctx.m_vecSphericalZ[v] = sphericalZ * (AY1 / ctx.m_numDirectRays);
			}

			ctx.m_vecCurr[v] = E * (1.0f / ctx.m_numDirectRays);
			ctx.m_vecResult[v] = ctx.m_vecCurr[v];
			ctx.m_vecIntermediateResult[ctx.m_currentBounce][v] = ctx.m_vecCurr[v];
		}
		else
		{
			Mat3f B = formBasis(n);

			Vec3f E(0.0f);
			int baseX = 2;
			int baseY = 3;
			int samplingIndex = 0;

			Vec3f sphericalC(0);
			Vec3f sphericalX(0);
			Vec3f sphericalY(0);
			Vec3f sphericalZ(0);

			for (int r = 0; r < ctx.m_numHemisphereRays; ++r)
			{
				// Inefficient QMC based on Halton sequence (copy-paste from area light)
				// Using different base for x and y to avoid correlation
				float x = 0;
				float y = 0;
				int i = samplingIndex;
				float fX = 1.0 / baseX;
				float fY = 1.0 / baseY;

				while (i > 0) {
					x += fX * (i % baseX);
					i /= baseX;
					fX /= baseX;
				}
				i = samplingIndex;
				while (i > 0) {
					y += fY * (i % baseY);
					i /= baseY;
					fY /= baseY;
				}
				samplingIndex++;

				// Above QMC should produce (x, y) in [0,1) unit square
				// so I use the concentric mapping from Shirley & Chiu 97 to map it to a disk.
				float phi, rr, u, v;

				float a = 2 * x - 1;
				float b = 2 * y - 1;

				if (a > -b) {
					if (a > b) {
						rr = a;
						phi = (FW_PI / 4) * (b / a);
					}
					else {
						rr = b;
						phi = (FW_PI / 4) * (2 - (a / b));
					}
				}
				else {
					if (a < b) {
						rr = -a;
						phi = (FW_PI / 4) * (4 + (b / a));
					}
					else {
						rr = -b;
						if (b != 0) {
							phi = (FW_PI / 4) * (6 - (a / b));
						}
						else {
							phi = 0;
						}
					}
				}

				u = rr * cos(phi);
				v = rr * sin(phi);

				Vec3f d;

				// If using spherical harmonics, distribute uniformly on the hemisphere
				// Otherwise, distribute by cosine distribution
				// Again from Shirley & Chiu 97
				if (ctx.m_useSpherical) {
					float z = 1 - rr * rr;
					float c = sqrt(1 - z * z) / rr;
					d = Vec3f(u * c, v * c, z);
				}
				else {
					d = Vec3f(u, v, sqrt(1 - rr * rr));
				}

				d = 100.0f * (B * d);

				const RaycastResult result = ctx.m_rt->raycast(o, d);
				if (result.tri != nullptr)
				{

					const Vec3i& indices = result.tri->m_data.vertex_indices;

					if ((result.tri->normal().dot(d)) > 0.0f)
						continue;

					float alpha = result.u;
					float beta = result.v;
					float gamma = 1.0f - alpha - beta;

					Vec3f Ei = gamma * ctx.m_vecPrevBounce[indices[0]] + alpha * ctx.m_vecPrevBounce[indices[1]] + beta * ctx.m_vecPrevBounce[indices[2]];

					Ei *= (1.0f / FW_PI);

					const auto mat = result.tri->m_material;
					if (mat->textures[MeshBase::TextureType_Diffuse].exists())
					{
						const Texture& tex = mat->textures[MeshBase::TextureType_Diffuse];
						const Image& teximg = *tex.getImage();
						Vec2f uv = gamma * result.tri->m_vertices[0].t + alpha * result.tri->m_vertices[1].t + beta * result.tri->m_vertices[2].t;
						Vec2i texeCoords = getTexelCoords(uv, teximg.getSize());

						Vec3f diffuse = teximg.getVec4f(texeCoords).getXYZ();
						Ei *= diffuse;
					}
					else
					{
						Ei *= mat->diffuse.getXYZ();
					}

					if (ctx.m_useSpherical) {
						sphericalC += Ei;
						sphericalX += Ei * d.normalized().x;
						sphericalY += Ei * d.normalized().y;
						sphericalZ += Ei * d.normalized().z;
					}

					E += Ei;
				}
			}
			ctx.m_vecCurr[v] = E * (FW_PI / ctx.m_numHemisphereRays);

			if (ctx.m_useSpherical) {
				ctx.m_vecSphericalC[v] += sphericalC * (AY0 * FW_PI / ctx.m_numHemisphereRays);
				ctx.m_vecSphericalX[v] += sphericalX * (AY1 * FW_PI / ctx.m_numHemisphereRays);
				ctx.m_vecSphericalY[v] += sphericalY * (AY1 * FW_PI / ctx.m_numHemisphereRays);
				ctx.m_vecSphericalZ[v] += sphericalZ * (AY1 * FW_PI / ctx.m_numHemisphereRays);
			}

			ctx.m_vecResult[v] = ctx.m_vecResult[v] + ctx.m_vecCurr[v];
			ctx.m_vecIntermediateResult[ctx.m_currentBounce][v] = ctx.m_vecCurr[v];

			if (ctx.m_visualizeLastBounce) {
				ctx.m_vecResult[v] = ctx.m_vecCurr[v];
			}
		}
		return;
	}
	// --------------------------------------------------------------------------

	void Radiosity::startRadiosityProcess(MeshWithColors* scene, AreaLight* light, RayTracer* rt, int numBounces, int numDirectRays, int numHemisphereRays, bool useSpherical, int bounceToVisualize, bool visualizeLastBounce)
	{
		// put stuff the asyncronous processor needs 
		m_context.m_scene = scene;
		m_context.m_rt = rt;
		m_context.m_light = light;
		m_context.m_currentBounce = 0;
		m_context.m_numBounces = max(numBounces, bounceToVisualize - 1);
		m_context.m_numDirectRays = numDirectRays;
		m_context.m_numHemisphereRays = numHemisphereRays;
		m_context.m_visualizeLastBounce = visualizeLastBounce;
		m_context.m_useSpherical = useSpherical;

		// resize all the buffers according to how many vertices we have in the scene
		m_context.m_vecResult.resize(scene->numVertices());
		m_context.m_vecIntermediateResult.resize(m_context.m_numBounces + 1);
		for (int i = 0; i <= m_context.m_numBounces; i++) {
			m_context.m_vecIntermediateResult[i].resize(scene->numVertices());
		}
		m_context.m_vecCurr.resize(scene->numVertices());
		m_context.m_vecPrevBounce.resize(scene->numVertices());
		m_context.m_vecResult.assign(scene->numVertices(), Vec3f(0, 0, 0));

		m_context.m_vecSphericalC.resize(scene->numVertices());
		m_context.m_vecSphericalX.resize(scene->numVertices());
		m_context.m_vecSphericalY.resize(scene->numVertices());
		m_context.m_vecSphericalZ.resize(scene->numVertices());

		m_context.m_vecSphericalC.assign(scene->numVertices(), Vec3f(0, 0, 0));
		m_context.m_vecSphericalX.assign(scene->numVertices(), Vec3f(0, 0, 0));
		m_context.m_vecSphericalY.assign(scene->numVertices(), Vec3f(0, 0, 0));
		m_context.m_vecSphericalZ.assign(scene->numVertices(), Vec3f(0, 0, 0));

		// fire away!
		m_launcher.setNumThreads(m_launcher.getNumCores());	// the solution exe is multithreaded
		m_launcher.popAll();
		m_launcher.push(vertexTaskFunc, &m_context, 0, scene->numVertices());
	}
	// --------------------------------------------------------------------------

	bool Radiosity::updateMeshColors(std::vector<Vec4f>& spherical1, std::vector<Vec4f>& spherical2, std::vector<float>& spherical3, bool spherical, int bounceToVisualize)
	{
		if (!m_context.m_scene || m_context.m_vecResult.size() == 0) return false;
		// Print progress.
		printf("%.2f%% done     \r", 100.0f * m_launcher.getNumFinished() / m_context.m_scene->numVertices());

		// Copy irradiance over to the display mesh.
		// Because we want outgoing radiosity in the end, we divide by PI here
		// and let the shader multiply the final diffuse reflectance in. See App::setupShaders() for details.
		for (int i = 0; i < m_context.m_scene->numVertices(); ++i) {

			Vec3f vecResult = m_context.m_vecResult[i];
			Vec3f vecSphericalC = m_context.m_vecSphericalC[i];
			Vec3f vecSphericalX = m_context.m_vecSphericalX[i];
			Vec3f vecSphericalY = m_context.m_vecSphericalY[i];
			Vec3f vecSphericalZ = m_context.m_vecSphericalZ[i];

			if (bounceToVisualize > 0) {
				if (bounceToVisualize > m_context.m_numBounces + 1) break;
				vecResult = m_context.m_vecIntermediateResult[bounceToVisualize - 1][i];
			}

			// Packing data for the spherical harmonic extra.
			// In order to manage with fewer vertex attributes in the shader, the third component is stored as the w components of other actually three-dimensional vectors.
			if (spherical) {
				m_context.m_scene->mutableVertex(i).c = vecSphericalC * (1.0f / FW_PI);
				spherical3[i] = vecSphericalZ.x * (1.0f / FW_PI);
				spherical1[i] = Vec4f(vecSphericalX, vecSphericalZ.y) * (1.0f / FW_PI);
				spherical2[i] = Vec4f(vecSphericalY, vecSphericalZ.z) * (1.0f / FW_PI);
			}
			else {
				m_context.m_scene->mutableVertex(i).c = vecResult * (1.0f / FW_PI);
			}
		}
		return true;
	}
	// --------------------------------------------------------------------------

	void Radiosity::checkFinish()
	{
		// have all the vertices from current bounce finished computing?
		if (m_launcher.getNumTasks() == m_launcher.getNumFinished())
		{
			// yes, remove from task list
			m_launcher.popAll();

			// more bounces desired?
			if (m_context.m_currentBounce < m_context.m_numBounces)
			{
				// move current bounce to prev
				m_context.m_vecPrevBounce = m_context.m_vecCurr;
				++m_context.m_currentBounce;
				// start new tasks for all vertices
				m_launcher.push(vertexTaskFunc, &m_context, 0, m_context.m_scene->numVertices());
				printf("\nStarting bounce %d\n", m_context.m_currentBounce);
			}
			else printf("\n DONE!\n");
		}
	}
	// --------------------------------------------------------------------------

} // namespace FW
