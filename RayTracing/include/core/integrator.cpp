#include <fstream>
#include <iostream>
#include <thread>
#include <mutex>
#include <atomic>
#include "rtweekend.h"
#include "integrator.h"
#include "interaction.h"
#include "reflection.h"
#include "light.h"
#include "core/mipmap.h"
#include "core/spectrum.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using std::cout;
using std::endl;
using std::thread;
using std::mutex;
static long long totalPaths = 0;
static long long zeroRadiancePaths = 0;

Integrator::~Integrator() {}

Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene, Sampler &sampler,
	const std::vector<int> &nLightSamples,
	bool handleMedia) {
	//ProfilePhase p(Prof::DirectLighting);
	Spectrum L(0.f);
	for (size_t j = 0; j < scene.lights.size(); ++j) {
		// Accumulate contribution of _j_th light to _L_
		const std::shared_ptr<Light> &light = scene.lights[j];
		int nSamples = nLightSamples[j];
		const Point2f *uLightArray = sampler.Get2DArray(nSamples);
		const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);
		if (!uLightArray || !uScatteringArray) {
			// Use a single sample for illumination from _light_
			Point2f uLight = sampler.Get2D();
			Point2f uScattering = sampler.Get2D();
			L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
				 handleMedia);
		}
		else {
			// Estimate direct lighting using sample arrays
			Spectrum Ld(0.f);
			for (int k = 0; k < nSamples; ++k)
				Ld += EstimateDirect(it, uScatteringArray[k], *light,
					uLightArray[k], scene, sampler,
					handleMedia);
			L += Ld / nSamples;
		}
	}
	return L;
}

Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene, Sampler &sampler,
	bool handleMedia, const Distribution1D *lightDistrib) {
	//ProfilePhase p(Prof::DirectLighting);
	// Randomly choose a single light to sample, _light_
	int nLights = int(scene.lights.size());
	if (nLights == 0) return Spectrum(0.f);
	int lightNum;
	Float lightPdf;
	if (lightDistrib) {
		lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
		if (lightPdf == 0) return Spectrum(0.f);
	}
	else {
		lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
		lightPdf = Float(1) / nLights;
	}
	const std::shared_ptr<Light> &light = scene.lights[lightNum];
	Point2f uLight = sampler.Get2D();
	Point2f uScattering = sampler.Get2D();
	return EstimateDirect(it, uScattering, *light, uLight,
		scene, sampler,  handleMedia) / lightPdf;
}

Spectrum EstimateDirect(const Interaction &it, const Point2f &uScattering,
	const Light &light, const Point2f &uLight,
	const Scene &scene, Sampler &sampler, bool handleMedia, bool specular) {
	BxDFType bsdfFlags =
		specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
	Spectrum Ld(0.f);
	// Sample light source with multiple importance sampling
	Vector3f wi;
	Float lightPdf = 0, scatteringPdf = 0;
	VisibilityTester visibility;
	Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
	//VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li << ", wi: "
	//	<< wi << ", pdf: " << lightPdf;
	if (lightPdf > 0 && !Li.IsBlack()) {
		// Compute BSDF or phase function's value for light sample
		Spectrum f;
		if (it.IsSurfaceInteraction()) {
			// Evaluate BSDF for light sampling strategy
			const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
			f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *
				AbsDot(wi, isect.shading.n);
			//scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
			//VLOG(2) << "  surf f*dot :" << f << ", scatteringPdf: " << scatteringPdf;
		}
		//else {
		//	// Evaluate phase function for light sampling strategy
		//	const MediumInteraction &mi = (const MediumInteraction &)it;
		//	Float p = mi.phase->p(mi.wo, wi);
		//	f = Spectrum(p);
		//	scatteringPdf = p;
		//	VLOG(2) << "  medium p: " << p;
		//}
		if (!f.IsBlack()) {
			// Compute effect of visibility for light source sample
			if (handleMedia) {
				Li *= visibility.Tr(scene, sampler);
				//VLOG(2) << "  after Tr, Li: " << Li;
			}
			else {
				if (!visibility.Unoccluded(scene)) {
					//VLOG(2) << "  shadow ray blocked";
					Li = Spectrum(0.f);
				}
				//else
					//VLOG(2) << "  shadow ray unoccluded";
			}

			// Add light's contribution to reflected radiance
			if (!Li.IsBlack()) {
				//if (IsDeltaLight(light.flags))
					Ld += f * Li / lightPdf;
					//cout << f << ' ' << Li << ' ' << lightPdf << endl;
				//else {
				//	Float weight =
				//		PowerHeuristic(1, lightPdf, 1, scatteringPdf);
				//	Ld += f * Li * weight / lightPdf;
				//	}
			}
		}
	}

	//// Sample BSDF with multiple importance sampling
	//if (!IsDeltaLight(light.flags)) {
	//	Spectrum f;
	//	bool sampledSpecular = false;
	//	if (it.IsSurfaceInteraction()) {
	//		// Sample scattered direction for surface interactions
	//		BxDFType sampledType;
	//		const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
	//		f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
	//			bsdfFlags, &sampledType);
	//		f *= AbsDot(wi, isect.shading.n);
	//		sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
	//	}
	//	//else {
	//	//	// Sample scattered direction for medium interactions
	//	//	const MediumInteraction &mi = (const MediumInteraction &)it;
	//	//	Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
	//	//	f = Spectrum(p);
	//	//	scatteringPdf = p;
	//	//}
	//	//VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " << scatteringPdf;
	//	if (!f.IsBlack() && scatteringPdf > 0) {
	//		// Account for light contributions along sampled direction _wi_
	//		Float weight = 1;
	//		if (!sampledSpecular) {
	//			lightPdf = light.Pdf_Li(it, wi);
	//			if (lightPdf == 0) return Ld;
	//			weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
	//		}

	//		// Find intersection and compute transmittance
	//		SurfaceInteraction lightIsect;
	//		Ray ray = it.SpawnRay(wi);
	//		Spectrum Tr(1.f);
	//		bool foundSurfaceInteraction =
	//			//handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr): 
	//			scene.Intersect(ray, &lightIsect);

	//		// Add light contribution from material sampling
	//		Spectrum Li(0.f);
	//		if (foundSurfaceInteraction) {
	//			if (lightIsect.primitive->GetAreaLight() == &light)
	//				Li = lightIsect.Le(-wi);
	//		}
	//		else
	//			Li = light.Le(ray);
	//		if (!Li.IsBlack()) Ld += f * Li * Tr * weight / scatteringPdf;
	//	}
	//}
	return Ld;
}


Spectrum SamplerIntegrator::SpecularReflect(
	const Ray &ray, const SurfaceInteraction &isect,
	const Scene &scene, Sampler &sampler, int depth) const {
	// Compute specular reflection direction _wi_ and BSDF value
	Vector3f wo = isect.wo, wi;
	Float pdf;
	BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

	// Return contribution of specular reflection
	const Normal3f &ns = isect.shading.n;
	if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
		// Compute ray differential _rd_ for specular reflection
		Ray rd = isect.SpawnRay(wi);
		//if (ray.hasDifferentials) {
		//	rd.hasDifferentials = true;
		//	rd.rxOrigin = isect.p + isect.dpdx;
		//	rd.ryOrigin = isect.p + isect.dpdy;
		//	// Compute differential reflected directions
		//	Normal3f dndx = isect.shading.dndu * isect.dudx +
		//		isect.shading.dndv * isect.dvdx;
		//	Normal3f dndy = isect.shading.dndu * isect.dudy +
		//		isect.shading.dndv * isect.dvdy;
		//	Vector3f dwodx = -ray.rxDirection - wo,
		//		dwody = -ray.ryDirection - wo;
		//	Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
		//	Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);
		//	rd.rxDirection =
		//		wi - dwodx + 2.f * Vector3f(Dot(wo, ns) * dndx + dDNdx * ns);
		//	rd.ryDirection =
		//		wi - dwody + 2.f * Vector3f(Dot(wo, ns) * dndy + dDNdy * ns);
		//}
		return f * Li(rd, scene, sampler, depth + 1) * AbsDot(wi, ns) /
			pdf;
	}
	else
		return Spectrum(0.f);
}

Spectrum SamplerIntegrator::SpecularTransmit(
	const Ray &ray, const SurfaceInteraction &isect,
	const Scene &scene, Sampler &sampler, int depth) const {
	Vector3f wo = isect.wo, wi;
	Float pdf;
	const Point3f &p = isect.p;
	const BSDF &bsdf = *isect.bsdf;
	Spectrum f = bsdf.Sample_f(wo, &wi, sampler.Get2D(), &pdf,
		BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
	Spectrum L = Spectrum(0.f);
	Normal3f ns = isect.shading.n;
	if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
		// Compute ray differential _rd_ for specular transmission
		Ray rd = isect.SpawnRay(wi);
		//if (ray.hasDifferentials) {
		//	rd.hasDifferentials = true;
		//	rd.rxOrigin = p + isect.dpdx;
		//	rd.ryOrigin = p + isect.dpdy;

		//	Normal3f dndx = isect.shading.dndu * isect.dudx +
		//		isect.shading.dndv * isect.dvdx;
		//	Normal3f dndy = isect.shading.dndu * isect.dudy +
		//		isect.shading.dndv * isect.dvdy;

		//	// The BSDF stores the IOR of the interior of the object being
		//	// intersected.  Compute the relative IOR by first out by
		//	// assuming that the ray is entering the object.
		//	Float eta = 1 / bsdf.eta;
		//	if (Dot(wo, ns) < 0) {
		//		// If the ray isn't entering, then we need to invert the
		//		// relative IOR and negate the normal and its derivatives.
		//		eta = 1 / eta;
		//		ns = -ns;
		//		dndx = -dndx;
		//		dndy = -dndy;
		//	}
		//	Vector3f dwodx = -ray.rxDirection - wo,
		//		dwody = -ray.ryDirection - wo;
		//	Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
		//	Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

		//	Float mu = eta * Dot(wo, ns) - AbsDot(wi, ns);
		//	Float dmudx =
		//		(eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNdx;
		//	Float dmudy =
		//		(eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNdy;

		//	rd.rxDirection =
		//		wi - eta * dwodx + Vector3f(mu * dndx + dmudx * ns);
		//	rd.ryDirection =
		//		wi - eta * dwody + Vector3f(mu * dndy + dmudy * ns);
		//}
		L = f * Li(rd, scene, sampler, depth + 1) * AbsDot(wi, ns) / pdf;
	}
	return L;
}

Spectrum WhittedIntegrator::Li(const Ray &ray, const Scene &scene,
	Sampler &sampler,int depth) const {
	//cout << depth << endl;
	Spectrum L(0.);
	// Find closest ray intersection or return background radiance
	SurfaceInteraction isect;
	if (!scene.Intersect(ray, &isect)) {
		for (const auto &light : scene.lights) L += light->Le(ray);
		//if (depth == 0)
		//{
		//	Spectrum LL;
		//	LL[2] = 0.8;
		//	return LL;
		//}
		return L;
	}
	//cout << 1 << endl;
	// Compute emitted and reflected light at ray intersection point

	// Initialize common variables for Whitted integrator
	const Normal3f &n = isect.shading.n;
	Vector3f wo = isect.wo;

	// Compute scattering functions for surface interaction
	isect.ComputeScatteringFunctions(ray);
	if (!isect.bsdf)
		return Li(isect.SpawnRay(ray.d), scene, sampler, depth);

	// Compute emitted light if ray hit an area light source
	L += isect.Le(wo);

	// Add contribution of each light source
	for (const auto &light : scene.lights) {
		Vector3f wi;
		Float pdf;
		VisibilityTester visibility;
		Spectrum Li =
			light->Sample_Li(isect, sampler.Get2D(), &wi, &pdf, &visibility);
		if (Li.IsBlack() || pdf == 0) continue;
		Spectrum f = isect.bsdf->f(wo, wi);
		if (!f.IsBlack() && visibility.Unoccluded(scene))
			L += f * Li * AbsDot(wi, n) / pdf;
	}
	if (depth + 1 < maxDepth) {
		// Trace rays for specular reflection and refraction
		L += SpecularReflect(ray, isect, scene, sampler,  depth);
		L += SpecularTransmit(ray, isect, scene, sampler,  depth);
	}
	return L;
}

void PathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
	lightDistribution =
		CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum PathIntegrator::Li(const Ray &r, const Scene &scene,
	Sampler &sampler,int depth) const {
	//ProfilePhase p(Prof::SamplerIntegratorLi);
	Spectrum L(0.f), beta(1.f);
	Ray ray(r);
	bool specularBounce = false;
	int bounces;

	Float etaScale = 1;

	for (bounces = 0;; ++bounces) {
		//cout << L<<' '<<bounces << endl;
		// Find next path vertex and accumulate contribution
		//VLOG(2) << "Path tracer bounce " << bounces << ", current L = " << L
		//	<< ", beta = " << beta;

		// Intersect _ray_ with scene and store intersection in _isect_
		SurfaceInteraction isect;
		//cout << 's' << endl;
		bool foundIntersection = scene.Intersect(ray, &isect);

		// Possibly add emitted light at intersection
		if (bounces == 0 || specularBounce) {
			// Add emitted light at path vertex or from the environment
			if (foundIntersection) {
				L += beta * isect.Le(-ray.d);
				//VLOG(2) << "Added Le -> L = " << L;
			}
			else {
				for (const auto &light : scene.infiniteLights)
					L += beta * light->Le(ray);
				//VLOG(2) << "Added infinite area lights -> L = " << L;
			}
		}

		//if (L.HasNaNs())
		//{
		//	cout << L << endl;
		//}

		// Terminate path if ray escaped or _maxDepth_ was reached
		if (!foundIntersection || bounces >= maxDepth)
		{
			break;
		}

		// Compute scattering functions and skip over medium boundaries
		isect.ComputeScatteringFunctions(ray, true);
		if (!isect.bsdf) {
			//VLOG(2) << "Skipping intersection due to null bsdf";
			ray = isect.SpawnRay(ray.d);
			bounces--;
			continue;
		}

		//if (L.HasNaNs())
		//{
		//	cout << L << endl;
		//}
		const Distribution1D *distrib = lightDistribution->Lookup(isect.p);

		// Sample illumination from lights to find path contribution.
		// (But skip this for perfectly specular BSDFs.)
		if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >
			0) {
			++totalPaths;
			Spectrum Ld = beta * UniformSampleOneLight(isect, scene,
				sampler, false, distrib);
			//cout << beta << endl;
			//VLOG(2) << "Sampled direct lighting Ld = " << Ld;
			if (Ld.IsBlack()) ++zeroRadiancePaths;
			//cout << Ld << endl;
			CHECK_GE(Ld.y(), 0.f);
			L += Ld;
		}

		//if (L.HasNaNs())
		//{
		//	cout << L << endl;
		//}

		// Sample BSDF to get new path direction
		Vector3f wo = -ray.d, wi;
		Float pdf;
		BxDFType flags;
		Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
			BSDF_ALL, &flags);
		//VLOG(2) << "Sampled BSDF, f = " << f << ", pdf = " << pdf;
		if (f.IsBlack() || pdf == 0.f) break;
		beta *= f * AbsDot(wi, isect.shading.n) / pdf;
		//VLOG(2) << "Updated beta = " << beta;
		CHECK_GE(beta.y(), 0.f);
		DCHECK(!std::isinf(beta.y()));
		specularBounce = (flags & BSDF_SPECULAR) != 0;
		if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
			Float eta = isect.bsdf->eta;
			// Update the term that tracks radiance scaling for refraction
			// depending on whether the ray is entering or leaving the
			// medium.
			etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
		}
		ray = isect.SpawnRay(wi);

		//cout << L << endl;
		// Account for subsurface scattering, if applicable

		// Possibly terminate the path with Russian roulette.
		// Factor out radiance scaling due to refraction in rrBeta.
		Spectrum rrBeta = beta * etaScale;
		if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
			Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
			if (sampler.Get1D() < q) break;
			beta /= 1 - q;
			DCHECK(!std::isinf(beta.y()));
		}
		//break;
	}
	//if (L.HasNaNs())
	//{
	//	cout << L << endl;
	//}
	//ReportValue(pathLength, bounces);
	//cout << L << endl;
	return L;
}
Spectrum SamplerIntegrator::Li_re(const Ray &r, const Scene &scene,
	Sampler &sampler, int depth) const {
	return Spectrum(0.f);
}

bool perfectspecular = false;
Spectrum PathIntegrator::Li_re(const Ray &r, const Scene &scene,
	Sampler &sampler, int depth) const {
	Spectrum L(0.f), beta(1.f);
	Ray ray(r);
	SurfaceInteraction isect;
	bool foundIntersection = scene.Intersect(ray, &isect);
	// Possibly add emitted light at intersection
	if (depth == 0 || perfectspecular|| scene.lights.size()==0) 
	{
		// Add emitted light at path vertex or from the environment
		if (foundIntersection) {
			L += isect.Le(-ray.d) ;
			//L += isect.Le(-ray.d) / std::max(1.f, DistanceSquared(ray.o, isect.p));
			if (!L.IsBlack())
				return L;
		}
		else if (scene.lights.size() == 0)
		{
			Vector3f unit_direction = Normalize(r.d);
			auto t = 0.5*(unit_direction[1] + 1.0);
			L+= (1.0 - t)*Spectrum(1.0, 1.0, 1.0) + t * Spectrum(0.5, 0.7, 1.0);
			//L += Spectrum(0.8f);
		}
	}
	//if (L.HasNaNs())
	//	cout << L << endl;

	// Terminate path if ray escaped or _maxDepth_ was reached
	if (!foundIntersection || depth >= maxDepth)
		return L;

	// Compute scattering functions and skip over medium boundaries
	isect.ComputeScatteringFunctions(ray, true);
	//if (L.HasNaNs())
	//	cout << L << endl;
	//const Distribution1D *distrib = lightDistribution->Lookup(isect.p);
	
	// Sample illumination from lights to find path contribution.
	// (But skip this for perfectly specular BSDFs.)
	if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >0) 
	{
		Spectrum Ld = UniformSampleOneLight(isect, scene,sampler, false, nullptr);
		//cout << beta << endl;
		if (Ld.IsBlack()) 
			++zeroRadiancePaths;
		//cout << Ld << endl;
		L += Ld;
	}
	//if (L.HasNaNs())
	//	cout << L << endl;
	Vector3f wo = -ray.d, wi;
	Float pdf;
	BxDFType flags;
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &flags);
	if (f.IsBlack() || pdf == 0.f) 
		return L ;

	beta *= f * AbsDot(wi, isect.shading.n) / pdf;
	perfectspecular = (flags & BSDF_SPECULAR) != 0;
	//Float dis = DistanceSquared(ray.o, isect.p);
	ray = isect.SpawnRay(wi);
	//cout << L << endl;
	if (depth > 3) {
		Float q = std::max((Float).05, 1 - beta.MaxComponentValue());
		if (sampler.Get1D() < q) 
			return L;
		beta /= 1 - q;
	}
	//cout << DistanceSquared(ray.o, isect.p) << endl;
	//cout << depth << endl;
	return L + beta*Li_re(ray, scene, sampler, depth + 1);
	//if (L.HasNaNs())
	//	cout << L << endl;
}

inline Float clamp(Float x, Float min, Float max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

inline void write_color(Spectrum pixel_color, int samples_per_pixel, unsigned char *data,int j,int i,int w) {
	auto r = pixel_color[0];
	auto g = pixel_color[1];
	auto b = pixel_color[2];

	if (r != r) r = 0.0;
	if (g != g) g = 0.0;
	if (b != b) b = 0.0;

	// Divide the color by the number of samples.

	auto scale = 1.0 / samples_per_pixel;
	r = sqrt(scale * r);
	g = sqrt(scale * g);
	b = sqrt(scale * b);

	// Write the translated [0,255] value of each color component.
	//out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
	//	<< static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
	//	<< static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';

	data[w*j * 3 + i * 3 + 0] = static_cast<int>(256 * clamp(r, 0.0, 0.999));
	data[w*j * 3 + i * 3 + 1] = static_cast<int>(256 * clamp(g, 0.0, 0.999));
	data[w*j * 3 + i * 3 + 2] = static_cast<int>(256 * clamp(b, 0.0, 0.999));
}
unsigned char *data;
vector<thread>pool;
//mutex mt;
//std::condition_variable cv;
int thread_count=8;
int iw, ih;
int maxpos;
int now;
mutex out_mt;
const Scene *sc;
void SamplerIntegrator::render_pixel(int nowpos)
{
	int j = nowpos / iw, i = nowpos - j * iw;
	if (i == iw - 1 && j % 10 == 0)
		cout << j << endl;
	//cout << j << ' ' << i << endl;
	//i = 90; j = 90;
	std::unique_ptr<Sampler>pixel_sampler = sampler->Clone(i + j * iw);
	Point2i pixel = Point2i(i, j);
	pixel_sampler->StartPixel(pixel);
	Spectrum colObj(0.0);
	do {
		Ray r;
		CameraSample cs;
		cs = pixel_sampler->GetCameraSample(pixel);
		//cs.pFilm[0] = i; cs.pFilm[1] = j;
		//cout << i << ' ' << j << ' ' << ' '<<cs.pFilm[0] << ' '<< cs.pFilm[1]<<endl;
		camera->GenerateRay(cs, &r);
		//float tHit;
		SurfaceInteraction  isect;
		//colObj += Li(r, *sc, *pixel_sampler, 0);
		perfectspecular = false;
		colObj += Li_re(r, *sc, *pixel_sampler, 0);

	} while (pixel_sampler->StartNextSample());
	//cout << i<<' '<<j<<' '<<colObj << endl;
	if (colObj.HasNaNs())
	{
		cout << i << ' ' << j << ' ' << colObj << endl;
		//colObj[0] = colObj[1] = colObj[2] = 1.0f;
	}
	colObj = colObj / pixel_sampler->samplesPerPixel;
	//cout << i << ' ' << j << ' ' << colObj << endl;
	write_color(colObj, 1, data, j, i, iw);
}
void SamplerIntegrator::Render(const Scene &scene,int image_width,int image_height) {
	iw = image_width;ih = image_height;
	sc = &scene;
	//Preprocess(scene, *sampler);
	data = new unsigned char[image_height*image_width*3];
	maxpos = image_height * image_width;
	if (thread_count == 1)
	{
		for (int j = 0; j < image_height; j++) {
			if (j % 10 == 0)
				cout << j << endl;
			for (int i = 0; i < image_width; i++) {
				//cout << j << ' ' << i << endl;
				i = 445; j = 400;
				//if (j < 85 || j>90 || i < 86 || i>90)
				//{
				//	Spectrum colObj(1.0f,0.0f,0.0f);
				//	write_color(colObj, 1, data, j, i, iw);
				//	continue;
				//}
				std::unique_ptr<Sampler>pixel_sampler = sampler->Clone(i + j * iw);
				Point2i pixel = Point2i(i, j);
				pixel_sampler->StartPixel(pixel);
				Spectrum colObj(0.0);
				do {
					Ray r;
					CameraSample cs;
					cs = pixel_sampler->GetCameraSample(pixel);
					cs.pFilm[0] = i; cs.pFilm[1] = j;
					//cout << i << ' ' << j << ' ' << ' '<<cs.pFilm[0] << ' '<< cs.pFilm[1]<<endl;
					camera->GenerateRay(cs, &r);
					//float tHit;
					SurfaceInteraction  isect;
					colObj += Li_re(r, *sc, *pixel_sampler, 0);

				} while (pixel_sampler->StartNextSample());
				cout << i<<' '<<j<<' '<<colObj << ' '<< pixel_sampler->samplesPerPixel<<endl;
				if (colObj.HasNaNs())
				{
					cout << 'h'<<' '<<i << ' ' << j << ' ' << colObj << endl;
					//colObj[0] = colObj[1] = colObj[2] = 1.0f;
				}
				else if (colObj.IsBlack())
				{
					cout << 'b' << ' ' << i << ' ' << j << ' ' << colObj << endl;
				}
				colObj = colObj / pixel_sampler->samplesPerPixel;
				cout << i << ' ' << j << ' ' << colObj << endl;
				write_color(colObj, 1, data, j, i, iw);
				return;
			}
		}
		stbi_write_png("image.png", image_width, image_height, 3, data, 0);
		std::cerr << "\nDone.\n";
		return;
	}
	for (int i = 1; i <= thread_count-1; ++i)
	{
		pool.push_back(thread([&] {
			{
				std::unique_lock<mutex> loc(out_mt);
				while (now < maxpos)
				{
					int tmp = now;
					now++;
					loc.unlock();
					render_pixel(tmp);
					loc.lock();
				}
			}
		}));
	}
	{
		std::unique_lock<mutex> loc(out_mt);
		while (now < maxpos)
		{
			int tmp = now;
			now++;
			loc.unlock();
			render_pixel(tmp);
			loc.lock();
		}
	}
	for (int i = 0; i < pool.size(); ++i)
		pool[i].join();
	stbi_write_png("image.png", image_width, image_height, 3, data, 0);
	std::cerr << "\nDone.\n";
	return;
}
