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
		Spectrum L(0.f);
	for (size_t j = 0; j < scene.lights.size(); ++j) {
				const std::shared_ptr<Light> &light = scene.lights[j];
		int nSamples = nLightSamples[j];
		const Point2f *uLightArray = sampler.Get2DArray(nSamples);
		const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);
		if (!uLightArray || !uScatteringArray) {
			Point2f uLight = sampler.Get2D();
			Point2f uScattering = sampler.Get2D();
			L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
				 handleMedia);
		}
		else {
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
	Vector3f wi;
	Float lightPdf = 0, scatteringPdf = 0;
	VisibilityTester visibility;
	Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
	if (lightPdf > 0 && !Li.IsBlack()) {
		Spectrum f;
		if (it.IsSurfaceInteraction()) {
			const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
			f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *AbsDot(wi, isect.shading.n);
		}
		if (!f.IsBlack()) {
			if (handleMedia) 
				Li *= visibility.Tr(scene, sampler);
			else {
				if (!visibility.Unoccluded(scene))
					Li = Spectrum(0.f);
				}
			if (!Li.IsBlack()) 
				Ld += f * Li / lightPdf;
		}
	}																																		
	return Ld;
}

Spectrum SamplerIntegrator::SpecularReflect(
	const Ray &ray, const SurfaceInteraction &isect,
	const Scene &scene, Sampler &sampler, int depth) const {
	Vector3f wo = isect.wo, wi;
	Float pdf;
	BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);
	const Normal3f &ns = isect.shading.n;
	if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
		Ray rd = isect.SpawnRay(wi);
		return f * Li(rd, scene, sampler, depth + 1) * AbsDot(wi, ns) /pdf;
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
	Spectrum f = bsdf.Sample_f(wo, &wi, sampler.Get2D(), &pdf,BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
	Spectrum L = Spectrum(0.f);
	Normal3f ns = isect.shading.n;
	if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
	Ray rd = isect.SpawnRay(wi);																																																						
	L = f * Li(rd, scene, sampler, depth + 1) * AbsDot(wi, ns) / pdf;
	}
	return L;
}

Spectrum WhittedIntegrator::Li(const Ray &ray, const Scene &scene,
	Sampler &sampler,int depth) const {
	Spectrum L(0.);
	SurfaceInteraction isect;
	if (!scene.Intersect(ray, &isect)) {
		for (const auto &light : scene.lights) L += light->Le(ray);
		return L;
	}	
	const Normal3f &n = isect.shading.n;
	Vector3f wo = isect.wo;
	isect.ComputeScatteringFunctions(ray);
	if (!isect.bsdf)
		return Li(isect.SpawnRay(ray.d), scene, sampler, depth);
		L += isect.Le(wo);
		for (const auto &light : scene.lights) 
		{
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
	Spectrum L(0.f), beta(1.f);
	Ray ray(r);
	bool specularBounce = false;
	int bounces;
	Float etaScale = 1;
	for (bounces = 0;; ++bounces) {					
		SurfaceInteraction isect;
		bool foundIntersection = scene.Intersect(ray, &isect);
		if (bounces == 0 || specularBounce) {
			if (foundIntersection) 
			L += beta * isect.Le(-ray.d);
		else {
			for (const auto &light : scene.infiniteLights)
				L += beta * light->Le(ray);
			}
		}
		if (!foundIntersection || bounces >= maxDepth)
			break;
		isect.ComputeScatteringFunctions(ray, true);
		if (!isect.bsdf) {
			ray = isect.SpawnRay(ray.d);
			bounces--;
			continue;
		}
		const Distribution1D *distrib = lightDistribution->Lookup(isect.p);
		if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >0) {
			++totalPaths;
			Spectrum Ld = beta * UniformSampleOneLight(isect, scene,sampler, false, distrib);
			if (Ld.IsBlack()) ++zeroRadiancePaths;
			L += Ld;
		}							
		Vector3f wo = -ray.d, wi;
		Float pdf;
		BxDFType flags;
		Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
			BSDF_ALL, &flags);
				if (f.IsBlack() || pdf == 0.f) break;
		beta *= f * AbsDot(wi, isect.shading.n) / pdf;
		specularBounce = (flags & BSDF_SPECULAR) != 0;
		if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
			Float eta = isect.bsdf->eta;
			etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
		}
		ray = isect.SpawnRay(wi);		
		Spectrum rrBeta = beta * etaScale;
		if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
			Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
			if (sampler.Get1D() < q) break;
			beta /= 1 - q;
		}
	}
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
	if (depth == 0 || perfectspecular|| scene.lights.size()==0) 
	{
		if (foundIntersection) {
			L += isect.Le(-ray.d) ;
			if (!L.IsBlack())
			return L;
		}
		else if (scene.lights.size() == 0)
		{
			Vector3f unit_direction = Normalize(r.d);
			auto t = 0.5*(unit_direction[1] + 1.0);
			L+= (1.0 - t)*Spectrum(1.0, 1.0, 1.0) + t * Spectrum(0.5, 0.7, 1.0);
		}
	}
	if (!foundIntersection || depth >= maxDepth)
		return L;
	isect.ComputeScatteringFunctions(ray, true);	
	if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >0) {
		Spectrum Ld = UniformSampleOneLight(isect, scene,sampler, false, nullptr);
		if (Ld.IsBlack()) 
			++zeroRadiancePaths;
		L += Ld;
	}
	Vector3f wo = -ray.d, wi;
	Float pdf;
	BxDFType flags;
	Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &flags);
	if (f.IsBlack() || pdf == 0.f) 
		return L ;
	beta *= f * AbsDot(wi, isect.shading.n) / pdf;
	perfectspecular = (flags & BSDF_SPECULAR) != 0;
	ray = isect.SpawnRay(wi);
	if (depth > 3) {
		Float q = std::max((Float).05, 1 - beta.MaxComponentValue());
		if (sampler.Get1D() < q) 
			return L;
		beta /= 1 - q;
	}
	return L + beta*Li_re(ray, scene, sampler, depth + 1);
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

	
	auto scale = 1.0 / samples_per_pixel;
	r = sqrt(scale * r);
	g = sqrt(scale * g);
	b = sqrt(scale * b);

				
	data[w*j * 3 + i * 3 + 0] = static_cast<int>(256 * clamp(r, 0.0, 0.999));
	data[w*j * 3 + i * 3 + 1] = static_cast<int>(256 * clamp(g, 0.0, 0.999));
	data[w*j * 3 + i * 3 + 2] = static_cast<int>(256 * clamp(b, 0.0, 0.999));
}
unsigned char *data;
vector<thread>pool;
int thread_count=8;
int iw, ih;
int maxpos;
int now;
mutex out_mt;
const Scene *sc;
void SamplerIntegrator::render_pixel(int nowpos){
	int j = nowpos / iw, i = nowpos - j * iw;
	if (i == iw - 1 && j % 10 == 0)
		cout << j << endl;
	std::unique_ptr<Sampler>pixel_sampler = sampler->Clone(i + j * iw);
	Point2i pixel = Point2i(i, j);
	pixel_sampler->StartPixel(pixel);
	Spectrum colObj(0.0);
	do {
		Ray r;
		CameraSample cs;
		cs = pixel_sampler->GetCameraSample(pixel);
		camera->GenerateRay(cs, &r);
		SurfaceInteraction  isect;
		perfectspecular = false;
		colObj += Li_re(r, *sc, *pixel_sampler, 0);
	} while (pixel_sampler->StartNextSample());
	if (colObj.HasNaNs()){
		cout << i << ' ' << j << ' ' << colObj << endl;
	}
	colObj = colObj / pixel_sampler->samplesPerPixel;
		write_color(colObj, 1, data, j, i, iw);
}
void SamplerIntegrator::Render(const Scene &scene,int image_width,int image_height) {
	iw = image_width;ih = image_height;
	sc = &scene;
	data = new unsigned char[image_height*image_width*3];
	maxpos = image_height * image_width;
	clock_t start = clock(), end;
	if (thread_count == 1){
		for (int j = 0; j < image_height; j++) {
			if (j % 10 == 0)
				cout << j << endl;
			for (int i = 0; i < image_width; i++) {
				std::unique_ptr<Sampler>pixel_sampler = sampler->Clone(i + j * iw);
				Point2i pixel = Point2i(i, j);
				pixel_sampler->StartPixel(pixel);
				Spectrum colObj(0.0);
				do {
					Ray r;
					CameraSample cs;
					cs = pixel_sampler->GetCameraSample(pixel);
					cs.pFilm[0] = i; cs.pFilm[1] = j;
					camera->GenerateRay(cs, &r);
					SurfaceInteraction  isect;
					colObj += Li_re(r, *sc, *pixel_sampler, 0);
				} while (pixel_sampler->StartNextSample());
				if (colObj.HasNaNs()){
					cout << 'h'<<' '<<i << ' ' << j << ' ' << colObj << endl;
				}
				colObj = colObj / pixel_sampler->samplesPerPixel;
				write_color(colObj, 1, data, j, i, iw);
			}
		}
		end = clock();
		double t = (double)(end - start) / CLOCKS_PER_SEC;
		cout << t << endl;
		stbi_write_png("image.png", image_width, image_height, 3, data, 0);
		std::cerr << "\nDone.\n";
		return;
	}
	for (int i = 1; i <= thread_count-1; ++i){
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
	end = clock();
	double t = (double)(end - start) / CLOCKS_PER_SEC;
	cout << t << endl;
	stbi_write_png("image.png", image_width, image_height, 3, data, 0);
	std::cerr << "\nDone.\n";
	return;
}
