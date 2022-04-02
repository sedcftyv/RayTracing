#include <fstream>
#include <iostream>
#include "rtweekend.h"
#include "integrator.h"
#include "interaction.h"
#include "reflection.h"
#include "light.h"

using std::cout;
using std::endl;

Integrator::~Integrator() {}


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
		if (depth == 0)
		{
			Spectrum LL;
			LL[2] = 0.8;
			return LL;
		}
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


inline void write_color(std::ostream &out, Spectrum pixel_color, int samples_per_pixel) {
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
	out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
		<< static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
		<< static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

void SamplerIntegrator::Render(const Scene &scene) {

	std::ofstream fout("image.ppm");
	int image_height = 300;
	int image_width = 300;
	fout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
	Vector3f Light(-2.0,4.0,-3.0);
	Point3f LightPosition(1.0, 4, -3.0);
	for (int j = 0; j < image_height; j++) {
		for (int i = 0; i < image_width; i++) {
			std::unique_ptr<Sampler>pixel_sampler = sampler->Clone(i + j * image_width);
			Point2i pixel = Point2i(i, j);
			pixel_sampler->StartPixel(pixel);
			Spectrum colObj(0.0);
			do {
				Ray r;
				CameraSample cs;
				cs = pixel_sampler->GetCameraSample(pixel);
				camera->GenerateRay(cs, &r);
				//float tHit;
				SurfaceInteraction  isect;
				//r.d = Vector3f(0.663268, 0.273197, -0.696735);
				colObj += Li(r,scene,*pixel_sampler,0);
				//cout << r.d << ' ' << colObj << endl;
				//if (scene.Intersect(r, &isect))
				//{
				//	//cout << isect.p << endl;
				//	Interaction p1;
				//	Vector3f wi;
				//	VisibilityTester vist;
				//	Float pdf_light;
				//	Spectrum colTmp(0.0);
				//	for (int count = 0; count < scene.lights.size(); count++)
				//	{
				//		Spectrum Li = scene.lights[count]->Sample_Li(isect, pixel_sampler->Get2D(), &wi, &pdf_light, &vist);
				//		if (vist.Unoccluded(scene))
				//		{
				//			//cout << Li << endl;
				//			isect.ComputeScatteringFunctions(r);
				//			Vector3f wo = isect.wo;

				//			Spectrum f = isect.bsdf->f(wo, wi);
				//			Float pdf_scattering = isect.bsdf->Pdf(wo, wi);
				//			//cout << f << ' ' << pdf << endl;
				//			colTmp += Li * pdf_scattering * f*3.0f / pdf_light;
				//		}
				//	}
				//	colTmp /= scene.lights.size();
				//	colObj += colTmp;
				//}
			} while (pixel_sampler->StartNextSample());
		//	cout << colObj << endl;
			//if (colObj.HasNaNs())
			//{
			//	cout << 1 << endl;
			//	colObj[0] = colObj[1] = colObj[2] = 1.0f;
			//}
			colObj = colObj / pixel_sampler->samplesPerPixel;
			//colObj = colObj / success;
			write_color(fout, colObj, 1);
		}
	}
	fout.close();
	std::cerr << "\nDone.\n";

}