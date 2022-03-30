#include <fstream>
#include <iostream>
#include "rtweekend.h"
#include "integrator.h"
#include "interaction.h"
#include "reflection.h"

Integrator::~Integrator() {}



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
	for (int j = 0; j < image_height; j++) {
		for (int i = 0; i < image_width; i++) {
			//float v = float(i + 0.5) / float(image_height); //random ( )
			//float u = float(j + 0.5) / float(image_width); // 0.5
			//cout << (lower_left_corner + u * horizontal + v * vertical) << endl;
			//Ray r(Vector3f(origin), (lower_left_corner + u * horizontal + v* vertical) - Vector3f(origin));
			std::unique_ptr<Sampler>pixel_sampler = sampler->Clone(i + j * image_width);
			Point2i pixel = Point2i(i, j);
			pixel_sampler->StartPixel(pixel);
			Spectrum colObj(0.0);
			do {
				Ray r;
				CameraSample cs;
				cs = pixel_sampler->GetCameraSample(pixel);
				camera->GenerateRay(cs, &r);
				float tHit;
				SurfaceInteraction  isect;
				//colObj[0] = 1.0f;
				//colObj[1] = 1.0f;

				if (scene.Intersect(r, &isect))
				{
					isect.ComputeScatteringFunctions(r);
					Vector3f wo = isect.wo;

					Vector3f LightNorm = Light - isect.p;
					LightNorm = Normalize(LightNorm);
					Vector3f wi = LightNorm;
					Spectrum f = isect.bsdf->f(wo, wi);
					float pdf = isect.bsdf->Pdf(wo, wi);
					//colObj += pdf * f*4.0f;

					Vector3f viewInv = -r.d;
					Vector3f H = Normalize(viewInv+LightNorm);
					float Ls = Dot(H, isect.n);
					Ls = (Ls > 0.0f) ? Ls : 0.0f;
					Ls = pow(Ls, 32);
					float Ld = Dot(LightNorm, isect.n);
					Ld = (Ld > 0.0f) ? Ld : 0.0f;
					float Li = 1.5*(0.2 + 0.2*Ld + 0.7*Ls);
					colObj += Li* pdf * f*4.0f;
					//float Li = Dot(Light, isect.n);
					//colObj[1] = std::abs(Li);
					//cout << r.d[1] << endl;
				}
			} while (pixel_sampler->StartNextSample());
			colObj = colObj / pixel_sampler->samplesPerPixel;
			write_color(fout, colObj, 1);
		}
	}

	fout.close();
	std::cerr << "\nDone.\n";

}