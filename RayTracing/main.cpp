#include <iostream>
#include <fstream>
#include "rtweekend.h"
//#include "color.h"
#include "ray.h"
#include "hittable_list.h"
//#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "aarect.h"
#include "box.h"
#include "bvh.h"
#include "pdf.h"
#include "core/geometry.h"
#include "core/transform.h"
#include "core/sphere.h"
#include "core/triangle.h"
#include "core/primitive.h"
#include "core/bvh.h"
#include "core/kdtreeaccel.h"
#include "core/spectrum.h"
#include "core/camera.h"
#include "core/objload.h"
using std::cout;
using std::endl;
using std::make_shared;

void write_color(std::ostream &out, Spectrum pixel_color, int samples_per_pixel) {
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

color ray_color(const ray& r, const color& background, const hittable& world, shared_ptr<hittable>& lights, int depth) {
	hit_record rec;
	if (depth <= 0)
		return color(0, 0, 0);

	if (!world.hit(r, 0.001, infinity, rec))
		return background;

	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
	scatter_record srec;

	if (!rec.mat_ptr->scatter(r, rec, srec))
		return emitted;

	if (srec.is_specular) {
		return srec.attenuation
			* ray_color(srec.specular_ray, background, world, lights, depth - 1);
	}
	auto p0 = make_shared<hittable_pdf>(lights, rec.p);
	auto p1 = make_shared<cosine_pdf>(rec.normal);
	mixture_pdf mixed_pdf(p0, p1);

	auto scattered = ray(rec.p, mixed_pdf.generate(), r.time());
	auto pdf_val = mixed_pdf.value(scattered.direction());
	//auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
	//mixture_pdf p(light_ptr, srec.pdf_ptr);

	//ray scattered = ray(rec.p, p.generate(), r.time());
	//double pdf_val = p.value(scattered.direction());

	return emitted
		+ srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
		* ray_color(scattered, background, world, lights, depth - 1) / pdf_val;
}


int main()
{
	std::ofstream fout("image.ppm");

	// Image
	auto aspect_ratio = 1;
	int image_width = 300;
	int samples_per_pixel = 100;
	int max_depth = 50;

	// World
	hittable_list world;
	//Vector3f lookfrom;
	//Vector3f lookat;
	auto vfov = 40.0;
	auto aperture = 0.0;
	color background(0, 0, 0);

	// Render
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	fout << "P3\n" << image_width << ' ' << image_height << "\n255\n";


	//triangle
	Transform tri_Object2World, tri_World2Object;
	std::vector<Point3f>p;
	std::vector<int>vi;
	int nTriangles = int(vi.size()) / 3;
	int nVertices = p.size();
	int *vertexIndices = new int[nTriangles * 3];
	int tmp = p.size();
	Point3f* P = new Point3f[tmp];
	for (int i = 0; i < vi.size(); ++i)
		vertexIndices[i] = vi[i];
	for (int j = 0; j < nVertices; ++j)
		P[j] = p[j];
	std::vector<std::shared_ptr<Shape>> tris = CreateTriangleMesh(&tri_Object2World, &tri_World2Object,false,nTriangles,vertexIndices,nVertices,P,nullptr,nullptr,nullptr,nullptr);
	std::vector<std::shared_ptr<Primitive>> prims;
	for (int i = 0; i < nTriangles; ++i)
		prims.push_back(make_shared<GeometricPrimitive>(tris[i]));
	Aggregate *agg;
	agg = new BVHAccel(prims);
}