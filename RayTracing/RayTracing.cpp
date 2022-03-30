#include <iostream>
#include <fstream>
#include "rtweekend.h"
//#include "color.h"
//#include "ray.h"
//#include "hittable_list.h"
//#include "sphere.h"
//#include "camera.h"
//#include "material.h"
//#include "aarect.h"
//#include "box.h"
//#include "bvh.h"
//#include "pdf.h"
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
#include "core/stratified.h"
#include "core/scene.h"
#include "core/integrator.h"
#include "core/texture.h"
using std::cout;
using std::endl;
using std::make_shared;
using std::unique_ptr;
using std::make_unique;
using std::shared_ptr;

//color ray_color(const ray& r, const color& background, const hittable& world, shared_ptr<hittable>& lights, int depth) {
//	hit_record rec;
//	if (depth <= 0)
//		return color(0, 0, 0);
//
//	if (!world.hit(r, 0.001, infinity, rec))
//		return background;
//
//	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
//	scatter_record srec;
//
//	if (!rec.mat_ptr->scatter(r, rec, srec))
//		return emitted;
//
//	if (srec.is_specular) {
//		return srec.attenuation
//			* ray_color(srec.specular_ray, background, world, lights, depth - 1);
//	}
//	auto p0 = make_shared<hittable_pdf>(lights, rec.p);
//	auto p1 = make_shared<cosine_pdf>(rec.normal);
//	mixture_pdf mixed_pdf(p0, p1);
//
//	auto scattered = ray(rec.p, mixed_pdf.generate(), r.time());
//	auto pdf_val = mixed_pdf.value(scattered.direction());
//	//auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
//	//mixture_pdf p(light_ptr, srec.pdf_ptr);
//
//	//ray scattered = ray(rec.p, p.generate(), r.time());
//	//double pdf_val = p.value(scattered.direction());
//
//	return emitted
//		+ srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
//		* ray_color(scattered, background, world, lights, depth - 1) / pdf_val;
//}


int main()
{
	//std::ofstream fout("image.ppm");

	// Image
	auto aspect_ratio = 1;
	int image_width = 300;
	int samples_per_pixel = 100;
	int max_depth = 50;

	// World
	//hittable_list world;
	//Vector3f lookfrom;
	//Vector3f lookat;
	auto vfov = 40.0;
	auto aperture = 0.0;
	color background(0, 0, 0);

	// Render
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	//fout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

	//triangle
	Transform tri_Object2World, tri_World2Object;
	std::vector<Point3f>p;
	std::vector<int>vi;
	string path("C:\\VS\\5_ARAP\\project\\bin\\Balls.obj");
	objload(path, p, vi);
	int nTriangles = vi.size() / 3;
	int nVertices = p.size();
	int *vertexIndices = new int[nTriangles * 3];
	int tmp = p.size();
	Point3f* P = new Point3f[tmp];
	for (int i = 0; i < vi.size(); ++i)
		vertexIndices[i] = vi[i];
	//Point3f ct;
	for (int j = 0; j < nVertices; ++j)
	{
		P[j] = p[j];
		P[j][2] += 10.0;
		//ct += P[j];
	}
	//ct /= nVertices;
	//cout << ct << endl;
	std::vector<std::shared_ptr<Shape>> tris = CreateTriangleMesh(&tri_Object2World, &tri_World2Object, false, nTriangles, vertexIndices, nVertices, P, nullptr, nullptr, nullptr, nullptr);
	int nTrianglesFloor = 2;
	int vertexIndicesFloor[6] = { 0,1,2,3,4,5 };
	int nVerticesFloor = 6;
	const float yPos_Floor = -7;
	Point3f P_Floor[6] = {
			Point3f(-20.0,yPos_Floor,20.0),Point3f(20.0,yPos_Floor,20.0),Point3f(-20.0,yPos_Floor,-20.0),
			Point3f(20.0,yPos_Floor,20.0),Point3f(20.0,yPos_Floor,-20.0),Point3f(-20.0,yPos_Floor,-20.0)};
	std::vector<std::shared_ptr<Shape>> meshFloor = CreateTriangleMesh(&tri_Object2World, &tri_World2Object, false, nTrianglesFloor, vertexIndicesFloor, nVerticesFloor, P_Floor, nullptr, nullptr, nullptr, nullptr);
	
	Spectrum floorColor(0.2, 0.3, 0.9);
	Spectrum Color(0.8, 0.1, 0.2);
	shared_ptr<Texture<Spectrum>> Kd1 = make_shared<ConstantTexture<Spectrum>>(floorColor);
	shared_ptr<Texture<Spectrum>> Kd2 = make_shared<ConstantTexture<Spectrum>>(Color);
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(Kd1, sigma, bump);
	shared_ptr<Material> m2 = make_shared<MatteMaterial>(Kd2, sigma, bump);

	std::vector<std::shared_ptr<Primitive>> prims;
	for (int i = 0; i < nTriangles; ++i)
		prims.push_back(make_shared<GeometricPrimitive>(tris[i], m2));
	for (int i = 0; i < nTrianglesFloor; ++i)
		prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i],m1));

	std::unique_ptr<Scene>worldScene;
	worldScene = std::make_unique<Scene>(make_unique<BVHAccel>(prims));

	shared_ptr<Camera> cam;
	Point3f eye(4.0f, 4.f, 2.f), look(0.0f, 0.0f, 10.0f);
	Vector3f up(0.0f, 1.0f, 0.0f);
	Transform lookat = LookAt(Vector3f(eye), Vector3f(look), up);
	Transform Camera2World = Inverse(lookat);
	cam = shared_ptr<Camera>(CreatePerspectiveCamera(Camera2World));
	shared_ptr<Sampler> ss = make_unique<StratifiedSampler>(2,2,true,1);
	Bounds2i pixelBounds;
	shared_ptr<Integrator>SI=make_shared<SamplerIntegrator>(cam,ss, pixelBounds);

	SI->Render(*worldScene);
}
