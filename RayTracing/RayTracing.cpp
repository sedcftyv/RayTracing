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
#include "core/light.h"
using std::cout;
using std::endl;
using std::make_shared;
using std::unique_ptr;
using std::make_unique;
using std::shared_ptr;
using std::vector;

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
void Cornellbox(shared_ptr<Scene> &scene)
{
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	Spectrum green(0.0f, 1.0f,0.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor) , sigma, bump);
	shared_ptr<Material> mgreen = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(green), sigma, bump);

	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	Transform tri_Object2WorldTri, tri_World2ObjectTri;
	std::vector<Point3f>p;
	std::vector<int>vi;
	string path("C:\\VS\\5_ARAP\\project\\bin\\Beetle_ABF.obj");
	objload(path, p, vi);
	int nTriangles = vi.size() / 3;
	int nVertices = p.size();
	int *vertexIndices = new int[nTriangles * 3];
	int tmp = p.size();
	Point3f* P = new Point3f[tmp];
	for (int i = 0; i < vi.size(); ++i)
		vertexIndices[i] = vi[i];
	for (int j = 0; j < nVertices; ++j)
		P[j] = p[j];
	tri_Object2WorldTri = Translate(Vector3f(2.5f, 1.0f, 1.5f))*Scale(2.0f, 2.0f, 2.0f)*RotateY(90)*tri_Object2WorldTri;
	tri_World2ObjectTri = Inverse(tri_Object2WorldTri);
	std::vector<std::shared_ptr<Shape>> tris = CreateTriangleMesh(&tri_Object2WorldTri, &tri_World2ObjectTri, false, nTriangles, vertexIndices, nVertices, P, nullptr, nullptr, nullptr, nullptr);
	for (int i = 0; i < nTriangles; ++i)
		prims.push_back(make_shared<GeometricPrimitive>(tris[i], mgreen, area));
	delete []P;
	delete []vertexIndices;

	float length_Floor = 5;
	vector<shared_ptr<Light>>lights;
	const float yPos_AreaLight = length_Floor-0.01;
	float len = 0.5;
	int nTrianglesAreaLight = 2;
	int vertexIndicesAreaLight[6] = { 0,1,2,3,4,5 };
	int nVerticesAreaLight = 6;
	Point3f P_AreaLight[6] = {
			Point3f(-len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,len),
			Point3f(len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,-len) };
	Transform tri_Object2World_AreaLight;
	tri_Object2World_AreaLight = Translate(Vector3f(2.5f, 0.0f, 2.5f))*tri_Object2World_AreaLight;
	Transform tri_World2Object_AreaLight = Inverse(tri_Object2World_AreaLight);
	std::vector<std::shared_ptr<Shape>> meshAreaLight = CreateTriangleMesh(&tri_Object2World_AreaLight, &tri_World2Object_AreaLight, false, nTrianglesAreaLight, vertexIndicesAreaLight, nVerticesAreaLight, P_AreaLight, nullptr, nullptr, nullptr, nullptr);
	for (int i = 0; i < nTrianglesAreaLight; ++i)
	{
		area = make_shared<DiffuseAreaLight>(tri_Object2World_AreaLight, Spectrum(15.f), 5, meshAreaLight[i], false);
		lights.push_back(area);
		prims.push_back(make_shared<GeometricPrimitive>(meshAreaLight[i], m1, area));
	}


	Transform tri_Object2World, tri_World2Object;
	const int nTrianglesFloor = 10;
	const int nVerticesFloor = nTrianglesFloor*3;
	int vertexIndicesFloor[nVerticesFloor];
	for (int i = 0; i < nVerticesFloor; ++i)
		vertexIndicesFloor[i] = i;
	Point3f P_Floor[nVerticesFloor] = {
		//底座
		Point3f(0.f,0.f,length_Floor),Point3f(length_Floor,0.f,length_Floor),Point3f(0.f,0.f,0.f),
		Point3f(length_Floor,0.f,length_Floor),Point3f(length_Floor,0.f,0.f),Point3f(0.f,0.f,0.f),
		//天花板
		Point3f(0.f,length_Floor,length_Floor),Point3f(0.f,length_Floor,0.f),Point3f(length_Floor,length_Floor,length_Floor),
		Point3f(length_Floor,length_Floor,length_Floor),Point3f(0.f,length_Floor,0.f),Point3f(length_Floor,length_Floor,0.f),
		//后墙
		Point3f(0.f,0.f,0.f),Point3f(length_Floor,0.f,0.f),Point3f(length_Floor,length_Floor,0.f),
		Point3f(0.f,0.f,0.f),Point3f(length_Floor,length_Floor,0.f),Point3f(0.f,length_Floor,0.f),
		//右墙
		Point3f(0.f,0.f,0.f),Point3f(0.f,length_Floor,length_Floor),Point3f(0.f,0.f,length_Floor),
		Point3f(0.f,0.f,0.f),Point3f(0.f,length_Floor,0.f),Point3f(0.f,length_Floor,length_Floor),
		//左墙
		Point3f(length_Floor,0.f,0.f),Point3f(length_Floor,length_Floor,length_Floor),Point3f(length_Floor,0.f,length_Floor),
		Point3f(length_Floor,0.f,0.f),Point3f(length_Floor,length_Floor,0.f),Point3f(length_Floor,length_Floor,length_Floor)};
	std::vector<std::shared_ptr<Shape>> meshFloor = CreateTriangleMesh(&tri_Object2World, &tri_World2Object, false, nTrianglesFloor, vertexIndicesFloor, nVerticesFloor, P_Floor, nullptr, nullptr, nullptr, nullptr);
	Spectrum blue(0.0, 0.0, 1.0);
	Spectrum red(1.0, 0.0, 0.0);
	Spectrum black(0.0, 0.0, 0.0);
	shared_ptr<Material> mblue = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(blue), sigma, bump);
	shared_ptr<Material> mred = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(red), sigma, bump);
	shared_ptr<Material> mblack = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(black), sigma, bump);

	area.reset();
	for (int i = 0; i < nTrianglesFloor; ++i)
	{
		if(i<2)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, area));
		else if(i<4)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, area));
		else if (i < 6)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, area));
		else if (i < 8)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mred, area));
		else if (i < 10)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mblue, area));
	}


	scene=make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

}
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
	//Transform tri_Object2World, tri_World2Object;
	//std::vector<Point3f>p;
	//std::vector<int>vi;
	//string path("C:\\VS\\5_ARAP\\project\\bin\\Beetle_ABF.obj");
	//objload(path, p, vi);
	//int nTriangles = vi.size() / 3;
	//int nVertices = p.size();
	//int *vertexIndices = new int[nTriangles * 3];
	//int tmp = p.size();
	//Point3f* P = new Point3f[tmp];
	//for (int i = 0; i < vi.size(); ++i)
	//	vertexIndices[i] = vi[i];
	//Point3f ct;
	//for (int j = 0; j < nVertices; ++j)
	//{
	//	P[j] = p[j];
	//	P[j][2] += 0.0;
	//	//cout << P[j] << endl;
	//	//ct += P[j];
	//}
	////ct /= nVertices;
	////cout << ct << endl;
	//std::vector<std::shared_ptr<Shape>> tris = CreateTriangleMesh(&tri_Object2World, &tri_World2Object, false, nTriangles, vertexIndices, nVertices, P, nullptr, nullptr, nullptr, nullptr);
	//int nTrianglesFloor = 2;
	//int vertexIndicesFloor[6] = { 0,1,2,3,4,5 };
	//int nVerticesFloor = 6;
	//const float yPos_Floor = -0.75;
	//Point3f P_Floor[6] = {
	//		Point3f(-20.0,yPos_Floor,20.0),Point3f(20.0,yPos_Floor,20.0),Point3f(-20.0,yPos_Floor,-20.0),
	//		Point3f(20.0,yPos_Floor,20.0),Point3f(20.0,yPos_Floor,-20.0),Point3f(-20.0,yPos_Floor,-20.0)};
	//std::vector<std::shared_ptr<Shape>> meshFloor = CreateTriangleMesh(&tri_Object2World, &tri_World2Object, false, nTrianglesFloor, vertexIndicesFloor, nVerticesFloor, P_Floor, nullptr, nullptr, nullptr, nullptr);
	//
	//Spectrum floorColor(1.0, 1.0, 1.0);
	//Spectrum Color(0.8, 0.1, 0.2);
	//Spectrum mirrorColor(1.f);
	//shared_ptr<Texture<Spectrum>> Kd1 = make_shared<ConstantTexture<Spectrum>>(floorColor);
	//shared_ptr<Texture<Spectrum>> Kd2 = make_shared<ConstantTexture<Spectrum>>(Color);
	//shared_ptr<Texture<Spectrum>> KrMirror = make_shared<ConstantTexture<Spectrum>>(mirrorColor);
	//shared_ptr<Texture<Spectrum>> glassKr = make_shared<ConstantTexture<Spectrum>>(1.0f);
	//shared_ptr<Texture<Spectrum>> glassKt = make_shared<ConstantTexture<Spectrum>>(1.0f);
	//shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	//shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	//shared_ptr<Texture<Float>> glassEta = make_shared<ConstantTexture<Float>>(1.2f);
	//shared_ptr<Material> m1 = make_shared<MatteMaterial>(Kd1, sigma, bump);
	//shared_ptr<Material> m2 = make_shared<MatteMaterial>(Kd2, sigma, bump);
	//shared_ptr<Material> mirrorMaterial = make_shared<MirrorMaterial>(KrMirror, bump);
	//shared_ptr<Material> glassMaterial = make_shared<GlassMaterial>(glassKr,glassKt,glassEta, bump);

	//shared_ptr<AreaLight>area;
	//std::vector<std::shared_ptr<Primitive>> prims;
	//for (int i = 0; i < nTriangles; ++i)
	//	prims.push_back(make_shared<GeometricPrimitive>(tris[i], m2,area));
	//for (int i = 0; i < nTrianglesFloor; ++i)
	//	prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1,area));

	//Spectrum LightI(30.f);
	//Transform LightToWorld;
	//LightToWorld  = Translate(Vector3f(1.0f, 4.5f, -6.0f))*LightToWorld;
	//shared_ptr<Light>pointlight = make_shared<PointLight>(LightToWorld, LightI);
	//vector<shared_ptr<Light>>lights;
	////lights.push_back(pointlight);

	//const float yPos_AreaLight = 1.3;
	//float len = 1;
	//Point3f P_AreaLight[6] = {
	//		Point3f(-len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,len),
	//		Point3f(len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,-len) };
	//Transform tri_Object2World_AreaLight;
	////tri_Object2World_AreaLight = Translate(Vector3f(0.7f, 0.0f, -2.0f))*tri_Object2World_AreaLight;
	//Transform tri_World2Object_AreaLight=Inverse(tri_Object2World_AreaLight);
	//std::vector<std::shared_ptr<Shape>> meshAreaLight = CreateTriangleMesh(&tri_Object2World_AreaLight, &tri_World2Object_AreaLight, false, nTrianglesFloor, vertexIndicesFloor, nVerticesFloor, P_AreaLight, nullptr, nullptr, nullptr, nullptr);
	//for (int i = 0; i < nTrianglesFloor; ++i)
	//{
	//	area = make_shared<DiffuseAreaLight>(tri_Object2World_AreaLight, Spectrum(5.f), 5, meshAreaLight[i], false);
	//	lights.push_back(area);
	//	prims.push_back(make_shared<GeometricPrimitive>(meshAreaLight[i], m1,area));
	//}


	std::shared_ptr<Scene>worldScene;
	Cornellbox(worldScene);
	//worldScene = std::make_unique<Scene>(make_unique<BVHAccel>(prims),lights);

	shared_ptr<Camera> cam;
	Point3f eye(2.5f, 2.5f, 6.0f);
	Point3f look(2.5f, 2.5f, 0.0f);
	Vector3f up(0.0f, 1.0f, 0.0f);
	Transform lookat = LookAt(Vector3f(eye), Vector3f(look), up);
	Transform Camera2World = Inverse(lookat);
	cam = shared_ptr<Camera>(CreatePerspectiveCamera(Camera2World));
	shared_ptr<Sampler> ss = make_unique<StratifiedSampler>(8,8,true,1);
	Bounds2i pixelBounds;
	shared_ptr<Integrator>wSI=make_shared<WhittedIntegrator>(5,cam,ss, pixelBounds);
	shared_ptr<Integrator>pSI = make_shared<PathIntegrator>(5, cam, ss, pixelBounds,1.f);

	pSI->Render(*worldScene);
}
