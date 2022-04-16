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
#include "core/infinite.h"
#include "core/modelload.h"
using std::cout;
using std::endl;
using std::make_shared;
using std::unique_ptr;
using std::make_unique;
using std::shared_ptr;
using std::vector;

void Cornellbox(shared_ptr<Scene> &scene)
{
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	Spectrum green(0.0f, 1.0f,0.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor) , sigma, bump);
	shared_ptr<Material> mgreen = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(green), sigma, bump);
	shared_ptr<Material> mirror = make_shared<MirrorMaterial>(make_shared<ConstantTexture<Spectrum>>(green), bump);

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

	//Transform InfinityLightToWorld;
	//InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	//Spectrum power(1.0f);
	//string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	//shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld, power, 10, hdrFile);
	//lights.push_back(infinityLight);


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

void MetalRoughSphere(shared_ptr<Scene> &scene)
{
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor), sigma, bump);
	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	vector<shared_ptr<Light>>lights;

	Transform InfinityLightToWorld;
	InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	Spectrum power(1.0f);
	string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld,power,10,hdrFile);
	lights.push_back(infinityLight);

	//const float yPos_AreaLight = -16;
	//float len = 10;
	//int nTrianglesAreaLight = 2;
	//int vertexIndicesAreaLight[6] = { 0,1,2,3,4,5 };
	//int nVerticesAreaLight = 6;
	//Point3f P_AreaLight[6] = {
	//		Point3f(-len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,len),
	//		Point3f(len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,-len) };
	//Transform tri_Object2World_AreaLight;
	//tri_Object2World_AreaLight = Translate(Vector3f(2.5f, 0.0f, 2.5f))*tri_Object2World_AreaLight;
	//Transform tri_World2Object_AreaLight = Inverse(tri_Object2World_AreaLight);
	//std::vector<std::shared_ptr<Shape>> meshAreaLight = CreateTriangleMesh(&tri_Object2World_AreaLight, &tri_World2Object_AreaLight, false, nTrianglesAreaLight, vertexIndicesAreaLight, nVerticesAreaLight, P_AreaLight, nullptr, nullptr, nullptr, nullptr);
	//for (int i = 0; i < nTrianglesAreaLight; ++i)
	//{
	//	area = make_shared<DiffuseAreaLight>(tri_Object2World_AreaLight, Spectrum(7.f), 5, meshAreaLight[i], true);
	//	lights.push_back(area);
	//	prims.push_back(make_shared<GeometricPrimitive>(meshAreaLight[i], m1, area));
	//}
	Transform tri_ObjectWorldModel = RotateX(90);
	//Transform tri_ObjectWorldModel = Translate(Vector3f(0.f,2.f,0.f))*RotateY(-30)*Scale(3, 3, 3);
	ModelLoad ML;
	string metalroughsphere = "C:\\glTF-Sample-Models-master\\2.0\\MetalRoughSpheres\\glTF\\MetalRoughSpheres.gltf";
	string car = "C:\\VS\\5_ARAP\\project\\bin\\Beetle_ABF.obj";
	ML.loadModel(metalroughsphere, tri_ObjectWorldModel);
	Spectrum black(1.f, 1.f, 1.f);
	shared_ptr<Material> mblack = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(black), sigma, bump);
	ML.buildNoTextureModel(tri_ObjectWorldModel, prims, mblack);
	scene = make_shared<Scene>(make_shared<BVHAccel>(prims), lights);
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

	std::shared_ptr<Scene>worldScene;
	//Cornellbox(worldScene);
	MetalRoughSphere(worldScene);
	//worldScene = std::make_unique<Scene>(make_unique<BVHAccel>(prims),lights);

	shared_ptr<Camera> cam;
	Point3f eye(2.5f, 2.5f, 6.0f);
	Point3f look(2.5f, 2.5f, 0.0f);
	Vector3f up(0.0f, 1.0f, 0.0f);
	//Transform lookat = LookAt(Vector3f(eye), Vector3f(look), up);
	Transform lookat = LookAt(Vector3f(0.f, 0.f, 23.0f), Vector3f(0.f, 0.f, 0.0f), up);
	//Transform lookat = LookAt(Vector3f(2.f, 3.f, 3.0f), Vector3f(0.f, 2.f, 0.0f), up);
	Transform Camera2World = Inverse(lookat);

	cam = shared_ptr<Camera>(CreatePerspectiveCamera(Camera2World));
	shared_ptr<Sampler> ss = make_unique<StratifiedSampler>(8,8,true,1);
	Bounds2i pixelBounds;
	shared_ptr<Integrator>wSI=make_shared<WhittedIntegrator>(5,cam,ss, pixelBounds);
	shared_ptr<Integrator>pSI = make_shared<PathIntegrator>(10, cam, ss, pixelBounds,1.f);

	pSI->Render(*worldScene);
	return 0;
}
