#include <iostream>
#include <fstream>
#include <time.h>
#include "rtweekend.h"
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

inline double random_double() {
	// Returns a random real in [0,1).
	return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
	// Returns a random real in [min,max).
	return min + (max - min)*random_double();
}

inline static Vector3f random() {
	return Vector3f(random_double(), random_double(), random_double());
}

inline static Vector3f random(double min, double max) {
	return Vector3f(random_double(min, max), random_double(min, max), random_double(min, max));
}

void Cornellbox(shared_ptr<Scene> &scene, Transform &lookat)
{
	//shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	Spectrum purple; purple[0] = 0.35; purple[1] = 0.12; purple[2] = 0.48;
	shared_ptr<Texture<Spectrum>> plasticKd = make_shared<ConstantTexture<Spectrum>>(purple);
	shared_ptr<Texture<Spectrum>> plasticKr = make_shared<ConstantTexture<Spectrum>>(Spectrum(1.f)-purple);
	shared_ptr<Texture<Float>> plasticRoughness = make_shared<ConstantTexture<Float>>(0.1f);
	shared_ptr<Material> plastic = make_shared<PlasticMaterial>(plasticKd,plasticKr,plasticRoughness,true);
	
	shared_ptr<Texture<Spectrum>> mr = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.f,0.0f,1.0f));
	shared_ptr<Texture<Spectrum>> c = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5f,0.5f,0.5f));
	shared_ptr<Material> MR = make_shared<MetalRoughnessMaterial>(c, mr);

	Spectrum eta(0.18f,0.15f,0.81f);
	shared_ptr<Texture<Spectrum>> etaM = make_shared<ConstantTexture<Spectrum>>(eta);
	shared_ptr<Texture<Spectrum>> kM = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.11f,0.11f,0.11f));
	shared_ptr<Texture<Float>> Roughness = make_shared<ConstantTexture<Float>>(0.2f);
	shared_ptr<Texture<Float>> RoughnessU = make_shared<ConstantTexture<Float>>(0.2f);
	shared_ptr<Texture<Float>> RoughnessV = make_shared<ConstantTexture<Float>>(0.2f);
	shared_ptr<Material> metal = make_shared<MetalMaterial>(etaM, kM, Roughness, RoughnessU, RoughnessV,false);

	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor) , sigma);

	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	//Transform *tri_Object2WorldTri=new Transform, *tri_World2ObjectTri=new Transform;
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
	//for (int j = 0; j < nVertices; ++j)
	//	P[j] = p[j];
	//*tri_Object2WorldTri = Translate(Vector3f(2.5f, 1.0f, 1.5f))*Scale(2.0f, 2.0f, 2.0f)*RotateY(90)*(*tri_Object2WorldTri);
	//*tri_World2ObjectTri = Inverse(*tri_Object2WorldTri);
	//std::vector<std::shared_ptr<Shape>> tris = CreateTriangleMesh(tri_Object2WorldTri, tri_World2ObjectTri, false, nTriangles, vertexIndices, nVertices, P, nullptr, nullptr, nullptr, nullptr);
	//for (int i = 0; i < nTriangles; ++i)
	//	prims.push_back(make_shared<GeometricPrimitive>(tris[i], MR, area));
	//delete []P;
	//delete []vertexIndices;

	vector<shared_ptr<Light>>lights;
	float length_Floor = 5;
	const float yPos_AreaLight = length_Floor-0.01;
	float len = 0.5;
	int nTrianglesAreaLight = 2;
	int vertexIndicesAreaLight[6] = { 0,1,2,3,4,5 };
	int nVerticesAreaLight = 6;
	Point3f P_AreaLight[6] = {
			Point3f(-len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,len),
			Point3f(len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,-len) };
	Transform *tri_Object2World_AreaLight=new Transform,*tri_World2Object_AreaLight=new Transform;
	*tri_Object2World_AreaLight = Translate(Vector3f(2.5f, 0.0f, 2.5f))*(*tri_Object2World_AreaLight);
	*tri_World2Object_AreaLight = Inverse(*tri_Object2World_AreaLight);
	std::vector<std::shared_ptr<Shape>> meshAreaLight = CreateTriangleMesh(tri_Object2World_AreaLight, tri_World2Object_AreaLight, false, nTrianglesAreaLight, vertexIndicesAreaLight, nVerticesAreaLight, P_AreaLight, nullptr, nullptr, nullptr, nullptr);
	for (int i = 0; i < nTrianglesAreaLight; ++i)
	{
		area = make_shared<DiffuseAreaLight>(*tri_Object2World_AreaLight, Spectrum(20.f), 1, meshAreaLight[i], false);
		lights.push_back(area);
		prims.push_back(make_shared<GeometricPrimitive>(meshAreaLight[i], m1, area));
	}

	//Transform InfinityLightToWorld;
	//InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	//Spectrum power(1.0f);
	//string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	//shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld, power, 10, hdrFile);
	//lights.push_back(infinityLight);

	shared_ptr<Material> mirror = make_shared<MirrorMaterial>(make_shared<ConstantTexture<Spectrum>>(Spectrum(1.f)));
	Transform *tri_Object2World=new Transform, *tri_World2Object=new Transform;
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
	std::vector<std::shared_ptr<Shape>> meshFloor = CreateTriangleMesh(tri_Object2World, tri_World2Object, false, nTrianglesFloor, vertexIndicesFloor, nVerticesFloor, P_Floor, nullptr, nullptr, nullptr, nullptr);
	Spectrum green(.12f, .45f, .15f);
	Spectrum red(.65f, .05f, .05f);
	Spectrum white(.73f, .73f, .73f);
	shared_ptr<Material> mgreen = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(green), sigma);
	shared_ptr<Material> mred = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(red), sigma);
	shared_ptr<Material> mwhite = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(white), sigma);

	for (int i = 0; i < nTrianglesFloor; ++i)
	{
		if(i<2)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mwhite, nullptr));//底
		else if(i<4)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mwhite, nullptr));//上
		else if (i < 6)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mwhite, nullptr));//后
		else if (i < 8)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mred, nullptr));//右
		else if (i < 10)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mgreen, nullptr));//左
	}

	tri_Object2World = new Transform, tri_World2Object = new Transform;
	*tri_Object2World = Translate(Vector3f(0.7f, 0.f, 2.8f))*RotateY(18)*(*tri_Object2World);
	*tri_World2Object = Inverse(*tri_Object2World);
	const int nTrianglesCube = 12;
	const int nVerticesCube = nTrianglesCube * 3;
	int vertexIndicesCube[nVerticesCube];
	for (int i = 0; i < nVerticesCube; ++i)
		vertexIndicesCube[i] = i;
	Float height = 1.2f;
	height = 1.486f;
	Point3f P_Cube[nVerticesCube] = {
		//底座
		Point3f(0.f,0.f,height),Point3f(height,0.f,height),Point3f(0.f,0.f,0.f),
		Point3f(height,0.f,height),Point3f(height,0.f,0.f),Point3f(0.f,0.f,0.f),
		//天花板
		Point3f(0.f,height,height),Point3f(0.f,height,0.f),Point3f(height,height,height),
		Point3f(height,height,height),Point3f(0.f,height,0.f),Point3f(height,height,0.f),
		//后墙
		Point3f(0.f,0.f,0.f),Point3f(height,0.f,0.f),Point3f(height,height,0.f),
		Point3f(0.f,0.f,0.f),Point3f(height,height,0.f),Point3f(0.f,height,0.f),
		//右墙
		Point3f(0.f,0.f,0.f),Point3f(0.f,height,height),Point3f(0.f,0.f,height),
		Point3f(0.f,0.f,0.f),Point3f(0.f,height,0.f),Point3f(0.f,height,height),
		//左墙
		Point3f(height,0.f,0.f),Point3f(height,height,height),Point3f(height,0.f,height),
		Point3f(height,0.f,0.f),Point3f(height,height,0.f),Point3f(height,height,height),
		//前墙
		Point3f(0.f,0.f,height),Point3f(height,0.f,height),Point3f(height,height,height),
		Point3f(0.f,0.f,height),Point3f(height,height,height),Point3f(0.f,height,height),
	};
	std::vector<std::shared_ptr<Shape>> meshCube = CreateTriangleMesh(tri_Object2World, tri_World2Object, false, nTrianglesCube, vertexIndicesCube, nVerticesCube, P_Cube, nullptr, nullptr, nullptr, nullptr);
	for (int i = 0; i < nTrianglesCube; ++i)
			prims.push_back(make_shared<GeometricPrimitive>(meshCube[i], mwhite, nullptr));

	tri_Object2World = new Transform, tri_World2Object = new Transform;
	*tri_Object2World = Translate(Vector3f(2.7f, 0.f, 1.0f))*RotateY(-15)*(*tri_Object2World);
	*tri_World2Object = Inverse(*tri_Object2World);
	const int nTrianglesRect = 12;
	const int nVerticesRect = nTrianglesRect * 3;
	int vertexIndicesRect[nVerticesRect];
	for (int i = 0; i < nVerticesRect; ++i)
		vertexIndicesRect[i] = i;
	Float width = height * 2;
	Point3f P_Rect[nVerticesRect] = {
		//底座
		Point3f(0.f,0.f,height),Point3f(height,0.f,height),Point3f(0.f,0.f,0.f),
		Point3f(height,0.f,height),Point3f(height,0.f,0.f),Point3f(0.f,0.f,0.f),
		//天花板
		Point3f(0.f,height,height),Point3f(0.f,height,0.f),Point3f(height,height,height),
		Point3f(height,height,height),Point3f(0.f,height,0.f),Point3f(height,height,0.f),
		//后墙
		Point3f(0.f,0.f,0.f),Point3f(height,0.f,0.f),Point3f(height,width,0.f),
		Point3f(0.f,0.f,0.f),Point3f(height,width,0.f),Point3f(0.f,width,0.f),
		//右墙
		Point3f(0.f,0.f,0.f),Point3f(0.f,width,height),Point3f(0.f,0.f,height),
		Point3f(0.f,0.f,0.f),Point3f(0.f,width,0.f),Point3f(0.f,width,height),
		//左墙
		Point3f(height,0.f,0.f),Point3f(height,width,height),Point3f(height,0.f,height),
		Point3f(height,0.f,0.f),Point3f(height,width,0.f),Point3f(height,width,height),
		//前墙
		Point3f(0.f,0.f,height),Point3f(height,0.f,height),Point3f(height,width,height),
		Point3f(0.f,0.f,height),Point3f(height,width,height),Point3f(0.f,width,height),
	};
	std::vector<std::shared_ptr<Shape>> meshRect = CreateTriangleMesh(tri_Object2World, tri_World2Object, false, nTrianglesRect, vertexIndicesRect, nVerticesRect, P_Rect, nullptr, nullptr, nullptr, nullptr);
	for (int i = 0; i < nTrianglesRect; ++i)
	{
		if(i<10)
			prims.push_back(make_shared<GeometricPrimitive>(meshRect[i], mwhite, nullptr));
		else
			prims.push_back(make_shared<GeometricPrimitive>(meshRect[i], mwhite, nullptr));
	}

	//Transform *tri_ObjectWorldModel = new Transform, *tri_WorldObjectModel = new Transform;
	//*tri_ObjectWorldModel = Translate(Vector3f(3.5f, 1.2f, 1.8f))*(*tri_ObjectWorldModel);
	//*tri_WorldObjectModel = Inverse(*tri_ObjectWorldModel);
	//std::shared_ptr<Shape> sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_WorldObjectModel, false, 1.2f);
	//prims.push_back(make_shared<GeometricPrimitive>(sphere, mirror, nullptr));

	scene=make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

	Point3f eye(2.5f, 2.5f, 11.8f);
	Point3f look(2.5f, 2.5f, 0.0f);
	Vector3f up(0.0f, 1.0f, 0.0f);
	lookat = LookAt(Point3f(eye), Point3f(look), up);

}

void MetalRoughSphere(shared_ptr<Scene> &scene,Transform &lookat)
{
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	//shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	//shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor), sigma);
	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	vector<shared_ptr<Light>>lights;

	//Transform InfinityLightToWorld;
	//InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	//Spectrum power(1.0f);
	//string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	//shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld,power,10,hdrFile);
	//lights.push_back(infinityLight);

	Transform *tri_ObjectWorldModel = new Transform;
	*tri_ObjectWorldModel=RotateZ(180)*RotateY(180)*RotateX(-90)*(*tri_ObjectWorldModel);
	//Transform tri_ObjectWorldModel = Translate(Vector3f(0.f,2.f,0.f))*RotateY(-30)*Scale(3, 3, 3);
	ModelLoad ML;
	string metalroughsphere = "C:\\glTF-Sample-Models-master\\2.0\\MetalRoughSpheres\\glTF\\MetalRoughSpheres.gltf";
	//string car = "C:\\VS\\5_ARAP\\project\\bin\\Beetle_ABF.obj";
	ML.loadModel(metalroughsphere, *tri_ObjectWorldModel);
	Spectrum black(1.f, 1.f, 1.f);
	shared_ptr<Material> mblack = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(black), sigma);
	//ML.buildNoTextureModel(tri_ObjectWorldModel, prims, mblack);
	ML.buildTextureModel(*tri_ObjectWorldModel, prims);
	scene = make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

	Vector3f up(0.0f, 1.0f, 0.0f);
	//lookat = LookAt(Vector3f(2.f, 3.f, 3.0f), Vector3f(0.f, 2.f, 0.0f), up);
	lookat = LookAt(Point3f(0.f, 0.f, 50.0f), Point3f(0.f, 0.f, 0.0f), up);
	//lookat = LookAt(Vector3f(0.f, 0.f, -24.0f), Vector3f(0.f, 0.f, 0.0f), up);
}

void DamagedHelmet(shared_ptr<Scene> &scene, Transform &lookat)
{
	Spectrum purple; purple[0] = 0.35; purple[1] = 0.12; purple[2] = 0.48;
	shared_ptr<Texture<Spectrum>> plasticKd = make_shared<ConstantTexture<Spectrum>>(purple);
	shared_ptr<Texture<Spectrum>> plasticKr = make_shared<ConstantTexture<Spectrum>>(Spectrum(1.f) - purple);
	shared_ptr<Texture<Float>> plasticRoughness = make_shared<ConstantTexture<Float>>(0.1f);
	shared_ptr<Material> plastic = make_shared<PlasticMaterial>(plasticKd, plasticKr, plasticRoughness, true);

	shared_ptr<Texture<Spectrum>> mr = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.f, 0.0f, 1.0f));
	shared_ptr<Texture<Spectrum>> c = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5f, 0.5f, 0.5f));
	shared_ptr<Material> MR = make_shared<MetalRoughnessMaterial>(c, mr);

	Spectrum eta(0.18f, 0.15f, 0.81f);
	shared_ptr<Texture<Spectrum>> etaM = make_shared<ConstantTexture<Spectrum>>(eta);
	shared_ptr<Texture<Spectrum>> kM = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.11f, 0.11f, 0.11f));
	shared_ptr<Texture<Float>> Roughness = make_shared<ConstantTexture<Float>>(0.2f);
	shared_ptr<Texture<Float>> RoughnessU = make_shared<ConstantTexture<Float>>(0.2f);
	shared_ptr<Texture<Float>> RoughnessV = make_shared<ConstantTexture<Float>>(0.2f);
	shared_ptr<Material> metal = make_shared<MetalMaterial>(etaM, kM, Roughness, RoughnessU, RoughnessV, false);

	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	Spectrum green(0.0f, 1.0f, 0.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor), sigma);
	shared_ptr<Material> mgreen = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(green), sigma);

	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	vector<shared_ptr<Light>>lights;

	float length_Floor = 5;
	const float yPos_AreaLight = length_Floor - 0.01;
	float len = 0.5;
	int nTrianglesAreaLight = 2;
	int vertexIndicesAreaLight[6] = { 0,1,2,3,4,5 };
	int nVerticesAreaLight = 6;
	Point3f P_AreaLight[6] = {
			Point3f(-len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,len),
			Point3f(len,yPos_AreaLight,len),Point3f(-len,yPos_AreaLight,-len),Point3f(len,yPos_AreaLight,-len) };
	Transform *tri_Object2World_AreaLight = new Transform, *tri_World2Object_AreaLight = new Transform;
	*tri_Object2World_AreaLight = Translate(Vector3f(2.5f, 0.0f, 2.5f))*(*tri_Object2World_AreaLight);
	*tri_World2Object_AreaLight = Inverse(*tri_Object2World_AreaLight);
	std::vector<std::shared_ptr<Shape>> meshAreaLight = CreateTriangleMesh(tri_Object2World_AreaLight, tri_World2Object_AreaLight, false, nTrianglesAreaLight, vertexIndicesAreaLight, nVerticesAreaLight, P_AreaLight, nullptr, nullptr, nullptr, nullptr);
	for (int i = 0; i < nTrianglesAreaLight; ++i)
	{
		area = make_shared<DiffuseAreaLight>(*tri_Object2World_AreaLight, Spectrum(25.f), 1, meshAreaLight[i], false);
		lights.push_back(area);
		prims.push_back(make_shared<GeometricPrimitive>(meshAreaLight[i], m1, area));
	}

	//Transform InfinityLightToWorld;
	//InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	//Spectrum power(1.0f);
	//string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	//shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld, power, 10, hdrFile);
	//lights.push_back(infinityLight);

	Transform *tri_Object2World = new Transform, *tri_World2Object = new Transform;
	const int nTrianglesFloor = 10;
	const int nVerticesFloor = nTrianglesFloor * 3;
	int vertexIndicesFloor[nVerticesFloor];
	for (int i = 0; i < nVerticesFloor; ++i)
		vertexIndicesFloor[i] = i;
	Float height = 0.f;
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
		Point3f(length_Floor,0.f,0.f),Point3f(length_Floor,length_Floor,0.f),Point3f(length_Floor,length_Floor,length_Floor) };
	std::vector<std::shared_ptr<Shape>> meshFloor = CreateTriangleMesh(tri_Object2World, tri_World2Object, false, nTrianglesFloor, vertexIndicesFloor, nVerticesFloor, P_Floor, nullptr, nullptr, nullptr, nullptr);
	Spectrum blue(0.0, 0.0, 1.0);
	Spectrum red(1.0, 0.0, 0.0);
	Spectrum black(0.0, 0.0, 0.0);
	shared_ptr<Material> mblue = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(blue), sigma);
	shared_ptr<Material> mred = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(red), sigma);
	shared_ptr<Material> mblack = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(black), sigma);

	for (int i = 0; i < nTrianglesFloor; ++i)
	{
		if (i < 2)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, nullptr));//底
		else if (i < 4)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, nullptr));//上
		else if (i < 6)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, nullptr));//后
		else if (i < 8)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, nullptr));//右
		else if (i < 10)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, nullptr));//左
	}

	Float ts = 2.5f;
	Transform tri_ObjectWorldModel = Translate(Vector3f(ts, ts, ts))*RotateY(-90)*RotateX(-90);
	ModelLoad ML;
	string damagedhelmet = "C:\\glTF-Sample-Models-master\\2.0\\DamagedHelmet\\glTF\\DamagedHelmet.gltf";
	ML.loadModel(damagedhelmet, tri_ObjectWorldModel);
	ML.buildTextureModel(tri_ObjectWorldModel, prims);

	scene = make_shared<Scene>(make_shared<BVHAccel>(prims), lights);
	Vector3f up(0.0f, 1.0f, 0.0f);
	//lookat = LookAt(Vector3f(2.f, 3.f, 3.0f), Vector3f(0.f, 2.f, 0.0f), up);
	//左视图
	//lookat = LookAt(Vector3f(-0.2f, 0.f, -2.0f), Vector3f(-0.2f, 0.f, 0.0f), up);
	//俯视图
	lookat = LookAt(Point3f(1.f+ ts, 0.7f+ ts, -1.f+ ts), Point3f(-0.5f+ ts, -0.2f+ ts, 0.3f+ ts), up);
}

void sphere(shared_ptr<Scene> &scene, Transform &lookat)
{
	shared_ptr<Texture<Spectrum>> mr = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.f, 0.0f, 1.0f));
	shared_ptr<Texture<Spectrum>> c = make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5f, 0.5f, 0.5f));
	shared_ptr<Material> MR = make_shared<MetalRoughnessMaterial>(c, mr);
	
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	shared_ptr<Material> matte = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor), sigma);


	shared_ptr<Material> mirror = make_shared<MirrorMaterial>(make_shared<ConstantTexture<Spectrum>>(Spectrum(1.0f)));

	vector<shared_ptr<Light>>lights;
	std::vector<std::shared_ptr<Primitive>> prims;

	//Transform InfinityLightToWorld;
	//InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	//Spectrum power(1.0f);
	//string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	//shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld, power, 10, hdrFile);
	//lights.push_back(infinityLight);

	Transform *tri_ObjectWorldModel=new Transform, *tri_World2Object=new Transform;
	*tri_ObjectWorldModel = Translate(Vector3f(-2.f, -2.f, -2.f))*(*tri_ObjectWorldModel);
	*tri_World2Object = Inverse(*tri_ObjectWorldModel);

	std::shared_ptr<Shape> sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_World2Object,false,2.0);
	prims.push_back(make_shared<GeometricPrimitive>(sphere, matte, nullptr));
	scene = make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

	Vector3f up(0.0f, 1.0f, 0.0f);
	lookat = LookAt(Point3f(3.f, 3.f, 3.f), Point3f(0.f, 0.f, 0.f), up);
}

void sphere1(shared_ptr<Scene> &scene, Transform &lookat)
{
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Material> ground_material = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5, 0.5, 0.5)), sigma);

	std::vector<std::shared_ptr<Primitive>> prims;
	Transform *tri_ObjectWorldModel = new Transform, *tri_World2Object = new Transform;
	*tri_ObjectWorldModel = Translate(Vector3f(0.f, -1000.f, 0.f))*(*tri_ObjectWorldModel);
	*tri_World2Object = Inverse(*tri_ObjectWorldModel);

	std::shared_ptr<Shape> sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_World2Object, false, 1000.0);
	prims.push_back(make_shared<GeometricPrimitive>(sphere, ground_material, nullptr));

	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			auto choose_mat = random_double();
			Point3f center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

			if ((center - Point3f(4, 0.2, 0)).Length() > 0.9) {
				shared_ptr<Material> sphere_material;

				if (choose_mat < 0.8) {
					// diffuse
					Spectrum albedo(random() * random());
					sphere_material = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(albedo), sigma);
					
					tri_ObjectWorldModel = new Transform, tri_World2Object = new Transform;
					*tri_ObjectWorldModel = Translate(Vector3f(center))*(*tri_ObjectWorldModel);
					*tri_World2Object = Inverse(*tri_ObjectWorldModel);

					std::shared_ptr<Shape> sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_World2Object, false, 0.2);
					prims.push_back(make_shared<GeometricPrimitive>(sphere, sphere_material, nullptr));
				}
				else{
					// metal
					Spectrum albedo = random(0.5, 1);
					sphere_material = make_shared<MirrorMaterial>(make_shared<ConstantTexture<Spectrum>>(albedo));
					tri_ObjectWorldModel = new Transform, tri_World2Object = new Transform;
					*tri_ObjectWorldModel = Translate(Vector3f(center))*(*tri_ObjectWorldModel);
					*tri_World2Object = Inverse(*tri_ObjectWorldModel);
					std::shared_ptr<Shape> sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_World2Object, false, 0.2);
					prims.push_back(make_shared<GeometricPrimitive>(sphere, sphere_material, nullptr));
				}
			}
		}
	}

	shared_ptr<Material> sphere_material = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(Spectrum(0.4, 0.2, 0.1)), sigma);
	tri_ObjectWorldModel = new Transform, tri_World2Object = new Transform;
	*tri_ObjectWorldModel = Translate(Vector3f(-4, 1, 0))*(*tri_ObjectWorldModel);
	*tri_World2Object = Inverse(*tri_ObjectWorldModel);
	sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_World2Object, false, 1.0);
	prims.push_back(make_shared<GeometricPrimitive>(sphere, sphere_material, nullptr));

	sphere_material = make_shared<MirrorMaterial>(make_shared<ConstantTexture<Spectrum>>(Spectrum(0.7, 0.6, 0.5)));
	tri_ObjectWorldModel = new Transform, tri_World2Object = new Transform;
	*tri_ObjectWorldModel = Translate(Vector3f(4, 1, 0))*(*tri_ObjectWorldModel);
	*tri_World2Object = Inverse(*tri_ObjectWorldModel);
	sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_World2Object, false, 1.0);
	prims.push_back(make_shared<GeometricPrimitive>(sphere, sphere_material, nullptr));

	vector<shared_ptr<Light>>lights;
	scene = make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

	Vector3f up(0.0f, 1.0f, 0.0f);
	lookat = LookAt(Point3f(13, 2, 3), Point3f(0.f, 0.f, 0.f), up);

}

void Sponza(shared_ptr<Scene> &scene, Transform &lookat)
{
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	//shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	//shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor), sigma);
	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	vector<shared_ptr<Light>>lights;

	//Transform InfinityLightToWorld;
	//InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	//Spectrum power(1.0f);
	//string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	//shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld,power,10,hdrFile);
	//lights.push_back(infinityLight);

	Transform *tri_ObjectWorldModel=new Transform;
	//Transform tri_ObjectWorldModel = Translate(Vector3f(0.f,2.f,0.f))*RotateY(-30)*Scale(3, 3, 3);
	ModelLoad ML;
	string Sponza = "C:\\glTF-Sample-Models-master\\2.0\\Sponza\\glTF\\Sponza.gltf";
	//string car = "C:\\VS\\5_ARAP\\project\\bin\\Beetle_ABF.obj";
	ML.loadModel(Sponza, *tri_ObjectWorldModel);
	Spectrum black(1.f, 1.f, 1.f);
	shared_ptr<Material> mblack = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(black), sigma);
	//ML.buildNoTextureModel(tri_ObjectWorldModel, prims, mblack);
	ML.buildTextureModel(*tri_ObjectWorldModel, prims);
	scene = make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

	Vector3f up(0.0f, 1.0f, 0.0f);
	//lookat = LookAt(Vector3f(2.f, 3.f, 3.0f), Vector3f(0.f, 2.f, 0.0f), up);
	lookat = LookAt(Point3f(0.f, 50.f, -10.0f), Point3f(0.f, 0.f, 0.0f), up);
	//lookat = LookAt(Vector3f(0.f, 0.f, -24.0f), Vector3f(0.f, 0.f, 0.0f), up);
}
int main()
{
	std::shared_ptr<Scene>worldScene;
	Transform lookat;
	Cornellbox(worldScene, lookat);
	//MetalRoughSphere(worldScene,lookat);
	//DamagedHelmet(worldScene, lookat);
	//sphere(worldScene, lookat);
	//sphere1(worldScene, lookat);
	//Sponza(worldScene, lookat);
	Transform Camera2World = Inverse(lookat);

	int resolution = 1440;
	resolution = 800;
	int image_width= resolution, image_height= resolution;
	Float focaldistance = 0.01f;
	Float fov = 40.0f;
	//fov = 90.0f;
	shared_ptr<Camera> cam = shared_ptr<Camera>(CreatePerspectiveCamera(Camera2World, focaldistance,fov,image_width, image_height));
	
	int spp = 3;
	shared_ptr<Sampler> ss = make_unique<StratifiedSampler>(spp, spp,true,1);
	Bounds2i pixelBounds;
	//shared_ptr<Integrator>wSI=make_shared<WhittedIntegrator>(5,cam,ss, pixelBounds);
	shared_ptr<Integrator>pSI = make_shared<PathIntegrator>(50, cam, ss, pixelBounds,1.f);

	//clock_t start=clock(), end;
	pSI->Render(*worldScene, image_width, image_height);
	//end = clock();
	//double t = (double)(end - start) / CLOCKS_PER_SEC;
	//cout << t << endl;
	//std::cerr << "\nDone.\n";
	//cout << worldScene.use_count() << endl;
	return 0;
}
