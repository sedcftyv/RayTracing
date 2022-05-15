﻿#include <iostream>
#include <fstream>
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
	Spectrum green(0.0f, 1.0f,0.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor) , sigma);
	shared_ptr<Material> mgreen = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(green), sigma);
	//shared_ptr<Material> mirror = make_shared<MirrorMaterial>(make_shared<ConstantTexture<Spectrum>>(green));

	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	Transform *tri_Object2WorldTri=new Transform, *tri_World2ObjectTri=new Transform;
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
	*tri_Object2WorldTri = Translate(Vector3f(2.5f, 1.0f, 1.5f))*Scale(2.0f, 2.0f, 2.0f)*RotateY(90)*(*tri_Object2WorldTri);
	*tri_World2ObjectTri = Inverse(*tri_Object2WorldTri);
	std::vector<std::shared_ptr<Shape>> tris = CreateTriangleMesh(tri_Object2WorldTri, tri_World2ObjectTri, false, nTriangles, vertexIndices, nVertices, P, nullptr, nullptr, nullptr, nullptr);
	//for (int i = 0; i < nTriangles; ++i)
	//	prims.push_back(make_shared<GeometricPrimitive>(tris[i], MR, area));
	delete []P;
	delete []vertexIndices;

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

	Transform *tri_Object2World=new Transform, *tri_World2Object=new Transform;
	const int nTrianglesFloor = 10;
	const int nVerticesFloor = nTrianglesFloor*3;
	int vertexIndicesFloor[nVerticesFloor];
	for (int i = 0; i < nVerticesFloor; ++i)
		vertexIndicesFloor[i] = i;
	Float height=0.f;
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
	Spectrum blue(0.0, 0.0, 1.0);
	Spectrum red(1.0, 0.0, 0.0);
	Spectrum black(0.0, 0.0, 0.0);
	shared_ptr<Material> mblue = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(blue), sigma);
	shared_ptr<Material> mred = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(red), sigma);
	shared_ptr<Material> mblack = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(black), sigma);

	shared_ptr<Material> mirror = make_shared<MirrorMaterial>(make_shared<ConstantTexture<Spectrum>>(Spectrum(1.0f)));
	Transform *tri_ObjectWorldModel = new Transform, *tri_WorldObjectModel = new Transform;
	*tri_ObjectWorldModel = Translate(Vector3f(1.5f, 1.5f, 1.5f))*(*tri_ObjectWorldModel);
	*tri_WorldObjectModel = Inverse(*tri_ObjectWorldModel);
	std::shared_ptr<Shape> sphere = make_shared<Sphere>(tri_ObjectWorldModel, tri_WorldObjectModel, false, 1.25f);
	prims.push_back(make_shared<GeometricPrimitive>(sphere, mirror, nullptr));


	area.reset();
	for (int i = 0; i < nTrianglesFloor; ++i)
	{
		if(i<2)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, area));//底
		else if(i<4)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, area));//上
		else if (i < 6)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], m1, area));//后
		else if (i < 8)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mred, area));//右
		else if (i < 10)
			prims.push_back(make_shared<GeometricPrimitive>(meshFloor[i], mblue, area));//左
	}

	scene=make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

	Point3f eye(2.5f, 2.5f, 6.0f);
	Point3f look(2.5f, 2.5f, 0.0f);
	Vector3f up(0.0f, 1.0f, 0.0f);
	lookat = LookAt(Point3f(eye), Point3f(look), up);

}

void MetalRoughSphere(shared_ptr<Scene> &scene,Transform &lookat)
{
	Spectrum LightColor(1.0f, 1.0f, 1.0f);
	shared_ptr<Texture<Float>> sigma = make_shared<ConstantTexture<Float>>(0.0f);
	//shared_ptr<Texture<Float>> bump = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Material> m1 = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(LightColor), sigma);
	shared_ptr<AreaLight>area;
	std::vector<std::shared_ptr<Primitive>> prims;
	vector<shared_ptr<Light>>lights;

	//Transform InfinityLightToWorld;
	//InfinityLightToWorld = RotateX(90)*InfinityLightToWorld;
	//Spectrum power(1.0f);
	//string hdrFile = "C:\\glTF-Sample-Models-master\\Walk_Of_Fame\\Mans_Outside_2k.hdr";
	//shared_ptr<Light>infinityLight = make_shared<InfiniteAreaLight>(InfinityLightToWorld,power,10,hdrFile);
	//lights.push_back(infinityLight);

	Transform tri_ObjectWorldModel = RotateZ(180)*RotateY(180)*RotateX(-90);
	//Transform tri_ObjectWorldModel = Translate(Vector3f(0.f,2.f,0.f))*RotateY(-30)*Scale(3, 3, 3);
	ModelLoad ML;
	string metalroughsphere = "C:\\glTF-Sample-Models-master\\2.0\\MetalRoughSpheres\\glTF\\MetalRoughSpheres.gltf";
	string car = "C:\\VS\\5_ARAP\\project\\bin\\Beetle_ABF.obj";
	ML.loadModel(metalroughsphere, tri_ObjectWorldModel);
	Spectrum black(1.f, 1.f, 1.f);
	shared_ptr<Material> mblack = make_shared<MatteMaterial>(make_shared<ConstantTexture<Spectrum>>(black), sigma);
	//ML.buildNoTextureModel(tri_ObjectWorldModel, prims, mblack);
	ML.buildTextureModel(tri_ObjectWorldModel, prims);
	scene = make_shared<Scene>(make_shared<BVHAccel>(prims), lights);

	Vector3f up(0.0f, 1.0f, 0.0f);
	//lookat = LookAt(Vector3f(2.f, 3.f, 3.0f), Vector3f(0.f, 2.f, 0.0f), up);
	lookat = LookAt(Point3f(0.f, 0.f, 24.0f), Point3f(0.f, 0.f, 0.0f), up);
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

int main()
{
	std::shared_ptr<Scene>worldScene;
	Transform lookat;
	//Cornellbox(worldScene, lookat);
	//MetalRoughSphere(worldScene,lookat);
	DamagedHelmet(worldScene, lookat);
	//sphere(worldScene, lookat);
	Transform Camera2World = Inverse(lookat);

	int resolution = 1440;
	//resolution = 800;
	int image_width= resolution, image_height= resolution;
	shared_ptr<Camera> cam = shared_ptr<Camera>(CreatePerspectiveCamera(Camera2World, image_width, image_height));
	
	int spp = 16;
	shared_ptr<Sampler> ss = make_unique<StratifiedSampler>(spp, spp,true,1);
	Bounds2i pixelBounds;
	shared_ptr<Integrator>wSI=make_shared<WhittedIntegrator>(5,cam,ss, pixelBounds);
	shared_ptr<Integrator>pSI = make_shared<PathIntegrator>(4, cam, ss, pixelBounds,1.f);

	pSI->Render(*worldScene, image_width, image_height);
	//cout << worldScene.use_count() << endl;
	return 0;
}
