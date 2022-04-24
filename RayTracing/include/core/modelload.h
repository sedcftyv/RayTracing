#pragma once

#ifndef PBRT_CORE_MODELLOAD_H
#define PBRT_CORE_MODELLOAD_H

#include<assimp/Importer.hpp>
#include<assimp/scene.h>
#include<assimp/postprocess.h>
#include<core/pbrt.h>
#include<core/mipmap.h>
#include<core/material.h>
#include<core/texture.h>
#include<core/spectrum.h>
#include<core/imagemap.h>

class ModelLoad {
public:
	vector<shared_ptr<TriangleMesh>> meshes;
	string directory;
	vector<string>texName;
	vector<string>diffTexName;
	vector<string>specTexName;

	void loadModel(string path, const Transform &ObjectToWorld);
	void processNode(aiNode *node, const aiScene *scene, const Transform &ObjectToWorld);
	shared_ptr<TriangleMesh> processMesh(aiMesh *mesh, const aiScene *scene, const Transform &ObjectToWorld);
	void buildNoTextureModel(Transform& tri_Object2World, vector<shared_ptr<Primitive>>&prims, shared_ptr<Material>material);
	void buildTextureModel(Transform& tri_Object2World, vector<shared_ptr<Primitive>>&prims);
};

 inline shared_ptr<Material> getDiffuseMaterial(string filename)
{
	unique_ptr<TextureMapping2D> map = make_unique<UVMapping2D>(1.f,1.f,0.f,0.f);
	ImageWrap wrapMode = ImageWrap::Repeat;
	bool trilerp = false;
	Float maxAniso = 8.f;
	Float scale = 1.f;
	bool gamma = false;
	shared_ptr<Texture<Spectrum>> Kt = make_shared<ImageTexture<RGBSpectrum, Spectrum>>(std::move(map), filename, trilerp, maxAniso, wrapMode, scale, gamma);
	shared_ptr<Texture<Float>> sigmaRed = make_shared<ConstantTexture<Float>>(0.0f);
	shared_ptr<Texture<Float>> bumpMap = make_shared<ConstantTexture<Float>>(0.0f);
	return make_shared<MatteMaterial>(Kt, sigmaRed, bumpMap);
}

 inline shared_ptr<Material> getPlasticMaterial(string diffFilename,string specFilename)
 {
	 unique_ptr<TextureMapping2D> map1 = make_unique<UVMapping2D>(1.f, 1.f, 0.f, 0.f);
	 unique_ptr<TextureMapping2D> map2 = make_unique<UVMapping2D>(1.f, 1.f, 0.f, 0.f);
	 ImageWrap wrapMode = ImageWrap::Repeat;
	 bool trilerp = false;
	 Float maxAniso = 8.f;
	 Float scale = 1.f;
	 bool gamma = false;
	 shared_ptr<Texture<Spectrum>> plasticKd = make_shared<ImageTexture<RGBSpectrum, Spectrum>>(std::move(map1), diffFilename, trilerp, maxAniso, wrapMode, scale, gamma);
	 shared_ptr<Texture<Spectrum>> plasticKr = make_shared<ImageTexture<RGBSpectrum, Spectrum>>(std::move(map2), specFilename, trilerp, maxAniso, wrapMode, scale, gamma);
	 shared_ptr<Texture<Float>> plasticRoughness = make_shared<ConstantTexture<Float>>(1.f);
	 shared_ptr<Texture<Float>> bumpMap = make_shared<ConstantTexture<Float>>(0.0f);
	 return make_shared<PlasticMaterial>(plasticKd, plasticKr,plasticRoughness, bumpMap,true);
 }
#endif