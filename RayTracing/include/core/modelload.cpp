#include<core/modelload.h>
#include<core/geometry.h>
#include<core/triangle.h>
#include<core/primitive.h>
#include<core/imagemap.h>

void ModelLoad::processNode(aiNode *node, const aiScene *scene, const Transform &ObjectToWorld)
{
	for (unsigned int i = 0; i < node->mNumMeshes; ++i)
	{
		aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
		meshes.push_back(processMesh(mesh, scene, ObjectToWorld));

		aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
		for (unsigned int i = 0; i < material->GetTextureCount(aiTextureType_DIFFUSE); ++i)
		{
			aiString str;
			material->GetTexture(aiTextureType_DIFFUSE, i, &str);
			string filename = str.C_Str();
			texName.push_back(filename);
		}
	}
	for (unsigned int i = 0; i < node->mNumChildren; ++i)
		processNode(node->mChildren[i], scene, ObjectToWorld);

}

void ModelLoad::loadModel(std::string path, const Transform &ObjectToWorld)
{
	Assimp::Importer import;
	const aiScene *scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs);
	if (!scene || scene->mFlags&AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
		return;
	directory = path.substr(0, path.find_last_of('/'));
	processNode(scene->mRootNode, scene, ObjectToWorld);
}

shared_ptr<TriangleMesh> ModelLoad::processMesh(aiMesh *mesh, const aiScene *scene, const Transform &ObjectToWorld)
{
	size_t nVertices = mesh->mNumVertices;
	size_t nTriangles = mesh->mNumFaces;
	int *vertexIndices = new int[nTriangles * 3];
	Point3f *P = new Point3f[nVertices];
	Vector3f *S = nullptr;
	Normal3f *N = new Normal3f[nVertices];
	Point2f *uv = new Point2f[nVertices];
	int *faceIndices = nullptr;
	for (unsigned int i = 0; i < mesh->mNumVertices; ++i)
	{
		P[i].x = mesh->mVertices[i].x;
		P[i].y = mesh->mVertices[i].y;
		P[i].z = mesh->mVertices[i].z;
		if (mesh->HasNormals())
		{
			N[i].x = mesh->mNormals[i].x;
			N[i].y = mesh->mNormals[i].y;
			N[i].z = mesh->mNormals[i].z;
		}

		if (mesh->mTextureCoords[0])
		{
			uv[i].x = mesh->mTextureCoords[0][i].x;
			uv[i].y = mesh->mTextureCoords[0][i].y;
		}
	}
	for (unsigned int i = 0; i < mesh->mNumFaces; ++i)
	{
		aiFace face = mesh->mFaces[i];
		for (unsigned int j = 0; j < face.mNumIndices; ++j)
			vertexIndices[3 * i + j] = face.mIndices[j];
	}
	if (!mesh->HasNormals())
	{
		delete[] N;
		N = nullptr;
	}
	if (!mesh->mTextureCoords[0])
	{
		delete[] uv;
		uv = nullptr;
	}
	shared_ptr<TriangleMesh> a = make_shared<TriangleMesh>(ObjectToWorld,nTriangles,vertexIndices,nVertices,P,S,N,uv,faceIndices);
	
	delete[] P;
	delete[] S;
	delete[] N;
	delete[] uv;
	return a;
}

void ModelLoad::buildNoTextureModel(Transform& tri_Object2World, vector<shared_ptr<Primitive>>&prims, shared_ptr<Material>material)
{
	vector<shared_ptr<Shape>>trisObj;
	Transform tri_World2Object = Inverse(tri_Object2World);
	for(int i=0;i<meshes.size();++i)
		for (int j = 0; j < meshes[i]->nTriangles; ++j)
		{
			shared_ptr<TriangleMesh> meshPtr = meshes[i];
			trisObj.push_back(make_shared<Triangle>(&tri_Object2World, &tri_World2Object, false, meshPtr, j));
		}
	for (int i = 0; i < trisObj.size(); ++i)
	{
		prims.push_back(make_shared<GeometricPrimitive>(trisObj[i],material,nullptr));
	}
}

//inline shared_ptr<Material> ModelLoad::getDiffuseMaterial(string filename)
//{
//	unique_ptr<TextureMapping2D> Map = make_unique<UVMapping2D>(1.f, 1.f, 0.f, 0.f);
//	ImageWrap wrapMode = ImageWrap::Repeat;
//	bool trilerp = false;
//	Float maxAniso = 8.f;
//	Float scale = 1.f;
//	bool gamma = false;
//	shared_ptr<Texture<Spectrum>> Kt = make_shared<ImageTexture<RGBSpectrum, Spectrum>>(std::move(Map), filename, trilerp, maxAniso, wrapMode, scale, gamma);
//	shared_ptr<Texture<Float>> sigmaRed = make_shared<ConstantTexture<Float>>(0.0f);
//	shared_ptr<Texture<Float>> bumpMap = make_shared<ConstantTexture<Float>>(0.0f);
//	return make_shared<MatteMaterial>(Kt, sigmaRed, bumpMap);
//}

//void ModelLoad::buildTextureModel(Transform& tri_Object2World, vector<shared_ptr<Primitive>>&prims)
//{
//	vector<shared_ptr<Shape>> trisObj;
//	Transform tri_World2Object = Inverse(tri_Object2World);
//	for (int i = 0; i < meshes.size(); ++i)
//	{
//		string filename = directory + "/" + texName[i];
//		shared_ptr<Material> material = getDiffuseMaterial(filename);
//		for (int j = 0; j < meshes[i]->nTriangles; ++j)
//		{
//			shared_ptr<TriangleMesh> meshPtr = meshes[i];
//			prims.push_back(make_shared<GeometricPrimitive>(make_shared<Triangle>(&tri_Object2World, &tri_World2Object, false, meshPtr, j), material, nullptr));
//		}
//	}
//}