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

//hittable_list two_spheres() {
//	hittable_list objects;
//
//	auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
//
//	objects.add(make_shared<sphere>(Vector3f(0, -10, 0), 10, make_shared<lambertian>(checker)));
//	objects.add(make_shared<sphere>(Vector3f(0, 10, 0), 10, make_shared<lambertian>(checker)));
//
//	return objects;
//}

//hittable_list two_perlin_spheres() {
//	hittable_list objects;
//
//	auto pertext = make_shared<noise_texture>(4);
//
//	objects.add(make_shared<sphere>(Vector3f(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
//	objects.add(make_shared<sphere>(Vector3f(0, 2, 0), 2, make_shared<lambertian>(pertext)));
//
//	return objects;
//}

//hittable_list earth() {
//	auto earth_texture = make_shared<image_texture>("earthmap.jpg");
//	auto earth_surface = make_shared<lambertian>(earth_texture);
//	auto globe = make_shared<sphere>(Vector3f(0, 0, 0), 2, earth_surface);
//
//	return hittable_list(globe);
//}

//hittable_list simple_light() {
//	hittable_list objects;
//
//	auto pertext = make_shared<noise_texture>(4);
//	objects.add(make_shared<sphere>(Vector3f(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
//	objects.add(make_shared<sphere>(Vector3f(0, 2, 0), 2, make_shared<lambertian>(pertext)));
//
//	auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
//	objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));
//
//	return objects;
//}

//hittable_list cornell_box() {
//	hittable_list objects;
//
//	auto red = make_shared<lambertian>(color(.65, .05, .05));
//	auto white = make_shared<lambertian>(color(.73, .73, .73));
//	auto green = make_shared<lambertian>(color(.12, .45, .15));
//	auto light = make_shared<diffuse_light>(color(15, 15, 15));
//
//	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
//	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
//	objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
//	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
//	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
//	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
//
//	//shared_ptr<material> aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
//
//	//shared_ptr<hittable> box1 = make_shared<box>(Vector3f(0, 0, 0), Vector3f(165, 330, 165), white);
//	//box1 = make_shared<rotate_y>(box1, 15);
//	//box1 = make_shared<translate>(box1, Vector3f(265, 0, 295));
//	//objects.add(box1);
//	//auto glass = make_shared<dielectric>(1.5);
//	//objects.add(make_shared<sphere>(Vector3f(190, 90, 190), 90, glass));
//
//	Transform t1 = Translate(Vector3f(390, 90, 190));
//	Point3f c1 = t1(Point3f(0, 0, 0));
//	Transform t2 = Translate(Vector3f(400, 150, 295));
//	Point3f c2 = t2(Point3f(0, 0, 0));
//	auto glass = make_shared<dielectric>(1.5);
//	auto aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
//	objects.add(make_shared<sphere>(Vector3f(c1), 90, glass));
//	//objects.add(make_shared<sphere>(Vector3f(c2), 90, aluminum));
//
//	//shared_ptr<hittable> box2 = make_shared<box>(Vector3f(0, 0, 0), Vector3f(165, 165, 165), white);
//	//box2 = make_shared<rotate_y>(box2, -18);
//	//box2 = make_shared<translate>(box2, Vector3f(130, 0, 65));
//	//objects.add(box2);
//
//	return objects;
//}

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

	//switch (6) {

	//case 2:
	//	world = two_spheres();
	//	background = color(0.70, 0.80, 1.00);
	//	lookfrom = Vector3f(13, 2, 3);
	//	lookat = Vector3f(0, 0, 0);
	//	vfov = 20.0;
	//	break;

	//case 3:
	//	world = two_perlin_spheres();
	//	background = color(0.70, 0.80, 1.00);
	//	lookfrom = Vector3f(13, 2, 3);
	//	lookat = Vector3f(0, 0, 0);
	//	vfov = 20.0;
	//	break;

	//case 4:
	//	world = earth();
	//	background = color(0.70, 0.80, 1.00);
	//	lookfrom = Vector3f(13, 2, 3);
	//	lookat = Vector3f(0, 0, 0);
	//	vfov = 20.0;
	//	break;

	//case 5:
	//	world = simple_light();
	//	samples_per_pixel = 400;
	//	background = color(0, 0, 0);
	//	lookfrom = Vector3f(26, 3, 6);
	//	lookat = Vector3f(0, 2, 0);
	//	vfov = 20.0;
	//	break;

	//case 6:
	//	world = cornell_box();
	//	aspect_ratio = 1.0;
	//	image_width = 300;
	//	samples_per_pixel = 50;
	//	background = color(0, 0, 0);
	//	lookfrom = Vector3f(278, 278, -800);
	//	lookat = Vector3f(278, 278, 0);
	//	vfov = 40.0;
	//	aperture = 0;
	//	break;
	//}
	// Camera
	//Vector3f vup(0, 1, 0);
	//auto dist_to_focus = 10.0;
	//camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

	// Render
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	fout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

	//shared_ptr<hittable> lights =
	//	  make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>());
	//	//make_shared<sphere>(Vector3f(190, 90, 190), 90, shared_ptr<material>());

	//hittable_list lights;
	//lights.add(make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr <material>(0)));
	//lights.add(make_shared<sphere>(Vector3f(190, 90, 190), 90, shared_ptr < material>(0)));

	//shared_ptr<hittable>lt = make_shared<hittable_list>(lights);

	//sphere
	Transform sphereT_Object2World, sphereT_World2Object;
	Shape *s = new Sphere(&sphereT_Object2World, &sphereT_World2Object, false,1.0);

	//triangle
	Transform tri_Object2World, tri_World2Object;

	std::vector<Point3f>p;
	std::vector<int>vi;
	string path("C:\\VS\\5_ARAP\\project\\bin\\Balls.obj");
	objload(path, p, vi);
	int nTriangles = vi.size()/3;	int nVertices = p.size();	int *vertexIndices = new int[vi.size()];	Point3f* P = new Point3f[p.size()];	for(int i = 0; i < vi.size(); ++i)		vertexIndices[i] = vi[i];	for(int i = 0; i < nVertices; ++i)		P[i] = p[i];	//int vertexIndices[6] = { 0 ,1 ,2 ,3 ,4 ,5 };	//Point3f P[6] = {Point3f(-1.0 ,1.0 ,0) , Point3f(-1.0 , -1.0 ,0) , Point3f(0.0 , 1.0 , 0) ,Point3f(1.0 , 1.0 , 0) , Point3f(1.0 , -1.0 ,0) , Point3f(0.0 , -1.0 ,0)};
	std::vector<std::shared_ptr<Shape>> tris = CreateTriangleMesh(&tri_Object2World, &tri_World2Object,false,nTriangles,vertexIndices,nVertices,P,nullptr,nullptr,nullptr,nullptr);
	std::vector<std::shared_ptr<Primitive>> prims;
	for (int i = 0; i < nTriangles; ++i)
		prims.push_back(make_shared<GeometricPrimitive>(tris[i]));
	Aggregate *agg;
	agg = new BVHAccel(prims);
	
	//Camera* cam;
	//Point3f eye(0.f, 0.f, -20.f), look(0.0, 0.0, 0.0f);
	//Vector3f up(0.0f, 1.0f, 0.0f);
	//Transform lookat = LookAt(Vector3f(eye), Vector3f(look), up);
	//Transform Camera2World = Inverse(lookat);
	//cam = CreatePerspectiveCamera(Camera2World);
	//Vector3f Light(1.0, 1.0, 1.0);
	////Vector3f lower_left_corner(-2.0, -2.0, -1.0);
	////Vector3f horizontal(4.0, 0.0, 0.0);
	////Vector3f vertical(0.0, 4.0,0.0);
	////Point3f origin(0.0, 0.0, -3.0);
	//for (int j = 0; j < image_height; j++) {
	//	for (int i = 0; i < image_width; i++) {
	//		//float v = float(i + 0.5) / float(image_height); //random ( )
	//		//float u = float(j + 0.5) / float(image_width); // 0.5
	//		//cout << (lower_left_corner + u * horizontal + v * vertical) << endl;
	//		//Ray r(Vector3f(origin), (lower_left_corner + u * horizontal + v* vertical) - Vector3f(origin));
	//		Ray r;
	//		CameraSample cs; 
	//		cs.pFilm = Point2f(i + random_double(), j + random_double());
	//		cam->GenerateRay(cs, &r);

	//		float tHit;
	//		SurfaceInteraction * isect = nullptr;
	//		Spectrum colObj(0.0);
	//		colObj[0] = 1.0f;
	//		colObj[1] = 1.0f;
	//		bool hit = agg->Intersect(r, isect);
	//		if(hit)
	//		{
	//			float Li = Dot(Light, isect->n);
	//			colObj[1] = std::abs(Li);
	//			//cout << r.d[1] << endl;
	//		}
	//		write_color(fout, colObj, 1);
	//	}
	//}


	////for (int j = image_height - 1; j >= 0; --j) {
	////	std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
	////	for (int i = 0; i < image_width; ++i) {
	////		color pixel_color(0, 0, 0);
	////		for (int s = 0; s < samples_per_pixel; ++s) {
	////			auto u = (i + random_double()) / (image_width - 1);
	////			auto v = (j + random_double()) / (image_height - 1);
	////			ray r = cam.get_ray(u, v);
	////			color tmp = ray_color(r, background, world, lt,max_depth);
	////			//cout << tmp << endl;
	////			pixel_color += tmp;
	////		}
	////		//cout << pixel_color << endl;
	////		write_color(fout, pixel_color, samples_per_pixel);
	////	}
	////}
	//fout.close();
	//std::cerr << "\nDone.\n";
}
