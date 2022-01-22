#include <iostream>
#include <fstream>
#include "rtweekend.h"
#include "color.h"
#include "vec3.h"
#include "ray.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "moving_sphere.h"
#include "aarect.h"
#include "box.h"
//#include "constant_medium.h"
#include "bvh.h"
#include "pdf.h"

using std::cout;
using std::endl;


hittable_list two_spheres() {
	hittable_list objects;

	auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));

	objects.add(make_shared<sphere>(point3(0, -10, 0), 10, make_shared<lambertian>(checker)));
	objects.add(make_shared<sphere>(point3(0, 10, 0), 10, make_shared<lambertian>(checker)));

	return objects;
}

hittable_list two_perlin_spheres() {
	hittable_list objects;

	auto pertext = make_shared<noise_texture>(4);

	objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
	objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

	return objects;
}

hittable_list earth() {
	auto earth_texture = make_shared<image_texture>("earthmap.jpg");
	auto earth_surface = make_shared<lambertian>(earth_texture);
	auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);

	return hittable_list(globe);
}

hittable_list simple_light() {
	hittable_list objects;

	auto pertext = make_shared<noise_texture>(4);
	objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
	objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

	auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
	objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));

	return objects;
}

hittable_list cornell_box() {
	hittable_list objects;

	auto red = make_shared<lambertian>(color(.65, .05, .05));
	auto white = make_shared<lambertian>(color(.73, .73, .73));
	auto green = make_shared<lambertian>(color(.12, .45, .15));
	auto light = make_shared<diffuse_light>(color(15, 15, 15));

	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
	objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

	//shared_ptr<material> aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
	shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, vec3(265, 0, 295));
	objects.add(box1);

	auto glass = make_shared<dielectric>(1.5);
	objects.add(make_shared<sphere>(point3(190, 90, 190), 90, glass));

	//shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
	//box2 = make_shared<rotate_y>(box2, -18);
	//box2 = make_shared<translate>(box2, vec3(130, 0, 65));
	//objects.add(box2);

	return objects;
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
	auto aspect_ratio = 16.0 / 9.0;
	int image_width = 600;
	int samples_per_pixel = 100;
	int max_depth = 50;

	// World
	hittable_list world;
	point3 lookfrom;
	point3 lookat;
	auto vfov = 40.0;
	auto aperture = 0.0;
	color background(0, 0, 0);

	switch (6) {

	case 2:
		world = two_spheres();
		background = color(0.70, 0.80, 1.00);
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;

	case 3:
		world = two_perlin_spheres();
		background = color(0.70, 0.80, 1.00);
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;

	case 4:
		world = earth();
		background = color(0.70, 0.80, 1.00);
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;

	case 5:
		world = simple_light();
		samples_per_pixel = 400;
		background = color(0, 0, 0);
		lookfrom = point3(26, 3, 6);
		lookat = point3(0, 2, 0);
		vfov = 20.0;
		break;

	case 6:
		world = cornell_box();
		aspect_ratio = 1.0;
		image_width = 600;
		samples_per_pixel = 1000;
		background = color(0, 0, 0);
		lookfrom = point3(278, 278, -800);
		lookat = point3(278, 278, 0);
		vfov = 40.0;
		aperture = 0;
		break;
	}
	// Camera
	vec3 vup(0, 1, 0);
	auto dist_to_focus = 10.0;
	camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

	// Render
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	fout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

	//shared_ptr<hittable> lights =
	//	  make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>());
	//	//make_shared<sphere>(point3(190, 90, 190), 90, shared_ptr<material>());

	hittable_list lights;
	lights.add(make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr <material>(0)));
	lights.add(make_shared<sphere>(point3(190, 90, 190), 90, shared_ptr < material>(0)));

	shared_ptr<hittable>lt = make_shared<hittable_list>(lights);

	for (int j = image_height - 1; j >= 0; --j) {
		std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
		for (int i = 0; i < image_width; ++i) {
			color pixel_color(0, 0, 0);
			for (int s = 0; s < samples_per_pixel; ++s) {
				auto u = (i + random_double()) / (image_width - 1);
				auto v = (j + random_double()) / (image_height - 1);
				ray r = cam.get_ray(u, v);
				color tmp = ray_color(r, background, world, lt,max_depth);
				//cout << tmp << endl;
				pixel_color += tmp;
			}
			//cout << pixel_color << endl;
			write_color(fout, pixel_color, samples_per_pixel);
		}
	}
	
	fout.close();
	std::cerr << "\nDone.\n";
}
