#pragma once

#include "glm/glm.hpp"
#include "Types.h"
#include "Scene.h"

#include "string"
#include "atlimage.h"

using namespace glm;

class CTracer
{
public:
  SRay MakeRay(glm::uvec2 pixelPos);  // Create ray for specified pixel
  glm::vec3 TraceRay(SRay ray); // Trace ray, compute its color
  void RenderImage(int xRes, int yRes, double blackhole_mass, int k, float x_cam, float y_cam, float z_cam);
  void SaveImageToFile(std::string fileName);
  CImage* LoadImageFromFile(std::string fileName);

  void CameraSetUp(glm::vec3 pos, double xRes, double yRes);
 bool RaySphereIntersection(glm::vec3 ray_pos, glm::vec3 ray_dir, glm::vec3 sphere_pos, double r);
  bool RayDiskIntersection(glm::vec3 ray_pos1, glm::vec3 ray_pos2, double r);

  vec3 DiskColor(vec3 x, vec3 y, float &alpha_disk);
  vec3 CTracer::StarsColor(vec3 dir);
  vector<glm::vec3> stars_pixels;
  int background_height, background_width;

public:
  SCamera m_camera;
  CScene* m_pScene;
  Sphere blackhole;
  Disk disk;
};


/*
#include "Tracer.h"
#include <cmath>
#include <climits>
#include <iostream>
#include <omp.h>
#include "glm\gtx\perpendicular.hpp"
#include "glm\gtx\rotate_vector.hpp"
#include <vector>

using namespace glm;
using namespace std;

const float Gravity = 0.06674f;
const float Mass = 8.57f;
const float Radius = 12.7f;
const double M_PI = 3.14159265f;

SRay CTracer::MakeRay(uvec2 pixelPos)
{
	SRay ray;
	ray.m_start = m_camera.m_pos;
	ray.m_dir = normalize(m_camera.m_forward + m_camera.m_right * ((pixelPos.x + 0.5f) / m_camera.m_resolution.x - 0.5f)
		+ m_camera.m_up * ((pixelPos.y + 0.5f) / m_camera.m_resolution.y - 0.5f));
	return ray;
}

vec3 CTracer::DiskColor(vec3 x, vec3 y){
	vec3 intersect = x + (y - x) * fabs(x.z / (y.z - x.z));
	float x_1 = (intersect.x / disk.r_max + 1) / 2, y_1 = (intersect.y / disk.r_max + 1) / 2;
	//cout << x_1 << ' ' << y_1 << endl;
	x_1 = fabs(x_1 * disk.width - 0.5f), y_1 = fabs(y_1 * disk.height - 0.5f);
	if (x_1 >= disk.width)
		x_1 = disk.width - 1;
	if (y_1 >= disk.height)
		y_1 = disk.height - 1;

	if (int(disk.width * round(x_1 - 0.5) + round(y_1 - 0.5)) >= disk.width * disk.height || (int)(disk.width * round(x_1 - 0.5) + round(y_1 - 0.5)) <= 0) {
		return disk.color[disk.width * disk.height - 1];
	}

	//cout << int(disk.width * round(x_1 - 0.5) + round(y_1 - 0.5)) << endl;
	return disk.color[int(disk.width * round(x_1 - 0.5) + round(y_1 - 0.5))];
}

vec3 CTracer::StarsColor(vec3 dir){
	float y = (atan2(dir.x, dir.y) + M_PI) / (2 * M_PI);
	float x = (asin(dir.z) + M_PI / 2) / M_PI;

	y *= background_width, x *= background_height;
	if (x < 0)
		x += background_height;
	if (y < 0)
		y += background_width;
	if (y >= background_width)
		y -= background_width;
	if (x >= background_height)
		x -= background_height;
	return stars_pixels[int(round(x - 0.5) * background_width + round(y - 0.5))];
}

vec3 CTracer::TraceRay(SRay ray)
{
	int iterations_limit = 5000;
	int cur_iteration = 0;

	float distance = length(ray.m_start);
	float base_dt = 0.5f, dt;								   	 //промежуток времени

	vec3 a;		//полное ускорение, нормальна€ coставл€юща€ ускорени

	vec3 prev_pos;
	vec3 cur_pos;
	double result_sphere;
	double result_disk;

	while (cur_iteration++ < iterations_limit) {
		float distance = length(ray.m_start);
		vec3 r = -ray.m_start;
		a = r * (float)(Gravity * Mass / (distance * distance * distance)); //ускорение в текущий момент времени
		a = perp(a, ray.m_dir);
		prev_pos = ray.m_start;

		dt = base_dt;
		ray.m_start += dt * ray.m_dir + a * dt * dt / 2.0f;

		cur_pos = ray.m_start;
		ray.m_dir += a * dt;														//изменение скорости за счЄт нормального ускорени€
		ray.m_dir = normalize(ray.m_dir);

		if (RaySphereIntersection(prev_pos, ray.m_dir, blackhole.pos, blackhole.r)) {
			return vec3();
		}
		if (RayDiskIntersection(prev_pos, cur_pos, disk.r_max)) {
			return DiskColor(prev_pos, cur_pos);
		}
	}
	return StarsColor(ray.m_dir);
}

bool CTracer::RayDiskIntersection(vec3 ray_pos1, vec3 ray_pos2, double r) {
	return (ray_pos1.z * ray_pos2.z < 0) && (length(vec2(ray_pos1.x, ray_pos1.y)) <= r);
}

bool CTracer::RaySphereIntersection(vec3 ray_pos, vec3 ray_dir, vec3 sphere_pos, double r) {
	return (length(ray_pos) < r);
}

void CTracer::CameraSetUp(vec3 pos, double xRes, double yRes)
{
	m_camera.m_pos = pos;
	m_camera.m_forward = -normalize(m_camera.m_pos);
	m_camera.m_right = normalize(vec3(-m_camera.m_pos.y, !(m_camera.m_pos.x || m_camera.m_pos.y) ? 1 : m_camera.m_pos.x, 0.0f));
	m_camera.m_up = normalize(cross(m_camera.m_right, m_camera.m_forward)) * (float)(yRes / xRes);
}

void GetImage(CImage* pImage, vector<vec3> &pixels, vector<float> &alphas = vector<float>(0)) {
	int step;
	if (pImage->GetBPP() == 32) {
		step = 4;
	}
	else if (pImage->GetBPP() == 24) {
		step = 3;
	}
	else {
		cout << "Invalid image!\n";
		_exit(-1);
	}
	auto pData = (unsigned char*)pImage->GetBits();
	auto pCurrentLine = pData;
	int pitch = pImage->GetPitch();

	int index = 0;
	for (int i = 0; i < pImage->GetHeight(); ++i) // Image lines
	{
		for (int j = 0; j < pImage->GetWidth(); ++j) // Pixels in line
		{
			unsigned char b = pCurrentLine[i * pitch + j * step];
			unsigned char g = pCurrentLine[i * pitch + j * step + 1];
			unsigned char r = pCurrentLine[i * pitch + j * step + 2];
			pixels[index] = vec3(r, g, b) / 255.0f;

			if (step == 4) {
				unsigned char alpha = pCurrentLine[i * pitch + j * step + 3];
				alphas[index] = (float)alpha / 255.0f;
			}
			index++;
		}
	}
}

void CTracer::RenderImage(int xRes, int yRes)
{
	// Reading input texture sample
	CImage* pImageDisk = LoadImageFromFile("data/disk_32.png");
	CImage* pImageStars = LoadImageFromFile("data/stars.jpg");

	disk.height = pImageDisk->GetHeight();
	disk.width = pImageDisk->GetWidth();
	background_height = pImageStars->GetHeight();
	background_width = pImageStars->GetWidth();

	stars_pixels.resize(background_height * background_width);
	disk.color.resize(disk.height * disk.width);
	disk.alpha.resize(disk.height * disk.width);

	GetImage(pImageDisk, disk.color, disk.alpha);
	GetImage(pImageStars, stars_pixels);

	// Rendering
	m_camera.m_resolution = uvec2(xRes, yRes);
	m_camera.m_pixels.resize(xRes * yRes);

	float pos_x = -130.0f;
	float pos_y = -250.0f;
	float pos_z = 100.0f;

	CameraSetUp(vec3(pos_x, pos_y, pos_z), xRes, yRes);

	blackhole.pos = vec3(0.0f, 0.0f, 0.0f);
	blackhole.r = Radius;

	disk.pos = vec3(0.0f, 0.0f, 0.0f);
	disk.r_max = 8 * Radius;

	int row = 0;
#pragma omp parallel for
	for (int i = 0; i < yRes; i++){
		for (int j = 0; j < xRes; j++)
		{
			SRay ray = MakeRay(uvec2(j, i));
			m_camera.m_pixels[i * xRes + j] = TraceRay(ray); //присваиваим посчитанные цвета
		}
		row++;
		cout << "\r" << 100 * row / yRes << "%\r";
	}
	cout << endl;
}

void CTracer::SaveImageToFile(std::string fileName)
{
	CImage image;

	int width = m_camera.m_resolution.x;
	int height = m_camera.m_resolution.y;

	image.Create(width, height, 24);

	int pitch = image.GetPitch();
	unsigned char* imageBuffer = (unsigned char*)image.GetBits();

	if (pitch < 0)
	{
		imageBuffer += pitch * (height - 1);
		pitch = -pitch;
	}

	int i, j;
	int imageDisplacement = 0;
	int textureDisplacement = 0;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			vec3 color = m_camera.m_pixels[textureDisplacement + j];

			imageBuffer[imageDisplacement + j * 3] = clamp(color.b, 0.0f, 1.0f) * 255.0f;
			imageBuffer[imageDisplacement + j * 3 + 1] = clamp(color.g, 0.0f, 1.0f) * 255.0f;
			imageBuffer[imageDisplacement + j * 3 + 2] = clamp(color.r, 0.0f, 1.0f) * 255.0f;
		}

		imageDisplacement += pitch;
		textureDisplacement += width;
	}

	image.Save(fileName.c_str());
	image.Destroy();
}

CImage* CTracer::LoadImageFromFile(std::string fileName)
{
	CImage* pImage = new CImage;

	if (SUCCEEDED(pImage->Load(fileName.c_str())))
		return pImage;
	else
	{
		delete pImage;
		return NULL;
	}
}*/