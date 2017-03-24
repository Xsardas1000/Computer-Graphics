#pragma once

#include "glm/glm.hpp"
#include <vector>

using namespace std;

struct SRay
{
  glm::vec3 m_start;
  glm::vec3 m_dir;
};

struct SCamera
{
  glm::vec3 m_pos;          // Camera position and orientation
  glm::vec3 m_forward;      // Orthonormal basis
  glm::vec3 m_right;
  glm::vec3 m_up;

  glm::vec2 m_viewAngle;    // View angles, rad
  glm::uvec2 m_resolution;  // Image resolution: w, h

  std::vector<glm::vec3> m_pixels;  // Pixel array
};

struct SMesh
{
  std::vector<glm::vec3> m_vertices;  // vertex positions
  std::vector<glm::uvec3> m_triangles;  // vetrex indices
};

class Sphere {
public:
	glm::vec3 pos;   //Координаты центра сферы
	double r;        //Радиус сферы
	double mass;

	Sphere(double xc = 0.0f, double yc = 0.0f, double zc = 0.0f, double rc = 1.0f)
	{
		pos = glm::vec3(xc, yc, zc);
		r = rc;
	}
	~Sphere() {  }
};

class Disk{
public:
	double r_max;    //внешний радиус
	glm::vec3 pos;   //центр диска
	float height, width;
	vector<glm::vec3> color;
	vector<float> alpha;
	Disk(double xc = 0, double yc = 0, double zc = 0, double r1 = 0.0f, double r2 = 1.0f) {
		r_max = r2, pos = glm::vec3(xc, yc, zc);
	}
	~Disk() {}
};
