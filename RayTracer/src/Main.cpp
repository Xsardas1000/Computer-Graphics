#include "Tracer.h"
#include "stdio.h"
#include <iostream>

using namespace std;

void main(int argc, char** argv)
{
  CTracer tracer;
  CScene scene;

  int xRes = 1280;  // Default parameters
  int yRes = 720;
  float x_cam = -250;
  float y_cam = -250;
  float z_cam = 20;
  double blackhole_mass = 8.57e36f;
  int k = 8;


  if(argc == 2) // There is input file in parameters
  {
    FILE* file = fopen(argv[1], "r");
    if(file)
    {
      int xResFromFile = 0;
      int yResFromFile = 0;
	  
      if(fscanf(file, "%d %d %f %d %f %f %f", &xResFromFile, 
											  &yResFromFile, 
											  &blackhole_mass, 
											  &k, 
											  &x_cam, 
											  &y_cam,
											  &z_cam) == 7)
      {
        xRes = xResFromFile;
        yRes = yResFromFile;
      }
      else
        printf("Invalid config format! Using default parameters.\r\n");

      fclose(file);
    }
    else
      printf("Invalid config path! Using default parameters.\r\n");
  }
  else
    printf("No config! Using default parameters.\r\n");

  tracer.m_pScene = &scene;
  tracer.RenderImage(xRes, yRes, blackhole_mass, k, x_cam, y_cam, z_cam);
  tracer.SaveImageToFile("Result.png");
}
