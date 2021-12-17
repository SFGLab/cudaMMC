/*
 * Density.h
 *
 *  Created on: Jan 15, 2015
 *      Author: psz
 */

#ifndef DENSITY_H_
#define DENSITY_H_

#include "../src/lib/common.h"

class Density {
public:
  Density();

  void fromFile(string filename); // read 3d density array
  void toFile(string filename);   // save 3d density array

  void toSegmentationFile(string filename); // write segmentation file (x, y, z,
                                            // density; line by line)
  void fromSegmentationFile(string filename); // read segmentation file (x, y,
                                              // z, density; line by line)

  void init(int size_x, int size_y, int size_z, bool is_static);
  void init(int size_x, int size_y, int size_z, int start_x, int start_y,
            int start_z, bool is_static);

  void normalize();

  void clear();

  void print();

  Density scaleDown(int factor); // used to create smaller maps by averaging
                                 // multiple cells into one

  void addPointMass(int x, int y, int z, float val);

  vector<vector<vector<float>>> t;

  int size_x, size_y, size_z;

  int range_x_start, range_x_end;
  int range_y_start, range_y_end;
  int range_z_start, range_z_end;
  float range_d_start, range_d_end;

  vector3 center;
  vector3 origin;

  float scale; // scale from settings (used for loading density data in
               // visualizer)

private:
  struct key_3d {
    int x, y, z;
    float value;
  };

  bool isInside(int x, int y, int z);

  void clearOdw();

  vector<vector<vector<bool>>> odw;

  bool is_static; // whether the map is static (true for input density, false
                  // for current structure based one)

  const int dx[6] = {-1, 1, 0, 0, 0, 0};
  const int dy[6] = {0, 0, 1, -1, 0, 0};
  const int dz[6] = {0, 0, 0, 0, 1, -1};
};

#endif /* DENSITY_H_ */
