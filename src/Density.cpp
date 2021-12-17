/*
 * Density.cpp
 *
 *  Created on: Jan 15, 2015
 *      Author: psz
 */

#include "../include/Density.h"
#include "../include/Settings.h"
#include "lib/common.h"
#include "lib/mtxlib.h"
#include <queue>

Density::Density() { size_x = size_y = size_z = 0; }

void Density::print() {
  printf("Density (%d x %d x %d)\n", size_x, size_y, size_z);
  printf("range X = %d..%d\n", range_x_start, range_x_end);
  printf("range Y = %d..%d\n", range_y_start, range_y_end);
  printf("range Z = %d..%d\n", range_z_start, range_z_end);
  printf("range d = %f..%f\n", range_d_start, range_d_end);
  print_vector(origin, "origin");
  print_vector(center, "center");
}

void Density::clear() {
  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      for (int k = 0; k < size_z; ++k) {
        t[i][j][k] = 0.0f;
        if (!is_static)
          odw[i][j][k] = false;
      }
    }
  }
}

void Density::clearOdw() {
  if (is_static)
    return;
  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      for (int k = 0; k < size_z; ++k)
        odw[i][j][k] = false;
    }
  }
}

void Density::init(int size_x, int size_y, int size_z, int start_x, int start_y,
                   int start_z, bool is_static) {
  this->init(size_x, size_y, size_z, is_static);
  range_x_start = start_x;
  range_y_start = start_y;
  range_z_start = start_z;
  range_x_end = start_x + size_x;
  range_y_end = start_y + size_y;
  range_z_end = start_z + size_z;

  if (size_x > 0)
    range_x_end--;
  if (size_y > 0)
    range_y_end--;
  if (size_z > 0)
    range_z_end--;
}

void Density::init(int size_x, int size_y, int size_z, bool is_static = true) {
  this->size_x = size_x;
  this->size_y = size_y;
  this->size_z = size_z;

  this->is_static = is_static;

  t.resize(size_x);
  for (int i = 0; i < size_x; ++i) {
    t[i].resize(size_y);
    for (int j = 0; j < size_y; ++j)
      t[i][j].resize(size_z);
  }

  if (!is_static) {
    odw.resize(size_x);
    for (int i = 0; i < size_x; ++i) {
      odw[i].resize(size_y);
      for (int j = 0; j < size_y; ++j)
        odw[i][j].resize(size_z);
    }
  }

  // save the center
  center = vector3((range_x_start + range_x_end) * 0.5f,
                   (range_y_start + range_y_end) * 0.5f,
                   (range_z_start + range_z_end) * 0.5f);
  origin = vector3(range_x_start, range_y_start, range_z_start);
}

void Density::normalize() {
  float total = 0.0f;
  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      for (int k = 0; k < size_z; ++k)
        total += t[i][j][k];
    }
  }

  float avg = total / ((float)size_x * size_y * size_z);

  if (avg < epsilon)
    return;

  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      for (int k = 0; k < size_z; ++k)
        t[i][j][k] /= avg;
    }
  }
}

void Density::addPointMass(int x, int y, int z, float val) {

  if (x < 0 || y < 0 || z < 0 || x >= size_x || y >= size_y || z >= size_z) {
    printf("mass outside: (%d %d %d)\n", x, y, z);
    return;
  }

  clearOdw();

  queue<key_3d> q;
  key_3d p, np;
  key_3d first = {x, y, z, val};

  q.push(first);
  odw[x][y][z] = true;

  while (q.size()) {
    p = q.front();
    q.pop();

    t[p.x][p.y][p.z] = max(t[p.x][p.y][p.z], p.value);
    // printf("%d %d %d %f\n", p.x, p.y, p.z, val);

    for (int i = 0; i < 6; i++) {
      np.x = p.x + dx[i];
      np.y = p.y + dy[i];
      np.z = p.z + dz[i];
      np.value = p.value * Settings::densityInfluence;

      if (np.x >= 0 && np.y >= 0 && np.z >= 0 && np.x < size_x &&
          np.y < size_y && np.z < size_z && np.value > 0.001f) {
        if (!odw[np.x][np.y][np.z]) {
          q.push(np);
          odw[np.x][np.y][np.z] = true;
        }
      }
    }
  }
}

Density Density::scaleDown(int factor) {
  Density d;

  int new_size_x = ceil(size_x / (float)factor);
  int new_size_y = ceil(size_y / (float)factor);
  int new_size_z = ceil(size_z / (float)factor);

  d.init(new_size_x, new_size_y, new_size_z, range_x_start, range_y_start,
         range_z_start, true);

  for (int i = 0; i < new_size_x; ++i) {
    for (int j = 0; j < new_size_y; ++j) {
      for (int k = 0; k < new_size_z; ++k) {

        double avg = 0.0;
        int cnt = 0; // count how many cells are averaged (cells on boundaries
                     // might have fewer)

        for (int ii = 0; ii < factor; ++ii) {
          for (int jj = 0; jj < factor; ++jj) {
            for (int kk = 0; kk < factor; ++kk) {
              if (i * factor + ii < size_x && j * factor + jj < size_y &&
                  k * factor + kk < size_z) {
                avg += t[i * factor + ii][j * factor + jj][k * factor + kk];
                cnt++;
              }
            }
          }
        }

        avg /= cnt;

        d.t[i][j][k] = avg;
      }
    }
  }

  d.print();
  return d;
}

bool Density::isInside(int x, int y, int z) { return t[x][y][z]; }

void Density::fromSegmentationFile(string filename) {

  FILE *f = open(filename, "r");
  if (f == NULL)
    return;

  int x, y, z;
  float d;

  range_x_start = 9999999;
  range_y_start = 9999999;
  range_z_start = 9999999;
  range_d_start = 9999999.0f;
  range_x_end = 0;
  range_y_end = 0;
  range_z_end = 0;
  range_d_end = 0.0f;

  while (!feof(f)) {
    if (fscanf(f, "%d %d %d %f", &x, &y, &z, &d) != 4)
      break;
    range_x_start = min(x, range_x_start);
    range_x_end = max(x, range_x_end);
    range_y_start = min(y, range_y_start);
    range_y_end = max(y, range_y_end);
    range_z_start = min(z, range_z_start);
    range_z_end = max(z, range_z_end);
    range_d_start = min(d, range_d_start);
    range_d_end = max(d, range_d_end);
  }

  int sizex = range_x_end - range_x_start + 1;
  int sizey = range_y_end - range_y_start + 1;
  int sizez = range_z_end - range_z_start + 1;

  init(sizex, sizey, sizez, true);

  print();

  fseek(f, 0, SEEK_SET);

  while (!feof(f)) {
    if (fscanf(f, "%d %d %d %f", &x, &y, &z, &d) != 4)
      break;

    x -= range_x_start;
    y -= range_y_start;
    z -= range_z_start;

    if (x >= sizex || y >= sizey || z >= sizez) {
      printf("out of range: %d %d %d (%d %d %d)\n", x, y, z, sizex, sizey,
             sizez);
      continue;
    }

    t[x][y][z] = d;
  }

  fclose(f);
}

void Density::toSegmentationFile(string filename) {

  FILE *f = open(filename, "w");
  if (f == NULL)
    return;

  for (int i = range_x_start; i <= range_x_end; ++i) {
    for (int j = range_y_start; j <= range_y_end; ++j) {
      for (int k = range_z_start; k <= range_z_end; ++k) {
        float val = t[i - range_x_start][j - range_y_start][k - range_z_start];
        if (val > epsilon)
          fprintf(f, "%d %d %d %lf\n", i, j, k, val);
      }
    }
  }

  fclose(f);
}

void Density::fromFile(string filename) {

  FILE *f = open(filename, "r");
  if (f == NULL)
    return;

  fscanf(f, "%d %d %d %d %d %d %f", &size_x, &size_y, &size_z, &range_x_start,
         &range_y_start, &range_z_start, &scale);
  range_x_end = range_x_start + size_x - 1;
  range_y_end = range_y_start + size_y - 1;
  range_z_end = range_z_start + size_z - 1;

  init(size_x, size_y, size_z, true);

  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      for (int k = 0; k < size_z; ++k)
        fscanf(f, "%f", &t[i][j][k]);
    }
  }

  fclose(f);
}

void Density::toFile(string filename) {

  FILE *f = open(filename, "w");
  if (f == NULL)
    return;

  fprintf(f, "%d %d %d %d %d %d %f\n", size_x, size_y, size_z, range_x_start,
          range_y_start, range_z_start, Settings::densityScale);

  for (int i = 0; i < size_x; ++i) {
    for (int j = 0; j < size_y; ++j) {
      for (int k = 0; k < size_z; ++k)
        fprintf(f, "%f ", t[i][j][k]);
      fprintf(f, "\n");
    }
    fprintf(f, "\n");
  }

  fclose(f);
}
