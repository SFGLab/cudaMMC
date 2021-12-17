/*
 * Cluster.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: psz
 */

#include "../include/Cluster.h"

Cluster::Cluster() { init(); }

Cluster::Cluster(int start, int end) {
  init();
  this->start = start;
  this->end = end;
  genomic_pos =
      (start + end) / 2; // Pozycja genomiczna jako położenie środka klastra
}

void Cluster::init() {
  start = 0;
  end = 0;
  genomic_pos = 0;
  orientation = 'N';
  base_start = 0;
  base_end = 0;
  parent = -1;
  is_fixed = false;
  level = 0;
  dist_to_next = 0.0;
  pos.set(0.0f, 0.0f, 0.0f);
}

void Cluster::print() {
  printf("pos=%d (%d-%d, %d-%d), lvl=%d, par=%d, next=%lf, pts=(%f %f %f)",
         genomic_pos, start, end, base_start, base_end, level, parent,
         dist_to_next, pos.x, pos.y, pos.z);
  if (orientation != 'N')
    printf(" %c", orientation);

  if (is_fixed)
    printf(" fixed");

  if (siblings.size() > 0) {
    printf(", sibl=[");
    for (unsigned int i = 0; i < siblings.size(); ++i)
      printf("%d ", siblings[i]);
    printf("]");
  }

  if (arcs.size() > 0) {
    printf(", arcs=[");
    for (unsigned int i = 0; i < arcs.size(); ++i)
      printf("%d ", arcs[i]);
    printf("]");
  }

  if (children.size() > 0) {
    printf(", children=[");
    printv(children, false);
    printf("]");
  }

  printf("\n");
}

bool Cluster::contains(int genomic_pos) {
  return genomic_pos >= start && genomic_pos <= end;
}

void Cluster::toFile(FILE *file) {
  fprintf(file, "%d %d %c %f %f %f %lu", start, end, orientation, pos.x, pos.y,
          pos.z, children.size());
  for (unsigned int i = 0; i < children.size(); ++i)
    fprintf(file, " %d", children[i]);
  fprintf(file, "\n");
}

void Cluster::toFilePreviousFormat(FILE *file) {
  fprintf(file, "%d %d %d %f %f %f %lu", (start + end) / 2, start, end, pos.x,
          pos.y, pos.z, children.size());
  for (unsigned int i = 0; i < children.size(); ++i)
    fprintf(file, " %d", children[i]);
  fprintf(file, "\n");
}

void Cluster::fromFile(FILE *file) {
  int st, end;
  float x, y, z;
  int children_cnt, tmp;
  char c;

  fscanf(file, "%d %d %c %f %f %f %d", &st, &end, &c, &x, &y, &z,
         &children_cnt);

  for (int i = 0; i < children_cnt; ++i) {
    fscanf(file, "%d", &tmp);
    children.push_back(tmp);
  }

  this->start = st;
  this->end = end;
  this->genomic_pos = (st + end) / 2;
  this->orientation = c;
  pos.set(x, y, z);
}
