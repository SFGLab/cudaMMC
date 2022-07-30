/*
 * Cluster.h
 *
 *  Created on: Jun 1, 2014
 *      Author: psz
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <stdio.h>
#include <vector>

#include <common.h>
#include <InteractionArc.h>

class Cluster {
public:
  Cluster();
  Cluster(int start, int end);

  void init();
  void print();

  void toFile(FILE *file);
  void toFilePreviousFormat(FILE *file);
  void fromFile(FILE *file);

  bool contains(int genomic_pos); // check if a genomic position is contained in
                                  // a given cluster

  vector3 pos;     // 3D position
  int genomic_pos; // genomic position
  int start, end;  // genomic start and end position
  char orientation;

  int parent;
  int level;

  int base_start, base_end;

  std::vector<int> arcs;
  std::vector<int> siblings;

  std::vector<int> children;

  bool is_fixed;
  double dist_to_next;
};

#endif /* CLUSTER_H_ */
