/*
 * HierarchicalChromosome.h
 *
 *  Created on: Jun 2, 2014
 *      Author: psz
 */

#ifndef HIERARCHICALCHROMOSOME_H_
#define HIERARCHICALCHROMOSOME_H_

#include <map>
#include <stdio.h>
#include <vector>

#include "Chromosome.h"
#include "Cluster.h"
#include "InteractionArcs.h"
#include <optional>
#ifdef __CUDACC__
#include <vector_types.h>
#else
struct float3 {};
#endif

//

class HierarchicalChromosome {
public:
  HierarchicalChromosome();

  void print(int level = 1);
  void printRegionsTree();
  void printRegionsTree(string chr);

  void toFile(string filename);
  void toFile(FILE *file);
  void toFilePreviousFormat(string filename);
  void fromFile(string filename);
  void fromFilePreviousFormat(FILE *file);
  void fromFile(FILE *file);
  bool fromTxt(string filename);

  bool fromHiCEvo(
      string filename); // create chr based on hic-evo hierarchical chromosome

  void useTopLevel();
  void useLowestLevel();

  void setLevel(int level);
  void levelDown();
  void expandRegion(int start, int end, bool include_external = true);

  vector3 getCenter(bool current = true);
  void
  center(bool current = true); // translate so that mass center is in the origin
  map<string, float> getAvgAnchorDistance();

  // find a 3D point corresponding to a given genomic position
  vector3 find3DPosition(string chr, int pos); // (use structure 'as it is')
  vector3 find3DSmoothPosition(string chr, int pos); // (use smoothed structure)

  void scale(float factor);
  // void align(const HierarchicalChromosomeMixed &hc);

  // int findClosestCluster(int genomic_position);	// given genomic position
  // finds closest cluster

  void createCurrentLevelStructure();

  Heatmap
  createStructuralHeatmap(string chr,
                          int level); // create structural heatmap (ie. with
                                      // pairwise distances between beads)
  float
  calcDistance(HierarchicalChromosome hc,
               int level); // find structural distance to hc on a given level

  void smoothSpline(int n);
  int clusterIndexToSmoothedIndex(int ind);

  // create a flat structure with beads covering the same genomic region
  // (exactly 'resolution_bp' bp)
  Chromosome createEqidistantModel(int resolution_bp, string chr = "");

  HierarchicalChromosome
  extractFragment(int start, int end); // extract fragment of structure.

  // calculate some basic stats for chromosomes at different scales, eg. mean &
  // max distance, chromosome diameter etc.
  Heatmap getDistancesHeatmap();
  void getSpatialDistributionStats();

  // returns indices of anchors closest to specified position
  // vector<int> findFlankingAnchors(int genomic_position);
  vector<int> findFlankingAnchors(string chr, int genomic_position_start,
                                  int genomic_position_end);
  // vector<int> findOuterLoop(string chr, int genomic_position_start, int
  // genomic_position_end);

  std::map<std::string, Chromosome> chr;
  std::map<std::string, Chromosome> chr_smooth;

  std::vector<string> chrs; // list of chromosomes

  std::vector<Cluster> clusters;
  std::map<std::string, std::vector<int>> current_level;
  std::map<std::string, int> chr_root; // index of chromosome's root

  InteractionArcs arcs;

  std::optional<std::vector<float3>> gpuHelper(const int, const int, const int,
                                               std::string);

private:
  void printRegionsTreeRecursive(int region_index_curr, int margin = 0);

  std::vector<int> expandRegion(int region_ind, int start, int end,
                                bool include_external);

  int smooth_factor;
};

#endif /* HIERARCHICALCHROMOSOME_H_ */
