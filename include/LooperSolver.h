/*
 * LooperSolver.h
 *
 *  Created on: Jun 1, 2014
 *      Author: psz
 */

#ifndef LOOPERSOLVER_H_
#define LOOPERSOLVER_H_

#include <fstream>
#include <queue>
#include <stdarg.h>
#include <string>
#include <vector>

#include "BedRegion.h"
#include "BedRegions.h"
#include "ChromosomesSet.h"
#include "Cluster.h"
#include "Density.h"
#include "Heatmap.h"
#include "HierarchicalChromosome.h"
#include "InteractionArc.h"
#include "InteractionArcs.h"

class LooperSolver {
public:
  LooperSolver(string label, string outdir);

  void output(int level, const char *format, ...);

  void print();
  void printActiveRegion();

  void setContactData(std::vector<string> chrs_list,
                      BedRegion region_of_interest, string anchors,
                      std::vector<string> factors,
                      std::vector<string> arcs_clusters,
                      std::vector<string> arcs_singletons,
                      std::vector<string> arcs_singletons_inter);

  void initDensityData();

  void createTreeGenome(); // create clusters tree
  int createTreeChromosome(string chr);

  void reconstructClustersArcsDistances();
  void reconstructClusterArcsDistances(int cluster, int current_ind,
                                       bool smooth = false,
                                       bool use_subanchor_heatmap = false);

  void reconstructClustersHeatmap();
  void reconstructClustersHeatmapSingleLevel(int heatmap_ind);

  Heatmap normalizeHeatmap(const Heatmap &heat);
  void normalizeHeatmapInter(Heatmap &heat, float scale = 2.0f);
  void normalizeHeatmapDiagonalTotal(
      Heatmap &heat,
      float val = 1.0f); // normalize for the total sum of near-diagonal cells

  // void findBinsBoundaries(std::map<std::string, std::vector<int> > &breaks,
  // std::map<std::string, int> &start_ind);

  float getHeatmapAvgBetweenNeighboringRegions(Heatmap &heat, int level);

  void createDistanceHeatmap(int heatmap_ind);
  Heatmap createIFHeatmapFromDistances(const Heatmap &h);

  void setTopLevel();
  void setLevel(int level);
  void levelDown();

  void reset(); // reset all needed variables between consecutive runs (eg. for
                // ensembles)
  void removeSubanchorBeads(); // remove the subanchor beads (remove the
                               // corresponding clusters, clean
                               // cluster::children etc)

  HierarchicalChromosome getModel();

  void repairStructureHeatmap();
  int findUnloopedChromatinLength(int cluster);

  float freqToDistanceHeatmap(float freq);
  float freqToDistanceHeatmapInter(float freq);

  float distToFreqHeatmap(float freq);
  float freqToDistance(int freq, bool memo = true);

  void densify(); // add empty clusters to split too long fragments

  void diagnose();
  void getSnapshot(string s = "<no_desc>");
  void addSnapshot(Chromosome chr, string s = "<no_desc>");

  void inputArbitraryHeatmap(Heatmap h);

  std::map<std::string, Chromosome> getCurrentChromosomes();

  std::vector<int> findGaps(string chr);
  std::vector<int> findGaps(int start_index);

  // find the optimal split of gaps, so that the bins have size ~'size'
  std::vector<int> findSplit(std::vector<int> gaps, int size, string chr);

  void setDebuggingOptions(int chr_number_limit = -1, int length_limit = -1);

  void showCoordinatesRange(int level);
  void showCoordinatesRange(bool active_region = false);
  void calcAnchorsStats();

  // calculate the length of region (st, end) but excluding the centromere
  int getGenomicLengthExcludingCentromere(string chr, int st, int end);

  // int root_index;

  // we keep all clusters (for all chromosomes) in a single array (it is easier
  // to represent structure tree that way) for arcs, we keep a separate list for
  // every chromosome (in 'arcs'). as a result, clusters[i].arcs contains
  // indices of arcs relative to the appropriate chromosome (so, multiple
  // clusters will have common arcs indices, but - if the clusters are on
  // different chrs - they refer to different arcs)
  InteractionArcs arcs;
  std::vector<Cluster> clusters;
  std::vector<Cluster> base_clusters;

  std::map<std::string, int> chr_root; // index of chromosome's root

  std::map<std::string, std::vector<int>> current_level;
  std::vector<int> active_region;

  std::vector<Heatmap> heatmap; // heatmaps for all singleton levels (eg.
                                // chromosome and segment levels)
  Heatmap heatmap_dist;         // heatmap with expected distances

  ChromosomesSet chr_set;
  ChromosomesSet chromosome_set;

  std::map<string, int> chr_first_cluster; // gives an index of first cluster

  std::map<string, int> chr_length;
  std::vector<double> chr_length_v;

  std::vector<string>
      chrs; // list of all chromosomes that are to be reconstructed

  // we can force only a specific chromosome region to be reconstructed (with -c
  // option) following variables keep track of it
  bool is_bed_region_provided; // if yes, then only a region described by
                               // 'region_of_interest' is reconstructed
  BedRegion selected_region;   // BED region that is to be reconstructed

  // find ib and segment containing the specified position and print short
  // summary
  void printStructuresForGenomicPosition(string chr, int pos);
  int findClusterForGenomicPosition(
      string chr,
      int pos); // find current_level-index of the specified positions
  int findCurrentLevelForGenomicPosition(
      string chr,
      int pos); // find current_level-index of the specified positions

  float ParallelMonteCarloHeatmap(float step_size);
  double MonteCarloHeatmap(float step_size);
  double MonteCarloHeatmapAndDensity(float step_size);

  // functions to convert coordinates between density and structure space
  vector3 densityCoordToStructure(vector3 coord);
  vector3 structureCoordToDensity(vector3 pos);

private:
  Heatmap createSingletonHeatmap(int diag = 1);
  void createSingletonSubanchorHeatmap(string chr,
                                       vector<int> &anchors_gap_len);
  void createExpectedDistSubanchorHeatmap(const Heatmap &subanchor_avg_dist);

  // heatmaps (segment and chromosome level)
  double calcScoreHeatmapActiveRegion(int moved = -1);

  double calcScoreDensity();

  // arcs (anchor level)
  double calcScoreDistancesActiveRegion();
  double calcScoreDistancesActiveRegion(int cluster_moved); // arcs

  // linker length (anchor level)
  double calcScoreStructureActiveRegionLengths();
  double calcScoreStructureActiveRegionLengths(int cluster_moved);

  // smoothing: linker lenght + angles (subanchor level)
  double calcScoreStructureSmooth(bool lengths, bool angles);
  double calcScoreStructureSmooth(int cluster_moved, bool lengths, bool angles);

  // score based on anchor orientation (eg. CTCF motif)
  double calcScoreOrientation(const vector<vector3> &orientation);
  double calcScoreOrientation(const vector<vector3> &orientation,
                              int anchor_index);

  double calcScoreSubanchorHeatmap(int cluster_moved = -1);

  // calc orientation of a cluster (index should point to an anchor in
  // active_region)
  vector3 calcOrientation(int cluster_index);

  float parallelMonteCarloArcs(float step_size);
  double MonteCarloArcs(
      float step_size); // simulation for anchor level, only clusters. score
                        // based on arcs and distance between neighbors
  double MonteCarloArcsSmooth(
      float step_size,
      bool use_subanchor_heatmap =
          false); // simulation for ~10kb level, densified clusters, base
                  // clusters fixed. score based on distances and angles

  void updateDensity();
  void repairDensityBoundary();

  void densifyActiveRegion(int cluster, bool fix = false);

  double genomicLengthToDistance(int dist);

  Heatmap downsampleSingletonHeatmap(
      Heatmap &heat); // create a lvl1 heatmap from lvl2 singletons heatmap

  Heatmap calculateRegionsDistances(
      const std::vector<int> clusters_ind); // create a heatmap with real 3D
                                            // distances between regions.

  int otherEnd(const InteractionArc &arc, int current_cluster);
  int otherEnd(int arc, int current_cluster);
  int otherEnd(string chr, int arc, int current_cluster);

  void interpolateChildrenPosition(
      std::vector<int> &regions); // interpolate lower level
  void interpolatePosition(std::vector<int> &regions,
                           std::vector<int> &clusters_template, int region_ind,
                           bool noise = false);

  void interpolateChildrenPositionSpline(
      std::vector<int> &regions,
      bool use_genomic_dist = false); // interpolate lower level with splines

  void positionInteractionBlocks(std::vector<int> &segments);

  float boundaryScore(string chr, vector3 pt);

  string activeRegionIndexToChr(int ind);

  std::set<int> getCurrentHeatmapBreaks();

  Heatmap calcTrueDistancesHeatmapForRegion(std::vector<int> &regions);

  void calcActiveAnchorsNeighbors();

  void calcAnchorExpectedDistancesHeatmap();

  void getChromosomeHeatmapBoundary(int p, int &start, int &end);

  string current_chr; // currently processed chromosome

  vector3 rw_pos; // helper for random walk

  bool zero_distance;

  float freq_to_distance[101];

  int tmp_var;

  std::vector<string> arcs_clusters;
  std::vector<string> arcs_singletons;
  std::vector<string> arcs_singletons_inter;

  bool split_singleton_files_by_chr;

  std::map<string, int>
      active_region_to_chr; // keep max cluster id for every chromosome (let us
                            // get chr name from cluster id)

  string output_dir;
  string label; // label, used to name heatmap files etc.

  // some debug options
  int debug_length_limit; // length of each chromosome to be reconstructured
                          // fragment of that length will be reconstructed
  int debug_chromosomes_limit; // number of chromosomes to be reconstructued

  BedRegions centromeres;        // centromere regions
  BedRegions interaction_blocks; // spli into interaction blocks

  BedRegions segments_predefined;

  Heatmap
      heatmap_anchor; // singletons heatmap for reconstruction of anchor level
  Heatmap heatmap_exp_dist_anchor;
  std::map<int, int> active_to_cluster_index;

  Heatmap heatmap_subanchor; // singletons heatmap for reconstruction of
                             // subanchor level
  Heatmap heatmap_dist_subanchor;

  std::map<int, vector<int>> active_anchors_neighbors;
  std::map<int, vector<float>> active_anchors_neighbors_weight;

  vector<int> heatmap_chromosome_boundaries; // keep indices of chromosome
                                             // boundaries for heatmaps

  Density density;      // 3D density map (from microscopy)
  Density density_curr; // density map of the current structure
};

#endif /* LOOPERSOLVER_H_ */
