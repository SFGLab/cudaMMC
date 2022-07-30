/*
 * Chromosome.h
 *
 *  Created on: Aug 3, 2013
 *      Author: psz
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <string.h>
#include <vector>

#include <locale>

#include <common.h>
#include <rmsd.h>
#include <Heatmap.h>
#include <Settings.h>

using namespace std;

class Chromosome {

public:
  Chromosome();

  void init();

  void print();

  void toFile(string filename);
  void toFile(FILE *file);
  void fromFile(string filename);
  void fromFile(FILE *file, int pts_cnt = 0);
  void fromFilePDB(string filename);

  void createRandom(int pts_cnt, float size = 0.1f, bool walk = true);
  void makeLine(int pts_cnt, float dist);
  void makeSpiral(int pts_cnt, float r, float angle, float spin_height);
  void makeSpiral(int pts_cnt, float r, float angle, vector3 p1,
                  vector3 p2); // make spiral starting in p1 and ending in p2

  void createFromSubchromosomes(const vector<Chromosome> &chrs);
  void randomizeSubchromosomes(float dispersion,
                               bool keep_subchr_structure = true);
  void setSubchromosomesIndices(const vector<int> &subchr_index);

  Chromosome clone();

  void center(); // move center to (0, 0, 0)
  void scale(float scale, bool center = false);

  Heatmap createHeatmap();
  Heatmap
  createInverseHeatmap(); // create simulated heatmap ("inverse HiC"-like)

  void align(const Chromosome &chr, bool resize = false,
             float max_angle = 2 * 3.14f);
  void adjustSize(const Chromosome &chr);
  // void alignWithResize(const Chromosome &chr);

  void updateSize();

  float getDistanceSqr(const Chromosome &chr);
  vector3 getCenter();

  void translate(const vector3 &v);
  void rotate(const matrix44 &mat);

  float getDiameter() const;
  std::vector<float> getConsecutiveBeadDistances();
  float
  getAvgConsecutiveBeadDistance(); // avg distance between consecutive beads

  float calcDistance(const Chromosome &chr);
  float calcRMSD(const Chromosome &chr);

  vector<vector<float>> calcBaseDensity(float sphere_r = 0.0f);

  void trim(int start, int end = 0);

  int size;

  vector<vector3> points;
  vector<int> genomic_position;
  vector<vector<vector3>> points_hierarchical;

  float score;

private:
  float findBestAlignmentRotation(const Chromosome &chr, int steps,
                                  vector3 &angle, float max_angle);

  int getMostDistantPointsPair(float &dist);
};

#endif /* CHROMOSOME_H_ */
