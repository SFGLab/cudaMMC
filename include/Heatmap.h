#ifndef HEATMAP_H_
#define HEATMAP_H_

#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include <locale>

#include <common.h>
#include <Settings.h>

using namespace std;

class Heatmap {

public:
  Heatmap();
  Heatmap(int size);

  void init(int size);
  void zero(); // zero all values
  void print();

  void toFile(string filename, bool total_count = true,
              bool zero_diagonal = false, bool as_integers = false);
  void fromFile(string filename);
  void fromFile(string filename, bool labels);
  void fromMDS(string filename, int size);

  bool isEmpty();
  vector<bool> getEmpty();
  Heatmap removeColumns(vector<bool> del);
  Heatmap removeEmptyColumns();

  void scale(float scale);
  void add(const Heatmap &heat);
  void add(float val);
  Heatmap diff(const Heatmap &heat,
               bool abs = false); // calc difference between two heatmaps
  void divide(const Heatmap &hmap);
  float calcDistance(const Heatmap &hmap);

  void smooth(float threshold, float factor);
  void smooth(float threshold, float factor, set<int> &breaks);

  int getDiagonalSize(); // gives size of diagonal (i.e. width of zero cells
                         // band on the diagonal)
  void getRange(float &min, float &max);
  float getAvg();             // calc avg value
  float getAvgNearDiagonal(); // calc avg value from the non-zero cells closest
                              // to diagonal
  void calcAvgValues(
      bool count_zeros = true); // calc avg values for k-distant cells, for all
                                // k. Save result to 'avg_value'

  void clearDiagonal(int diag = 1);

  vector<float> toVector(int diag = 0);

  size_t size;
  float **v;

  int start;
  int resolution;

  int diagonal_size;

  std::vector<float> avg_value;

private:
  void clear();
};

#endif /* HEATMAP_H_ */
