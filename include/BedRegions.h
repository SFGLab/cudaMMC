#ifndef BEDREGIONS_H_
#define BEDREGIONS_H_

#include <stdio.h>
#include <string>
#include <vector>

#include <BedRegion.h>

using namespace std;

class BedRegions {
public:
  BedRegions();

  void fromFile(std::string path);

  void print();
  void addNewIntervals(std::string chr, int start, int end, int step);

  std::vector<BedRegion> regions;
};

#endif /* BEDREGIONS_H_ */
