/*
 * BedRegions.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: psz
 */

#include "../include/BedRegions.h"

BedRegions::BedRegions() {}

void BedRegions::fromFile(std::string filename) {
  FILE *f;
  f = fopen(filename.c_str(), "r");
  if (f == NULL) {
    printf("Error opening file [%s]!\n", filename.c_str());
    return;
  }

  char chr[16], line[100];
  int start, end;

  while (!feof(f)) {
    if (fscanf(f, "%s %d %d", chr, &start, &end) != 3)
      continue;
    // printf("%s %d %d\n", chr, start, end);
    fgets(line, 100, f); // read to the end of line

    BedRegion b((std::string)chr, start, end);
    regions.push_back(b);
  }
  fclose(f);
}

void BedRegions::print() {
  printf("regions: %d\n", (int)regions.size());
  for (unsigned int i = 0; i < regions.size(); ++i) {
    printf("[%s] %d %d\n", regions[i].chr.c_str(), regions[i].start,
           regions[i].end);
  }
}

void BedRegions::addNewIntervals(std::string chr, int start, int end,
                                 int step) {

  if (step < 1) {
    printf("Step must be at least 1!\n");
    return;
  }

  while (start < end) {
    BedRegion reg(chr, start, start);
    regions.push_back(reg);
    start += step;
  }
}
