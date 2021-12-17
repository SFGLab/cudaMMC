/*
 * BedRegion.h
 *
 *  Created on: Jun 19, 2014
 *      Author: psz
 */

#ifndef BEDREGION_H_
#define BEDREGION_H_

#include <stdio.h>
#include <string>

class BedRegion {
public:
  BedRegion(std::string _chr, int _start, int _end);
  BedRegion();

  static bool tryParse(std::string str);
  bool parse(std::string str);

  void print();

  bool contains(int pos);

  std::string chr;
  int start;
  int end;
};

#endif /* BEDREGION_H_ */
