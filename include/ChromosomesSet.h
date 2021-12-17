/*
 * ChromosomesSet.h
 *
 *  Created on: May 25, 2014
 *      Author: psz
 */

#ifndef CHROMOSOMESSET_H_
#define CHROMOSOMESSET_H_

#include "Chromosome.h"
#include <string.h>
#include <vector>

class ChromosomesSet {
public:
  ChromosomesSet();

  void print();

  void add(std::map<std::string, Chromosome> chr);
  void add(std::map<std::string, Chromosome> chr, string desc);

  void toFile(string filename);
  void toFile(FILE *file);
  void fromFile(string filename);
  void fromFile(FILE *file);

  std::vector<std::map<std::string, Chromosome>> chromosome;
  std::vector<string> desc;
};

#endif /* CHROMOSOMESSET_H_ */
