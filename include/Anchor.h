/*
 * Anchor.h
 *
 *  Created on: Apr 27, 2014
 *      Author: psz
 */

#ifndef ANCHOR_H_
#define ANCHOR_H_

#include "InteractionArc.h"
#include <string>

class Anchor {
public:
  Anchor();
  // Anchor(std::string chr, int start, int end);
  Anchor(std::string chr, int start, int end, char orientation);

  int length();
  bool contains(int pos);

  std::string chr;
  int start;
  int end;
  int center;
  char orientation; // denotes the corresponding motif orientation [possible
                    // values: 'L', 'R', 'N' (stand for not available)]
};

#endif /* ANCHOR_H_ */
