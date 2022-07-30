#ifndef INTERACTIONARC_H_
#define INTERACTIONARC_H_

#include <stdio.h>
#include <string>

/*
 * class representing an interaction between two genomic regions.
 *
 * connected clusters are denoted by 'start' and 'end'
 * Initially these indices are local and refer to anchor index on a specific
 * chromosome, but later (in LooperSolver::CreateTree()) they are shifted to
 * reflect global, absolute indexing ('cluster' in LooperSolver).
 *
 * it may happen that there are multiple interactions between two specific
 * regions. if they are of the same factor, we can simply add their scores and
 * keep only one arc; if not, we want to keep all interactions for different
 * factors (for example to display them) but also to consider all of them as a
 * single, stronger interaction. So, we create an additional, "summary"
 * interaction, with factor=-1 and eff_score being the sum of scores for all
 * corresponding arcs. For those arcs we set eff_score=0. Thus, when creating
 * the structure we consider 'eff_score', and when displaying - 'score'.
 */

class InteractionArc {
public:
  InteractionArc();
  InteractionArc(int start, int end, int score, int factor = -1);

  void toFile(FILE *file);
  void fromFile(FILE *file);

  void init();
  void print();
  int length();

  bool operator<(const InteractionArc &arc) const {
    if (start < arc.start)
      return true;
    if (start == arc.start) {
      if (end < arc.end)
        return true;
      if (end == arc.end)
        return factor < arc.factor;
    }
    return false;
  }

  int start;
  int end;
  int genomic_start;
  int genomic_end;
  int score;
  int eff_score; // effective score
  int factor;
};

#endif /* INTERACTIONARC_H_ */
