#include <InteractionArc.h>

InteractionArc::InteractionArc() { init(); }

InteractionArc::InteractionArc(int _start, int _end, int _score, int _factor) {
  init();
  start = _start;
  end = _end;
  score = _score;
  eff_score = score;
  factor = _factor;
}

void InteractionArc::init() {
  start = -1;
  end = -1;
  genomic_start = -1;
  genomic_end = -1;
  score = 0;
  eff_score = 0;
  factor = -1;
}

void InteractionArc::toFile(FILE *file) {
  fprintf(file, "%d %d %d %d %d %d %d\n", start, end, genomic_start,
          genomic_end, score, eff_score, factor);
}

void InteractionArc::fromFile(FILE *file) {
  fscanf(file, "%d %d %d %d %d %d %d", &start, &end, &genomic_start,
         &genomic_end, &score, &eff_score, &factor);
}

int InteractionArc::length() { return end - start + 1; }

void InteractionArc::print() {
  printf("%d %d %d (eff_sc: %d, factor: %d)\n", start, end, score, eff_score,
         factor);
}
