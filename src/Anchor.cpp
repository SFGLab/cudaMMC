#include <Anchor.h>

Anchor::Anchor() {
  start = 0;
  end = 0;
  center = 0;
  chr = "";
  orientation = 'N';
}

Anchor::Anchor(std::string chr, int start, int end, char orientation) {
  this->start = start;
  this->end = end;
  this->chr = chr;
  this->center = (start + end) / 2;
  this->orientation = orientation;
}

int Anchor::length() {
  if (start == 0 && end == 0)
    return 0;
  return end - start + 1;
}

bool Anchor::contains(int pos) { return pos >= start && pos <= end; }
