#include <Heatmap.h>

Heatmap::Heatmap() {
  this->size = 0;
  v = 0;
  init(0);
}

Heatmap::Heatmap(int size) {
  this->size = 0;
  init(size);
}

bool Heatmap::isEmpty() { return size == 0; }

void Heatmap::init(int size) {
  start = -1;
  resolution = 0;
  diagonal_size = 0;

  if (this->size > 0) {
    for (size_t i = 0; i < this->size; i++)
      free(v[i]);
    free(v);
  }

  this->size = size;
  v = (float **)malloc(size * sizeof(float *));
  if (v == NULL)
    error("malloc failed!");
  for (int i = 0; i < size; i++) {
    v[i] = (float *)malloc(size * sizeof(float));
    if (v[i] == NULL)
      error("malloc failed!");
    for (int l = 0; l < size; l++)
      v[i][l] = 0.0f;
  }
}

void Heatmap::zero() {
  for (size_t i = 0; i < size; i++) {
    for (size_t l = 0; l < size; l++)
      v[i][l] = 0.0f;
  }
}

void Heatmap::clearDiagonal(int diag) {
  for (size_t i = 0; i < size; ++i) {
    for (int p = 0; p < diag && i + p < size; ++p) {
      v[i][i + p] = 0.0f;
      v[i + p][i] = 0.0f;
    }
  }

  diagonal_size = diag;
}

int Heatmap::getDiagonalSize() {
  // try different widths
  for (size_t w = 0; w < size; ++w) {

    // check if all cells with given width are zero
    for (size_t i = 0; i + w < size; ++i) {
      if (v[i][i + w] > 1e-6)
        return w;
    }
  }
  return 0;
}

void Heatmap::smooth(float threshold, float factor) {
  // printf("smooth heatmap\n");

  float avg = getAvg();
  threshold *= avg;
  // printf("threshold = %f\n", threshold);

  int diag = getDiagonalSize();
  float diff, adiff;
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = i + diag; j + 1 < size; ++j) {
      diff = v[i][j + 1] - v[i][j];
      adiff = fabs(diff);
      if (adiff > threshold) {
        // printf("%f %f, %f %f\n", diff, threshold, v[i][j], v[i][j+1]);
        diff /= 2.0f;
        v[i][j] += diff * factor;
        v[i][j + 1] -= diff * factor;
        v[j][i] = v[i][j];
        v[j + 1][i] = v[i][j + 1];
        // printf("new %f %f\n", v[i][j], v[i][j+1]);
      }
    }
  }
}

void Heatmap::smooth(float threshold, float factor, set<int> &breaks) {
  float avg = getAvg();
  threshold *= avg;
  // printf("threshold = %f\n", threshold);

  int diag = getDiagonalSize();
  float diff, adiff;
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j + 1 < size; ++j) {
      // for (int j = i+diag; j+1 < size; ++j) {

      if (j >= i - diag && j < i + diag)
        continue;
      if (breaks.find(j + 1) != breaks.end())
        continue;

      diff = v[i][j + 1] - v[i][j];
      adiff = fabs(diff);
      if (adiff > threshold) {
        // printf("%f %f, %f %f\n", diff, threshold, v[i][j], v[i][j+1]);
        diff /= 2.0f;
        v[i][j] += diff * factor;
        v[i][j + 1] -= diff * factor;
        v[j][i] = v[i][j];
        v[j + 1][i] = v[i][j + 1];
        // printf("new %f %f\n", v[i][j], v[i][j+1]);
      }
    }
  }
}

void Heatmap::getRange(float &min, float &max) {
  max = 0.0f;
  min = 1e9;
  int diag = getDiagonalSize();
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = i + diag; j < size; ++j) {
      if (v[i][j] < min)
        min = v[i][j];
      if (v[i][j] > max)
        max = v[i][j];
    }
  }
}

float Heatmap::getAvg() {
  float ret = 0.0f;
  int cnt = 0;
  int diag = getDiagonalSize();
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = i + diag; j < size; ++j) {
      ret += v[i][j];
      cnt++;
    }
  }
  return ret / cnt;
}

float Heatmap::getAvgNearDiagonal() {
  int diag = getDiagonalSize();
  float ret = 0.0f;
  for (size_t i = 0; i + diag < size; ++i)
    ret += v[i][i + diag];
  return ret / (size - diag);
}

void Heatmap::calcAvgValues(bool count_zeros) {
  avg_value.clear();
  float ret = 0.0f;
  int cnt = 0;

  // for all distances
  for (size_t k = 0; k < size; ++k) {
    ret = 0.0f;
    cnt = 0;
    for (size_t i = 0; i + k < size; ++i) {
      if (v[i][i + k] > 1e-6 || count_zeros) {
        ret += v[i][i + k];
        cnt++;
      }
    }

    // add value - either average, or, if there are no on-zero cells, value for
    // k-1
    avg_value.push_back(cnt > 0 ? (ret / cnt)
                                : avg_value[avg_value.size() - 1]);
  }
}

vector<bool> Heatmap::getEmpty() {
  vector<bool> r;
  for (size_t i = 0; i < size; i++) {
    bool ok = false;
    for (size_t l = 0; l < size; ++l)
      if (v[i][l] > epsilon)
        ok = true;
    r.push_back(!ok);
  }
  return r;
}

Heatmap Heatmap::removeColumns(vector<bool> del) {

  int del_cnt = 0;
  for (size_t i = 0; i < del.size(); ++i)
    if (del[i])
      del_cnt++;

  Heatmap h = Heatmap(size - del_cnt);

  // printf("1\n");
  int curr_row = 0, curr_col = 0;
  for (size_t i = 0; i < size; i++) {
    if (del[i])
      continue;

    curr_col = 0;
    for (size_t l = 0; l < size; ++l) {
      if (del[l])
        continue;
      //	printf("%d %d\n", curr_row, curr_col);
      h.v[curr_row][curr_col] = v[i][l];
      curr_col++;
    }

    curr_row++;
  }

  return h;
}

Heatmap Heatmap::removeEmptyColumns() {
  bool *empty = new bool[size];
  int empty_cnt = 0;

  for (size_t i = 0; i < size; i++) {

    empty[i] = false;

    bool ok = false;
    for (size_t l = 0; l < size; ++l)
      if (v[i][l] > epsilon)
        ok = true;

    if (!ok) {
      empty[i] = true;
      empty_cnt++;
    }
  }

  //	printf("%d\n", size - empty_cnt);
  Heatmap h = Heatmap(size - empty_cnt);

  // printf("1\n");
  int curr_row = 0, curr_col = 0;
  for (size_t i = 0; i < size; i++) {
    if (empty[i])
      continue;

    curr_col = 0;
    for (size_t l = 0; l < size; ++l) {
      if (empty[l])
        continue;
      //	printf("%d %d\n", curr_row, curr_col);
      h.v[curr_row][curr_col] = v[i][l];
      curr_col++;
    }

    curr_row++;
  }

  return h;
}

void Heatmap::clear() {
  /*
  if (v == NULL) return ;
  for (int i=0; i<size; i++) {
          if (v[i] != NULL) delete v[i];
  }
  delete v;*/
}

void Heatmap::print() {
  for (size_t i = 0; i < size; ++i) {
    for (size_t l = 0; l < size; ++l) {
      printf("%f ", v[i][l]);
    }
    printf("\n");
  }
}

void Heatmap::scale(float scale) {
  for (size_t i = 0; i < size; ++i) {
    for (size_t l = 0; l < size; ++l)
      v[i][l] *= scale;
  }
}

void Heatmap::add(float val) {
  int diag = getDiagonalSize();
  for (size_t i = 0; i < size; ++i) {
    for (size_t l = i + diag; l < size; ++l) {
      v[i][l] += val;
      if (i != l)
        v[l][i] += val;
    }
  }
}

void Heatmap::add(const Heatmap &heat) {
  if (this->size != heat.size) {
    printf("heatmaps have different sizes (%lu, %lu)!\n", this->size,
           heat.size);
    exit(0);
  }
  for (size_t i = 0; i < size; ++i) {
    for (size_t l = 0; l < size; ++l) {
      this->v[i][l] += heat.v[i][l];
    }
  }
}

Heatmap Heatmap::diff(const Heatmap &heat, bool abs) {
  Heatmap h;

  if (this->size != heat.size) {
    printf("heatmaps have different sizes (%lu, %lu)!\n", this->size,
           heat.size);
    return h;
  }

  float val;

  h.init(this->size);
  for (size_t i = 0; i < size; ++i) {
    for (size_t l = 0; l < size; ++l) {
      val = this->v[i][l] - heat.v[i][l];
      h.v[i][l] = abs ? fabs(val) : val;
    }
  }
  return h;
}

void Heatmap::divide(const Heatmap &hmap) {
  if (size != hmap.size) {
    printf("heatmap divide: size mismatch (%lu vs. %lu)\n", size, hmap.size);
    return;
  }

  float d;
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      d = max(0.000001f, hmap.v[i][j]);
      v[i][j] /= d;
    }
  }
}

float Heatmap::calcDistance(const Heatmap &hmap) {
  float d = 0.0f, ret = 0.0f;
  assert(size == hmap.size);
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      d = v[i][j] - hmap.v[i][j];
      ret += d * d;
    }
  }
  return ret;
}

vector<float> Heatmap::toVector(int diag) {
  vector<float> r;
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = i + diag; j < size; ++j) {
      r.push_back(v[i][j]);
    }
  }
  return r;
}

void Heatmap::toFile(string filename, bool total_count, bool zero_diagonal,
                     bool as_integers) {
  FILE *f = open(filename, "w");

  if (total_count)
    fprintf(f, "%lu\n", size);
  float val;
  for (size_t i = 0; i < size; ++i) {
    for (size_t l = 0; l < size; ++l) {
      val = v[i][l];
      if (i == l && zero_diagonal)
        val = 0.0f;
      if (as_integers)
        fprintf(f, "%d ", (int)val);
      else
        fprintf(f, "%f ", val);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}

void Heatmap::fromFile(string filename) {
  FILE *f = open(filename, "r");

  char line[64];
  if (fgets(line, 64, f) == NULL)
    return;
  int args = countWords(line);

  fseek(f, 0, 0);

  int n;

  if (args == 1) {
    fscanf(f, "%d", &n);
  } else if (args == 3) {
    fscanf(f, "%d %d %d", &n, &start, &resolution);
  } else {
    printf("Unrecognized heatmap header: [%s]\n", line);
    exit(0);
  }

  init(n);
  for (size_t i = 0; i < size; ++i) {
    for (size_t l = 0; l < size; ++l) {
      fscanf(f, "%f", &v[i][l]);
    }
  }
  fclose(f);

  diagonal_size = getDiagonalSize();
}

void Heatmap::fromFile(string filename, bool labels) {
  FILE *f = open(filename, "r");
  if (f == NULL)
    return;

  int word_cnt = 0;
  size_t ir = 0, ic = 0;

  char c;
  char line[4096];
  char word[128];
  int p = 0;
  float val;

  while (1) {

    if (fgets(line, 4096, f) == NULL)
      break;

    // printf("\nline = %d %d %s\n", ir, ic, line);
    p = 0;
    ic = 0;

    for (size_t i = 0; i < strlen(line); ++i) {
      c = line[i];
      // printf("[%d, %c]\n", i, c);

      if (c == ' ' || c == '\t' || c == '\n') {
        if (p > 0) {
          // mamy slowo
          val = -1.0f;

          if (word[0] >= '0' && word[0] <= '9')
            val = atof(word);
          else if (word[0] == 'N')
            val = 0.0f;

          if (val > -0.5f)
            v[ir - 1][ic++] = val;

          while (p > 0)
            word[p--] = '\0'; // czyscimy stringa
          word_cnt++;
        }
      } else
        word[p++] = c;
    }

    if (p > 0) {
      val = -1.0f;
      if (word[0] >= '0' && word[0] <= '9')
        val = atof(word);
      else if (word[0] == 'N')
        val = 0.0f;

      if (val > -0.5f)
        v[ir - 1][ic++] = val;

      while (p > 0)
        word[p--] = '\0'; // czyscimy stringa
      word_cnt++;
    }

    if (ir == 0) {
      clear();
      init(word_cnt);
    }

    ir++;
    if (ir > size)
      break;
  }

  fclose(f);
}

void Heatmap::fromMDS(string filename, int size) {
  FILE *f = open(filename, "r");
  if (f == NULL)
    return;

  init(size);

  // ignore 3 first lines
  char line[64];
  for (int i = 0; i < 4; ++i)
    fgets(line, 64, f);

  float val;
  int num = (size - 1) * size / 2;
  int row = 0, col = 1;
  for (int i = 0; i < num; ++i) {
    fscanf(f, "%*d %*d %*d %*d %f", &val);
    // printf("%d %f    %d %d\n", i, val, row, col);
    v[row][col] = v[col][row] = val;
    col++;
    if (col >= size) {
      row++;
      col = row + 1;
      if (row >= size)
        break; // just to make sure
    }
  }
}
