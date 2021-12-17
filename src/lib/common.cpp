/*
 * common.cpp
 *
 *  Created on: Aug 4, 2013
 *      Author: psz
 */

#include "common.h"

int random(int range) { return rand() % range; }

float random_uniform() { return rand() / (float)RAND_MAX; }

float random(float range, bool negative) {
  if (negative)
    return (2.0f * random_uniform() - 1.0f) * range;
  return range * random_uniform();
}

vector3 random_vector(float max_size, bool in2D) {
  if (in2D)
    return vector3(random(max_size, true), random(max_size, true), 0.0);
  return vector3(random(max_size, true), random(max_size, true),
                 random(max_size, true));
}

matrix44 random_rot_matrix(float max_angle) {
  vector3 axis = random_vector(1.0f);
  axis = axis.normalize();
  matrix44 rot = RotateRadMatrix44(axis, DegToRad(random(max_angle)));
  return rot;
}

float dist(const vector3 &v1, const vector3 &v2) {
  return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) +
              (v1.z - v2.z) * (v1.z - v2.z));
}

float angle(vector3 v1, vector3 v2) {
  v1.normalize();
  v2.normalize();
  float dot = DotProduct(v1, v2);
  // assert(dot >= -1.0 && dot <= 1.0);
  return 1.0f - (dot + 1.0f) / 2.0f;
  // return acos((double)dot);
}

float angle_norm(vector3 v1, vector3 v2) {
  float dot = DotProduct(v1, v2);
  // assert(dot >= -1.0 && dot <= 1.0);
  return 1.0f - (dot + 1.0f) / 2.0f;
  // return acos((double)dot);
}

vector<int> interpolateInterval(int posl, int posr, int dl, int dr, float mn) {
  vector<int> v;
  //	int min_dist = 100;

  while (1) {
    int d = posr - posl;
    // printf("%d (%d) - %d (%d)  %d\n", posl, dl, posr, dr, d);

    if (2 * dl >= d && 2 * dr >= d) {
      v.push_back((posl + posr) /
                  2); // the hole is small now, we can just put finishing point
      break;
    }

    if (dl < dr) {
      dl *= mn;
      // if (dl >= min_dist) {
      posl += dl;
      v.push_back(posl);
      //}
    } else {
      dr *= mn;
      // if (dr >= min_dist) {
      posr -= dr;
      v.push_back(posr);
      //}
    }
  }

  sort(v.begin(), v.end());
  return v;
}

vector<int> interpolatePoints(vector<int> points, float mn) {
  if (points.size() < 3)
    return points;
  vector<int> v, tmp;
  int dl, dr;

  for (size_t i = 0; i < points.size() - 1; ++i) {
    dl = (i > 0) ? points[i] - points[i - 1] : points[i + 1] - points[i];
    dr = (i + 2 < points.size()) ? points[i + 2] - points[i + 1]
                                 : points[i + 1] - points[i];
    tmp = interpolateInterval(points[i], points[i + 1], dl, dr, mn);

    v.push_back(points[i]);
    for (size_t j = 0; j < tmp.size(); ++j)
      v.push_back(tmp[j]);
  }
  v.push_back(points[points.size() - 1]);
  return v;
}

/*
vector<float> interpolateInterval(float posl, float posr, float dl, float dr,
float mn) { vector<float> v; float cl = 0.0f; while (1) {

                float d = posr - posl;
                printf("%f (%f) - %f (%f)  %f\n", posl, dl, posr, dr, d);

                if (2.0f * dl >= d && 2.0f * dr >= d) {
                        v.push_back((posl+posr)/2.0f); // the hole is small now,
we can just put finishing point break;
                }

                if (dl < dr) {
                        dl *= mn;
                        posl += dl;
                        v.push_back(posl);
                }
                else {
                        dr *= mn;
                        posr -= dr;
                        v.push_back(posr);
                }
        }

        sort(v.begin(), v.end());
        return v;
}

vector<float> interpolatePoints(vector<float> points, float mn) {
        if (points.size() < 3) return points;
        vector<float> v, tmp;
        float dl, dr;


        for (int i = 0; i < points.size()-1; ++i) {
                dl = (i > 0) ? points[i] - points[i-1] : points[i+1] -
points[i]; dr = (i+2 < points.size()) ? points[i+2] - points[i+1] : points[i+1]
- points[i]; tmp = interpolateInterval(points[i], points[i+1], dl, dr, mn);

                v.push_back(points[i]);
                for (int j = 0; j < tmp.size(); ++j) v.push_back(tmp[j]);
        }
        v.push_back(points[points.size()-1]);
        return v;
}
 */

vector3 displace(vector3 v, float displacement, bool in2D) {
  if (in2D)
    return vector3(v.x + random(displacement, true),
                   v.y + random(displacement, true), v.z);
  return vector3(v.x + random(displacement, true),
                 v.y + random(displacement, true),
                 v.z + random(displacement, true));
}

void displace_ref(vector3 &v, float displacement) {
  v.x += random(displacement, true);
  v.y += random(displacement, true);
  v.z += random(displacement, true);
}

float distanceToInteractionFrequency(float distance) {
  return 1.0f / pow(distance, 1.0f);
}

float interpolate(float a, float b, float p) { return a + (b - a) * p; }

vector3 interpolate(vector3 a, vector3 b, float p) { return a + (b - a) * p; }

void random_shuffle(float *arr, int size) {
  int p;
  float tmp;
  for (int j = 0; j < 10; ++j) {
    for (int i = 0; i < size; ++i) {
      p = random(size);
      tmp = arr[i];
      arr[i] = arr[p];
      arr[p] = tmp;
    }
  }
}
//
// char* ftext(const char *fmt, ...) {
//	char *text = new char[256];
//	if (fmt == NULL) return "";
//	va_list ap;
//	va_start(ap, fmt);
//	vsprintf(text, fmt, ap);
//	va_end(ap);
//	return text;
//}

// char* concat(const char *a, const char *b) {
//	int len = strlen(a) + strlen(b) + 1;
//	printf("%d\n", len);
//	char *text = new char[len];
//	strcpy(text, a);
//	strcat(text, b);
//	text[len-1] = '\0';
//	return text;
// }

void printv(vector<bool> v) {
  for (size_t i = 0; i < v.size(); ++i) {
    printf("%c", v[i] ? 'T' : '.');
  }
  printf("\n");
}

void printv(const vector<int> &v, bool count, bool newline, const char *label) {
  if (strlen(label) > 0)
    printf("%s: ", label);
  if (count)
    printf("[total: %lu] ", v.size());
  for (size_t i = 0; i < v.size(); ++i)
    printf(" %d", v[i]);
  if (newline)
    printf("\n");
}

void printv(const vector<float> &v, bool count, bool newline,
            const char *label) {
  if (strlen(label) > 0)
    printf("%s: ", label);
  if (count)
    printf("[total: %lu] ", v.size());
  for (size_t i = 0; i < v.size(); ++i)
    printf("%f ", v[i]);
  if (newline)
    printf("\n");
}

void printv(const vector<string> &v, bool count, bool newline,
            const char *label) {
  if (strlen(label) > 0)
    printf("%s: ", label);
  if (count)
    printf("[total: %lu] ", v.size());
  for (size_t i = 0; i < v.size(); ++i)
    printf("%s ", v[i].c_str());
  if (newline)
    printf("\n");
}

void printm(const std::map<int, int> &map, const char *label) {
  if (strlen(label) > 0)
    printf("%s:\n", label);
  for (auto el : map) {
    printf("[%d] -> %d\n", el.first, el.second);
  }
}

void printm(const std::map<int, vector<int>> &map, const char *label) {
  if (strlen(label) > 0)
    printf("%s:\n", label);
  for (auto el : map) {
    printf("[%d] -> ", el.first);
    for (int val : el.second)
      printf("%d, ", val);
    printf("\n");
  }
}

void printm(const std::map<int, vector<float>> &map, const char *label) {
  if (strlen(label) > 0)
    printf("%s:\n", label);
  for (auto el : map) {
    printf("[%d] -> ", el.first);
    for (float val : el.second)
      printf("%f, ", val);
    printf("\n");
  }
}

void printm(const std::map<std::string, int> &map, bool flat, bool newline,
            const char *label) {
  if (strlen(label) > 0)
    printf("%s: ", label);
  for (auto el : map) {
    if (flat)
      printf(" %d", el.second);
    else
      printf("[%s] -> %d\n", el.first.c_str(), el.second);
  }
  if (newline)
    printf("\n");
}

void printm(std::map<std::string, int> &map,
            const std::vector<std::string> &keys, bool flat, bool newline,
            const char *label) {
  if (strlen(label) > 0)
    printf("%s: ", label);
  for (string key : keys) {
    if (flat)
      printf(" %d", map[key]);
    else
      printf("[%s] -> %d\n", key.c_str(), map[key]);
  }
  if (newline)
    printf("\n");
}

void printm(const std::map<std::string, std::vector<int>> &map, bool flat,
            const char *label) {
  if (strlen(label) > 0)
    printf("%s:\n", label);
  if (flat)
    for (auto el : map) {
      printv(el.second, false, false);
      printf("\n");
    }
  else
    for (auto el : map)
      printv(el.second, true, true, el.first.c_str());
}

bool withChance(float chance) { return random_uniform() < chance; }

bool file_exists(const std::string &name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  }
  return false;
}

FILE *open(string path, string mode) {
  FILE *f = NULL;
  f = fopen(path.c_str(), mode.c_str());
  if (f == NULL) {
    printf("Error opening file [%s]!\n", path.c_str());
    return NULL;
  }
  return f;
}

void saveToCSV(string path, vector<int> v, bool add_index_column) {
  FILE *f = open(path.c_str(), "w");
  if (!f)
    return;

  for (size_t i = 0; i < v.size(); ++i) {
    if (add_index_column)
      fprintf(f, "%lu\t%d\n", i, v[i]);
    else
      fprintf(f, "%d\n", v[i]);
  }
  fclose(f);
}

void saveToCSV(string path, vector<float> v, bool add_index_column) {
  FILE *f = open(path.c_str(), "w");
  if (!f)
    return;

  for (size_t i = 0; i < v.size(); ++i) {
    if (add_index_column)
      fprintf(f, "%lu\t%f\n", i, v[i]);
    else
      fprintf(f, "%f\n", v[i]);
  }
  fclose(f);
}

void saveToCSV(string path, vector<bool> v, bool add_index_column) {
  FILE *f = open(path.c_str(), "w");
  if (!f)
    return;

  for (size_t i = 0; i < v.size(); ++i) {
    if (add_index_column)
      fprintf(f, "%lu\t%d\n", i, v[i] ? 1 : 0);
    else
      fprintf(f, "%d\n", v[i] ? 1 : 0);
  }
  fclose(f);
}

void saveToCSV(string path, vector<vector<float>> v, bool add_index_column) {
  FILE *f = open(path.c_str(), "w");
  if (!f)
    return;

  if (v.size() == 0)
    return;

  for (size_t j = 0; j < v[0].size(); ++j) { // rows
    if (add_index_column)
      fprintf(f, "%lu\t", j);
    for (size_t i = 0; i < v.size() - 1; ++i) {
      fprintf(f, "%f\t", v[i][j]);
    }
    fprintf(f, "%f\n", v[v.size() - 1][j]);
  }
  fclose(f);
}

void saveToCSV(string path, vector<vector<int>> v, bool add_index_column) {
  FILE *f = open(path.c_str(), "w");
  if (!f)
    return;

  if (v.size() == 0)
    return;

  for (size_t j = 0; j < v[0].size(); ++j) { // rows
    if (add_index_column)
      fprintf(f, "%lu\t", j);
    for (size_t i = 0; i < v.size() - 1; ++i) {
      fprintf(f, "%d\t", v[i][j]);
    }
    fprintf(f, "%d\n", v[v.size() - 1][j]);
  }
  fclose(f);
}

int countWords(char *line) {
  int c = 1, p = 0;
  while (line[p] != '\0') {
    if (line[p] == ' ' || line[p] == '\t')
      c++;
    p++;
  }
  return c;
}

bool vector_contains(const std::vector<int> &v, int val) {
  return std::find(v.begin(), v.end(), val) != v.end();
  // for (int i = 0; i < v.size(); ++i) if (v[i] == val) return true;
  // return false;
}

void vector_insert_unique(std::vector<int> &v, int val) {
  if (!vector_contains(v, val))
    v.push_back(val);
}

void error(string msg) {
  printf("\n%s\n", msg.c_str());
  exit(0);
}

// Catmull–Rom interpolating spline
vector3 interpolateSpline(float t, const vector3 &p1, const vector3 &p2,
                          const vector3 &p3, const vector3 &p4) {
  float t2 = t * t;
  float t3 = t2 * t;
  float b1 = .5 * (-t3 + 2 * t2 - t);
  float b2 = .5 * (3 * t3 - 5 * t2 + 2);
  float b3 = .5 * (-3 * t3 + 4 * t2 + t);
  float b4 = .5 * (t3 - t2);
  return (p1 * b1 + p2 * b2 + p3 * b3 + p4 * b4);
}

std::vector<vector3> interpolateSpline(const std::vector<vector3> &points,
                                       int splits) {
  std::vector<vector3> v;
  vector3 pt;
  float t, dt;
  int k, n = points.size();

  vector3 pt_begin = points[0] - (points[1] - points[0]);
  vector3 pt_end = points[n - 1] + (points[n - 1] - points[n - 2]);

  // dt = 1.0f / (splits+1);
  dt = 1.0f / (splits);
  for (int i = 0; i + 1 < n; ++i) {
    t = 0.0f;
    // printf("change points\n");
    for (t = 0.0f, k = 0; t < 1.0f && k <= splits; t += dt, k++) {

      if (i == 0)
        pt = interpolateSpline(t, pt_begin, points[0], points[1], points[2]);
      else if (i == n - 2)
        pt = interpolateSpline(t, points[i - 1], points[i], points[i + 1],
                               pt_end);
      else
        pt = interpolateSpline(t, points[i - 1], points[i], points[i + 1],
                               points[i + 2]);

      // printf("i=%d t=%f, (%f %f %f)\n", i, t, pt.x, pt.y, pt.z);
      v.push_back(pt);
    }
  }
  return v;
}

////////////////////////////////////
////////////////////////////////////
////////////////////////////////////

std::vector<vector3> interpolateSplineCentripetal(std::vector<vector3> &points,
                                                  int splits) {
  std::vector<vector3> v;
  vector3 pt;
  int n = points.size();

  vector3 pt_begin = points[0] - (points[1] - points[0]);
  vector3 pt_end = points[n - 1] + (points[n - 1] - points[n - 2]);

  points.insert(points.begin(), pt_begin);
  points.insert(points.end(), pt_end);

  double time[4];

  for (int i = 0; i + 3 < n; ++i) {
    // printf("change points\n");

    for (int j = 0; j < 4; j++)
      time[j] = j;

    double tstart = 1;
    double tend = 2;

    // if (!curveType.equals(CatmullRomType.Uniform)) {
    double total = 0;
    for (int j = 1; j < 4; j++) {
      double tt = (points[i + j] - points[i + j - 1]).lengthSqr();
      total += pow(tt, 0.25);
      time[j] = total;
    }
    tstart = time[1];
    tend = time[2];

    int segments = splits - 1;
    v.push_back(points[i + 1]);
    for (int j = 1; j < segments; j++) {
      vector3 rr =
          interpolate(points[i], points[i + 1], points[i + 2], points[i + 3],
                      time, tstart + (j * (tend - tstart)) / segments);
      v.push_back(rr);
    }
  }

  points.erase(points.end() - 1);
  points.erase(points.begin());
  return v;
}

vector3 interpolate(vector3 &p0, vector3 &p1, vector3 &p2, vector3 &p3,
                    double time[], double t) {
  vector3 L01 = p0 * (time[1] - t) / (time[1] - time[0]) +
                p1 * (t - time[0]) / (time[1] - time[0]);
  vector3 L12 = p1 * (time[2] - t) / (time[2] - time[1]) +
                p2 * (t - time[1]) / (time[2] - time[1]);
  vector3 L23 = p2 * (time[3] - t) / (time[3] - time[2]) +
                p3 * (t - time[2]) / (time[3] - time[2]);
  vector3 L012 = L01 * (time[2] - t) / (time[2] - time[0]) +
                 L12 * (t - time[0]) / (time[2] - time[0]);
  vector3 L123 = L12 * (time[3] - t) / (time[3] - time[1]) +
                 L23 * (t - time[1]) / (time[3] - time[1]);
  vector3 C12 = L012 * (time[2] - t) / (time[2] - time[1]) +
                L123 * (t - time[1]) / (time[2] - time[1]);
  return C12;
}

double interpolate(double p[], double time[], double t) {
  double L01 = p[0] * (time[1] - t) / (time[1] - time[0]) +
               p[1] * (t - time[0]) / (time[1] - time[0]);
  double L12 = p[1] * (time[2] - t) / (time[2] - time[1]) +
               p[2] * (t - time[1]) / (time[2] - time[1]);
  double L23 = p[2] * (time[3] - t) / (time[3] - time[2]) +
               p[3] * (t - time[2]) / (time[3] - time[2]);
  double L012 = L01 * (time[2] - t) / (time[2] - time[0]) +
                L12 * (t - time[0]) / (time[2] - time[0]);
  double L123 = L12 * (time[3] - t) / (time[3] - time[1]) +
                L23 * (t - time[1]) / (time[3] - time[1]);
  double C12 = L012 * (time[2] - t) / (time[2] - time[1]) +
               L123 * (t - time[1]) / (time[2] - time[1]);
  return C12;
}

vector3 mirrorPoint(vector3 fixed, vector3 pt) { return 2.0f * fixed - pt; }

void print_vector(const vector3 &v, const char *label) {
  if (strcmp(label, "") == 0)
    printf("%f %f %f", v.x, v.y, v.z);
  else
    printf("%s: (%f %f %f)\n", label, v.x, v.y, v.z);
}

string ftext(const string fmt_str, ...) {
  int final_n, n = ((int)fmt_str.size()) *
                   2; /* reserve 2 times as much as the length of the fmt_str */
  string str;
  std::unique_ptr<char[]> formatted;
  va_list ap;
  while (1) {
    formatted.reset(
        new char[n]); /* wrap the plain char array into the unique_ptr */
    strcpy(&formatted[0], fmt_str.c_str());
    va_start(ap, fmt_str);
    final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
    va_end(ap);
    if (final_n < 0 || final_n >= n)
      n += abs(final_n - n + 1);
    else
      break;
  }
  return string(formatted.get());
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

// attempts to split string by an unknown delimiter
std::vector<std::string> split_unknown(const std::string &s) {
  const string delims = " \t;,";
  std::vector<std::string> elems;
  for (size_t i = 0; i < delims.size(); i++) {
    split(s, delims[i], elems);
    if (elems.size() > 1)
      break;
    else
      elems.clear();
  }
  return elems;
}

void split_file_path(const std::string &s, std::string &dir,
                     std::string &filename) {
  size_t split_pos = s.rfind("/");
  if (split_pos == string::npos)
    error("couldn't split file path");
  dir = s.substr(0, split_pos + 1);
  filename = s.substr(split_pos + 1);
}
