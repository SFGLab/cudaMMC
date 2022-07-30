#include <Chromosome.h>

Chromosome::Chromosome() {
  size = 0;
  score = 0.0f;
}

void Chromosome::init() {
  points.clear();
  genomic_position.clear();
  size = 0;
}

void Chromosome::print() {
  for (int i = 0; i < size; i++)
    printf("(%f %f %f)\n", points[i].x, points[i].y, points[i].z);
  printf("\n");
}

void Chromosome::updateSize() { size = points.size(); }

void Chromosome::scale(float scale, bool center) {
  if (center)
    this->center();
  for (int i = 0; i < size; ++i) {
    points[i] *= scale;
  }
}

void Chromosome::fromFile(string filename) {
  FILE *f = open(filename, "r");
  if (f == NULL)
    return;

  int pts_cnt = 0;

  char line[32];
  if (fgets(line, 32, f) == NULL)
    return;
  fseek(f, 0, 0);

  int args = countWords(line);
  if (args == 1) {
    fscanf(f, "%d", &pts_cnt);
    this->fromFile(f, pts_cnt);
  } else if (args == 3) {
    this->fromFile(f, -1);
  } else {
    fclose(f);
    error("Unrecognized format. First line should contain the size or 3D "
          "coordinates");
  }

  fclose(f);
}

float Chromosome::findBestAlignmentRotation(const Chromosome &chr, int steps,
                                            vector3 &angle, float max_angle) {
  float step = max_angle / (steps - 1);
  float d_ang = max_angle / -2.0f;

  matrix44 m;
  float min_dist = 1e12f, dist;
  for (int ix = 0; ix < steps; ++ix) {
    for (int iy = 0; iy < steps; ++iy) {
      for (int iz = 0; iz < steps; ++iz) {
        m = RotateRadMatrix44('x', step * ix - d_ang);
        m = m * RotateRadMatrix44('y', step * iy - d_ang);
        m = m * RotateRadMatrix44('z', step * iz - d_ang);

        Chromosome ch = clone();
        ch.rotate(m);

        dist = ch.getDistanceSqr(chr);
        // printf("%d %d %d %f\n", ix, iy, iz, dist);
        if (dist < min_dist) {
          // printf("  ! %f %f      %f %d %f\n", dist, min_dist, step, ix,
          // step*ix);
          min_dist = dist;
          angle.x = step * ix;
          angle.y = step * iy;
          angle.z = step * iz;
        }
      }
    }
  }

  return min_dist;
}

void Chromosome::adjustSize(const Chromosome &chr) {

  float score;
  Chromosome chr_tmp;

  float low = 1e-6, high = 1e2, dscale;
  float sc, prev_sc;
  bool sign_change;

  while (high - low > 1e-6) {

    if (low < 0.0f)
      low = 1e-6;
    dscale = (high - low) / 7;

    sign_change = false;
    score = chr_tmp.calcDistance(chr);
    prev_sc = score;
    for (int i = 0; i < 8; ++i) {
      sc = low + dscale * i;
      chr_tmp = this->clone();
      chr_tmp.scale(sc);
      score = chr_tmp.calcDistance(chr);
      if (score > prev_sc) { // start to raise
        // modify range
        low = low + dscale * (i - 2);
        high = sc;
        sign_change = true;
        break;
      }
      prev_sc = score;
    }

    if (!sign_change) {
      // values are decreasing, we need to expand range to the right
      low = high - dscale;
      high = high + dscale;
    }
  }

  scale(low);
}

void Chromosome::align(const Chromosome &chr, bool resize, float max_angle) {
  int steps = 50;

  // center();

  vector3 angle(0.0f, 0.0f, 0.0f);

  if (resize)
    adjustSize(chr);

  findBestAlignmentRotation(chr, steps, angle, max_angle);

  matrix44 m = RotateRadMatrix44('x', angle.x);
  m = m * RotateRadMatrix44('y', angle.y);
  m = m * RotateRadMatrix44('z', angle.z);
  rotate(m);
}

float Chromosome::getDistanceSqr(const Chromosome &chr) {
  float dist = 0.0f;
  for (int i = 0; i < size; i++) {
    dist += (points[i].x - chr.points[i].x) * (points[i].x - chr.points[i].x);
    dist += (points[i].y - chr.points[i].y) * (points[i].y - chr.points[i].y);
    dist += (points[i].z - chr.points[i].z) * (points[i].z - chr.points[i].z);
  }
  return dist / size;
}

void Chromosome::translate(const vector3 &v) {
  for (int i = 0; i < size; i++)
    points[i] += v;
}

void Chromosome::rotate(const matrix44 &mat) {
  for (int i = 0; i < size; i++)
    points[i] = points[i] * mat;
}

vector3 Chromosome::getCenter() {
  vector3 center(0.0f, 0.0f, 0.0f);
  for (int i = 0; i < size; i++)
    center += points[i];
  return center / size;
}

void Chromosome::center() {
  vector3 center = getCenter();
  print_vector(center, "center");
  for (int i = 0; i < size; i++)
    points[i] -= center;
}

Chromosome Chromosome::clone() {
  Chromosome ch;
  ch.points = points;
  ch.points_hierarchical = points_hierarchical;
  ch.score = score;
  ch.size = size;
  return ch;
}

Heatmap Chromosome::createHeatmap() {
  Heatmap h(size);
  float d;
  for (int i = 0; i < size; ++i) {
    for (int l = i + 1; l < size; ++l) {
      d = dist(points[i], points[l]);
      h.v[i][l] = d;
      h.v[l][i] = h.v[i][l];
    }
  }
  return h;
}

Heatmap Chromosome::createInverseHeatmap() {
  Heatmap h(size);
  float d;
  for (int i = 0; i < size; ++i) {
    for (int l = i + 1; l < size; ++l) {
      d = dist(points[i], points[l]);
      // h.v[i][l] = 1.0f / pow(d, Settings::;
      h.v[i][l] = 1.0f / d;
      h.v[l][i] = h.v[i][l];
    }
  }
  return h;
}

void Chromosome::createRandom(int pts_cnt, float size, bool walk) {
  vector3 v(0.0, 0.0, 0.0);
  points.clear();
  for (int i = 0; i < pts_cnt; i++) {
    points.push_back(v);
    if (walk)
      v = displace(v, size);
    else
      v = random_vector(size);
  }
  this->size = pts_cnt;
}

void Chromosome::makeLine(int pts_cnt, float dist) {
  vector3 v(0.0f, 0.0f, 0.0f);
  for (int i = 0; i < pts_cnt; ++i) {
    v.x += dist;
    points.push_back(v);
  }
  center();
  size = pts_cnt;
}

void Chromosome::makeSpiral(int pts_cnt, float r, float angle,
                            float spin_height) {

  init();
  vector3 v(0.0f, 0.0f, 0.0f);

  float alpha = 0.0f;

  for (int i = 0; i < pts_cnt; i++) {
    v.x = r * cos(alpha);
    v.z = r * sin(alpha);
    v.y += spin_height;
    alpha += angle;
    points.push_back(v);
  }
  size = pts_cnt;
}

void Chromosome::makeSpiral(int pts_cnt, float r, float angle, vector3 p1,
                            vector3 p2) {

  // we do as in 'normal', axis-oriented spring, but we use vectors
  // perpendicular to p1-p2 axis

  vector3 dir = p2 - p1; // direction of vector
  vector3 axis1(-dir.y, dir.x, 0.0f);
  if (fabs(axis1.x) < 1e-6 && fabs(axis1.y) < 1e-6)
    axis1.x = -dir.z;
  axis1 = Normalized(axis1);
  vector3 axis2 = Normalized(CrossProduct(Normalized(dir), axis1));

  // print_vector(dir, "dir");
  // print_vector(axis1, "a1");
  // print_vector(axis2, "a2");

  init();

  vector3 v(0.0f, 0.0f, 0.0f);

  float alpha = 0.0f;
  float t = 0.0f;
  float dt = 1.0f / (pts_cnt - 1);

  // correction
  vector3 d1 = r * cos(0.0f) * axis1 + r * sin(0.0f) * axis2;
  float ang = alpha + angle * (pts_cnt - 1);
  vector3 d2 = r * cos(ang) * axis1 + r * sin(ang) * axis2;

  // print_vector(d1, "corr1");
  // print_vector(d2, "corr2");

  for (int i = 0; i < pts_cnt; i++) {
    v = p1 + dir * t + r * cos(alpha) * axis1 + r * sin(alpha) * axis2;
    v = v - (1.0f - t) * d1 - t * d2;
    // print_vector(v, " ");
    points.push_back(v);

    alpha += angle;
    t += dt;
  }
  size = pts_cnt;
}

float Chromosome::getDiameter() const {
  float diameter = 0.0f, d;
  for (int i = 0; i < size; ++i) {
    for (int l = i + 1; l < size; ++l) {
      d = (points[i] - points[l]).lengthSqr();
      diameter = max(diameter, d);
    }
  }
  return sqrt(diameter);
}

float Chromosome::calcDistance(const Chromosome &chr) {
  float d1, d2, dist = 0.0f;
  for (int i = 0; i < size; ++i) {
    for (int l = i + 1; l < size; ++l) {
      d1 = (points[i] - points[l]).length();
      d2 = (chr.points[i] - chr.points[l]).length();
      // printf("= %f %f %f\n", d1, d2, fabs(d1- d2));
      dist += fabs(d1 - d2);
    }
  }
  return (dist / ((size * (size - 1) / 2)));
}

float Chromosome::calcRMSD(const Chromosome &chr) {

  double str1[size][3], str2[size][3];
  double res = -1.0;

  for (int i = 0; i < size; ++i) {
    str1[i][0] = points[i].x;
    str1[i][1] = points[i].y;
    str1[i][2] = points[i].z;
    str2[i][0] = chr.points[i].x;
    str2[i][1] = chr.points[i].y;
    str2[i][2] = chr.points[i].z;
  }

  fast_rmsd(str1, str2, size, &res);
  // printf("rmsd: %lf\n", res);
  // fast_rmsd(double ref_xlist[][3], double mov_xlist[][3], int n_list, double*
  // rmsd);
  return res;
}

vector<vector<float>> Chromosome::calcBaseDensity(float sphere_r) {
  // we approximate linear structure by additional points added to the structure
  // beetwen every two beads we add N points

  vector<vector3> pts;
  int N = 5, pts_cnt;
  float dist;
  vector3 curr, diff;
  vector<int> intersection;
  bool is_contained = false;
  vector<vector<float>> out;
  vector<float> tmp;
  out.push_back(tmp);
  out.push_back(tmp);
  out.push_back(tmp);

  if (sphere_r < 1e-6)
    sphere_r = this->getDiameter() / 5.0f;

  // create a list of points
  for (int i = 0; i < size - 1; ++i) {
    pts.push_back(points[i]);
    diff = (points[i + 1] - points[i]) / (N + 1);
    curr = points[i];
    for (int j = 0; j < N; ++j) {
      curr += diff;
      pts.push_back(curr);
    }
  }
  pts.push_back(points[size - 1]); // add last point

  sphere_r = sphere_r * sphere_r; // for optimization, we calculate distance
                                  // squared, so we need to square radius too
  pts_cnt = pts.size();
  for (int i = 0; i < pts_cnt; i++) {

    is_contained = false;
    if (i % 100 == 0)
      printf("%d/%d\n", i, size);
    printf("%f %f %f\n", pts[i].x, pts[i].y, pts[i].z);
    intersection.clear();
    for (int j = 0; j < size; ++j) {
      dist = (pts[i] - pts[j]).lengthSqr();
      if ((dist <= sphere_r && !is_contained) ||
          (dist > sphere_r && is_contained)) {
        is_contained = !is_contained;
        intersection.push_back(j);
      }
    }

    // we need to "close" the intersection, if the end point is contained in the
    // last sphere
    if (is_contained)
      intersection.push_back(size - 1);

    int base_density = 0, base_looping = 0;

    // intersections come in pairs (a[k], b[k]).
    // for density, we sum b[k]-a[k] for k such that i belong to (a[k], b[k])
    // for looping, we sum intervals for all other k
    for (size_t j = 0; j < intersection.size(); j += 2) {
      if (intersection[j] <= i && intersection[j + 1] >= i)
        base_density += intersection[j + 1] - intersection[j] + 1;
      else
        base_looping += intersection[j + 1] - intersection[j] + 1;
    }

    // printf("%d %d %d   %d\n", i, base_density, base_looping,
    // intersection.size());
    out[0].push_back(base_density);
    out[1].push_back(base_looping);
    out[2].push_back(base_density + base_looping);
  }
  return out;
}

int Chromosome::getMostDistantPointsPair(float &dist) {
  float max = 0.0f, curr;
  int maxi = 0;
  for (int i = 0; i < size - 1; ++i) {
    curr = (points[i] - points[i + 1]).lengthSqr();
    if (curr > max) {
      max = curr;
      maxi = i;
    }
  }
  dist = sqrt(max);
  return maxi;
}

std::vector<float> Chromosome::getConsecutiveBeadDistances() {
  std::vector<float> v;
  for (int i = 0; i < size - 1; ++i)
    v.push_back((points[i] - points[i + 1]).length());
  return v;
}

float Chromosome::getAvgConsecutiveBeadDistance() {
  float v = 0.0f;
  for (int i = 0; i < size - 1; ++i) {
    v += (points[i] - points[i + 1]).length();
  }
  return v / (size - 1);
}

void Chromosome::trim(int start, int end) {
  if (size >= start) {
    if (end == 0)
      points.erase(points.begin() + start, points.end());
    else
      points.erase(points.begin() + start, points.begin() + end);
  }
}

void Chromosome::toFile(string filename) {
  FILE *f;
  f = open(filename, "w");
  if (f == NULL)
    return;

  if (size > 0) {
    fprintf(f, "%d\n", size);
    this->toFile(f);
  } else {
    fprintf(f, "0");
  }

  fclose(f);
}

void Chromosome::toFile(FILE *file) {
  for (int i = 0; i < size; i++) {
    if (genomic_position.size() > static_cast<std::size_t>(i))
      fprintf(file, "%f %f %f %d\n", points[i].x, points[i].y, points[i].z,
              genomic_position[i]);
    else
      fprintf(file, "%f %f %f\n", points[i].x, points[i].y, points[i].z);
  }
}

void Chromosome::fromFilePDB(string filename) {
  FILE *f = open(filename, "r");
  if (f == NULL)
    return;

  int i;
  init();

  char line[4096];
  char word[8];

  while (1) {

    if (fgets(line, 4096, f) == NULL)
      break;

    if (strncmp(line, "ATOM", 4) != 0)
      continue;
    // printf("line = [%s]", line);

    float x, y, z;
    for (i = 0; i < 5; ++i)
      word[i] = line[33 + i];
    x = atof(word);
    for (i = 0; i < 5; ++i)
      word[i] = line[41 + i];
    y = atof(word);
    for (i = 0; i < 5; ++i)
      word[i] = line[49 + i];
    z = atof(word);

    points.push_back(vector3(x, y, z));
  }

  fclose(f);
}

void Chromosome::fromFile(FILE *file, int pts_cnt) {

  init();
  if (pts_cnt < 0) {
    // unknown length
    double x, y, z;
    while (!feof(file)) {
      if (fscanf(file, "%lf %lf %lf", &x, &y, &z) == 3)
        points.push_back(vector3(x, y, z));
    }
    size = points.size();
  } else {
    if (pts_cnt == 0)
      fscanf(file, "%d", &pts_cnt);

    float x, y, z;
    for (int i = 0; i < pts_cnt; i++) {
      fscanf(file, "%f %f %f", &x, &y, &z);
      points.push_back(vector3(x, y, z));
    }
    size = pts_cnt;
  }
}
