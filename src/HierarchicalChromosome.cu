/*
 * HierarchicalChromosomeMixed.cpp
 *
 *  Created on: Jun 2, 2014
 *      Author: psz
 */

#include "../include/HierarchicalChromosome.h"
#include <iostream>
#include <optional>

#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

HierarchicalChromosome::HierarchicalChromosome() { smooth_factor = -1; }

void HierarchicalChromosome::print(int level) {
  printf("*** clusters: %d\n", (int)clusters.size());

  for (string chr : chrs) {
    printf(" [%s] root: %d, arcs: %d, curr_level size: %d\n", chr.c_str(),
           chr_root[chr], (int)arcs.arcs[chr].size(),
           (int)current_level[chr].size());
  }

  if (level >= 1) {
    printf("clusters\n");
    for (unsigned int i = 0; i < clusters.size(); ++i) {
      printf("[%u] ", i);
      clusters[i].print();
    }

    printf("current regions\n");
    for (string chr : chrs) {
      printf("%s\n", chr.c_str());
      for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
        printf("[%u] %d\n", i, current_level[chr][i]);
      }
    }

    // printf("arcs: (cnt=%d)\n", arcs.arcs.size());
    // arcs.print();
  }
}

vector3 HierarchicalChromosome::getCenter(bool current) {
  vector3 center(0.0f, 0.0f, 0.0f);
  if (current) {
    int cnt = 0;
    for (string chr : chrs) {
      for (unsigned int i = 0; i < current_level[chr].size(); i++)
        center += clusters[current_level[chr][i]].pos;
      cnt += current_level[chr].size();
    }
    center /= cnt;
  } else {
    for (unsigned int i = 0; i < clusters.size(); ++i)
      center += clusters[i].pos;
    center /= clusters.size();
  }
  return center;
}

void HierarchicalChromosome::center(bool current) {
  vector3 center = getCenter(current);

  if (current) {
    for (string chr : chrs) {
      for (unsigned int i = 0; i < current_level[chr].size(); ++i)
        clusters[current_level[chr][i]].pos -= center;
    }
  } else {
    for (unsigned int i = 0; i < clusters.size(); ++i)
      clusters[i].pos -= center;
  }
}

map<string, float> HierarchicalChromosome::getAvgAnchorDistance() {
  map<string, float> map;
  float d, tot_dist = 0.0f;
  int tot_cnt;
  for (string chr : chrs) {
    d = 0.0f;
    for (unsigned int i = 1; i < current_level[chr].size(); ++i) {
      d += (clusters[current_level[chr][i]].pos -
            clusters[current_level[chr][i - 1]].pos)
               .length();
    }
    tot_dist += d;
    tot_cnt += current_level[chr].size();

    d /= current_level[chr].size();
    map[chr] = d;
  }
  map["genome"] = tot_dist / tot_cnt;
  return map;
}

void HierarchicalChromosome::scale(float factor) {
  vector3 center = getCenter();
  for (unsigned int i = 0; i < clusters.size(); ++i) {
    clusters[i].pos = center + (clusters[i].pos - center) * factor;
  }
  createCurrentLevelStructure();
}

// int HierarchicalChromosomeMixed::findClosestCluster(int genomic_position) {
//
//	levelDown();
//	levelDown();
//
//	int n = current_level.size();
//	if (n < 2) return n-1;	// -1 if no clusters, 0 if there is only one
//	if (genomic_position < clusters[current_level[0]].base_start ||
//genomic_position > clusters[current_level[0]].base_start);
//
//
//	//printv(current_level)
//
//	int l = 0, r = n-1;
//	int mid = 0;
//	while (r>l+1) {
//		mid = (r+l) / 2;
//
//		printf("%d %d, mid=%d\n", l, r, mid);
//		if (genomic_position < mid) r = mid;
//		else l = mid;
//
//	}
// }

// void HierarchicalChromosomeMixed::align(const HierarchicalChromosomeMixed
// &hc) { 	int i, j, n = current_level.size(); 	int steps = 20; 	float ax, ay, az;
//	float score, err, best = 10e9;
//	matrix44 m, mbest;
//	map<int, int> map;	// map index from this->clusters to hc->clusters
//
//	center();
//
//
//	std::vector<vector3> pos;	// save original position
//	for (i = 0; i < n; ++i) pos.push_back(clusters[current_level[i]].pos);
//
//	for (i = 0; i < steps; ++i) {
//		ax = random(2.0f * 3.14f, false);
//		ay = random(2.0f * 3.14f, false);
//		az = random(2.0f * 3.14f, false);
//
//		// create rot matrix
//		m = RotateRadMatrix44('x', ax);
//		m = m * RotateRadMatrix44('y', ay);
//		m = m * RotateRadMatrix44('z', az);
//
//		// rotate
//		for (j = 0; j < n; ++j) clusters[current_level[i]].pos =
//clusters[current_level[i]].pos * m;
//
//		// calc score
//		score = 0.0f;
//		for (j = 0; j < n; ++j) {
//			//err = (clusters[current_level[i]].pos -
//hc.clusters[current_level[i]].pos).length();
//			//clusters[current_level[i]].pos =
//clusters[current_level[i]].pos * m;
//		}
//
//		// restore previous
//		for (j = 0; j < n; ++j) clusters[current_level[i]].pos =
//clusters[current_level[i]].pos * m;
//	}
//
//
// }

int HierarchicalChromosome::clusterIndexToSmoothedIndex(int ind) {
  if (smooth_factor == -1)
    return -1;
  return ind * (smooth_factor + 1);
}

void HierarchicalChromosome::printRegionsTree() {
  for (string chr : chrs) {
    printRegionsTree(chr);
  }
}

void HierarchicalChromosome::printRegionsTree(string chr) {
  printf("%s\n", chr.c_str());
  for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
    printRegionsTreeRecursive(current_level[chr][i]);
  }
}

void HierarchicalChromosome::printRegionsTreeRecursive(int region_index_curr,
                                                       int margin) {

  for (int i = 0; i < margin; ++i)
    printf("   ");
  clusters[region_index_curr].print(); // print current

  // print subregions
  for (unsigned int i = 0; i < clusters[region_index_curr].children.size();
       ++i) {
    printRegionsTreeRecursive(clusters[region_index_curr].children[i],
                              margin + 1);
  }
}

void HierarchicalChromosome::useTopLevel() {
  if (clusters.size() > 0) {
    current_level.clear();
    for (string chr : chrs)
      current_level[chr].push_back(chr_root[chr]);
    createCurrentLevelStructure();
  }
}

void HierarchicalChromosome::useLowestLevel() {
  for (int i = 0; i < 5; ++i)
    levelDown();
  createCurrentLevelStructure();
}

void HierarchicalChromosome::setLevel(int level) {
  useTopLevel();
  while (level--)
    levelDown();
}

void HierarchicalChromosome::levelDown() {
  std::vector<int> tmp;
  for (string chr : chrs) {

    for (unsigned int i = 0; i < current_level[chr].size(); ++i) {

      // if cluster has no children, then keep it
      if (clusters[current_level[chr][i]].children.size() == 0) {
        tmp.push_back(current_level[chr][i]);
      } else {
        for (unsigned int j = 0;
             j < clusters[current_level[chr][i]].children.size(); ++j) {
          tmp.push_back(clusters[current_level[chr][i]].children[j]);
        }
      }
    }

    current_level[chr] = tmp;
    tmp.clear();
  }
  createCurrentLevelStructure();
}

void HierarchicalChromosome::expandRegion(int start, int end,
                                          bool include_external) {
  //	std::vector<int> v;
  //
  //	for (int i = 0; i < current_level.size(); ++i) {
  //		std::vector<int> vc = expandRegion(current_level[i], start, end,
  //include_external); 		if (vc.size() > 0) v.insert(v.end(), vc.begin(),
  //vc.end()); 		else v.insert(v.end(), current_level[i]);
  //	}
  //
  //	current_level = v;
  //	createCurrentLevelStructure();
}

std::vector<int> HierarchicalChromosome::expandRegion(int region_ind, int start,
                                                      int end,
                                                      bool include_external) {
  std::vector<int> v;
  //	int p;
  //	//printf("expand region %d\n", region_ind);
  //
  //	if (end < clusters[region_ind].start || start >
  //clusters[region_ind].end) return v;  // not included
  //
  //	for (int i = 0; i < clusters[region_ind].children.size(); ++i) {
  //		p = clusters[region_ind].children[i];
  //		//printf("%d %d %d %d\n", end , clusters[p].start , start ,
  //clusters[p].end); 		if (end < clusters[p].start || start > clusters[p].end) {
  //			v.insert(v.end(), p);
  //			//printf("ngh\n");
  //			continue;
  //		}
  //
  //		if (clusters[p].children.size() == 0) {
  //			v.insert(v.end(), p);
  //			//printf("no child\n");
  //			continue;
  //		}
  //
  //		std::vector<int> vc = expandRegion(p, start, end,
  //include_external); 		if (vc.size() > 0) v.insert(v.end(), vc.begin(),
  //vc.end()); 		else v.insert(v.end(), current_level[i]);
  //	}
  //
  //	/*
  //	for (int i = 0; i < clusters[region_ind].children.size(); ++i) {
  //		p = clusters[region_ind].children[i];
  //		printf("%d %d %d %d\n", end , clusters[p].start , start ,
  //clusters[p].end); 		if (end < clusters[p].start || start > clusters[p].end) {
  //			if (include_external) v.insert(v.end(), p);
  //			printf("ig\n");
  //			continue;
  //		}
  //
  //		if (clusters[p].children.size() == 0) {
  //			v.insert(v.end(), p);
  //			printf("no child\n");
  //			continue;
  //		}
  //
  //		printf("ok [%d] ", p);
  //
  //		std::vector<int> vc = expandRegion(p, start, end,
  //include_external); 		printv(vc, true, true, "expanded children");
  //		v.insert(v.end(), vc.begin(), vc.end());
  //	}
  //	 */
  //	//printv(v, true, true, "v");
  return v;

  /*
  std::vector<int> tmp;
  for (int i = 0; i<current_level.size(); ++i) {
          if (clusters[current_level[i]].children.size() == 0) {
                  tmp.push_back(current_level[i]);
          }
          else {
                  for (int j = 0; j <
  clusters[current_level[i]].children.size(); ++j) {
                          tmp.push_back(clusters[current_level[i]].children[j]);
                  }
          }
  }

  current_level = tmp;
  createCurrentLevelStructure();
   */
}

void HierarchicalChromosome::createCurrentLevelStructure() {
  chr.clear();
  for (string c : chrs) {
    chr[c].init();
    for (unsigned int i = 0; i < current_level[c].size(); ++i) {
      chr[c].points.push_back(clusters[current_level[c][i]].pos);
      chr[c].genomic_position.push_back(
          clusters[current_level[c][i]].genomic_pos);
    }
    chr[c].updateSize();
  }
}

void HierarchicalChromosome::smoothSpline(int n) {

  createCurrentLevelStructure();

  for (string c : chrs) {
    std::vector<vector3> v;
    for (unsigned int i = 0; i < chr[c].points.size(); ++i)
      v.push_back(chr[c].points[i]);

    std::vector<vector3> vi = interpolateSpline(v, n);
    // std::vector<vector3> vi = interpolateSplineCentripetal(v, n);

    chr_smooth[c].init();
    for (unsigned int i = 0; i < vi.size(); ++i)
      chr_smooth[c].points.push_back(vi[i]);
    chr_smooth[c].updateSize();
  }

  smooth_factor = n;
}

Chromosome HierarchicalChromosome::createEqidistantModel(int resolution_bp,
                                                         string chr_ind) {
  createCurrentLevelStructure();

  Chromosome ch;

  if (chrs.size() == 0)
    error("No chromosomes in the model!");
  if (chr_ind == "")
    chr_ind = chrs[0];

  // safety checks
  if (std::find(chrs.begin(), chrs.end(), chr_ind) == chrs.end())
    error(ftext("Chromosome [%s] not found!", chr_ind.c_str()));
  int size = current_level[chr_ind].size();
  if (size == 0)
    error("Chromosome has size==0");

  if (size == 1) {
    ch.points.push_back(clusters[current_level[chr_ind][0]].pos);
    return ch;
  }

  int pos_st = clusters[current_level[chr_ind][0]].genomic_pos;
  int pos_end = clusters[current_level[chr_ind][size - 1]].genomic_pos;

  vector3 pt;
#ifdef CUMMC_SLOW_PLEASE_GO_SLOW_NEVER_SPEED_UP_I_HATE_FAST_PROGRAMS_PERFORMANCE_SUCKS
  // CPU
  for (int pos = pos_st; pos <= pos_end; pos += resolution_bp) {
    pt = find3DSmoothPosition(chr_ind, pos);
    ch.points.push_back(pt);
    ch.genomic_position.push_back(pos);
  }

#else
  // GPU
  auto res = gpuHelper(pos_st, pos_end, resolution_bp, chr_ind);
  if (res.has_value()) {
    auto pos = pos_st;
    for (auto &elem : res.value()) {
      vector3 v;
      v[0] = elem.x;
      v[1] = elem.y;
      v[2] = elem.z;
      ch.points.push_back(v);
      ch.genomic_position.push_back(pos);
      pos += resolution_bp;
    }
  }
#endif
  ch.updateSize();

  return ch;
}

void HierarchicalChromosome::toFile(string filename) {
  FILE *f = open(filename, "w");
  if (f == NULL)
    return;
  toFile(f);
  fclose(f);
}

void HierarchicalChromosome::toFile(FILE *file) {
  // fprintf(file, "%d %d %d %d\n", clusters.size(), arcs.arcs_cnt, root,
  // arcs.factors.size());
  fprintf(file, "%d %d %d %d\n", (int)clusters.size(), (int)arcs.factors.size(),
          -1, (int)chrs.size());
  for (unsigned int i = 0; i < arcs.factors.size(); ++i)
    fprintf(file, "%s ", arcs.factors[i].c_str());
  fprintf(file, "\n");
  for (unsigned int i = 0; i < clusters.size(); ++i)
    clusters[i].toFile(file);
  for (string c : chrs) {
    fprintf(file, "%s %d %d\n", c.c_str(), chr_root[c], arcs.arcs_cnt[c]);
    for (int i = 0; i < arcs.arcs_cnt[c]; ++i)
      arcs.arcs[c][i].toFile(file);
  }

  // for (int i = 0; i < arcs.arcs_cnt; ++i) {
  // arcs.arcs[i].toFile(file);
  // arcs.raw_arcs[i].toFile(file);
  // }
  // for (int i = 0; i < top_level_regions.size(); ++i) fprintf(file, "%d ",
  // top_level_regions[i]); global_chr.toFile(file);
}

void HierarchicalChromosome::toFilePreviousFormat(string filename) {
  FILE *f = open(filename, "w");
  if (f == NULL)
    return;

  string chr = chrs[0];
  fprintf(f, "%lu %lu %d %lu\n", clusters.size(), arcs.arcs[chr].size(),
          chr_root[chr], arcs.factors.size());
  // fprintf(f, "%lu %lu %d %lu\n", clusters.size(), arcs.arcs[chr].size(),
  // root_index, arcs.factors.size());
  for (unsigned int i = 0; i < arcs.factors.size(); ++i)
    fprintf(f, "%s ", arcs.factors[i].c_str());
  fprintf(f, "\n");
  for (unsigned int i = 0; i < clusters.size(); ++i)
    clusters[i].toFilePreviousFormat(f);

  int n = arcs.arcs[chr].size();
  InteractionArc *arc;
  for (int i = 0; i < n; ++i) {
    arc = &arcs.arcs[chr][i];
    fprintf(f, "%d %d %d %d %d %d 0 0\n", arc->start, arc->end, arc->score,
            arc->factor, arc->genomic_start, arc->genomic_end);
  }

  fclose(f);

  // there were no chromosomes then!
  //	for (string c: chrs) {
  //		fprintf(file, "%s %d %d\n", c.c_str(), chr_root[c],
  //arcs.arcs_cnt[c]); 		for (int i = 0; i < arcs.arcs_cnt[c]; ++i)
  //arcs.arcs[c][i].toFile(file);
  //	}

  /*
   * fprintf(file, "%d %d %d %d\n", clusters.size(), arcs.arcs_cnt, root,
  arcs.factors.size()); for (int i = 0; i < arcs.factors.size(); ++i)
  fprintf(file, "%s ", arcs.factors[i].c_str()); fprintf(file, "\n"); for (int i
  = 0; i < clusters.size(); ++i) clusters[i].toFile(file);
   */
}

void HierarchicalChromosome::fromFile(string filename) {
  FILE *f = open(filename, "r");
  if (f == NULL)
    return;

  // try to decide which format are we dealing with
  // second value now is n_factors (which is small), before it was n_arcs (which
  // is big). UPDATE 161212: the number of arcs can be small for domain models,
  // so this method is not reliable. Now checking both the number of clusters
  // and arcs (if many clusters (>100) but only a few loops (<5) -> previous
  // mode, otherwise the 2nd value is the number of factors)
  int n, n_factors, chrs_cnt, tmp;
  fscanf(f, "%d %d %d %d", &n, &n_factors, &tmp, &chrs_cnt);
  fseek(f, 0, SEEK_SET);
  if (n > 100 && n_factors < 5) {
    printf("Reading with fromFile function\n");
    this->fromFile(f);
  } else {
    printf("Reading with previousFormat function\n");
    this->fromFilePreviousFormat(f);
  }
  fclose(f);
}

void HierarchicalChromosome::fromFile(FILE *file) {

  printf("read chromosome\n");
  clusters.clear();
  current_level.clear();

  int n, n_factors, chrs_cnt, tmp;
  fscanf(file, "%d %d %d %d", &n, &n_factors, &tmp, &chrs_cnt);

  for (int i = 0; i < n_factors; ++i) {
    char str[20];
    fscanf(file, "%s", str);
    arcs.factors.push_back(str);
  }

  printf("clusters...\n");
  for (int i = 0; i < n; i++) {
    Cluster c;
    c.fromFile(file);
    clusters.push_back(c);
  }

  printf("chromosomes...\n");
  int n_arcs, root_ind;
  for (int i = 0; i < chrs_cnt; ++i) {
    char str[8];
    fscanf(file, "%s %d %d", str, &root_ind, &n_arcs);
    printf("%s...\n", str);
    chrs.push_back(str);
    chr_root[str] = root_ind;

    for (int j = 0; j < n_arcs; ++j) {
      InteractionArc arc;
      arc.fromFile(file);
      arcs.arcs[str].push_back(arc);

      if (arc.start != -1 && arc.end != -1) {
        vector_insert_unique(clusters[arc.start].arcs, j);
        vector_insert_unique(clusters[arc.end].arcs, j);
      }
    }

    arcs.arcs_cnt[str] = n_arcs;
  }

  printf("done!\n");

  /*
  for (int i = 0; i < clusters.size(); ++i) {
          for (int j = 0; j < clusters[i].children.size(); ++j) {
                  clusters[clusters[i].children[j]].parent = i;
          }
  }
   */
  useTopLevel();
}

void HierarchicalChromosome::fromFilePreviousFormat(FILE *file) {

  int total_read = 0;
  printf("old format detected\n");
  clusters.clear();
  current_level.clear();

  int n, n_arcs, n_factors, root;
  total_read += fscanf(file, "%d %d %d %d", &n, &n_arcs, &root, &n_factors);

  chrs.push_back("chr");
  chr_root["chr"] = root;

  printf("factors\n");
  for (int i = 0; i < n_factors; ++i) {
    char str[20];
    total_read += fscanf(file, "%s", str);
    arcs.factors.push_back(str);
  }

  printf("clusters\n");
  int st, end, gpos, tmp;
  float x, y, z;
  int children_cnt;
  for (int i = 0; i < n; i++) {
    Cluster c;
    // c.fromFile(file);
    total_read += fscanf(file, "%d %d %d %f %f %f %d", &gpos, &st, &end, &x, &y,
                         &z, &children_cnt);

    for (int i = 0; i < children_cnt; ++i) {
      total_read += fscanf(file, "%d", &tmp);
      c.children.push_back(tmp);
    }

    c.start = st;
    c.end = end;
    c.genomic_pos = gpos;
    c.pos.set(x, y, z);

    clusters.push_back(c);
  }

  // int st, end, score, factor;
  printf("arcs\n");
  int score_raw, factor_raw;
  for (int i = 0; i < n_arcs; ++i) {
    InteractionArc arc;
    total_read += fscanf(file, "%d %d %d %d", &arc.start, &arc.end, &arc.score,
                         &arc.factor);
    total_read += fscanf(file, "%d %d %d %d", &arc.genomic_start,
                         &arc.genomic_end, &score_raw, &factor_raw);

    arc.eff_score = arc.score;

    arcs.arcs["chr"].push_back(arc);

    if (arc.start != -1 && arc.end != -1) {
      vector_insert_unique(clusters[arc.start].arcs, i);
      vector_insert_unique(clusters[arc.end].arcs, i);
    }
  }
  arcs.arcs_cnt["chr"] = n_arcs;
  printf("READ NUMS %d\n", total_read);
}

bool HierarchicalChromosome::fromHiCEvo(string filename) {
  printf("read chromosome from HiC-evo\n");
  FILE *f = open(filename, "r");
  if (f == NULL)
    return false;

  clusters.clear();
  current_level.clear();

  arcs.clear(); // clear anchors, arcs, factors, chrs list etc.

  string chr_name = "chr";
  chrs.push_back(chr_name);
  chr_root[chr_name] = 0;
  // this->

  arcs.arcs_cnt[chr_name] = 0;

  int n;
  // fscanf(file, "%d %d %d %d", &n, &n_factors, &root_index, &chrs_cnt);
  fscanf(f, "%d", &n); // # of regions

  int tmp, start, end, par, n_child;
  float x, y, z;
  for (int i = 0; i < n; i++) {

    fscanf(f, "%d %d %d %f %f %f %d", &start, &end, &par, &x, &y, &z, &n_child);
    Cluster c(start, end);
    c.pos.x = x;
    c.pos.y = y;
    c.pos.z = z;
    c.parent = par;

    for (int j = 0; j < n_child; ++j) {
      fscanf(f, "%d", &tmp);
      c.children.push_back(tmp);
    }

    // c.print();
    clusters.push_back(c);
  }

  fclose(f);

  //	for (int i = 0; i < clusters.size(); ++i) {
  //		for (int j = 0; j < clusters[i].children.size(); ++j)
  //clusters[clusters[i].children[j]].parent = i;
  //	}

  current_level.clear();
  current_level[chr_name].push_back(0);
  return true;
}

bool HierarchicalChromosome::fromTxt(string filename) {
  printf("read chromosome from txt\n");
  FILE *f = open(filename, "r");
  if (f == NULL)
    return false;

  clusters.clear();
  current_level.clear();

  arcs.clear(); // clear anchors, arcs, factors, chrs list etc.

  string chr_name = "chrN";
  chrs.push_back(chr_name);
  chr_root[chr_name] = 0;

  arcs.arcs_cnt[chr_name] = 0;

  int n;
  float x, y, z;

  fscanf(f, "%d", &n); // # of regions

  Cluster root(1000, 1000 * n + 100);
  root.children.push_back(1);
  clusters.push_back(root);

  Cluster root_chr(1000, 1000 * n + 100);
  for (int i = 2; i <= n + 1; i++)
    root_chr.children.push_back(i);
  root_chr.parent = 0;
  clusters.push_back(root_chr);

  for (int i = 0; i < n; i++) {
    fscanf(f, "%f %f %f", &x, &y, &z);
    Cluster c(1000 * (i + 1), 1000 * (i + 1) + 100);
    c.pos.x = x;
    c.pos.y = y;
    c.pos.z = z;
    c.parent = 1;
    clusters.push_back(c);
  }

  fclose(f);

  current_level.clear();
  current_level[chr_name].push_back(0);
  return true;
}

// vector<int> HierarchicalChromosomeMixed::findFlankingAnchors(int
// genomic_position) { 	vector<int> r; 	for (int i = 0; i < current_level.size();
//++i) {
//		//		if (clusters[current_level[i]].end < start ||
//clusters[current_level[i]].start > end) continue; // outside the range
//		//
//		//
//hc.clusters.push_back(clusters[current_level[i]]);
//		//	}
//	}
//	return r;
// }

vector<int> HierarchicalChromosome::findFlankingAnchors(
    string chr, int genomic_position_start, int genomic_position_end) {
  vector<int> r;
  int p;

  int flank_left = -1, flank_right = -1;
  for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
    // printf("i=%d ", i);
    p = current_level[chr][i];
    // printf("p=%d\n", p);

    if (flank_left == -1) {
      if (clusters[p].end > genomic_position_start) {
        flank_left = p;
        r.push_back(p);
        r.push_back(i);
      }
    } else if (flank_right == -1) {
      if (clusters[p].start > genomic_position_end) {
        flank_right = p;
        r.push_back(p);
        r.push_back(i);
        break;
      }
    }
  }
  if (r.size() < 4)
    r.clear(); // something went wrong
  return r;
}

HierarchicalChromosome HierarchicalChromosome::extractFragment(int start,
                                                               int end) {

  // TODO: update to genome
  HierarchicalChromosome hc;

  //	useLowestLevel();
  //
  //	for (int i = 0; i < current_level.size(); ++i) {
  //		if (clusters[current_level[i]].end < start ||
  //clusters[current_level[i]].start > end) continue; // outside the range
  //
  //		hc.clusters.push_back(clusters[current_level[i]]);
  //	}
  //
  //	Cluster r(hc.clusters[0].start, hc.clusters[hc.clusters.size()-1].end);
  //	for (int i = 0; i < hc.clusters.size(); ++i) {
  //		r.children.push_back(i);
  //		hc.clusters[i].parent = hc.clusters.size();
  //		hc.clusters[i].level = 1;
  //	}
  //
  //
  //	// update arcs
  //	for (int i = 0; i < arcs.raw_arcs.size(); ++i) {
  //
  //		//printf("raw %d %d\n", raw_arcs[i].start, raw_arcs[i].end);
  //		int st = -1, end = -1;
  //		for (int j = 0; j < hc.clusters.size() && (st==-1||end==-1); ++j)
  //{ 			if (hc.clusters[j].contains(arcs.raw_arcs[i].start)) st = j; 			if
  //(hc.clusters[j].contains(arcs.raw_arcs[i].end)) end = j;
  //		}
  //
  //		if (st != -1 && end != -1) {
  //			//printf("new arc: %d %d\n", st, end);
  //			hc.arcs.raw_arcs.push_back(arcs.raw_arcs[i]);
  //
  //			InteractionArc arc(st, end, arcs.raw_arcs[i].score,
  //arcs.raw_arcs[i].factor); 			hc.arcs.arcs.push_back(arc);
  //		}
  //	}
  //	hc.arcs.arcs_cnt = hc.arcs.arcs.size();
  //	hc.arcs.factors = arcs.factors;
  //
  //	// add root
  //	hc.root = hc.clusters.size();
  //	hc.clusters.push_back(r);

  return hc;
}

Heatmap HierarchicalChromosome::getDistancesHeatmap() {

  // calc size
  int n = 0;
  for (string chr : chrs)
    n += current_level[chr].size();

  Heatmap h(n);
  int x = 0, y;
  for (string chr_id : chrs) {
    for (unsigned int i = 0; i < current_level[chr_id].size(); ++i) {
      y = 0;
      for (string chr_id2 : chrs) {
        for (unsigned int j = 0; j < current_level[chr_id2].size(); ++j) {
          h.v[x][y] = (clusters[current_level[chr_id][i]].pos -
                       clusters[current_level[chr_id2][j]].pos)
                          .length();
          y++;
        }
      }

      x++;
    }
  }
  return h;
}

// calculate some basic stats for chromosomes at different scales, eg. mean &
// max distance, chromosome diameter etc.
void HierarchicalChromosome::getSpatialDistributionStats() {
  // calculate, at different scales (chromosomal, subchromosomal and cluster):
  // * diameter of chromosomes
  // * distances between chromosomes

  int i, j;
  float d;
  float min, max, avg;
  for (int lvl = 0; lvl <= 0; lvl++) {
    setLevel(lvl);
    printf("\n* level %d\n", lvl);

    printf("diameters:\n");
    for (string chr_id : chrs) {
      d = chr[chr_id].getDiameter();
      printf("%s: %f\n", chr_id.c_str(), d);
    }

    printf("distance from nucleus center:\n");
    vector3 center = getCenter(true);
    for (string chr_id : chrs) {
      d = 0.0f;
      for (unsigned int i = 0; i < current_level[chr_id].size(); ++i) {
        d += (clusters[current_level[chr_id][i]].pos - center).length();
      }
      d /= current_level[chr_id].size();
      printf("%s: %f\n", chr_id.c_str(), d);
    }

    printf("\ninter chromosomal distances:\n");
    for (string chr_id : chrs) {
      for (string chr_id2 : chrs) {
        min = 1e9;
        max = 0.0f;
        avg = 0.0f;

        for (unsigned int i = 0; i < current_level[chr_id].size(); ++i) {
          for (unsigned int j = 0; j < current_level[chr_id2].size(); ++j) {
            d = (clusters[current_level[chr_id][i]].pos -
                 clusters[current_level[chr_id2][j]].pos)
                    .length();
            if (d < min)
              min = d;
            if (d > max)
              max = d;
            avg += d;
          }
        }

        avg /= (float)(current_level[chr_id].size() *
                       current_level[chr_id2].size());
        printf("%f ", avg);
      }
      printf("\n");
    }
  }
}

// create heatmap with pairwise distances between beads for a given chromosome
// level denotes on which hierarchy level we calculate distances:
//   0 - chromosome (ie. calc distances between chromosomal-beads)
//   1 - segment
//   2 - subanchor
// 'chr' is ignored for level==0, otherwise it denote the chromosome we
// calculate the distances we assume that the structure is at the correct level
Heatmap HierarchicalChromosome::createStructuralHeatmap(string chr, int level) {
  Heatmap h;
  float d = 0.0f;
  if (level == 0) {
    int n = chrs.size();
    h.init(n);
    // i,j - indices of chromosomes
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        d = (clusters[current_level[chrs[i]][0]].pos -
             clusters[current_level[chrs[j]][0]].pos)
                .length();
        h.v[i][j] = h.v[j][i] = d;
      }
    }
  } else if (level == 1) {
    int n = current_level[chr].size();
    h.init(n);
    // i,j - indices of segment beads
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        d = (clusters[current_level[chr][i]].pos -
             clusters[current_level[chr][j]].pos)
                .length();
        h.v[i][j] = h.v[j][i] = d;
      }
    }
  } else if (level == 2) {
    // only consider subanchor beads within the same segment
    //// first, create a map subanchor_id -> segment_id
    //// then for each pair

    // find the total heatmap size
    int n = 0;
    for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
      n += clusters[current_level[chr][i]].children.size();
    }

    h.init(n);

    // now create the heatmap, segment by segment

    int shift = 0;
    for (unsigned int i = 0; i < current_level[chr].size(); ++i) {
      int m = clusters[current_level[chr][i]]
                  .children.size(); // size of the current segment
      // for each segment iterate over pairs of subanchor beads
      for (int j = 0; j < m; ++j) {
        int ind1 = clusters[current_level[chr][i]].children[j];
        for (int k = j + 1; k < m; ++k) {
          int ind2 = clusters[current_level[chr][i]].children[k];
          d = (clusters[ind1].pos - clusters[ind2].pos).length();
          h.v[shift + j][shift + k] = h.v[shift + k][shift + j] = d;
        }
      }
      shift += m;
    }
  }
  return h;
}

float HierarchicalChromosome::calcDistance(HierarchicalChromosome hc,
                                           int level) {
  float ret = 0.0f;
  if (level == 0) {
    // level==0 means chromosome level.
  } else {

    for (string c : chrs) {
    }
  }
  return ret;
}

vector3 HierarchicalChromosome::find3DPosition(string chr_ind, int pos) {
  // printf("find %s %d\n", chr_ind.c_str(), pos);
  if (std::find(chrs.begin(), chrs.end(), chr_ind) == chrs.end()) {
    printf("error: no such chromosome! [%s]\n", chr_ind.c_str());
    return vector3(0.0, 0.0, 0.0);
  }
  int size = current_level[chr_ind].size();
  if (size == 0) {
    printf("error: current level has size 0!\n");
    return vector3(0.0f, 0.0f, 0.0f);
  }

  int p;

  // border cases ('pos' is before the first or after the last bead)
  // we can extend it to check if 'pos' is before the *end* of the first, or
  // after the *start* of the end
  p = current_level[chr_ind][0];
  if (pos < clusters[p].end)
    return clusters[p].pos;
  p = current_level[chr_ind][size - 1];
  if (pos >= clusters[p].start)
    return clusters[p].pos;

  int prev_p = -1;
  for (int i = 0; i < size; ++i) {
    p = current_level[chr_ind][i]; // index of i-th cluster
    // check if 'pos' is contained within
    if (pos >= clusters[p].start && pos <= clusters[p].end)
      return clusters[p].pos;

    // if 'pos' is within 2 beads interpolate between them
    if (pos < clusters[p].start) {
      float st = (pos - clusters[prev_p].end) /
                 (float)(clusters[p].start - clusters[prev_p].end);
      return interpolate(clusters[prev_p].pos, clusters[p].pos, st);
    }

    prev_p = p;
  }
  return vector3(
      0.0f, 0.0f,
      0.0f); // It will never happen - just a compilation warninng suppression
}

void parallelFind3DSmoothPositionHelper(const int chr_ind, const int pos_st,
                                        const int pos_end,
                                        const int resolution_bp) {
  // Prepare data in buffers
  // Run
  // Copy data back to output vector
}

vector3 HierarchicalChromosome::find3DSmoothPosition(string chr_ind, int pos) {

  // safety checks
  if (std::find(chrs.begin(), chrs.end(), chr_ind) == chrs.end())
    return vector3();

  int size = current_level[chr_ind].size();
  if (size == 0)
    return vector3();

  int p;

  // border cases ('pos' is before the first or after the last bead)
  // we can extend it to check if 'pos' is before the *end* of the first, or
  // after the *start* of the end
  p = current_level[chr_ind][0];
  if (pos < clusters[p].end)
    return clusters[p].pos;
  p = current_level[chr_ind][size - 1];
  if (pos >= clusters[p].start)
    return clusters[p].pos;

  int prev_p = -1;
  for (int i = 0; i < size; ++i) {
    p = current_level[chr_ind][i]; // index of i-th cluster

    // check if 'pos' is contained within
    // if (pos >= clusters[p].start && pos <= clusters[p].end) return
    // clusters[p].pos;

    // if 'pos' is within 2 beads interpolate between them
    if (pos < clusters[p].genomic_pos) {

      // p1, p2 - points we are currently between which
      // p0, p3 - flanking points
      vector3 p1 = clusters[prev_p].pos;
      vector3 p2 = clusters[p].pos;
      vector3 p0 = i > 1 ? clusters[current_level[chr_ind][i - 2]].pos
                         : mirrorPoint(p1, p2);
      vector3 p3 = i + 1 < current_level[chr_ind].size()
                       ? clusters[current_level[chr_ind][i + 1]].pos
                       : mirrorPoint(p2, p1);

      float t = (pos - clusters[prev_p].genomic_pos) /
                (float)(clusters[p].genomic_pos - clusters[prev_p].genomic_pos);

      // printf("%d (%d %d) p=%d t=%f \n", pos, clusters[prev_p].genomic_pos,
      // clusters[p].genomic_pos, p, t);
      vector3 r = interpolateSpline(t, p0, p1, p2, p3);
      return r;

      // float st = (pos - clusters[prev_p].end) / (float)(clusters[p].start -
      // clusters[prev_p].end); return interpolate(clusters[prev_p].pos,
      // clusters[p].pos, st);
    }

    prev_p = p;
  }
  return vector3(); // It will never happen - just a compilation warninng
                    // suppression
}

__global__ void parallelFind3DSmoothPosition(
    const int pos_start, const int pos_stop, const int resolution_bp,
    const size_t size, const int _p1, const int _p2, const int p1_end,
    const int p2_start, const int *current_levels,
    const size_t current_levels_size, const int *clusters_genomic_positions,
    const float3 *clusters_positions, float3 *out);

std::optional<std::vector<float3>>
HierarchicalChromosome::gpuHelper(const int pos_start, const int pos_stop,
                                  const int resolution_bp, string chr_ind) {
  std::optional<std::vector<float3>> ret;
  // safety checks
  if (std::find(chrs.begin(), chrs.end(), chr_ind) == chrs.end())
    return ret;
  int size = current_level[chr_ind].size();
  if (size == 0)
    return ret;
  // figure something out with p
  const auto p1 = current_level[chr_ind][0];
  const auto p2 = current_level[chr_ind][size - 1];

  const auto p1_end = clusters[p1].end;
  const auto p2_start = clusters[p2].start;

  // move to gpu
  std::vector<int> current_levels_h =
      current_level[chr_ind]; // index of i-th cluster
  for (const int e : current_level[chr_ind]) {
    if (e < 0) {
      std::cout << " UJEMNE DANE " << std::endl;
    }
  }
  const auto current_levels_size = current_levels_h.size();

  float3 *clusters_positions, *out;
  int *clusters_genomic_positions, *current_levels;
  cudaMalloc(&clusters_positions, clusters.size() * sizeof(float3));
  cudaMalloc(&clusters_genomic_positions, clusters.size() * sizeof(int));
  cudaMalloc((int **)&current_levels, current_levels_h.size() * sizeof(int));
  const auto out_size = ((pos_stop - pos_start + 1) / resolution_bp + 1);
  cudaMalloc(&out, out_size * sizeof(float3));

  std::vector<float3> clusters_positions_h(clusters.size());
  std::vector<int> clusters_genomic_positions_h(clusters.size());

  for (auto i = 0; i < clusters.size(); ++i) {
    clusters_genomic_positions_h[i] = clusters[i].genomic_pos;
    clusters_positions_h[i].x = clusters[i].pos[0];
    clusters_positions_h[i].y = clusters[i].pos[1];
    clusters_positions_h[i].z = clusters[i].pos[2];
  }

  cudaMemcpy(clusters_positions, &clusters_positions_h[0],
             clusters.size() * sizeof(float3), cudaMemcpyHostToDevice);
  cudaMemcpy(clusters_genomic_positions, &clusters_genomic_positions_h[0],
             clusters.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(current_levels, current_levels_h.data(),
             current_levels_h.size() * sizeof(int), cudaMemcpyHostToDevice);

  gpuErrchk(cudaDeviceSynchronize());

  parallelFind3DSmoothPosition<<<1024, 100>>>(
      pos_start, pos_stop, resolution_bp, size, p1, p2, p1_end, p2_start,
      current_levels, current_levels_size, clusters_genomic_positions,
      clusters_positions, out

  );

  gpuErrchk(cudaPeekAtLastError());

  gpuErrchk(cudaDeviceSynchronize());

  // copy results back
  ret = std::vector<float3>(out_size);
  cudaMemcpy(ret.value().data(), out, out_size * sizeof(float3),
             cudaMemcpyDeviceToHost);
  cudaFree(clusters_positions);
  cudaFree(clusters_genomic_positions);
  cudaFree((void *)out);
  cudaFree((void *)current_levels);

  return ret;
}

__forceinline__ float3 __device__ mirrorPointDevice(const float3 &a,
                                                    float3 &b) {
  return make_float3(a.x * 2.0f - b.x, a.y * 2.0f - b.y, a.z * 2.0f - b.z);
}

__forceinline__ __device__ float3 operator+(const float3 &a, const float3 &b) {

  return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__forceinline__ __device__ float3 operator-(const float3 &a, const float3 &b) {

  return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__forceinline__ __device__ float3 operator*(const float &a, const float3 &b) {

  return make_float3(a * b.x, a * b.y, a * b.z);
}

__forceinline__ __device__ float3 operator*(const float3 &a, const float &b) {

  return b * a;
}

__forceinline__ __device__ float3 interpolateSplineDevice(const float t,
                                                          const float3 &p1,
                                                          const float3 &p2,
                                                          const float3 &p3,
                                                          const float3 &p4) {
  const float t2 = t * t;
  const float t3 = t2 * t;
  const float b1 = .5f * (-t3 + 2.0f * t2 - t);
  const float b2 = .5f * (3.0f * t3 - 5.0f * t2 + 2.0f);
  const float b3 = .5f * (-3.0f * t3 + 4.0f * t2 + t);
  const float b4 = .5f * (t3 - t2);
  return (p1 * b1 + p2 * b2 + p3 * b3 + p4 * b4);
}

__global__ void parallelFind3DSmoothPosition(
    const int pos_start, const int pos_stop, const int resolution_bp,
    const size_t size, const int _p1, const int _p2, const int p1_end,
    const int p2_start, const int *current_levels,
    const size_t current_levels_size, const int *clusters_genomic_positions,
    const float3 *clusters_positions, float3 *out) {

  auto idx = threadIdx.x + blockIdx.x * blockDim.x;
  const auto gridSize = gridDim.x * blockDim.x;

  for (long long pos = pos_start + resolution_bp * idx; pos <= pos_stop;
       pos += gridSize * resolution_bp) {

    if (pos < p1_end) {
      out[idx] = clusters_positions[_p1];
      return;
    }
    if (pos >= p2_start) {
      out[idx] = clusters_positions[_p2];
      return;
    }

    int prev_p = -1;
    for (size_t i = 0; i < current_levels_size; ++i) {
      int p = current_levels[i]; // index of i-th cluster
      if (pos < clusters_genomic_positions[p]) {
        float3 p1 = clusters_positions[prev_p];
        float3 p2 = clusters_positions[p];
        float3 p0 = i > 1 ? clusters_positions[current_levels[i - 2]]
                          : mirrorPointDevice(p1, p2);
        float3 p3 = i + 1 < current_levels_size
                        ? clusters_positions[current_levels[i + 1]]
                        : mirrorPointDevice(p2, p1);
        const float t = (pos - clusters_genomic_positions[prev_p]) /
                        (float)(clusters_genomic_positions[p] -
                                clusters_genomic_positions[prev_p]);
        out[idx] = interpolateSplineDevice(t, p0, p1, p2, p3);
        return;
      }
      prev_p = p;
    }

    idx += gridSize;
  }
}
