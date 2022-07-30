#include <InteractionArcs.h>

InteractionArcs::InteractionArcs() {
  selected_region.chr = "";
  selected_region.start = 0;
  selected_region.end = 0; // this means that no region is selected (ie. we work
                           // on whole chromosomes)
}

void InteractionArcs::clear() {
  anchors.clear();
  arcs.clear();
  raw_arcs.clear();

  factors.clear();
  chrs.clear();
  arcs_cnt.clear();
  anchors_cnt.clear();

  expected_distances.clear();
}

void InteractionArcs::selectRegion(BedRegion region) {
  selected_region = region;
}

// we have anchors[], a sorted list of anchors, and raw_arcs[], a list of arcs
// (using genomic positions) we want to fill arc[] so that it contains a list of
// arcs with values referring to indices of 'anchors', and not genomic positions
void InteractionArcs::markArcs(bool ignore_missing) {

  // for every chromosome:
  // go through all raw arcs, find the corresponding anchor, create new arc
  // (anchor based) and add it to list

  int cnt = 0;
  int last_start = -1;
  std::map<int, std::vector<InteractionArc>> tmp_arcs;

  // normally we want to print warning about mismatching arcs (ones for which
  // anchors are missing), but if we have selected a region (either by providing
  // specific region or by limiting region in debug mode) then we want to
  // surpress them (because we can expect that there are going to be some
  // mismatches anyway)
  // TODO: we do filtering when reading arcs, so maybe there should be no
  // mismatches? bool mismatched_arcs_as_errors = ignore_missing ||
  // (selected_region.end == 0);

  for (string chr : chrs) {
    printf(" %s...\n", chr.c_str());

    arcs[chr].clear(); // we may run markArcs() multiple times. make sure we
                       // won't duplicate arcs

    std::sort(raw_arcs[chr].begin(), raw_arcs[chr].end());

    cnt = raw_arcs[chr].size();

    // we need i==cnt to process the remaining arcs
    for (int i = 0; i <= cnt; ++i) { // for every arc

      int st = -1, end = -1;
      if (i < cnt) {
        for (int j = 0; j < anchors_cnt[chr]; ++j) {
          if (anchors[chr][j].length() > 1) {
            if (anchors[chr][j].contains(raw_arcs[chr][i].start))
              st = j;
            if (anchors[chr][j].contains(raw_arcs[chr][i].end))
              end = j;
          }
        }

        if ((st == -1 || end == -1)) {
          printf("! error: non-matching arc\n");
          raw_arcs[chr][i].print();
          continue;
        }

        if (st == end) {
          // printf("st==end %d %d\n", st, end);
          // raw_arcs[chr][i].print();
          continue; // ignore looping arcs
        }
      }

      if (st != last_start || i == cnt) {

        // add all cached arcs
        for (auto el : tmp_arcs) {
          // here 'el.second' is a vector with arcs having common start and end
          // if there is only one arc, then we can simply add it to the list
          if (el.second.size() == 1) {
            arcs[chr].push_back(el.second[0]);
          } else {

            // sort arcs (so that they are ordered by factor)
            std::sort(el.second.begin(), el.second.end());

            // check how many different factors there are
            bool multiple_factors = false;
            for (size_t j = 1; j < el.second.size(); ++j)
              if (el.second[j].factor != el.second[j - 1].factor)
                multiple_factors = true;

            int total_score = 0;
            int factor_score = 0;
            int first_of_factor = 0;
            for (size_t j = 0; j <= el.second.size(); ++j) {

              // if factor is changing update the arc
              if (j == el.second.size() ||
                  (j > 0 && el.second[j].factor != el.second[j - 1].factor)) {
                el.second[first_of_factor].score = factor_score;
                el.second[first_of_factor].eff_score =
                    multiple_factors ? 0 : factor_score;
                arcs[chr].push_back(el.second[first_of_factor]);

                first_of_factor = j;
                total_score += factor_score;
                factor_score = 0;
              }

              if (j < el.second.size())
                factor_score += el.second[j].score;
            }

            // if we had multiple factors create a single, summary arc
            if (multiple_factors) {
              InteractionArc arc(el.second[0].start, el.first, 0, -1);
              arc.eff_score = total_score;
              arcs[chr].push_back(arc);
            }
          }

          el.second.clear();
        }
        tmp_arcs.clear();
        last_start = st;
      }

      // we will gather all arcs starting at anchor 'st', and process them
      // together (we do that because we need to merge arcs between the same
      // anchors but with different factors)
      InteractionArc arc(st, end, raw_arcs[chr][i].score,
                         raw_arcs[chr][i].factor);
      arc.genomic_start = raw_arcs[chr][i].start;
      arc.genomic_end = raw_arcs[chr][i].end;
      tmp_arcs[end].push_back(arc);
    }

    arcs_cnt[chr] = arcs[chr].size(); // update count
    tmp_arcs.clear();
  }
}

void InteractionArcs::removeEmptyAnchors() {
  printf("removing empty anchors\n");

  for (auto el : anchors) {
    std::vector<Anchor> correct_anchors;
    std::set<int> correct_anchor_positions;

    for (int i = 0; i < arcs_cnt[el.first]; ++i) {
      correct_anchor_positions.insert(arcs[el.first][i].start);
      correct_anchor_positions.insert(arcs[el.first][i].end);
    }

    for (const auto &pos : correct_anchor_positions) {
      correct_anchors.push_back(anchors[el.first][pos]);
    }

    this->anchors[el.first].clear();

    for (Anchor a : correct_anchors)
      this->anchors[el.first].push_back(a);

    int removed_cnt = anchors_cnt[el.first] - anchors[el.first].size();
    anchors_cnt[el.first] -= removed_cnt; // update number of anchors
    printf("   %s %d\n", el.first.c_str(), removed_cnt);
  }
}

void InteractionArcs::rewire() {
  // for each anchor find the total strength of anchored interactions (ie. the
  // sum of PET counts of arcs ending at the anchor) probabilistically remove
  // anchors based on this strength

  // printf("rewire!\n");
  for (auto el : anchors) { // iterate over chromosomes
    for (int i = anchors_cnt[el.first] - 1; i >= 0;
         i--) { // iterate over anchors
      // printf("anchor %d\n", i);
      //  find strength of interactions
      int sum = 0;
      for (int j = 0; j < arcs_cnt[el.first]; j++) {
        if (arcs[el.first][j].start == i || arcs[el.first][j].end == i) {
          sum += arcs[el.first][j].eff_score;
          // printf("   arc %d %d   %d %d\n", arcs[el.first][j].start,
          // arcs[el.first][j].end, arcs[el.first][j].eff_score,
          // arcs[el.first][j].score);
        }
      }
      // printf("total: %d\n", sum);

      if (withChance((float)sum / Settings::rewiringCertaintyThreshold)) {
        anchors[el.first].erase(anchors[el.first].begin() + i);
      }
    }

    anchors_cnt[el.first] = anchors[el.first].size();
  }
}

// void InteractionArcs::densify(int density, int margin) {
//	int range, d, p, cnt;
//	int min_d = 50;
//	std::vector<Anchor> new_anchors;
//
//	if (margin < 0) margin = density / 2;
//
//	range = anchors[1].center - anchors[0].center;
//	d = range / (density+1);
//	p = anchors[0].center - margin * d;
//	for (int i = 0; i < margin; ++i, p+=d) {
//		Anchor tmp("", p, p);
//		new_anchors.push_back(tmp);
//	}
//
//	for (int i = 1; i < anchors_cnt; ++i) {
//		range = anchors[i].center - anchors[i-1].center;
//		d = range / (density+1); // desired delta
//		cnt = density;
//		if (d < min_d) {
//			d = range / (1 + (range / min_d));
//			cnt = range / min_d;
//		}
//
//		new_anchors.push_back(anchors[i-1]);
//		for (int j=0, p=anchors[i-1].center+d; j<cnt && p <
//anchors[i].center; j++, p+=d) {
//			//printf("%d (%d %d) %d\n", i, anchors[i-1].center,
//anchors[i].center, d); 			Anchor tmp("", p, p); 			new_anchors.push_back(tmp);
//		}
//	}
//	new_anchors.push_back(anchors[anchors_cnt-1]);
//
//	range = anchors[anchors_cnt-1].center - anchors[anchors_cnt-2].center;
//	d = range / (density+1);
//	p = anchors[anchors_cnt-1].center + d;
//	for (int i = 0; i < margin; ++i, p+=d) {
//		Anchor tmp("", p, p);
//		new_anchors.push_back(tmp);
//	}
//
//	anchors = new_anchors;
//	anchors_cnt = anchors.size();
//	markArcs();
// }

// void InteractionArcs::printSummary() {
//
//	for (int i = 1; i < anchors.size(); ++i) {
//		printf("len=%d, dist=%d\n", anchors[i].length(),
//anchors[i].center-anchors[i-1].center);
//	}
//
//	//for (int i = 1; i < points.size() && anchor_start.size(); ++i) {
//	//	printf("len=%d, dist=%d\n", anchor_end[i]-anchor_start[i],
//anchor_start[i]-anchor_end[i-1]);
//	//}
// }

// note: anchor file must be sorted
void InteractionArcs::loadAnchorsData(string anchors_path) {

  printf("read anchors [%s]\n", anchors_path.c_str());

  FILE *f = open(anchors_path.c_str(), "r");
  if (f == NULL)
    exit(0);

  int posa, posb;
  char chr_tmp[10];
  char orientation = 'N';
  string chr;

  // create a set with all selected chromosomes (to speed up reading data)
  std::set<std::string> chrs_set;
  for (std::string chr : chrs)
    chrs_set.insert(chr);

  bool is_region_selected = selected_region.end > 0;

  // recognize format. allowed formats are "chr start end" and "chr start end
  // orientation"
  bool motif_info = false;

  char line[128];
  if (fgets(line, 64, f) != NULL) {
    int cnt = countWords(line);
    if (cnt != 3 && cnt != 4)
      error("unrecognized format of anchor file");
    if (cnt == 4) {
      motif_info = true;
      printf("motif orientation info detected\n");
    }
  }
  fseek(f, 0, SEEK_SET); // move cursor back to the file beginning

  // read anchor file line by line
  for (int i = 0; !feof(f); ++i) {

    int cols = 0; // number of values (columns) read from the file
    if (motif_info)
      cols = fscanf(f, "%s %d %d %c", chr_tmp, &posa, &posb, &orientation);
    else
      cols = fscanf(f, "%s %d %d", chr_tmp, &posa, &posb);

    if (cols < 3)
      continue; // probably just an empty line at the end

    chr = string(chr_tmp);

    // check if chromosome is ok
    // then check if the anchor fits into the selected region (at least one end
    // must be within)
    if (chrs_set.find(chr) == chrs_set.end())
      continue;
    if (is_region_selected &&
        (!selected_region.contains(posa) && !selected_region.contains(posb)))
      continue;

    Anchor anchor(chr, posa, posb, orientation);
    anchors[chr].push_back(anchor);
  }
  fclose(f);

  printf("anchors loaded:\n");
  for (auto el : anchors) {
    anchors_cnt[el.first] = el.second.size();
    printf("   %s: %d\n", el.first.c_str(), anchors_cnt[el.first]);
  }
}

// Reads clusters from file
// predefined_segments - bed regions (usually defining segments); clusters
// overlapping their ends will be ignored. Used to ensure there will be gaps
// where custom split is
void InteractionArcs::loadPetClustersData(string pet_clusters_path,
                                          string factor_name,
                                          BedRegions predefined_segments) {
  printf("read clusters [%s] (factor = %s)\n", pet_clusters_path.c_str(),
         factor_name.c_str());

  FILE *f = open(pet_clusters_path.c_str(), "r");
  if (f == NULL)
    exit(0);

  int ast, aend, bst, bend;
  int posa, posb;
  float sc;
  char buf_tmp[1024];
  char chr_a[10], chr_b[10];

  int added = 0;            // how many arcs were included
  int long_arcs_cnt = 0;    // how many arcs were filtered out as too long
  int arcs_over_splits = 0; // how many arcs were dropped because they stretched
                            // over the custom split sites

  int factor = factors.size();
  factors.push_back(factor_name);

  bool is_region_selected = selected_region.end > 0;

  // create a set with all selected chromosomes (to speed up reading data)
  std::set<std::string> chrs_set;
  for (std::string chr : chrs)
    chrs_set.insert(chr);

  for (int i = 1; !feof(f); ++i) {
    if (i % 1000000 == 0)
      printf(".");
    if (fscanf(f, "%s %d %d %s %d %d %f", chr_a, &ast, &aend, chr_b, &bst,
               &bend, &sc) < 7)
      continue;
    fgets(buf_tmp, 1024, f); // read to the end of line

    if (strcmp(chr_a, chr_b) != 0)
      continue; // only intra contacts

    // check if the arc belongs to the selected chromosomes
    if (chrs_set.find(chr_a) == chrs_set.end())
      continue;

    posa = (ast + aend) / 2;
    posb = (bst + bend) / 2;
    if (posa > posb) {
      int tmp = posa;
      posa = posb;
      posb = tmp;
    }

    if (is_region_selected &&
        (!selected_region.contains(posa) || !selected_region.contains(posb)))
      continue;

    // printf("[%s] %d %d [%s] %d %d %f\n", chr_a, ast, aend, chr_b, bst, bend,
    // sc);

    // if custom split is used, then ignore interactions spanning over the split
    // sites
    bool ok = true;
    for (std::size_t j = 0; j < predefined_segments.regions.size(); ++j) {
      if (strcmp(predefined_segments.regions[j].chr.c_str(), chr_a) == 0) {
        if ((posa <= predefined_segments.regions[j].start &&
             posb >= predefined_segments.regions[j].start) ||
            (posa <= predefined_segments.regions[j].end &&
             posb >= predefined_segments.regions[j].end)) {
          ok = false;
          break;
        }
      }
    }
    if (!ok) {
      arcs_over_splits++;
      continue;
    }
    //

    InteractionArc arc(posa, posb, sc, factor);

    // long arcs are used to refine the segment level heatmap
    if (posb - posa > Settings::maxPETClusterLength) {
      long_arcs_cnt++;
      long_arcs[chr_a].push_back(arc);
      continue;
    }

    // find a proper place for an arc (finding gaps requires arcs to be ordered)
    // (arcs are sorted in file by beginning, and we take center, so the last
    // few arcs can be unordered) (note that there are also different factors)
    int p = raw_arcs[chr_a].size();
    while (p > 0 && raw_arcs[chr_a][p - 1].start > arc.start)
      p--;

    raw_arcs[chr_a].insert(raw_arcs[chr_a].begin() + p, arc);
    added++;
  }
  fclose(f);

  printf("added %d, discarded long arcs: %d, over custom split: %d\n", added,
         long_arcs_cnt, arcs_over_splits);
}

void InteractionArcs::toFile(string filename) {
  FILE *f = open(filename, "w");
  if (f == NULL)
    return;

  int i;
  InteractionArc *arc;
  fprintf(f, "%d %d\n", (int)chrs.size(), (int)factors.size());

  for (string factor : factors)
    fprintf(f, "%s ", factor.c_str());
  fprintf(f, "\n");

  for (string chr : chrs) {
    fprintf(f, "%s %d %d %d\n", chr.c_str(), anchors_cnt[chr], arcs_cnt[chr],
            (int)long_arcs[chr].size());
    for (i = 0; i < anchors_cnt[chr]; ++i) {
      fprintf(f, "%d %d %c\n", anchors[chr][i].start, anchors[chr][i].end,
              anchors[chr][i].orientation);
    }

    for (int i = 0; i < arcs_cnt[chr]; ++i) {
      arc = &arcs[chr][i];
      fprintf(f, "%d %d %d %d %d %d %d\n", arc->start, arc->end,
              arc->genomic_start, arc->genomic_end, arc->score, arc->eff_score,
              arc->factor);
    }

    for (size_t i = 0; i < long_arcs[chr].size(); ++i) {
      arc = &long_arcs[chr][i];
      fprintf(f, "%d %d %d %d %d %d %d\n", arc->start, arc->end,
              arc->genomic_start, arc->genomic_end, arc->score, arc->eff_score,
              arc->factor);
    }
  }
  fclose(f);
}

bool InteractionArcs::fromFile(string filename) {

  FILE *f = open(filename, "r");
  if (f == NULL)
    return false;

  int i, j, cnt, cnt_arcs, cnt_anchors, cnt_factors, long_arcs_cnt;
  int st, end, st_g, end_g, sc, eff_sc, fact;
  char chr_tmp[10];
  char orientation;
  string chr;

  if (fscanf(f, "%d %d", &cnt, &cnt_factors) != 2)
    return false;
  for (i = 0; i < cnt_factors; i++) {
    fscanf(f, "%s", chr_tmp);
    factors.push_back(chr_tmp);
  }
  printf("chr cnt: %d\n", cnt);
  for (i = 0; i < cnt; ++i) {
    fscanf(f, "%s %d %d %d", chr_tmp, &cnt_anchors, &cnt_arcs, &long_arcs_cnt);

    chr = string(chr_tmp);
    printf("[%s] %d %d %d\n", chr.c_str(), cnt_anchors, cnt_arcs,
           long_arcs_cnt);

    chrs.push_back(chr);
    anchors_cnt[chr] = cnt_anchors;
    arcs_cnt[chr] = cnt_arcs;

    for (j = 0; j < cnt_anchors; ++j) {
      fscanf(f, "%d %d %c", &st, &end, &orientation);
      Anchor a(chr, st, end, orientation);
      anchors[chr].push_back(a);
    }

    for (j = 0; j < cnt_arcs; ++j) {
      fscanf(f, "%d %d %d %d %d %d %d", &st, &end, &st_g, &end_g, &sc, &eff_sc,
             &fact);
      InteractionArc arc(st, end, sc, fact);
      arc.genomic_start = st_g;
      arc.genomic_end = end_g;
      arc.eff_score = eff_sc;
      // arc.anchor_start = anch_st;
      // arc.anchor_end = anch_end;
      arcs[chr].push_back(arc);
    }

    for (j = 0; j < long_arcs_cnt; ++j) {
      fscanf(f, "%d %d %d %d %d %d %d", &st, &end, &st_g, &end_g, &sc, &eff_sc,
             &fact);
      InteractionArc arc(st, end, sc, fact);
      arc.genomic_start = st_g;
      arc.genomic_end = end_g;
      arc.eff_score = eff_sc;
      long_arcs[chr].push_back(arc);
    }
  }
  fclose(f);
  return true;
}

void InteractionArcs::print(int display_limit) {

  int i;
  int display;
  Anchor *anc;
  InteractionArc *arc;
  for (auto el : anchors) {

    printf(" [%s] anchors: %d, arcs: %d\n", el.first.c_str(),
           anchors_cnt[el.first], arcs_cnt[el.first]);

    printf(" [%s] anchors %d\n", el.first.c_str(), anchors_cnt[el.first]);
    display = min(anchors_cnt[el.first], display_limit);
    for (i = 0; i < display; ++i) {
      anc = &anchors[el.first][i];
      printf("   [%d] %d %d\n", i, anc->start, anc->end);
    }

    printf(" [%s] arcs %d\n", el.first.c_str(), arcs_cnt[el.first]);
    display = min(arcs_cnt[el.first], display_limit);
    for (i = 0; i < display; ++i) {
      arc = &arcs[el.first][i];
      printf("   [%d] (%d %d), sc=%d (%d), f=%d (%d %d)\n", i, arc->start,
             arc->end, arc->score, arc->eff_score, arc->factor,
             arc->genomic_start, arc->genomic_end);
    }

    //		printf(" [%s] raw arcs %d\n", el.first.c_str(),
    //raw_arcs[el.first].size()); 		display = min((int)raw_arcs[el.first].size(),
    //display_limit); 		for (i = 0; i < display; ++i) { 			arc =
    //&raw_arcs[el.first][i]; 			printf("   [%d] (%d %d), sc=%d, f=%d)\n", i,
    //arc->start, arc->end, arc->score, arc->factor);
    //		}
  }
}
