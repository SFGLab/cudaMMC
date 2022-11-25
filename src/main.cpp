#include <algorithm>
#include <chrono>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <time.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <BedRegion.h>
#include <HierarchicalChromosome.h>
#include <LooperSolver.h>
#include <common.h>

std::map<char, std::string> args;

void usage(const char *act = "", bool quit = true) {
  printf("Usage: -a <action> -s <path> -n <label>\n");
  printf("\t-a\taction to perform (available values: 'create')\n");
  printf("\t-s\tpath to settings file\n");
  printf("\t-n\tlabel to be used in names of output files\n");
  printf("\t-c\tchromosomes (or region) to be included in the simulation (eg "
         "'chr14:1:2500000', default: 'chr1-chr22,chrX')\n");
  printf("\t-o\toutput directory (with trailing '/'; must exists; default: "
         "'./')\n");
  printf("Debugging options\n");
  printf("\t-u\tmaximal number of chromosomes to reconstruct\n");
  printf("\t-l\tlength of chromosome fragment that will be reconstructed (in "
         "bp)\n");
  printf("Notes\n");
  printf("\tPET clusters should have filenames of the form "
         "'GM12878_CTCF_Rep1.withLinker.clusters.intra.PET<N>+.merged_anchors."
         "simple.txt', where <M> corresponds to the flag -m\n");

  if (quit) {
    args.clear();
    exit(0);
  }
}

void error_msg(string msg) {
  printf("%s\n", msg.c_str());
  args.clear();
  exit(0);
}

bool flag_set(char c) { return args.find(c) != args.end(); }

string get_arg(char c, string _default = "") {
  if (args.find(c) == args.end())
    return _default;
  return args[c];
}

void get3DPositions(string input_structure, BedRegions &regions,
                    string output_file, bool add_header = true) {
  // read 3D model
  HierarchicalChromosome chr;
  chr.fromFile(input_structure);
  chr.setLevel(3);

  printv(chr.chrs, "chromosomes found");

  FILE *f = open(output_file, "w");
  if (f == NULL)
    error("could not open the output file");
  if (add_header)
    fprintf(f, "%lu\n", regions.regions.size());

  for (uint i = 0; i < regions.regions.size(); i++) {
    int p = (regions.regions[i].start + regions.regions[i].end) / 2;
    vector3 pos = chr.find3DPosition(regions.regions[i].chr, p);
    printf("%d -> %f %f %f\n", p, pos.x, pos.y, pos.z);
    fprintf(f, "%f %f %f\n", pos.x, pos.y, pos.z);
  }

  fclose(f);
}

void preparePosition() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('o'))
    error_msg("No output file specified [-o]\n");
  if (!flag_set('b'))
    error_msg("No bed regions specified [-b]\n");

  string input = get_arg('i');
  string output = get_arg('o');
  string regions_path = get_arg('b');
  bool add_header = flag_set('h');

  // read regions to extract the position for
  BedRegions regions;
  regions.fromFile(regions_path);
  printf("%d regions loaded\n", (int)regions.regions.size());

  get3DPositions(input, regions, output, add_header);
}

void prepareSmooth() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('r'))
    error_msg("No resolution specified [-r]\n");

  string input = get_arg('i');
  int res = atoi(get_arg('r').c_str());

  string output = flag_set('o') ? get_arg('o') : input + ".smooth.txt";
  string regions_path = get_arg('b');
  int level = flag_set('l') ? atoi(get_arg('l').c_str()) : 3;
  string chr = flag_set('c') ? get_arg('c') : "";

  HierarchicalChromosome hc;
  hc.fromFile(input);
  hc.setLevel(level);
  // hc.createCurrentLevelStructure();

  Chromosome ch = hc.createEqidistantModel(res, chr);
  ch.toFile(output);

  // create artificial regions
  // BedRegions regions;
  // regions.addNewIntervals("chr", 1000, 2000, 50);
  // regions.print();
  // get3DPositions(input, regions, output, true);
}

// calculate structural distance matrix (i.e. matrix with distances between
// every pair of beads) input:
//    - hc: hierarchical chromosome model
//    - chr: chromosome to calculate the distances for (ignored if level=0)
//    - level: level of the structure (0-chr, 1-segment, 2-subanchor)
// output: filename where the output matrix is to be saved
void runCalcDistanceMatrix(HierarchicalChromosome hc, string chr, int level,
                           string output) {
  hc.setLevel(level);
  Heatmap h = hc.createStructuralHeatmap(chr, level);
  h.toFile(output);
}

void prepareCalcDistanceMatrix() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('l'))
    error_msg("No level specified [-l]\n");
  if (!flag_set('c'))
    error_msg("No chromosome specified (may be anything if level=0) [-c]\n");

  string input = get_arg('i');
  string output =
      flag_set('o') ? get_arg('o') : ftext("%s.dist.heat", input.c_str());
  string chr = get_arg('c');
  int level = atoi(get_arg('l').c_str());

  HierarchicalChromosome hc;
  hc.fromFile(input);
  runCalcDistanceMatrix(hc, chr, level, output);
}

// given a path to .hcm file create a flat representation of the specified
// chromosome on a given level if chr=="" then save all the chromosomes save the
// result to txt file output_path
void flatten(string path, int level, string chr = "", string output_path = "") {
  HierarchicalChromosome hc;
  hc.fromFile(path);
  hc.setLevel(level);
  hc.createCurrentLevelStructure();

  string out_path = output_path == "" ? path : output_path;
  if (chr == "") {
    for (auto const &item : hc.chr) {
      printf("do [%s]\n", item.first.c_str());
      if (hc.chr.size() > 1)
        hc.chr[item.first].toFile(
            ftext("%s.%s.txt", path.c_str(), item.first.c_str()));
      else
        hc.chr[item.first].toFile(out_path);
    }
  } else {
    if (hc.chr.find(chr) == hc.chr.end()) {
      printf("Could not find chromosome [%s]\nExisting chromosomes:",
             chr.c_str());
      for (auto const &item : hc.chr)
        printf(" %s", item.first.c_str());
      return;
    }
    hc.chr[chr].toFile(ftext("%s.txt", path.c_str()));
  }
}

void prepareFlatten() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('l'))
    error_msg("No level specified [-l]\n");

  string input = get_arg('i');
  int level = atoi(get_arg('l').c_str());

  string chr = flag_set('c') ? get_arg('c') : "";
  string output = flag_set('o') ? get_arg('o') : input + ".txt";
  flatten(input, level, chr, output);
}

void extract(string input_path, string output_path, int start, int end) {
  printf("extract fragment: %d %d\n", start, end);
  HierarchicalChromosome h;
  h.fromFile(input_path);
  printf("size = %d\n", (int)h.clusters.size());
  if (h.clusters.size() > 0) {
    HierarchicalChromosome extr = h.extractFragment(start, end);
    extr.toFile(output_path);
  }
}

void prepareExtract() {
  if (!flag_set('i'))
    error_msg("No input specified [-i]\n");
  if (!flag_set('o'))
    error_msg("No input specified [-o]\n");
  if (!flag_set('s'))
    error_msg("No start position specified [-s]\n");
  if (!flag_set('e'))
    error_msg("No end position specified [-e]\n");

  string input = get_arg('i');
  string output = get_arg('o');

  int start = atoi(get_arg('s').c_str());
  int end = atoi(get_arg('e').c_str());

  extract(input, output, start, end);
}

void ensembleAnalysis(string input_dir, string pattern, string output_dir) {
  // read all structures (hierarchical), and for each level (chr, segment,
  // subanchor) get the chr then calculate the matrix of pairwise distances
  // between structures
  printf("ensemble analysis:\n");

  // first read all hierarchical structures to the memory
  vector<HierarchicalChromosome> hchr; // keep all structures loaded
  string repl = "{N}";

  if (pattern.find(repl) == string::npos)
    error_msg(ftext("pattern does not contain replacement token (\"%s\")",
                    repl.c_str()));

  int p = 0;
  while (1) {
    string str = pattern;
    str.replace(str.find(repl), repl.size(), ftext("%d", p));
    string path = input_dir + str;
    if (!file_exists(path))
      break;
    printf("[%s]\n", path.c_str());
    p++;

    HierarchicalChromosome h;
    h.fromFile(path);
    hchr.push_back(h);

    if (p > 10000)
      break; // safety check to prevent infinite loop
  }

  int N = hchr.size(); // number of structures in the ensemble
  if (N < 2) {
    printf("Too few structures (%d) to conduct an analysis\n", N);
    return;
  }

  // now calculate the distances matrices
  vector<Heatmap>
      heat; // keep structural heatmaps for all structures on a given level
  Heatmap heat_pairs(N); // heatmap with pair-wise distances

  vector<string> chrs = hchr[0].chrs; // save chromosomes

  // first chr level
  if (chrs.size() < 2)
    printf("Only %d chromosome available, skip chr level\n", (int)chrs.size());
  else {
    printf("chr level\n");

    printf("create structural heatmaps\n");

    for (int j = 0; j < N; ++j) { // for all structures
      heat.push_back(hchr[j].createStructuralHeatmap("", 0));
    }

    printf("create distance heatmap\n");
    for (int i = 0; i < N; ++i) { // for all structures
      for (int j = i + 1; j < N; ++j) {
        heat_pairs.v[i][j] = heat_pairs.v[j][i] = heat[i].calcDistance(heat[j]);
      }
    }

    printf("save\n");
    heat_pairs.toFile(input_dir + "structural_distances_lvl0.heat", false);
    heat_pairs.zero();

    heat.clear();
  }

  // segment and subanchor levels
  int i, j;
  for (j = 0; j < N; ++j)
    hchr[j].levelDown(); // move to the segment level

  float d;
  for (int lvl = 1; lvl <= 2; ++lvl) {
    // i is current level (1-segment, 2-subanchor)
    printf("do level %d\n", lvl);

    // go through each chromosome one by one
    for (string chr : chrs) {
      for (j = 0; j < N; ++j) { // for each structure in the ensemble
        hchr[j].setLevel(lvl);
        heat.push_back(hchr[j].createStructuralHeatmap(chr, lvl));
      }

      for (i = 0; i < N; ++i) { // for all structures
        for (j = i + 1; j < N; ++j) {
          d = heat[i].calcDistance(heat[j]);
          heat_pairs.v[i][j] += d;
          heat_pairs.v[j][i] += d;
        }
      }
    }

    heat_pairs.toFile(input_dir + ftext("structural_distances_lvl%d.heat", lvl),
                      false);

    heat_pairs.zero();
    heat.clear();
  }
}

void prepareEnsembleAnalysis() {
  if (!flag_set('i'))
    error_msg("No input directory specified [-i]\n");
  if (!flag_set('p'))
    error_msg("No pattern specified [-p]\n");

  string input = get_arg('i');
  string pattern = get_arg('p');
  string output = flag_set('o') ? get_arg('o') : input;

  ensembleAnalysis(input, pattern, output);
}

void rewiringAnalysis(string input_dir, string pattern, int level = 1,
                      int resolution = 1000000, string output_dir = "") {
  // create structural heatmaps for all structures (hierarchical)
  printf("rewiring analysis:\n");

  string repl = "{N}";
  if (pattern.find(repl) == string::npos)
    error_msg(ftext("pattern does not contain the replacement token (\"%s\")",
                    repl.c_str()));

  int p = 0;
  while (1) {
    string str = pattern;
    str.replace(str.find(repl), repl.size(), ftext("%d", p));
    string path = input_dir + str;
    if (!file_exists(path)) {
      if (p == 0)
        error(ftext("No file (%s) found\n", path.c_str()));
      break;
    }
    printf("[%s]\n", path.c_str());
    p++;

    HierarchicalChromosome h;
    h.fromFile(path);
    h.setLevel(level);
    // hc.createCurrentLevelStructure();

    // go through each chromosome one by one
    for (string chr : h.chrs) {

      printf("   chr %s\n", chr.c_str());

      Chromosome ch = h.createEqidistantModel(resolution, chr);
      Heatmap h = ch.createInverseHeatmap();
      printf("save to %s\n",
             ftext("%s_%s_%d.heat", path.c_str(), chr.c_str(), p).c_str());
      h.toFile(ftext("%s_%s_%d.heat", path.c_str(), chr.c_str(), p));
    }

    if (p > 10000)
      break; // safety check to prevent infinite loop
  }
}

void prepareRewiring() {
  if (!flag_set('i'))
    error_msg("No input directory specified [-i]\n");
  if (!flag_set('p'))
    error_msg("No pattern specified [-p]\n");

  string input = get_arg('i');
  string pattern = get_arg('p');
  int level = flag_set('l') ? atoi(get_arg('l').c_str()) : 1;
  int resolution = flag_set('r') ? atoi(get_arg('r').c_str()) : 1e6;
  string output = flag_set('o') ? get_arg('o') : input;

  rewiringAnalysis(input, pattern, level, resolution, output);
}

// example formats: "chr1", "chr2,chr3", "chr8,chr10-chr15", "genome"
std::vector<std::string> parseChromosomeDescription(string desc) {
  std::vector<string> ret;

  if (desc == "genome") {
    for (int i = 1; i <= 22; i++)
      ret.push_back(ftext("chr%d", i));
    ret.push_back("chrX");
    return ret;
  }

  std::vector<string> tmp = split(desc);
  for (string token : tmp) {
    if (token.find('-') != string::npos) {
      std::vector<string> tmp2 = split(desc, '-');
      if (tmp2.size() == 2) {
        int st = atoi(tmp2[0].substr(3).c_str());
        int end = atoi(tmp2[1].substr(3).c_str());
        for (int i = st; i <= end; i++)
          ret.push_back(ftext("chr%d", i));
      }
    } else
      ret.push_back(token);
  }

  return ret;
}

// ensemble_mode - whether to create 1 structure or a whole ensemble. accepted
// values:
//	   -1 - single structure
//		0 - generate N independent structures
//		1 - generate structures in hierarchical way (start with N
//structures for chr level, then for each of these create N 			structures using the
//parent as a template, and so on.)
// ensemble_multiplicity - controls how many structures to create for an
// ensemble, corresponds to 'N' selected factors - if provided, uses only the
// cluster and singleton files related to specified clusters (order in stg file
// matters). If used, the number of factors must be equal to
//    		the number of cluster and singleton files. If empty, all input files
//    are used
void runLooper(std::vector<string> chrs, BedRegion region_of_interest,
               string label, string outdir, bool use_new_file_format,
               int max_level, int chr_number_limit = -1, int length_limit = -1,
               int ensemble_size = 1,
               const std::vector<string> selected_factors = vector<string>()) {

  printf("run [%s]\n", label.c_str());

// if on unix, attempt to create a directory
#ifdef __unix__
  int res = mkdir(outdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (res != 0 && errno != EEXIST)
    error("Couldn't create the output directory");
#endif

  string path_data = Settings::dataDirectory;
  string path_anchors = Settings::dataAnchors;
  std::vector<string> path_pets = split(Settings::dataPetClusters);
  std::vector<string> path_singletons = split(Settings::dataSingletons);
  std::vector<string> path_singletons_inter =
      split(Settings::dataSingletonsInter);
  std::vector<string> factors = split(Settings::dataFactors);

  if (selected_factors.size() > 0) {
    // check if there is a propoer number of cluster/singletons files
    if (path_pets.size() != path_singletons.size() ||
        factors.size() != path_pets.size()) {
      error("Invalid number of input files! There should be exactly one "
            "cluster and singleton file per factor when -g option is used");
      return;
    }

    // remove unwanted files
    for (uint i = factors.size() - 1; i >= 1; --i) {
      bool found = false;
      for (uint j = 0; j < selected_factors.size(); ++j) {
        if (factors[i] == selected_factors[j])
          found = true;
      }
      if (!found) {
        printf("ignore factor: %s\n", factors[i].c_str());
        path_pets.erase(path_pets.begin() + i);
        path_singletons.erase(path_singletons.begin() + i);
        factors.erase(factors.begin() + i);
      }
    }
  }

  string anchors = ftext("%s%s", path_data.c_str(), path_anchors.c_str());

  std::vector<string> arcs_clusters;
  std::vector<string> arcs_singletons;
  std::vector<string> arcs_singletons_inter;

  for (uint i = 0; i < path_pets.size(); ++i) {
    arcs_clusters.push_back(
        ftext("%s%s", path_data.c_str(), path_pets[i].c_str()));
  }

  for (uint i = 0; i < path_singletons.size(); ++i) {
    arcs_singletons.push_back(
        ftext("%s%s", path_data.c_str(), path_singletons[i].c_str()));
  }

  for (uint i = 0; i < path_singletons_inter.size(); ++i) {
    arcs_singletons_inter.push_back(
        ftext("%s%s", path_data.c_str(), path_singletons_inter[i].c_str()));
  }

  LooperSolver lsm(label, outdir);
  lsm.setDebuggingOptions(chr_number_limit, length_limit);

  lsm.setContactData(chrs, region_of_interest, anchors, factors, arcs_clusters,
                     arcs_singletons, arcs_singletons_inter);

  if (Settings::useDensity)
    lsm.initDensityData();

  // print_vector(lsm.densityCoordToStructure(vector3(17, 42, 28)), "telomer1");
  // print_vector(lsm.densityCoordToStructure(vector3(13, 62, 20)), "telomer2");

  printf("create tree\n");
  lsm.createTreeGenome();

  if (ensemble_size == 1) {
    // proceed the normal way
    printf("\nreconstruct heatmap\n");
    lsm.reconstructClustersHeatmap();

    if (max_level >= LVL_INTERACTION_BLOCK) {
      printf("reconstruct arcs\n");
      lsm.reconstructClustersArcsDistances();
    }

    printf("save\n");
    HierarchicalChromosome hc = lsm.getModel();

    // hc.center(false); // be careful! when using density this will shift the
    // structure from the density map

    lsm.showCoordinatesRange(1);

    if (use_new_file_format)
      hc.toFile(ftext("%sloops_new_%s.hcm", outdir.c_str(), label.c_str()));
    else
      hc.toFilePreviousFormat(
          ftext("%sloops_%s.hcm", outdir.c_str(), label.c_str()));
  } else if (ensemble_size > 1) {
    // generate N independent structures
    printf("\nreconstruct %d structures\n", ensemble_size);

    // note: subanchor beads are added during
    // reconstructClustersArcsDistances(), so we need to remove them between
    // runs otherwise all the subanchor beads from one run will be treated as
    // anchor beads in the next one

    for (int i = 0; i < ensemble_size; ++i) {
      printf("structure %d/%d\n", i + 1, ensemble_size);

      printf("reconstruct heatmap\n");
      lsm.reconstructClustersHeatmap();

      if (max_level >= LVL_INTERACTION_BLOCK) {
        printf("reconstruct arcs\n");
        lsm.reconstructClustersArcsDistances();
      }

      HierarchicalChromosome hc = lsm.getModel();
      // hc.center(false);

      if (use_new_file_format)
        hc.toFile(
            ftext("%sloops_new_%s_%d.hcm", outdir.c_str(), label.c_str(), i));
      else
        hc.toFilePreviousFormat(
            ftext("%sloops_%s_%d.hcm", outdir.c_str(), label.c_str(), i));

      // remove subanchor beads
      lsm.removeSubanchorBeads();
      lsm.reset();
    }
  } else
    error("ensemble size not correct");
}

void prepareLooper() {
  printf("prepare\n");

  // we need to know what regions we should reconstruct
  // we can either provide a number of chromosomes, or a single chromosome
  // region in the first case we need to put all chromosomes id to 'chrs' in the
  // second we put region info in 'region_of_interest', but we also need to
  // update 'chrs' to contain the selected chromosome (we need it because all
  // tree manipulations are based on 'chrs')

  string chromosomes =
      "genome"; // which chromosomes (eg. "chr2", "chr1-chr22,chrX")
  BedRegion region_of_interest;
  std::vector<string> chrs;
  if (flag_set('c')) {
    chromosomes = get_arg('c');
    if (BedRegion::tryParse(
            chromosomes)) { // check if the value is a BED region
      // it seems like it is, so parse the description and update 'chrs'
      region_of_interest.parse(chromosomes);
      chrs.push_back(region_of_interest.chr);
    } else
      chrs = parseChromosomeDescription(chromosomes);
  } else {
    // if there was no -c parameter use default parameters
    chrs = parseChromosomeDescription(chromosomes);
  }

  string name = "";
  if (flag_set('n'))
    name = get_arg('n');
  else {
    name = chromosomes;
    std::replace(name.begin(), name.end(), ':', '_');
    std::replace(name.begin(), name.end(), '-', '_');
  }

  if (flag_set('t') && flag_set('d'))
    error("-t and -d options are exclusive, you may only provide one of them");

  std::vector<string> selected_factors = vector<string>();
  if (flag_set('g'))
    selected_factors = split(get_arg('g'));

  // if there is a structural template provided
  if (flag_set('t')) {
    Settings::templateSegment = get_arg('t');
    if (flag_set('p'))
      Settings::templateScale = atof(get_arg('p').c_str());
  }

  // if there is a distance heatmap provided
  if (flag_set('d')) {
    Settings::distHeatmap = get_arg('d');
    if (flag_set('b'))
      Settings::distHeatmapScale = atof(get_arg('b').c_str());
  }
  if (flag_set('f')) {
    Settings::dataSegmentHeatmap = get_arg('f');
    if (flag_set('b'))
      Settings::distHeatmapScale = atof(get_arg('b').c_str());
  }

  // -x allows to override the values in the settings file
  if (flag_set('x')) {
    vector<string> arr = split(get_arg('x'), ',');
    for (std::size_t i = 0; i < arr.size(); ++i) {
      vector<string> arr_par = split(arr[i], ':');

      if (arr_par[0] == "ib_random_walk_jumps")
        Settings::ibRandomWalkJumps = atof(arr_par[1].c_str());
      else if (arr_par[0] == "freq_dist_scale")
        Settings::freqToDistHeatmapScale = atof(arr_par[1].c_str());
      else if (arr_par[0] == "freq_dist_power")
        Settings::freqToDistHeatmapPower = atof(arr_par[1].c_str());
      else if (arr_par[0] == "freq_dist_scale_inter")
        Settings::freqToDistHeatmapScaleInter = atof(arr_par[1].c_str());
      else if (arr_par[0] == "freq_dist_power_inter")
        Settings::freqToDistHeatmapPowerInter = atof(arr_par[1].c_str());
      else if (arr_par[0] == "count_dist_a")
        Settings::countToDistA = atof(arr_par[1].c_str());
      else if (arr_par[0] == "count_dist_scale")
        Settings::countToDistScale = atof(arr_par[1].c_str());
      else if (arr_par[0] == "count_dist_shift")
        Settings::countToDistShift = atof(arr_par[1].c_str());
      else if (arr_par[0] == "count_dist_base_level")
        Settings::countToDistBaseLevel = atof(arr_par[1].c_str());
      else if (arr_par[0] == "genomic_dist_power")
        Settings::genomicLengthToDistPower = atof(arr_par[1].c_str());
      else if (arr_par[0] == "genomic_dist_scale")
        Settings::genomicLengthToDistScale = atof(arr_par[1].c_str());
      else if (arr_par[0] == "genomic_dist_base")
        Settings::genomicLengthToDistBase = atof(arr_par[1].c_str());
    }
  }

  bool use_new_file_format = flag_set('z');
  string outdir = flag_set('o') ? get_arg('o') : "./";
  int length_limit = flag_set('l') ? atoi(get_arg('l').c_str()) : -1;
  int chr_number_limit = flag_set('u') ? atoi(get_arg('u').c_str()) : -1;
  int max_level = flag_set('v') ? atoi(get_arg('v').c_str()) : 5;

  int ensemble_size = flag_set('m') ? atoi(get_arg('m').c_str()) : 1;

  bool ignore_cache = flag_set('h');
  if (ignore_cache)
    Settings::useInputCache = false;

  args.clear(); // clean here, otherwise heap corruption occurs when program
                // abnormally stops (i.e. exit(0))

  runLooper(chrs, region_of_interest, name, outdir, use_new_file_format,
            max_level, chr_number_limit, length_limit, ensemble_size,
            selected_factors);
}

int main(int argc, char **argv) {
  printf("start\n");
  setbuf(stdout, NULL);

  // long int *tmp = NULL;
  auto tmp = chrono::high_resolution_clock::now().time_since_epoch().count();
  srand((unsigned int)tmp);
  printf("time = %u\n", (unsigned int)tmp);

  opterr = 0;
  int opt = 0;

  // Density d;
  // d.fromSegmentationFile("/home/pszemsza/Desktop/density_modeling/density_maps/intensity_chr1_D7_object1.csv");
  // Density dd = d.scaleDown(2);
  // dd.toSegmentationFile("/home/pszemsza/Desktop/density_modeling/density_maps/intensity_chr1_D7_object1_scaled2.csv");
  // exit(0);

  while ((opt = getopt(argc, argv,
                       "s:a:o:c:n:l:u:t:d:p:b:e:m:i:v:zf:hr:g:x:")) != EOF) {
    if (optarg == NULL)
      optarg = (char *)"";
    printf("opt = [%c %s]\n", opt, optarg);
    args[opt] = optarg;
  }

  if (!flag_set('a') || args['a'] == "c" || args['a'] == "create") {

    if (!flag_set('s')) {
      printf("no settings file specified!\n");
      usage();
    }

    Settings stg;
    string stg_path = get_arg('s');
    printf("Load settings from [%s]\n", stg_path.c_str());
    stg.loadFromINI(stg_path);
    // stg.print();
    // exit(0);

    prepareLooper();
  }

  else if (args['a'] == "e" || args['a'] == "ensemble")
    prepareEnsembleAnalysis();
  else if (args['a'] == "r" || args['a'] == "extract")
    prepareExtract();
  else if (args['a'] == "p" || args['a'] == "position")
    preparePosition();
  else if (args['a'] == "s" || args['a'] == "smooth")
    prepareSmooth();
  else if (args['a'] == "d" || args['a'] == "distance")
    prepareCalcDistanceMatrix();
  else if (args['a'] == "f" || args['a'] == "flatten")
    prepareFlatten();
  else if (args['a'] == "w" || args['a'] == "rewiring")
    prepareRewiring();
  else {
    printf("No matching action specified [%s]!\n", args['a'].c_str());
    usage();
  }

  args.clear();
  printf("end\n");
  return 0;
}
