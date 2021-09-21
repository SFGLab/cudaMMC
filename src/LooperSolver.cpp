/*
 * LooperSolverMixed.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: psz
 */

#include "../include/LooperSolver.h"

LooperSolver::LooperSolver(string label, string outdir) {

	this->label = label;
	this->output_dir = outdir;

	rw_pos.set(0.0f, 0.0f, 0.0f);
	zero_distance = true;
	tmp_var = 0;

	setDebuggingOptions();	// use default values (ie. disable debug options)

	for (int i = 2; i <= 100; ++i) freq_to_distance[i] = freqToDistance(i, false);

	if (Settings::outputLevel >= 7) {
		printf("* freq to distance:\n");
		for (int i = 2; i <= 100; i+=5) printf("%d -> %f\n", i, freq_to_distance[i]);
		printf("\n");

		// test genomic to 3D distance
		printf("* genomic to 3D distance\n");
		int dist[] = {1000, 2000, 5000, 10000, 15000, 20000, 50000, 100000, 200000};
		for (int i=0; i<9; i++) printf("%dkb -> %lf\n", dist[i]/1000, genomicLengthToDistance(dist[i]));

		printf("* heatmap frequency to 3D distance\n");
		float freq[] = {2.0f, 1.5f, 1.0f, 0.75f, 0.5f, 0.25f, 0.1f, 0.05f, 0.01f, 0.001f};
		for (int i=0; i<10; i++) printf("%f -> %lf\n", freq[i], freqToDistanceHeatmap(freq[i]));
	}

	if (Settings::dataCentromeres.size()) {
		centromeres.fromFile(Settings::dataCentromeres);
		printf("centromeres loaded\n");
	}

	if (Settings::dataSegmentsSplit.size()) {
		segments_predefined.fromFile(Settings::dataSegmentsSplit);
		printf("use predefined segments (cnt = %d)\n", (int)segments_predefined.regions.size());
	}
	else error("You need to provide the segment split");
}

void LooperSolver::print() {
	printf("*** clusters: %d\n", (int)clusters.size());
	printf("clusters\n");

	for (size_t i = 0; i < clusters.size(); ++i) {
		printf("[%lu] ", i);
		clusters[i].print();
	}
	printf("\n");
}

void LooperSolver::initDensityData() {
	printf("init density\n");
	density.fromSegmentationFile(Settings::densityMapFile);
	density.normalize();
	density.toFile(ftext("%sdensity_segm.txt", output_dir.c_str()));

	density_curr.init(density.size_x, density.size_y, density.size_z, density.range_x_start, density.range_y_start, density.range_z_start, false);
}


void LooperSolver::printActiveRegion() {
	printf("*** active region: size = %d\n", (int)active_region.size());
	for (size_t i = 0; i < active_region.size(); ++i) {
		printf("[%lu - %d] ", i, active_region[i]);
		if (i+1<active_region.size()) {
			printf("(dist exp=%lf, true=%lf) ", clusters[active_region[i]].dist_to_next, (clusters[active_region[i]].pos-clusters[active_region[i+1]].pos).length() );
		}
		clusters[active_region[i]].print();
	}
}


void LooperSolver::reconstructClustersHeatmap() {

	heatmap.clear();

	if (Settings::randomWalk) {

		// position segment-level points with a random walk
		printf("random walk for segment level\n");

		setLevel(LVL_SEGMENT);	// set segment level

		float size = genomicLengthToDistance(1000000);
		for (string chr: chrs) {
			printf("%s\n", chr.c_str());
			rw_pos.set(0.0f, 0.0f, 0.0f);
			rw_pos = displace(rw_pos, size, Settings::use2D);	// first point

			for (uint i = 0; i < current_level[chr].size(); ++i) {
				rw_pos = displace(rw_pos, 50.0f, Settings::use2D);
				clusters[current_level[chr][i]].pos = rw_pos;
			}

			// smoothly interpolate children (i.e. anchor level)
			interpolateChildrenPositionSpline(current_level[chr]);
		}
	}
	else {

		// if entire region is contained in a single segment then there is no heatmap
		// we can just position the corresponding bead in (0,0,0) and go to the anchor level
		setLevel(LVL_SEGMENT);	// segment level
		if (current_level.size() == 1 && current_level[chrs[0]].size() <= 1) {
			printf("current region size = 1, move to (0,0,0)\n");
			clusters[current_level[chrs[0]][0]].pos.set(0.0, 0.0, 0.0);
			return;
		}

		bool chromosome_level = chrs.size() > 1;	// do we need to create a heatmap on chromosome level?

		string heat_file;

		// ******* chromosome level
		if (chromosome_level) {
			printf("\n* chromosome level\n");

			setLevel(LVL_CHROMOSOME);
			//heat_file = ftext("%ssingletons_%s_chromosomes.heat", output_dir.c_str(), label.c_str());
			heat_file = ftext("%ssingletons_chromosomes_%s.heat", output_dir.c_str(), label.c_str());
			Heatmap hc, hnc;

			if (Settings::useInputCache && file_exists(heat_file)) {
				printf("load chromosome heatmap from cache file [%s]\n", heat_file.c_str());
				hc.fromFile(heat_file);
			}
			else {
				printf("create chromosome heatmap\n");
				hc = createSingletonHeatmap(1);
				hc.toFile(heat_file);
			}

			hnc = hc;
			normalizeHeatmapDiagonalTotal(hnc, 1.0f);

			if (Settings::outputLevel >= 4) {
				float max, min;
				hnc.getRange(min, max);
				printf("hmap range: %f %f\n", min, max);
			}

			//hnc.toFile(ftext("%ssingletons_%s_chromosomes_norm.heat", output_dir.c_str(), label.c_str()));
			hnc.toFile(ftext("%ssingletons_chromosomes_norm_%s.heat", output_dir.c_str(), label.c_str()));

			heatmap_chromosome_boundaries.clear(); // make sure it is empty for chr level (MC calc score)

			// note that order of pushing heatmaps is important, reconstructClustersHeatmapSingleLevel() take index of heatmap as a parameter
			heatmap.push_back(hnc);

			createDistanceHeatmap(0);
			heatmap_dist.toFile(ftext("%sheat_dist_chromosomes_%s.heat", output_dir.c_str(), label.c_str()), true);

			reconstructClustersHeatmapSingleLevel(0);

			if (Settings::outputLevel >= 3) {
				printf("structure dimensions (chromosomes)\n");
				showCoordinatesRange();
				printf("\n");
			}
		}
		else {
			// even if we don't have chromosome level we need to put a heatmap for we want the
			// indices of levels and heatmaps to be equal
			Heatmap tmp;
			heatmap.push_back(tmp);
		}


		setLevel(LVL_SEGMENT); // segments level

		// the structure on the segment level can be obtained from the template structure
		// it can be either .hcm model or .chr/.txt file with coordinates
		if (Settings::templateSegment.size() > 1) {
			printf("use structural template (%s)\n", Settings::templateSegment.c_str());
			Chromosome tpl;
			if (Settings::templateSegment.find(".hcm") != std::string::npos) {
				HierarchicalChromosome hc;
				hc.fromFile(Settings::templateSegment);
				hc.setLevel(LVL_INTERACTION_BLOCK);
				hc.createCurrentLevelStructure();
				tpl = hc.chr[hc.chrs[0]];
			}
			else tpl.fromFile(Settings::templateSegment.c_str());

			tpl.scale(Settings::templateScale);
			for (string chr: chrs) {
				for (uint i = 0; i < current_level[chr].size(); ++i) {
					clusters[current_level[chr][i]].pos = tpl.points[i];
				}
			}
		}
		else {
			// ***** create heatmap for segment level

			printf("\n\n* segment level\n");

			// if there is a distance heatmap provided
			if (Settings::distHeatmap.size() > 0) {
				setLevel(LVL_SEGMENT);
				int tot_size = 0;
				for (string chr: chrs) tot_size += current_level[chr].size();
				printf("dist heatmap size = %d\n", tot_size);
				heatmap_dist.fromMDS(Settings::distHeatmap, tot_size);
				heatmap_dist.scale(Settings::distHeatmapScale);
			}
			else {
				// if not, we use singletons segment heatmap to create one
				Heatmap h, hn;

				//heat_file = ftext("%ssingletons_segment.heat", output_dir.c_str(), label.c_str());
				heat_file = ftext("%ssingletons_segment_%s.heat", output_dir.c_str(), label.c_str());

				// if there is a frequency heatmap provided
				if (Settings::dataSegmentHeatmap.length() > 2) {
					if (file_exists(Settings::dataSegmentHeatmap)) {
						printf("load heatmap from file [%s]\n", Settings::dataSegmentHeatmap.c_str());
						h.fromFile(Settings::dataSegmentHeatmap);
					}
					else {
						printf("segment heatmap file does not exist [%s]\n", Settings::dataSegmentHeatmap.c_str());
						exit(0);
					}

				}
				else if (Settings::useInputCache && file_exists(heat_file)) {
					printf("load heatmap from cache file [%s]\n", heat_file.c_str());
					h.fromFile(heat_file);
				}
				else {
					printf("create segment heatmap\n");
					h = createSingletonHeatmap(1);
					h.toFile(heat_file);
				}

				hn = normalizeHeatmap(h);
				normalizeHeatmapDiagonalTotal(hn, 1.0f);
				normalizeHeatmapInter(hn, Settings::heatmapInterScaling);
				//hn.toFile(ftext("%ssingletons_%s_segment_norm.heat", output_dir.c_str(), label.c_str()));
				hn.toFile(ftext("%ssingletons_segment_norm_%s.heat", output_dir.c_str(), label.c_str()));
				heatmap.push_back(hn);
				createDistanceHeatmap(1);
				heatmap_dist.toFile(ftext("%sheat_dist_segment_%s.heat", output_dir.c_str(), label.c_str()), true);
			}

			heatmap_chromosome_boundaries.clear();
			if (chrs.size() > 1) {
				int p = 0;

				heatmap_chromosome_boundaries.push_back(p);
				for (string chr: chrs) {
					p += current_level[chr].size();
					heatmap_chromosome_boundaries.push_back(p);
				}
				printv(heatmap_chromosome_boundaries, true, true, "chromosome boundaries");
			}

			reconstructClustersHeatmapSingleLevel(1);
		}

		if (Settings::outputLevel >= 3) {
			printf("structure dimensions (segments)\n");
			showCoordinatesRange();
		}

		//Heatmap true_dist = getModel().getDistancesHeatmap();
		//true_dist.toFile(ftext("%strue_dist.heat", output_dir.c_str()));
	}
}

// level - which heatmap to use, 0 is the most general
void LooperSolver::reconstructClustersHeatmapSingleLevel(int level) {
	setLevel(level);  // set the corresponding level (here we index from 0, and levels are 1-based)

	// set initial structure and average distance between neighboring beads
	float avg_dist = 0.0f;
	std::vector<vector3> initial_structure;	// store initial position for all regions
	if (level == LVL_CHROMOSOME) avg_dist = heatmap_dist.getAvg() * Settings::noiseCoefficientLevelChr;  // position randomly around the origin (a point for chromosome)
	if (level == LVL_SEGMENT) {

		avg_dist = heatmap_dist.getAvg() * Settings::noiseCoefficientLevelSegment;
		for (string chr: chrs) {
			int par = clusters[current_level[chr][0]].parent;
			vector3 origin = clusters[par].pos;

			if (Settings::useDensity) origin = density.center;
			if (Settings::useTelomerePositions) origin = (Settings::telomere_1_position+Settings::telomere_2_position) / 2.0f;

			print_vector(origin, "set origin to");

			for (size_t i=0; i<current_level[chr].size(); i++) {
				clusters[current_level[chr][i]].pos = origin;
			}
		}
	}


	output(3, " avg distance = %f\n", avg_dist);

	// create active region (concatenate current_level for all chromosomes)
	// for every el of active_region we need to know which chromosome it belongs to, so
	// we keep track of how many segments there are for every chromosome (see activeRegionIndexToChr())
	active_region.clear();
	active_region_to_chr.clear();
	for (string chr: chrs) {
		active_region.insert(active_region.end(), current_level[chr].begin(), current_level[chr].end());
		active_region_to_chr[chr] = active_region.size();
	}

	printm(active_region_to_chr, false, true, "active_region_to_chr");
	printv(active_region, true, true, "active region");

	for (size_t i = 0; i < active_region.size(); ++i) initial_structure.push_back(clusters[active_region[i]].pos);

	int steps = level == LVL_CHROMOSOME ? Settings::simulationStepsLevelChr : Settings::simulationStepsLevelSegment;

	double score, best_score = -1.0;
	std::vector<vector3> best_structure;	// to save best structure
	for (int k = 0; k < steps; ++k) {

		// random initial position (use initial_structure and add some noise)
		for (size_t i = 0; i < active_region.size(); ++i) {
			clusters[active_region[i]].pos = initial_structure[i] + random_vector(avg_dist, Settings::use2D) ;
		}

		// if segment level and we have telomere positions
		// TODO: currently works only if a single chromosome is used
		if (level == LVL_SEGMENT && Settings::useTelomerePositions) {

//			print_vector(Settings::telomere_1_position, "tel 1 org");
//			print_vector(densityCoordToStructure(Settings::telomere_1_position), "tel 1 tr");
//
//			print_vector(Settings::telomere_2_position, "tel 2 org");
//			print_vector(densityCoordToStructure(Settings::telomere_2_position), "tel 2 tr");

			clusters[active_region[0]].pos = densityCoordToStructure(Settings::telomere_1_position - density.origin);
			clusters[active_region[0]].is_fixed = true;
			int last = active_region.size()-1;
			clusters[active_region[last]].pos = densityCoordToStructure(Settings::telomere_2_position - density.origin);
			clusters[active_region[last]].is_fixed = true;
		}

		output(2, "MC, heatmap, level %d, run %d/%d\n", level, k+1, steps);

		score = ParallelMonteCarloHeatmap(avg_dist);

		output(2, "score = %lf %lf\n", score, best_score);

		if (Settings::useDensity) {
			score = MonteCarloHeatmapAndDensity(avg_dist);
			output(2, "score (density) = %lf %lf\n", score, best_score);
		}

		// if current structure is better than the best one - save it
		if (score < best_score || best_score < 0.0) {
			best_structure.clear();
			for (size_t i = 0; i < active_region.size(); ++i) best_structure.push_back(clusters[active_region[i]].pos);
			best_score = score;
		}
	}

	// restore best structure
	for (size_t i = 0; i < active_region.size(); ++i) clusters[active_region[i]].pos = best_structure[i];


	if (Settings::useDensity) {
		updateDensity();
		density_curr.toFile(ftext("%sdensity_norm.txt", output_dir.c_str()));
	}

	//for (size_t i = 0; i < active_region.size(); ++i) print_vector(clusters[active_region[i]].pos, ftext("p_%d: ", i).c_str());
}

double LooperSolver::MonteCarloHeatmap(float step_size) {

	step_size *= 0.5;

	double T = Settings::maxTempHeatmap;		// set current temperature to max
	double tp = 0.0;	// transition probability
	vector3 displacement;
	int i, p, ind;
	bool ok;
	double score_prev, score_curr;		// global score (used to check if we can stop)
	double score_density;
	double milestone_score;		// score measured every N steps, used to see whether the solution significantly improves
	int success = 0;			// overall number of successes
	int milestone_success = 0;	// calc how many successes there were since last milestone
	string chr; 				// tmp variable to keep track to which chromosome a specific point belongs
	int size = active_region.size();

	if (size <= 1) return 0.0;	// there is nothing to do

	// calc initial score
	score_curr = calcScoreHeatmapActiveRegion();	// score from heatmap
	score_density = calcScoreDensity();

	score_prev = score_curr;
	milestone_score = score_curr;

	output(2, "initial score: %lf (density=%lf)\n", score_curr, score_density);

	int milestone_cnt = 0;

	i = 1;
	while (true) {

		// select point to mutate and find chromosome on which it is located
		p = random(size);	// index in 'active_region'
		ind = active_region[p];	// index in 'clusters'

		if (clusters[ind].is_fixed) continue;	// check if bead is fixed (for example telomere)

		displacement = random_vector(step_size, Settings::use2D);	// generate random displacement vector
		clusters[ind].pos += displacement;
		score_curr = calcScoreHeatmapActiveRegion();
		ok = score_curr <= score_prev;

		if (!ok && T > 0.0) {
			tp = Settings::tempJumpScaleHeatmap * exp(-Settings::tempJumpCoefHeatmap * (score_curr/score_prev) / T);
			ok = withChance(tp);
			//printf("%f %f %f %f   %d\n", score_curr, score_prev, T, tp, ok);
		}

		if (ok) {
			success++;
			milestone_success++;
		}
		else {
			clusters[ind].pos -= displacement;
			score_curr = score_prev;		// score doesn't change
		}

		T *= Settings::dtTempHeatmap;    // cooling
		// check if we should stop
		if (i % Settings::MCstopConditionStepsHeatmap == 0) {
			score_density = calcScoreDensity();
			//printf("milestone density score = %lf\n", score_density);
			output(2, "milestone [%d]: score = %lf, last = %lf (%f), T=%lf successes = %d, last = %d  (den=%lf), step=%f\n", milestone_cnt, score_curr,
					milestone_score, score_curr/milestone_score, T, success, milestone_success, score_density, step_size);

			// stop if the improvement since the last milestone is less than 0.5%, or if the score is too small (may happen when heatmap is small and points are perfectly arranged)
			if ((score_curr > Settings::MCstopConditionImprovementHeatmap * milestone_score &&
					milestone_success < Settings::MCstopConditionMinSuccessesHeatmap) || score_curr < 1e-6) {
				break;
			}
			milestone_score = score_curr;
			milestone_success = 0;

			milestone_cnt++;
		}

		score_prev = score_curr;
		i++;
	}

	return score_curr;
}

double LooperSolver::MonteCarloHeatmapAndDensity(float step_size) {

	//step_size *= 0.25;

	double T = Settings::maxTempHeatmapDensity;		// set current temperature to max
	double tp = 0.0;	// transition probability
	vector3 displacement;
	int i, p, ind;
	bool ok;
	double score_prev, score_curr;		// global score (used to check if we can stop)
	double score_density, score_heatmap;
	double score_density_prev, score_heatmap_prev;
	double milestone_score;		// score measured every N steps, used to see whether the solution significantly improves
	int success = 0;			// overall number of successes
	int milestone_success = 0;	// calc how many successes there were since last milestone
	string chr; 				// tmp variable to keep track to which chromosome a specific point belongs
	int size = active_region.size();

	if (size <= 1) return 0.0;	// there is nothing to do

	// calc initial score
	score_heatmap = calcScoreHeatmapActiveRegion();	// score from heatmap
	score_density = calcScoreDensity();
	score_curr = score_heatmap + score_density;
	//score_curr = score_density;


	score_prev = score_curr;
	milestone_score = score_curr;

	output(2, "initial score: %lf (heat=%lf, density=%lf); step=%f, steps=%d\n", score_curr, score_heatmap, score_density, step_size, Settings::MCstopConditionStepsHeatmapDensity);

	int milestone_cnt = 0;

	i = 1;
	while (true) {

		// select point to mutate and find chromosome on which it is located
		p = random(size);	// index in 'active_region'
		ind = active_region[p];	// index in 'clusters'

		if (clusters[ind].is_fixed) continue;	// check if bead is fixed (for example telomere)

		displacement = random_vector(step_size, Settings::use2D);	// generate random displacement vector
		clusters[ind].pos += displacement;


		score_heatmap = calcScoreHeatmapActiveRegion();
		score_density = calcScoreDensity();
		score_curr = score_heatmap + score_density;
		//score_curr = score_density;

		ok = score_curr <= score_prev;

		if (!ok && T > 0.0) {
			tp = Settings::tempJumpScaleHeatmapDensity * exp(-Settings::tempJumpCoefHeatmapDensity * (score_curr/score_prev) / T);
			//printf("!!!");
			ok = withChance(tp);
		}

		if (ok) {
			success++;
			milestone_success++;

			//output(2, "ok: %lf (heat=%lf, density=%lf)   %f\n", score_curr, score_heatmap, score_density, T);
			//if (success > 10) break;
		}
		else {
			clusters[ind].pos -= displacement;
			score_curr = score_prev;		// score doesn't change
			score_density = score_density_prev;
			score_heatmap = score_heatmap_prev;
		}

		T *= Settings::dtTempHeatmapDensity;    // cooling

		// check if we should stop
		if (i % Settings::MCstopConditionStepsHeatmapDensity == 0) {

			//score_density = calcScoreDensity();
			//score_curr = calcScoreHeatmapActiveRegion();

			output(2, "milestone [%d]: score = %lf, last = %lf (%f), T=%lf successes = %d, last = %d  (dens=%lf, heat=%lf), step=%f\n", milestone_cnt, score_curr,
					milestone_score, score_curr/milestone_score, T, success, milestone_success, score_density, score_heatmap, step_size);

			// stop if the improvement since the last milestone is less than 0.5%, or if the score is too small (may happen when heatmap is small and points are perfectly arranged)
			if ((score_curr > Settings::MCstopConditionImprovementHeatmapDensity * milestone_score &&
					milestone_success < Settings::MCstopConditionMinSuccessesHeatmapDensity) || score_curr < 1e-6) {
				break;
			}
			milestone_score = score_curr;
			milestone_success = 0;

			milestone_cnt++;

			repairDensityBoundary(); // NOTE
		}

		score_prev = score_curr;
		score_heatmap_prev = score_heatmap;
		score_density_prev = score_density;
		i++;
	}

	return score_curr;
}


void LooperSolver::setContactData(std::vector<string> chrs_list, BedRegion region_of_interest, string anchors, std::vector<string> factors,
		std::vector<string> arcs_clusters, std::vector<string> arcs_singletons, std::vector<string> arcs_singletons_inter) {

	printf("set contact data\n");

	if (factors.size()!=arcs_clusters.size()) error("Number of clusters is different then number of factors (see [data]->[clusters] and [factors])");

	this->arcs_clusters = arcs_clusters;
	this->arcs_singletons = arcs_singletons;
	this->arcs_singletons_inter = arcs_singletons_inter;
	this->selected_region = region_of_interest;
	this->is_bed_region_provided = region_of_interest.end > 0;

	if (is_bed_region_provided) {
		printf(" selected region: ");
		region_of_interest.print();
		arcs.selectRegion(region_of_interest);
	}
	else {
		printv(chrs_list, true, true, "selected chromosomes");
	}

	arcs.chrs = chrs_list;

	bool data_loaded = false;	// have we already loaded the data?

	printf(" load data\n");

	// first, try from cache file
	string input_file;
	if (Settings::useInputCache) {
		input_file = ftext("%sinput_cache.dat", output_dir.c_str());
		if (file_exists(input_file)) {
			printf(" read input data from cache file [%s]\n", input_file.c_str());
			if (arcs.fromFile(input_file)) data_loaded = true;
			else printf(" problem with cache file, recreate it...\n");
		}
	}
	else {
		printf("not using input cache\n");
	}

	if (!data_loaded) {

		arcs.loadAnchorsData(anchors);
		for (uint i = 0; i < factors.size(); ++i) arcs.loadPetClustersData(arcs_clusters[i], factors[i], segments_predefined);

		printf(" mark arcs\n");
		// arcs.markArcs(false);	// create cluster-indexed arcs
		arcs.parallelMarkArcs(false);

		arcs.removeEmptyAnchors();

		if (Settings::rewiringEnabled) arcs.rewire();

		printf(" mark arcs (again)\n");
		// arcs.markArcs(true);
		arcs.parallelMarkArcs(true);

		if (Settings::useInputCache) arcs.toFile(input_file);
	}

	chrs = chrs_list;
}


void LooperSolver::getSnapshot(string s) {
	return;
	chromosome_set.add(getCurrentChromosomes(), s);
}

void LooperSolver::addSnapshot(Chromosome chr, string s) {
	//chromosome_set.add(chr, s);
	//chr_set.add(getCurrentChromosome(), s);
}

std::map<std::string, Chromosome> LooperSolver::getCurrentChromosomes() {
	HierarchicalChromosome hc = getModel();
	hc.createCurrentLevelStructure();
	return hc.chr;
}

HierarchicalChromosome LooperSolver::getModel() {
	HierarchicalChromosome hc;
	hc.clusters = clusters;
	hc.arcs = arcs;
	hc.current_level = current_level;
	hc.chrs = chrs;
	//hc.root_index = root_index;
	hc.chr_root = chr_root;
	return hc;
}

int LooperSolver::findUnloopedChromatinLength(int cluster) {

	printf("find unlooped %d\n", cluster);
	int r = 0;
	int arcs = 0;	// arcs cnt
	for (int i = clusters[cluster].base_start; i <= clusters[cluster].base_end; ++i) {
		printf("base %d %d\n", i, arcs);
		for (size_t j = 0; j < base_clusters[i].siblings.size(); ++j) {
			//printf("arc %d\n", base_clusters[i].siblings[j]);
			if (base_clusters[i].siblings[j] > i) arcs++;
			else arcs--;
		}

		if (arcs == 0) {
			printf("gap %d %d\n", i, base_clusters[i].end);
		}
	}
	exit(0);
	return r;
}

// add empty clusters to densify large regions with no clusters
// assumes that the cluster tree is correctly built
void LooperSolver::densify() {

	// go through all leaves. if span between two is too large, add additional clusters every 'avg' bp

	//std::vector<int> vd;
	int lim = 600000;	// max allowed span
	int d;

	int base_level = clusters[0].level;	// find what is the level of leaves
	int base_clusters_cnt = 0;

	// find average span
	float avg = 0.0f;
	for (int i = 1; clusters[i].level == base_level; ++i) {	// go through all base level clusters
		d = clusters[i].start - clusters[i-1].end;
		avg += d;
		base_clusters_cnt++;
	}
	avg /= base_clusters_cnt;

	printf("average span: %f\n", avg);

	for (int i = 1; clusters[i].level == base_level; ++i) {	// go through all base level clusters
		//d = clusters[i].genomic_pos - clusters[i-1].genomic_pos;

		d = clusters[i].start - clusters[i-1].end;

		if (d > lim) {
			printf("large span: %d %d (%d %d)\n", i, d, clusters[i-1].genomic_pos, clusters[i].genomic_pos);

			// now add 'n' clusters between clusters 'i' and 'i-1'. set their parents to parents of 'i' and 'i-1', whichever is closer

			int add = floor(d / avg);		// how many clusters to add in between
			int shift = d / add;			// distance between succesive clusters
			int middle = clusters[i-1].end + d / 2;	 // find the middle point (used to decide to which parent to add)

			int p = clusters[i-1].end;
			int parent = clusters[i-1].parent;
			bool before_middle = true;
			//float dst = 1.0f / (add + 1);
			//float st = dst;

			//printf("dist=%d add=%d shift=%d mid=%d p=%d\n", d, add, shift, middle, p);
			std::vector<int> children;

			for (int j = 0; j < add; ++j) {
				p += shift;

				Cluster c(p, p);
				//c.pos = interpolate(clusters[active_region[i]].pos, clusters[active_region[i+1]].pos, st);
				c.level = base_level;
				c.parent = parent;
				clusters.push_back(c);

				// look for the switch of parent
				if (p > middle && before_middle) {
					printv(children, true, true, "insert to end");
					clusters[parent].children.insert(clusters[parent].children.end(), children.begin(), children.end());
					children.clear();
					parent = clusters[i].parent;
					before_middle = false;
				}

				children.push_back(clusters.size()-1);

				c.print();
			}

			printv(children, true, true, "insert to beg");
			clusters[parent].children.insert(clusters[parent].children.begin(), children.begin(), children.end());
		}
	}
}


std::vector<int> LooperSolver::findGaps(string chr) {
	std::vector<int> gaps;

	// for every anchor we want to count how many arcs are "above" it. if it's 0, than we have a gap

	int start = chr_first_cluster[chr];
	gaps.push_back(start);
	int arcs_cnt = 0;	// arcs cnt
	for (size_t i = start; i < clusters.size(); ++i) {
		for (size_t j = 0; j < clusters[i].arcs.size(); ++j) {		// check all arcs going in/out
			int other_end = otherEnd(chr, clusters[i].arcs[j], i);

			//printf("%d %d   (%d %d) %d\n", i, j, clusters[i].arcs[j], other_end, arcs_cnt);
			if ((int)i == other_end) continue;	// shouldn't happen

			if (other_end > (int)i) arcs_cnt++;	// check whether we start or end an arc
			else arcs_cnt--;
		}

		if (arcs_cnt == 0) {
			gaps.push_back(i);
			//printf("gap %d %d\n", i, (clusters[i+1].start+clusters[i].end)/2);
			//if (i+1 < base_clusters.size()) vector_insert_unique(gaps_pos, (base_clusters[i].end + base_clusters[i+1].start) / 2);
			//length.push_back(base_clusters[i].end - prev_pos);
			//prev_pos = base_clusters[i+1].start;
		}
	}

	vector_insert_unique(gaps, clusters.size()-1);
	//saveToCSV(ftext("%sloops_length_%s.csv", output_dir.c_str(), label.c_str()), length);
	return gaps;
}



// find the optimal split of gaps, so that the bins have approximately size of 'exp_size'
// 'gaps' should start with 0 (gaps_pos = leftmost position) and end with index of last gap (gap pos = rightmost pos)
std::vector<int> LooperSolver::findSplit(std::vector<int> gaps, int exp_size, string chr) {

	if (clusters.size() < gaps.size())
		error(ftext("find split: size mismatch, gaps=%d, clusters=%d", gaps.size(), clusters.size()));

	// if we do have predefined segment split use it
	// just iterate over all gaps, find the gap span and then check if there is a split defined in this region
	// we only consider the start coordinate of bed region
	if (segments_predefined.regions.size() > 0) {
		vector <int> splits;
		splits.push_back(gaps[0]);

		// index of next predefined segment (refers to segments_predefined)
		// we need to move the index to the correct position (important when there are multiple chromosomes)
		int curr_ind = -1;
		int last_seg_index = segments_predefined.regions.size() - 1;
		for (size_t i = 0; i < segments_predefined.regions.size(); ++i) {
			if (segments_predefined.regions[i].chr == chr) {
				if (curr_ind == -1) curr_ind = i;
			}
			else if (curr_ind != -1) {
				last_seg_index = i-1;
				break;
			}
		}

		//printf("curr ind = %d\n", curr_ind);
		if (curr_ind == -1) error(ftext("no segment split found for %s", chr.c_str()));

		for (size_t i = 1; i+1 < gaps.size(); ++i) {

			int gap_start = clusters[gaps[i]].end;
			int gap_end = clusters[gaps[i]+1].start;
			int seg_break = curr_ind <= last_seg_index ? segments_predefined.regions[curr_ind].start : -1;

			//printf("gap: %d %d   %d\n", gap_start, gap_end, seg_break);

			// MOD
			//if (seg_break>=0 && seg_break < gap_start) error("seg break > gap start");
			//if (seg_break>=0 && seg_break < gap_start) seg_break++;
			while (seg_break>=0 && seg_break < gap_start) {
				curr_ind++;
				seg_break = curr_ind <= last_seg_index ? segments_predefined.regions[curr_ind].start : -1;
			}

			if (seg_break >= gap_start && seg_break <= gap_end) {
				curr_ind++;
				splits.push_back(gaps[i]);
			}
		}

		vector_insert_unique(splits, gaps[gaps.size()-1]);
		return splits;
	}

	vector<int> L; // lengths of current regions (segments)
	vector<int> S; // lengths of spans between segments

	// init L and S
	for (size_t i = 1; i < gaps.size(); ++i) {
		int prev_cluster = i==1 ? 0 : gaps[i-1]+1;
		L.push_back(clusters[gaps[i]].end - clusters[prev_cluster].start);
		if (i+1<gaps.size()) S.push_back(clusters[gaps[i]+1].start - clusters[gaps[i]].end);
	}

	printv(L, true, true, "L");
	printv(S, true, true, "S");

	for (size_t i = 1; i < gaps.size(); ++i) {
		//printf("%d\n", i);
		int prev_cluster = i==1 ? 0 : gaps[i-1]+1;

		int next_cluster = i == gaps.size()-1 ? gaps[gaps.size()-1] : gaps[i]+1;
		int center = i+1 < gaps.size() ? (clusters[gaps[i]].end+clusters[next_cluster].start)/2 : 0;
		//printf("%d %d\n", clusters[gaps[i]].end, clusters[next_cluster].start);
		//printf("(%d - %d), gap span=%d, span center=%d\n", clusters[prev_cluster].start, clusters[gaps[i]].end, S[i-1], center);
	}

	return gaps;
}


void LooperSolver::createTreeGenome() {

	// create a true root cluster
	// root has no genomic position, as it includes (potentialy) different chromosomes
	// for the same reason it doesn't have a base start/end
	//Cluster root;
	//root.level = 0;

	// create trees for chromosomes (put them one by one on clusters[])
	for (string chr: chrs) {
		int chr_root_index = createTreeChromosome(chr);

		chr_root[chr] = chr_root_index;
		//root.children.push_back(chr_root_index);	// add root of chr tree as a child
	}

	//root_index = clusters.size(); // set true root index

	// update parent for root's children
	//for (int i = 0; i < root.children.size(); ++i) clusters[root.children[i]].parent = clusters.size();
	//clusters.push_back(root);
}

int LooperSolver::createTreeChromosome(string chr) {
	int cluster_start = clusters.size();	// index of first cluster for a given chromosome
	chr_first_cluster[chr] = cluster_start;

	// levels: root=0, chromosome=1, segment=2, interaction block=3, anchor=4
	int curr_level = 4;

	// create clusters for anchors
	for (int i = 0; i < arcs.anchors_cnt[chr]; ++i) {
		Cluster c(arcs.anchors[chr][i].start, arcs.anchors[chr][i].end);
		c.orientation = arcs.anchors[chr][i].orientation;
		c.base_start = cluster_start + i;
		c.base_end = cluster_start + i;
		c.level = curr_level;
		clusters.push_back(c);
	}

	for (int i = 0; i < arcs.arcs_cnt[chr]; ++i) {

		// arcs.[start,end] should refer to global indices of clusters
		// thus, we need to shift indices
		arcs.arcs[chr][i].start += cluster_start;
		arcs.arcs[chr][i].end += cluster_start;

		// clusters.arcs refer to local, chromosome based indices of arcs
		clusters[arcs.arcs[chr][i].start].arcs.push_back(i);
		clusters[arcs.arcs[chr][i].end].arcs.push_back(i);
	}


	Cluster rootc;	// create root node

	// find gaps - spots with no arcs above them
	printf("find gaps (%s)\n", chr.c_str());
	std::vector<int> gaps = findGaps(chr);
	printv(gaps, true, true, "gaps");

	printf("find splits\n");
	std::vector<int> splits = findSplit(gaps, Settings::segmentSize, chr);
	printv(splits, true, true, "splits");

	// we have something like:
	//    gaps: [total: 49]  18950 18956 18962 18968 18982 18991 19001 19005 19014 19031 19038 19046 19065 19073 19085 19087 19091 19102 19109 19126 19133 19134 19135 19139 19156 19165 19172 19181 19183 19191 19198 19205 19213 19230 19267 19274 19276 19295 19299 19302 19304 19309 19317 19334 19346 19354 19360 19363 19366
	//    splits: [total: 31]  18950 18962 18982 18991 19001 19014 19038 19046 19065 19073 19087 19091 19102 19109 19133 19139 19156 19181 19205 19213 19230 19267 19276 19299 19304 19309 19317 19334 19346 19354 19366
	// gaps refers to interaction blocks, splits - to segments
	// to create clusters for interaction blocks (ib) we simply include all clusters between consecutive gaps
	// to create clusters for segments we need to refer to ib indices.

	//int first_gap_index = clusters.size(); // index of first ib cluster (used when creating segment clusters)

	int next_split_ind = 1;
	int start_gap_index = clusters.size();


	// mark interactions blocks
	curr_level--;
	for (size_t i = 1; i < gaps.size(); ++i) {

		// create cluster contained within clusters (gaps[i-1], gaps[i])
		int prev_gap = i==1 ? gaps[i-1] : gaps[i-1] + 1;	// boundary gaps (both are inclusive)
		int curr_gap = gaps[i];

		int start_pos = clusters[prev_gap].start;		// genomic position of start/end
		int end_pos = clusters[curr_gap].end;

		//printf("new cluster [%d], gaps=(%d, %d), pos=(%d %d)\n", i, prev_gap, curr_gap, start_pos, end_pos);

		Cluster c(start_pos, end_pos);	// interaction block cluster (anywhere from 100kb to ~Mb)
		c.base_start = prev_gap;
		c.base_end = curr_gap;
		c.level = curr_level;

		// set as a parent for corresponding anchors
		for (int k = prev_gap; k <= curr_gap; ++k) {
			c.children.push_back(k);
			clusters[k].parent = clusters.size();
		}

		clusters.push_back(c);	// add it to the list

		if (gaps[i] == splits[next_split_ind]) {

			int end_gap_index = clusters.size()-1;

			Cluster cs(clusters[start_gap_index].start, clusters[end_gap_index].end);	// segment cluster
			cs.base_start = clusters[start_gap_index].base_start;
			cs.base_end = clusters[end_gap_index].base_end;
			cs.level = curr_level-1;	// we are 1 level higher now

			int segment_cluster_ind = clusters.size();	// cluster-index that will be given to segment cluster that is currently being created

			// set as a parent for corresponding ib
			for (int k = start_gap_index; k <= end_gap_index; ++k) {
				cs.children.push_back(k);
				clusters[k].parent = segment_cluster_ind;
			}

			rootc.children.push_back(segment_cluster_ind);
			clusters.push_back(cs);		// add it to the list

			start_gap_index = clusters.size();
			next_split_ind++;
		}
	}

	// make chromosome root
	curr_level -= 2;
	int root_ind = clusters.size();	// save root index

	rootc.level = curr_level;
	rootc.start = clusters[rootc.children[0]].start;
	rootc.end = clusters[rootc.children[rootc.children.size()-1]].end; // TODO: was start, should be?
	rootc.genomic_pos = (rootc.start + rootc.end) / 2;

	rootc.base_start = clusters[rootc.children[0]].base_start;
	rootc.base_end = clusters[rootc.children[rootc.children.size()-1]].base_end;

	for (size_t i = 0; i < rootc.children.size(); ++i) clusters[rootc.children[i]].parent = clusters.size();
	clusters.push_back(rootc);
	return root_ind;
}

// set top level of a structure - a root for every chromosome
void LooperSolver::setTopLevel() {
	current_level.clear();
	for (string chr: chrs) {
		current_level[chr].push_back(chr_root[chr]);
	}
}

void LooperSolver::setLevel(int level) {
	setTopLevel();
	while (level--) levelDown();
}

Heatmap LooperSolver::downsampleSingletonHeatmap(Heatmap &heat) {

	printf("downsampleSingletonHeatmap\n");

	// calculate size of the heatmap
	int n = 0;
	for (string chr: chrs) n += current_level[chr].size();

	printf("downsamples heatmap size: %d\n", n);
	Heatmap h(n);

	std::map<std::string, std::vector<int> > children_ind;	// index of first child for every chromosome
	for (string chr: chrs) {
		children_ind[chr].push_back(0);
		for (size_t i = 0; i < current_level[chr].size(); ++i) children_ind[chr].push_back(children_ind[chr][i] + clusters[current_level[chr][i]].children.size());
	}
	//for (int i = 0; i < n; ++i) children_ind.push_back(children_ind[i] + clusters[current_level[i]].children.size());
	printm(children_ind, false, "child starting index");

	// for all chromosome pairs
	int indA = 0, indB = 0;
	for (string chrA: chrs) {
		size_t cntA = current_level[chrA].size();
		indB = 0;
		for (string chrB: chrs) {
			size_t cntB = current_level[chrB].size();
			// go through all 16Mb regions (ie. current level)
			for (size_t i = 0; i < cntA; ++i) {
				for (size_t j = i; j < cntB; ++j) {
					// take average value from 2Mb level heatmap for all children
					// (cell (i,j) of current region correspond to a square subset of smaller heatmap: [a,b]x[p,q])
					size_t cnt = clusters[current_level[chrA][i]].children.size();
					size_t cnt2 = clusters[current_level[chrB][j]].children.size();
					float avg = 0.0f;
					for (size_t k = 0; k < cnt; ++k) {
						for (size_t l = 0; l < cnt2; ++l) 	{
							//printf("%s %s; %d %d; %d %d, ind=%d %d\n", chrA.c_str(), chrB.c_str(), i, j, k, l, indA, indB);
							avg += heat.v[children_ind[chrA][i]+k][children_ind[chrB][j]+l];
						}
					}
					avg /= cnt * cnt2;
					h.v[indB+j][indA+i] = h.v[indA+i][indB+j] = avg;
				}
			}
			indB += cntB;
		}
		indA += cntA;
	}
	h.clearDiagonal(1);
	return h;
}


// if selected region is small, it may happen that it is contained within one segment level bead. Then we can not
// create a heatmap
Heatmap LooperSolver::createSingletonHeatmap(int diag) {

	// find boundary positions for bins
	std::map<std::string, std::vector<int> > breaks;
	std::map<std::string, int> start_ind;	// starting index for chromosomes (index of first column in a heatmap)
	int pos;

	int curr_ind = 0;
	for (string chr: chrs) {
		breaks[chr].push_back(0);

		for (size_t i = 0; i+1 < current_level[chr].size(); i++) {
			pos = (clusters[current_level[chr][i]].end + clusters[current_level[chr][i+1]].start) / 2;
			breaks[chr].push_back(pos);
		}
		//breaks[chr].push_back()
		breaks[chr].push_back(1e9);

		start_ind[chr] = curr_ind;
		curr_ind += breaks[chr].size() - 1;
	}

	printm(breaks, true, "split positions");
	printm(start_ind, chrs, true, true, "start indices");

	Heatmap h;
	h.init(curr_ind);

	// create a set of chromosomes to speed up checking if chromosome is to be processed
	std::set<std::string> chrs_set;
	for (std::string chr: chrs) chrs_set.insert(chr);

	const int max_position = 2e9;

	// we need to keep min and max genomic position for every chromosome (to calculate lengths of first and last bin)
	std::map<std::string, int> chr_max_pos;
	std::map<std::string, int> chr_min_pos;
	for (string chr: chrs) {
		chr_max_pos[chr] = 0;
		chr_min_pos[chr] = max_position;
	}

	bool is_region_selected = selected_region.end > 0;

	char chr1[10], chr2[10];
	int sta, stb, enda, endb, sc;
	char buf_tmp[1024];

	FILE *f = NULL;

	vector<std::string> singleton_files;
	singleton_files.insert(singleton_files.end(), arcs_singletons.begin(), arcs_singletons.end());
	if (chrs.size() > 1) singleton_files.insert(singleton_files.end(), arcs_singletons_inter.begin(), arcs_singletons_inter.end());

	// if split by chr and single chromosome, then update names (update is a bit more problematic for multiple chromosomes)
	if (Settings::dataSplitSingletonFilesByChr && chrs.size()==1) {
		string path, fname, file_chr;
		for (uint i = 0; i < singleton_files.size(); ++i) {
			split_file_path(singleton_files[i], path, fname);
			singleton_files[i] = ftext("%schr/%s.%s", path.c_str(), fname.c_str(), chrs[0].c_str());
		}
	}

	// read data from files
	for (uint fi = 0; fi < singleton_files.size(); ++fi) {
		printf("\nread [%s]\n", singleton_files[fi].c_str());

		f = open(singleton_files[fi].c_str(), "r");
		if (!f) error("could not read the file!");

		int cnt_arcs_ok = 0;  // number of arcs added
		int cnt_arcs_diag = 0;  // number of arcs that fall on the diagonal
		int cnt_arcs_bad_chr = 0; // number of arcs on other chromosomes (not specified with -c or chrM)
		int cnt_lines_err = 0;  // number of lines that could not be read
		int k = 0;
		while (!feof(f)) {
			if (k++ % 1000000 == 0) printf(".");

			// * there are several possible formats of files, e.g.:
			// 1: chr9	77361535	77361615	chr9	78624184	78624253	1
			// 2: chr18	51650936	51651093	chr3	110168590	110168697	2	0	0
			// 3: chr6	34306237	34306388	chrX	109134352	109134456	1	+	+
			// we read first 7 column, and then fgets() to reach the end of line

			int args_read = fscanf(f, "%s %d %d %s %d %d %d", chr1, &sta, &stb, chr2, &enda, &endb, &sc);
			if (args_read < 7) {
				// this is usually true only for the last line in the file, if it is empty
				// reporting that as an error may be misleading, so only report non-cmoplete and non-empty lines
				if (args_read > 0) cnt_lines_err++;
				continue;
			}

			fgets(buf_tmp, 1024, f);

			// check whether we should process current chromosomes
			if (chrs_set.find(chr1) == chrs_set.end() || chrs_set.find(chr2) == chrs_set.end()) {
				cnt_arcs_bad_chr++;
				continue;
			}

			// check specified region
			if (is_region_selected && (!selected_region.contains(sta) || !selected_region.contains(endb))) continue;

			if (sta < 0 || stb < 0 || enda < 0 || endb < 0) {
				printf("negative positions! %d %d, %d %d)\n", sta, stb, enda, endb);
				if (sta < 0) sta = 0;
				if (enda < 0) enda = 0;
			}

			// find indices of start and end of the arc.
			// 1. find correct index within a chromosome by iterating over breaks for this chromosome
			// 2. add 'shift' for a given chromosomes using 'start_ind'
			int st = -1;
			int end = -1;
			for (size_t i = 0; i+1 < breaks[chr1].size(); ++i) {
				if (breaks[chr1][i] <= sta && sta <= breaks[chr1][i+1]) {
					st = start_ind[chr1] + i;
					break;
				}
			}
			for (size_t i = 0; i+1 < breaks[chr2].size(); ++i) {
				if (breaks[chr2][i] <= enda && enda <= breaks[chr2][i+1]) {
					end = start_ind[chr2] + i;
					break;
				}
			}

			if (st == -1 || end == -1) {
				printf("non matching arc! %d %d  (%d %d)\n", st, end, sta, enda);
				exit(0);
			}

			if (st == end) {
				cnt_arcs_diag++;
				continue;   // ignore diagonal
			}

			//printf("%d %d  %d %d  %s %s   %d %d\n", sta, stb, enda, endb, chr1, chr2, st, end);

			// update min/max position for chromosomes
			chr_max_pos[chr1] = max(chr_max_pos[chr1], stb);
			chr_min_pos[chr1] = min(chr_min_pos[chr1], sta);
			chr_max_pos[chr2] = max(chr_max_pos[chr2], endb);
			chr_min_pos[chr2] = min(chr_min_pos[chr2], enda);

			// TODO: we might want to use different scoring (something non-linear?)
			h.v[st][end] += sc;
			h.v[end][st] += sc;

			cnt_arcs_ok++;
		}

		fclose(f);

		printf("\n%d interactions added, ignored %d (bad chr) and %d (diagonal)", cnt_arcs_ok, cnt_arcs_bad_chr, cnt_arcs_diag);
		if (cnt_lines_err > 0) printf(", failed to read %d lines", cnt_lines_err);
		printf("\n");
	}

	// now add long-rang interactions effect
	printf("add long-range interactions\n");
	for (string chr: chrs) {
		printf("%s %d\n", chr.c_str(), (int)arcs.long_arcs[chr].size());
		for (InteractionArc arc: arcs.long_arcs[chr]) {
			// find bins for the current arc

			sta = arc.start;
			enda = arc.end;

			int st = -1;
			int end = -1;
			for (size_t i = 0; i+1 < breaks[chr].size(); ++i) {
				if (breaks[chr][i] <= sta && sta <= breaks[chr][i+1]) {
					st = start_ind[chr] + i;
					break;
				}
			}
			for (size_t i = 0; i+1 < breaks[chr].size(); ++i) {
				if (breaks[chr][i] <= enda && enda <= breaks[chr][i+1]) {
					end = start_ind[chr] + i;
					break;
				}
			}

			if (st == end) continue;

			// calc effect of the interaction
			float val = Settings::longPETClustersEffectScale * pow(arc.score, Settings::longPETClustersEffectPower);

			//arc.print();
			//printf("st/end = %d %d, val = %f\n", st, end, val);

			h.v[st][end] += val;
			h.v[st][end] += val;
		}
	}

	// update breaks with actual min/max position
	printf("\nupdate min/max positions\n");
	for (string chr: chrs) {
		printf("%s: %d-%dbp\n", chr.c_str(), chr_min_pos[chr], chr_max_pos[chr]);
		breaks[chr][0] = chr_min_pos[chr];
		breaks[chr][breaks[chr].size()-1] = chr_max_pos[chr];

		if (chr_min_pos[chr] == max_position) error("no singletons found for chromosome");
	}

	// calculate lengths of bins (ie. genomic span)
	std::vector<float> len;	// length of bins
	for (string chr: chrs) {
		for (size_t i = 0; i+1 < breaks[chr].size(); ++i) {
			int d = getGenomicLengthExcludingCentromere(chr, breaks[chr][i], breaks[chr][i+1]);
			len.push_back((float)d / 1e6);
		}
	}

	printv(len, true, true, "bins length");

	// normalize by bin sizes

	for (size_t i = 0; i < h.size; ++i) {
		for (size_t j = i+1; j < h.size; ++j) {
			float sc = (float)h.v[i][j] / (len[i] * len[j]);
			h.v[i][j] = sc;
			h.v[j][i] = sc;
		}
	}
	h.clearDiagonal(diag);
	return h;
}

void LooperSolver::createSingletonSubanchorHeatmap(string chr, vector<int> &anchors_gap_len) {
	int anchor_margin = 0;

	// find boundary positions for bins
	std::vector<int> breaks;

	// boundaries of the current region
	int region_start = clusters[active_region[0]].start - anchor_margin;
	int region_end = clusters[active_region[active_region.size()-1]].end + anchor_margin;

	printf("heatmap boundary: %d %d\n", region_start, region_end);

	float dt = 1.0 / (Settings::loopDensity);
	breaks.push_back(region_start);
	anchors_gap_len.push_back(clusters[active_region[0]].end - clusters[active_region[0]].start);

	for (size_t i = 1; i < active_region.size(); i++) {
		int span_start = clusters[active_region[i-1]].end + anchor_margin;
		int span_end = clusters[active_region[i]].start - anchor_margin;
		int len = span_end - span_start;

		int anchor_len = clusters[active_region[i]].end - clusters[active_region[i]].start;
		//printf("%d %d %d %d\n", span_start, span_end, len, anchor_len);
		anchors_gap_len.push_back(len);
		anchors_gap_len.push_back(anchor_len);

		breaks.push_back(span_start);
		for (int j=0, t=dt; j<Settings::loopDensity-1; j++, t+=dt) {
			breaks.push_back(span_start + len * t);
		}
		breaks.push_back(span_end);
	}
	breaks.push_back(region_end);

	heatmap_anchor.init(active_region.size());
	heatmap_subanchor.init(breaks.size()-1);

	// read data from files
	string path, fname;
	string file_chr;
	char chr1[10], chr2[10];
	int st, end;
	int sta, stb, enda, endb, st_pos, end_pos, sc;
	char buf_tmp[1024];
	FILE *f = NULL;

	for (size_t fi = 0; fi < arcs_singletons.size(); ++fi) {
		// ignore files with inter contacts (for anchor level we need intra only)
		//if (arcs_singletons[fi].find("inter") != std::string::npos) continue;
		if (Settings::dataSplitSingletonFilesByChr) {
			split_file_path(arcs_singletons[fi], path, fname);
			file_chr = ftext("%schr/%s.%s", path.c_str(), fname.c_str(), chr.c_str());
		}
		else file_chr = arcs_singletons[fi];
		printf("read [%s]\n", file_chr.c_str());

		f = open(file_chr.c_str(), "r");
		if (!f) {
			if (Settings::dataSplitSingletonFilesByChr) error("No chromosome-splitted file found!");
			error("could not read the file!");
		}

		int cnt_arcs_ok = 0;  // number of arcs added
		int cnt_arcs_diag = 0;  // number of arcs that fall on the diagonal
		int cnt_lines_err = 0;  // number of lines that could not be read

		while (!feof(f)) {

			// * there are several possible formats of files, e.g.:
			// 1: chr9	77361535	77361615	chr9	78624184	78624253	1
			// 2: chr18	51650936	51651093	chr3	110168590	110168697	2	0	0
			// 3: chr6	34306237	34306388	chrX	109134352	109134456	1	+	+
			// we read first 7 column, and then fgets() to reach the end of line

			if (fscanf(f, "%s %d %d %s %d %d %d", chr1, &sta, &stb, chr2, &enda, &endb, &sc) < 7) {
				cnt_lines_err++; // something went wrong (possibly empty last line in a file)
				continue;
			}

			fgets(buf_tmp, 1024, f);

			st_pos = (sta + stb) / 2;
			end_pos = (enda + endb) / 2;

			// check specified region
			if (st_pos < region_start || end_pos > region_end) continue;

			// find indices of start and end of the arc
			st = -1;
			end = -1;
			for (size_t i = 1; i < breaks.size() && (st==-1||end==-1); ++i) {
				if (st == -1 && st_pos <= breaks[i]) st = i - 1;
				if (end == -1 && end_pos <= breaks[i]) end = i - 1;
			}

			if (st == -1 || end == -1) {
				printf("non matching arc! %d %d   %d %d  (%d %d)\n", st, end, st_pos, end_pos, sta, enda);
				exit(0);
			}

			if (st == end) {
				cnt_arcs_diag++;
				continue;   // ignore diagonal
			}

			//printf("%d %d  %d %d  %s %s   %d %d\n", sta, stb, enda, endb, chr1, chr2, st, end);

			// TODO: we might want to use different scoring (something non-linear?)
			heatmap_subanchor.v[st][end] += sc;
			heatmap_subanchor.v[end][st] += sc;

			cnt_arcs_ok++;
		}

		fclose(f);
		printf("%d interactions added, %d ignored (diagonal); failed to read %d lines\n", cnt_arcs_ok, cnt_arcs_diag, cnt_lines_err);
	}

	int n = heatmap_subanchor.size;

	vector<bool> is_anchor(n, false);
	for (int i = 0; i < n; ++i) if (i % (Settings::loopDensity+1) == 0) is_anchor[i] = true;


	// create anchor heatmap
	float val;
	int n_anchor = heatmap_anchor.size;
	for (int i = 0; i < n_anchor; ++i) {
		for (int j = i+1; j < n_anchor; ++j) {
			val = heatmap_subanchor.v[i*(Settings::loopDensity+1)][j*(Settings::loopDensity+1)];
			val /= (anchors_gap_len[2*i] * anchors_gap_len[2*j]) / 1000000.0f;  // normalize by anchor length
			heatmap_anchor.v[j][i] = heatmap_anchor.v[i][j] = val;
		}
	}

	// normalize by anchor/loop counts
	// first, find count of singletons anchored in anchors and loops
	float end_anchor = 0.0f, end_loop = 0.0f, end_loop_anchor = 0.0f;
	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			if (is_anchor[i] && is_anchor[j]) end_anchor += heatmap_subanchor.v[i][j];
			else if (!is_anchor[i] && !is_anchor[j]) end_loop += heatmap_subanchor.v[i][j];
			else end_loop_anchor += heatmap_subanchor.v[i][j];
		}
	}

	float tot = end_anchor + end_loop_anchor + end_loop;
	printf("types [anchor/loop]: %.0f %.0f %.0f\n", end_anchor, end_loop_anchor, end_loop);

	if (tot < 1e-6) {
		printf("subanchor heatmap is empty\n");
		return;	// heatmap is empty
	}

	// normalize by average count
	float avg = heatmap_subanchor.getAvg();
	printf("avg count: %f\n", avg);
	if (avg < 1e-6) return;

	heatmap_subanchor.scale(1.0f / avg);

	// normalize by bin sizes
	printf("bin len norm: ");
	float d1, d2;
	int i1, i2;
	for (int i = 0; i < n; ++i) {
		i1 = i/(Settings::loopDensity+1);
		i1 *= 2;
		d1 = i % (Settings::loopDensity+1) == 0 ? anchors_gap_len[i1] : ((float)anchors_gap_len[i1+1] / Settings::loopDensity);
		d1 /= 1000.0f;
		printf("%f, ", d1);
		for (int j = i+1; j < n; ++j) {
			i2 = j/(Settings::loopDensity+1);
			i2 *= 2;
			d2 = j % (Settings::loopDensity+1) == 0 ? anchors_gap_len[i2] : ((float)anchors_gap_len[i2+1] / Settings::loopDensity);
			d2 /= 1000.0f;
			float sc = (float)heatmap_subanchor.v[i][j] / (d1 * d2);
			heatmap_subanchor.v[i][j] = sc;
			heatmap_subanchor.v[j][i] = sc;
		}
	}
}

Heatmap LooperSolver::normalizeHeatmap(const Heatmap &heat) {
	int i,j;
	int n = heat.size;
	Heatmap h;
	h.init(n);
	h.diagonal_size = heat.diagonal_size;

	// assumption: total number of reads should be equal for every region
	// find expected number of reads and multiply every row by R[i] so that the new sum is equal to the expected
	// as expected number of reads use the average

	float expected_sum = 0.0f;
	for (i = 0; i < n; ++i) {	// for every row...
		float sum = 0.0f;
		for (j = 0; j < n; ++j) sum += heat.v[i][j]; // ... find the sum of counts
		expected_sum += sum;
	}
	expected_sum /= n;	// take the average

	// now go through rows
	for (i = 0; i < n; ++i) {	// for every row...

		// ... find the sum of counts
		float sum = 0.0f;
		for (j = 0; j < n; ++j) sum += heat.v[i][j];

		float mn = expected_sum / sum;

		for (j = 0; j < n; ++j) h.v[i][j] = heat.v[i][j] * mn;
	}

	for (i = 0; i < n; ++i) {	// for every row...
		for (j = i+1; j < n; ++j) {
			h.v[i][j] = (h.v[i][j] + h.v[j][i]) / 2.0f;
			h.v[j][i] = h.v[i][j];
		}
	}
	return h;
}

void LooperSolver::createDistanceHeatmap(int heatmap_ind) {
	int n = heatmap[heatmap_ind].size;
	heatmap_dist.init(n);
	heatmap_dist.diagonal_size = heatmap[heatmap_ind].diagonal_size;

	printf(" heatmap size: %d, diag = %d\n", n, heatmap_dist.diagonal_size);

	float val = 0.0f;
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			val = heatmap[heatmap_ind].v[i][j];
			if (val < 1e-6) heatmap_dist.v[i][j] = 0.0f;
			else {
				heatmap_dist.v[i][j] = (abs(i-j) < heatmap_dist.diagonal_size) ? -1.0f :
						(heatmap_ind == 0 ? freqToDistanceHeatmapInter(val) : freqToDistanceHeatmap(val));
			}
			heatmap_dist.v[j][i] = heatmap_dist.v[i][j];
		}
	}

	// check heatmap values
	//float avg = heatmap[heatmap_ind].getAvg();
	float avg = heatmap_dist.getAvg();
	float max_dist = avg * Settings::heatmapDistanceHeatmapStretching;
	printf(" dist heatmap postprocesing. avg = %f, max = %f\n", avg, max_dist);
	printf(" limit large values...");
	int tmp = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (heatmap_dist.v[i][j] > max_dist) {
				printf(" %f ", heatmap_dist.v[i][j]);
				heatmap_dist.v[i][j] = max_dist;
				tmp++;
				if (tmp % 50 == 0) printf("\n");
			}
		}
	}
	printf("\n");
}

void LooperSolver::levelDown() {
	std::vector<int> tmp;
	for (string chr: chrs) {
		for (size_t i = 0; i<current_level[chr].size(); ++i) {
			// if cluster has no children, then keep it
			if (clusters[current_level[chr][i]].children.size() == 0) {
				tmp.push_back(current_level[chr][i]);
			}
			else {
				for (size_t j = 0; j < clusters[current_level[chr][i]].children.size(); ++j) {
					tmp.push_back(clusters[current_level[chr][i]].children[j]);
				}
			}
		}
		current_level[chr] = tmp;
		tmp.clear();
	}
}

void LooperSolver::normalizeHeatmapInter(Heatmap &heat, float scale) {
	if (chrs.size() <= 1) return;	// there is only one chr, scaling inter has no sense

	printf("normalize heatmap inter (scale = %f)\n", scale);

	int j, k;

	// find boundary positions for bins
	std::map<std::string, int> start_ind;	// starting index for chromosomes (index of first column in a heatmap)

	int curr_ind = 0;

	string chr;
	for (size_t j=0; j<chrs.size(); j++) {
		chr = chrs[j];
		start_ind[chr] = curr_ind;
		curr_ind += current_level[chr].size();
	}
	start_ind["end"] = curr_ind;

	std::vector<string> chrs_list = chrs;
	chrs_list.push_back("end");

	printf(" hmap size = %lu\n", heat.size);
	printm(start_ind, chrs, true, true, "start indices");
	printv(chrs_list, true, true, "chrs_list");

	heat.scale(scale);

	for (size_t i = 0; i+1 < chrs_list.size(); ++i) {
		for (k=start_ind[chrs_list[i]]; k<start_ind[chrs_list[i+1]]; k++) {
			for (j=start_ind[chrs_list[i]]; j<start_ind[chrs_list[i+1]]; j++) {
				heat.v[k][j] /= scale;
			}
		}
	}
}

void LooperSolver::normalizeHeatmapDiagonalTotal(Heatmap &heat, float val) {
	printf("normalize heatmap for near-diagonal total (%f)\n", val);

	size_t diag = heat.getDiagonalSize();
	printf(" diag size = %lu\n", diag);

	// calculate avg of first non-zero diagonal
	double avg = 0.0;
	for (size_t i = 0; i+diag < heat.size; ++i) avg += heat.v[i][i+diag];
	avg /= heat.size - diag;
	printf(" diagonal avg = %lf\n", avg);

	float mn = val / avg;
	for (size_t i = 0; i < heat.size; ++i) {
		for (size_t j = 0; j < heat.size; ++j) {
			heat.v[i][j] *= mn;
		}
	}
}


float LooperSolver::getHeatmapAvgBetweenNeighboringRegions(Heatmap &heat, int level) {
	int i,j,k;
	int n;
	float global_avg = 0.0;

	setLevel(level);
	for (string chr: chrs) {
		// for every intra-heatmap

		vector<int> v;
		v.push_back(0);
		for (size_t l = 0; l < current_level[chr].size(); ++l) v.push_back(v[l] + clusters[current_level[chr][l]].children.size());
		n = v.size();
		printv(v, true, true, "hmap avg v");

		if (n > 2) {
			float avg = 0.0f;
			for (i = 1; i+1 < n; ++i) {	// go through pairs of neighbors
				float curr_avg = 0.0f;
				for (j = v[i-1]; j < v[i]; ++j) for (k = v[i]; k < v[i+1]; ++k) curr_avg += heat.v[j][k];
				curr_avg /= (v[i] - v[i-1] + 1) * (v[i+1] - v[i] + 1);
				avg += curr_avg;
			}

			avg /= n-2;
			printf("local avg (%s) = %f\n", chr.c_str(), avg);

			global_avg += avg;
		}
	}

	global_avg /= chrs.size();

	printf("avg for neighboring regions = %f\n", global_avg);
	return global_avg;
}


double LooperSolver::calcScoreDistancesActiveRegion() {
	double sc = 0.0, diff;
	int i,j;
	int n = active_region.size();
	vector3 v;

	for (i = 0; i < n; ++i) {
		for (j = i+1; j < n; ++j) {

			v = clusters[active_region[i]].pos - clusters[active_region[j]].pos;

			// anchors not connected by an arc are denoted by -1
			// add penalty if they are close to each other
			if (heatmap_exp_dist_anchor.v[i][j] < 0.0f) {
				sc += 1.0f / v.length();
				continue;
			}

			if (heatmap_exp_dist_anchor.v[i][j] < 1e-6) continue; // ignore near-zero values (instability)


			diff = (v.length() - heatmap_exp_dist_anchor.v[i][j]) / heatmap_exp_dist_anchor.v[i][j];
			//printf("%d %d  -> (%lf %lf) %lf\n", i, j, v.length(), heatmap_exp_dist_anchor.v[i][j], diff);
			sc += diff * diff * (diff >= 0.0f ? Settings::springConstantStretchArcs : Settings::springConstantSqueezeArcs);
		}
	}
	return sc;
}

// calculate score based on distances between the selected bead and clusters connected with it by arcs
double LooperSolver::calcScoreDistancesActiveRegion(int cluster_moved) {
	double sc = 0.0, diff;
	int i, st;
	st = active_region[cluster_moved];
	int n = active_region.size();
	vector3 v;
	for (i = 0; i < n; ++i) {
		if (i == cluster_moved) continue;

		v = clusters[st].pos - clusters[active_region[i]].pos;

		if (heatmap_exp_dist_anchor.v[i][cluster_moved] < 0.0f) {
			//	sc += 1.0f / v.length();
			//sc += 1.0f;
			continue;
		}

		if (heatmap_exp_dist_anchor.v[i][cluster_moved] < 1e-6) continue;	// anchors not connected by an arc are denoted by -1

		diff = (v.length() - heatmap_exp_dist_anchor.v[i][cluster_moved]) / heatmap_exp_dist_anchor.v[i][cluster_moved];
		//printf("%d %d (%d %d)  -> (%lf %lf) %lf\n", i, cluster_moved, st, active_region[i], v.length(), heatmap_exp_dist_anchor.v[i][cluster_moved], diff);
		sc += diff * diff * (diff >= 0.0f ? Settings::springConstantStretchArcs : Settings::springConstantSqueezeArcs);
	}
	return sc;
}

double LooperSolver::calcScoreStructureActiveRegionLengths() {
	double sc = 0.0;
	double diff, dtn;
	vector3 v;
	for (size_t i = 0; i+1 < active_region.size(); i++) {
		vector3 v = clusters[active_region[i]].pos - clusters[active_region[i+1]].pos;
		dtn = clusters[active_region[i]].dist_to_next;
		if (dtn < 1e-6) dtn = 1e-6;
		diff = (v.length() - dtn) / dtn;
		sc += diff * diff * (diff >= 0.0f ? Settings::springConstantStretch : Settings::springConstantSqueeze);
	}
	return sc;
}

// calculate score based on distances to previous and next cluster
double LooperSolver::calcScoreStructureActiveRegionLengths(int cluster_moved) {
	double sc = 0.0;
	double diff, dtn;
	vector3 v;
	for (int i = cluster_moved-1; i <= cluster_moved; i++) {
		if (i<0 || i+1>=(int)active_region.size()) continue;

		v = clusters[active_region[i]].pos - clusters[active_region[i+1]].pos;
		dtn = clusters[active_region[i]].dist_to_next;
		if (dtn < 1e-5) dtn = 1e-5;
		diff = (v.length() - dtn) / dtn;
		sc += diff * diff * (diff >= 0.0f ? Settings::springConstantStretch : Settings::springConstantSqueeze);
	}
	return sc;
}

double LooperSolver::calcScoreStructureSmooth(bool lengths, bool angles) {
	double diff, ang, dtn;
	vector3 v, v2;

	double sca = 0.0, scb = 0.0;

	for (size_t i = 0; i+1 < active_region.size(); i++) {

		v = clusters[active_region[i]].pos - clusters[active_region[i+1]].pos;

		if (lengths) {
			dtn = clusters[active_region[i]].dist_to_next;
			if (dtn < 1e-6) dtn = 1e-6;
			diff = (v.length() - dtn) / dtn;
			sca += diff * diff * (diff >= 0.0f ? Settings::springConstantStretch : Settings::springConstantSqueeze);
		}

		if (angles && i>0) {
			ang = angle(v, v2);
			scb += ang * ang * ang *  Settings::springAngularConstant;

			//printf("sm %d %lf %f\n", i, ang * ang * Settings::springAngularConstant, ang);
		}

		//printf("diff = %lf, ang=%lf   %f %f %f  (%f %f %f)\n", diff, ang, v.x, v.y, v.z,
		//		clusters[active_region[i]].pos.x, clusters[active_region[i]].pos.y, clusters[active_region[i]].pos.z);
		v2 = v; 	// now, in the next loop, we will have two vectors and so we can calculate angles
	}

	//printf("%lf %lf\n", sca, scb);
	return sca * Settings::weightDistSmooth + scb * Settings::weightAngleSmooth;
}

double LooperSolver::calcScoreStructureSmooth(int cluster_moved, bool lengths, bool angles) {
	double diff, ang, dtn;
	vector3 v, v2;

	double sca = 0.0, scb = 0.0;

	for (int i = cluster_moved-2; i < cluster_moved+2; i++) {

		if (i < 0 || i+1 >= (int)active_region.size()) continue;

		v = clusters[active_region[i]].pos - clusters[active_region[i+1]].pos;

		// distance change only in two vectors - direct neighbours
		if (lengths && (i==cluster_moved-1 || i==cluster_moved)) {
			dtn = clusters[active_region[i]].dist_to_next;
			if (dtn < 1e-6) dtn = 1e-6;
			diff = (v.length() - dtn) / dtn;
			sca += diff * diff * (diff >= 0.0f ? Settings::springConstantStretch : Settings::springConstantSqueeze);
		}

		// three angles change
		if (angles && i > cluster_moved-2 && i > 0) {
			ang = angle(v, v2);
			scb += ang * ang * ang * Settings::springAngularConstant;
			//printf("sm p=%d %d %lf %f\n", cluster_moved, i, ang * ang * Settings::springAngularConstant, ang);
		}

		//printf("diff = %lf, ang=%lf   %f %f %f  (%f %f %f)\n", diff, ang, v.x, v.y, v.z,
		//		clusters[active_region[i]].pos.x, clusters[active_region[i]].pos.y, clusters[active_region[i]].pos.z);
		v2 = v; 	// now, in the next loop, we will have two vectors and so we can calculate angles
	}

	//printf("%lf %lf\n", sca, scb);
	return sca + scb;
}

double LooperSolver::calcScoreOrientation(const vector<vector3> &orientation) {
	double ang, err = 0.0;
	for (auto el: active_anchors_neighbors) {
		int n = el.second.size();
		for (int i = 0; i < n; ++i) {
			ang = angle_norm(orientation[el.first], (Settings::motifsSymmetric ? 1.0 : -1.0) * orientation[el.second[i]]);
			err += ang * ang * active_anchors_neighbors_weight[el.first][i];
		}

		//for (int end: el.second) {
		//	ang = angle_norm(orientation[el.first], orientation[end]);
		//	err += ang * ang;
		//}
	}

	// old version: pairwise orientation
	//	for (int i = 0; i+1 < orientation.size(); ++i) {
	//		ang = angle_norm(orientation[i], orientation[i+1]);
	//		err += ang * ang;
	//	}
	return err * Settings::motifOrientationWeigth;
}

double LooperSolver::calcScoreOrientation(const vector<vector3> &orientation, int anchor_index) {
	double err = 0.0, ang;
	if (active_anchors_neighbors.find(anchor_index) == active_anchors_neighbors.end()) {
		printm(active_anchors_neighbors, "active_anchors_neighbors");
		printf("anchor ind = %d\n", anchor_index);
		error("key not found");
	}

	// for each anchor's neighbor calc orientation score
	// consider whether the symmetric or assymetric paradigm is chosen (by possibly reversing motif orientation)
	for (int val: active_anchors_neighbors[anchor_index]) {
		ang = angle_norm(orientation[anchor_index], (Settings::motifsSymmetric ? 1.0 : -1.0) * orientation[val]);
		err += ang * ang;
	}


	//	if (anchor_index > 0) {
	//		ang = angle_norm(orientation[anchor_index], orientation[anchor_index-1]);
	//		err += ang * ang;
	//	}
	//	if (anchor_index+1 < orientation.size()) {
	//		ang = angle_norm(orientation[anchor_index], orientation[anchor_index+1]);
	//		err += ang * ang;
	//	}
	return err * Settings::motifOrientationWeigth;
}


// orientation may influence only direct neighbors
// the procedure is as follows:
// - check if clustered moved actually modified orientation (to do so, it must be an anchor or a direct neighbor)
// - find neighboring anchors (note: there may be only one)
// - calculate orientations for all of them
// - add score based on the angles between them
//double LooperSolver::calcScoreOrientation(int cluster_moved) {
//
//	// during subanchor level simulation anchors have 'is_fixed' flag set. we can use it to identify anchors
//
//	// first check if orientation changed after moving of 'cluster_moved'
//	bool is_changed = false;
//	for (int i=cluster_moved-1; i<=cluster_moved+1; i++) {
//		if (i<0 || i>= active_region.size()) continue;
//		if (clusters[active_region[cluster_moved]].is_fixed) {
//			is_changed = true;
//			break;
//		}
//	}
//
//	if (!is_changed) return 0.0;   // there was no change
//
//	return 0.0;
//}


// calc score for active_region based on heatmap_dist
// moved - index of a point last moved. if none is specified returns sum of score for all points
double LooperSolver::calcScoreHeatmapActiveRegion(int moved) {
	double err = 0.0, cerr;
	float d;

	size_t n = active_region.size();
	if (heatmap_dist.size != n) error(ftext("heatmap sizes mismatch, dist size=%d, active region=%d", heatmap_dist.size, n));

	if (moved == -1) {
		for (size_t i = 0; i < n; ++i) err += calcScoreHeatmapActiveRegion(i);
	}
	else {
		int st = 0;
		int end = n- 1;
		if (heatmap_chromosome_boundaries.size() > 0) getChromosomeHeatmapBoundary(moved, st, end);

		for (int i = st; i <= end; ++i) {
			if (abs(i-moved) >= heatmap_dist.diagonal_size) {

				if (heatmap_dist.v[i][moved] < 1e-6) continue;	// ignore 0 values

				d = (clusters[active_region[i]].pos - clusters[active_region[moved]].pos).length();
				cerr = (d - heatmap_dist.v[i][moved]) / heatmap_dist.v[i][moved];
				err += cerr * cerr;
				//printf("%d d=%f, exp_d=%f, %lf %lf\n", i, d, exp_len, cerr, err);
			}
		}
	}
	return err;
}

void LooperSolver::updateDensity() {
	// clear density_curr
	// for all points:
	//    add point masses, distribute influence
	// compare density and density_curr
	//    normalize them both

	density_curr.clear();

	//showCoordinatesRange(true);

	size_t n = active_region.size();
	vector3 pos;
	vector3 shift = density.center - density.origin;
	//vector3 origin = density.center - vector3(density. * Settings::densityScale;
	int px, py, pz;
	for (size_t i = 0; i < n; ++i) {

		// find 3d density-coordinates for the current bead
		//pos = (clusters[active_region[i]].pos - density.center) * Settings::densityScale + shift;
		pos = (clusters[active_region[i]].pos - density.center) * Settings::densityScale;
		pos.x /= 3.0f;
		pos += shift;
		px = (int)pos.x;
		py = (int)pos.y;
		pz = (int)pos.z;

		if (px < 0 || py < 0 || pz < 0 || px >= density.size_x || py >= density.size_y || pz >= density.size_z) {
			//printf("\npoint outside: %d %d %d\n", px, py, pz);
			//print_vector(clusters[active_region[i]].pos, "org pos");
			//print_vector(pos, "dens pos");
		}
		else {
			density_curr.addPointMass(px, py, pz, 1.0f);
		}
	}

	density_curr.normalize();
}

void LooperSolver::repairDensityBoundary() {

	density_curr.clear();

	size_t n = active_region.size();
	vector3 pos;
	vector3 shift = density.center - density.origin;
	int px, py, pz;
	for (size_t i = 0; i < n; ++i) {

		// find 3d density-coordinates for the current bead
		//pos = (clusters[active_region[i]].pos - density.center) * Settings::densityScale + shift;

		while (1) {

			pos = (clusters[active_region[i]].pos - density.center) * Settings::densityScale;
			pos.x /= 3.0f;
			pos += shift;
			px = (int)pos.x;
			py = (int)pos.y;
			pz = (int)pos.z;

			if (px < 0 || py < 0 || pz < 0 || px >= density.size_x || py >= density.size_y || pz >= density.size_z || density.t[px][py][pz] < epsilon) {
				//vector3 shift = density.center - clusters[active_region[i]].pos;

				vector3 cen = density.center + random_vector(20.0f);

				vector3 shift = cen - clusters[active_region[i]].pos;
				shift = shift.normalize();
				shift *= 2.0f;

				//				printf("\n\ni=%d\n", i);
				//				print_vector(density.center, "center");
				//				print_vector(clusters[active_region[i]].pos, "pos");
				//				print_vector(shift, "shift");
				//				printf("%d %d %d\n", px, py, pz);

				clusters[active_region[i]].pos += shift;
			}
			else break;
		}

		//if (found) exit(0);
	}

	density_curr.normalize();
}

double LooperSolver::calcScoreDensity() {

	double error = 0.0, e;

	density_curr.clear();

	//showCoordinatesRange(true);

	size_t n = active_region.size();
	vector3 pos;
	vector3 shift = density.center - density.origin;
	//vector3 origin = density.center - vector3(density. * Settings::densityScale;
	int px, py, pz;
	for (size_t i = 0; i < n; ++i) {

		// find 3d density-coordinates for the current bead
		//pos = (clusters[active_region[i]].pos - density.center) * Settings::densityScale + shift;
		pos = (clusters[active_region[i]].pos - density.center) * Settings::densityScale;
		pos.x /= 3.0f;
		pos += shift;
		px = (int)pos.x;
		py = (int)pos.y;
		pz = (int)pos.z;

		if (px < 0 || py < 0 || pz < 0 || px >= density.size_x || py >= density.size_y || pz >= density.size_z) {
			// find distance (use manhattan metric) to closest "proper" cell
			int d = 0;
			if (px < 0) d += -px;
			else if (px >= density.size_x) d += px - density.size_x + 1;
			if (py < 0) d += -py;
			else if (py >= density.size_y) d += py - density.size_y + 1;
			if (pz < 0) d += -pz;
			else if (pz >= density.size_z) d += pz - density.size_z + 1;

			error += d * 1000000;

			//printf("\npoint outside: %d %d %d\n", px, py, pz);
			//print_vector(clusters[active_region[i]].pos, "org pos");
			//print_vector(pos, "dens pos");
		}
		else if (density.t[px][py][pz] < epsilon) {
			//error += (clusters[active_region[i]].pos - density.center).length() * 1000.0f;
			//error += 100000;
			//error += 100.0f * abs(min(py-15, pz-15));
			error += 1000000.0f +  1000000.0f / (pos - density.center).length();
		}
		else {
			density_curr.addPointMass(px, py, pz, 1.0f);
		}
	}

	density_curr.normalize();

	for (int i = 0; i < density.size_x; ++i) {
		for (int j = 0; j < density.size_y; ++j) {
			for (int k = 0; k < density.size_z; ++k) {

				//e = (density.t[i][j][k] > epsilon) ? (density.t[i][j][k] - density_curr.t[i][j][k]) / density.t[i][j][k] : 1000.0f * density_curr.t[i][j][k];

				if (density.t[i][j][k] > epsilon) {
					e = (density.t[i][j][k] - density_curr.t[i][j][k]) / density.t[i][j][k];
					e = 0.0;
				}
				else {
					e = 5.0f * density_curr.t[i][j][k];
				}

				//e = (density.t[i][j][k] > epsilon) ? (density.t[i][j][k] - density_curr.t[i][j][k]) / density.t[i][j][k] : 1000.0f;
				//if (density.t[i][j][k] < epsilon) printf("near zero: %lf  %lf\n", density.t[i][j][k], density_curr.t[i][j][k]);
				error += e * e;
			}
		}
	}

	//density_curr.toFile(ftext("%sdensity_norm.txt", output_dir.c_str()));
	return error * Settings::densityWeight;
}

double LooperSolver::calcScoreSubanchorHeatmap(int moved) {
	double err = 0.0, cerr;
	size_t n = active_region.size();
	assert(heatmap_dist_subanchor.size == n);

	if (moved == -1) {
		for (size_t i = 0; i < n; ++i) {
			err += calcScoreSubanchorHeatmap(i);
			//printf("%lf, ", calcScoreSubanchorHeatmap(i));
		}
		//printf("\n");
	}
	else {
		double d;
		for (size_t i = 0; i < n; ++i) {
			if ((int)i==moved) continue;

			if (heatmap_dist_subanchor.v[i][moved] < 1e-6) continue;	// ignore 0 values

			d = (clusters[active_region[i]].pos - clusters[active_region[moved]].pos).length();
			cerr = (d - heatmap_dist_subanchor.v[i][moved]) / heatmap_dist_subanchor.v[i][moved];
			err += cerr * cerr;
			//printf("%d d=%f, exp_d=%f, %lf %lf\n", i, d, heatmap_dist_subanchor.v[i][moved], cerr, err);
		}
		return err * Settings::subanchorHeatmapDistWeight;
	}
	return err;
}

void LooperSolver::densifyActiveRegion(int cluster, bool fix) {
	// TODO: number of points based on length?

	//printActiveRegion();
	if (active_region.size() == 0) return;	// nothing to do

	int add = Settings::loopDensity;

	std::vector<int> tmp_active;
	std::vector<int> tmp_children;

	int level = clusters[active_region[0]].level + 1;	// level for newly added points

	for (size_t i = 0; i+1 < active_region.size(); ++i) {

		if (fix) clusters[active_region[i]].is_fixed = true;

		// get rightmost pos of left region and leftmost of right one
		int range = clusters[active_region[i+1]].start - clusters[active_region[i]].end;
		int d = range / (add+1);
		int p = clusters[active_region[i]].end;

		// update left anchor dist (because we insert some points that will be closer)
		//regions[active_region[i]].dist_to_next = d_genomic;

		tmp_active.push_back(active_region[i]);

		// add new points
		//regions[current_regions[0]].print();
		float dst = 1.0f / (add + 1);
		float st = dst;

		for (int j = 0; j < add; ++j) {
			p += d;

			Cluster c(p, p);
			c.pos = interpolate(clusters[active_region[i]].pos, clusters[active_region[i+1]].pos, st);
			c.level = level;

			clusters.push_back(c);

			tmp_active.push_back(clusters.size()-1);

			st += dst;
		}
	}

	tmp_active.push_back(active_region[active_region.size()-1]);
	if (fix) clusters[active_region[active_region.size()-1]].is_fixed = true;

	active_region = tmp_active;

	// we added new points thus modifying the structure tree, so we also need to update the parent info
	clusters[cluster].children = active_region;
}


double LooperSolver::genomicLengthToDistance(int length) {
	return Settings::genomicLengthToDistBase + Settings::genomicLengthToDistScale * pow(length/1000.0, Settings::genomicLengthToDistPower);
}

int LooperSolver::otherEnd(const InteractionArc &arc, int current_cluster) {
	if (arc.start == current_cluster) return arc.end;
	if (arc.end == current_cluster) return arc.start;
	printf("no other end found (%d)\n", current_cluster);
	exit(0);
	return -1;
}

// returns index of cluster connected to 'current_cluster' by arc with index 'arc'
// Note that 'current_chr' must be properly set
int LooperSolver::otherEnd(int arc, int current_cluster) {
	if (arcs.arcs[current_chr][arc].start == current_cluster) return arcs.arcs[current_chr][arc].end;
	if (arcs.arcs[current_chr][arc].end == current_cluster) return arcs.arcs[current_chr][arc].start;
	return -1;
}

int LooperSolver::otherEnd(string chr, int arc, int current_cluster) {
	if (arcs.arcs[chr][arc].start == current_cluster) return arcs.arcs[chr][arc].end;
	if (arcs.arcs[chr][arc].end == current_cluster) return arcs.arcs[chr][arc].start;
	return -1;
}

float LooperSolver::freqToDistanceHeatmap(float freq) {
	return Settings::freqToDistHeatmapScale * pow(freq, Settings::freqToDistHeatmapPower);
}

float LooperSolver::freqToDistanceHeatmapInter(float freq) {
	return Settings::freqToDistHeatmapScaleInter * pow(freq, Settings::freqToDistHeatmapPowerInter);
}

float LooperSolver::distToFreqHeatmap(float freq) {
	return pow(freq / Settings::freqToDistHeatmapScale, 1.0f / Settings::freqToDistHeatmapPower);
}

float LooperSolver::freqToDistance(int freq, bool memo) {
	if (memo && freq >= 1 && freq <= 100) return freq_to_distance[freq];
	return Settings::countToDistBaseLevel + Settings::countToDistScale / exp(Settings::countToDistA * (double)(freq+Settings::countToDistShift));
}

// reconstructs model on anchor level
// segments are reconstructed independently from each other, one by one
// for every segment:
// - get all anchors and set them as 'active_region'
// - reconstruct anchor-level model (ie. using a single point for every anchor) (to get their proper 3D position)
// - densify model (put additional points between every two consecutive anchors)
// - refine subanchor-level model, keeping anchors fixed and moving only these additional points (to obtain nice loops)
void LooperSolver::reconstructClustersArcsDistances() {

	setLevel(LVL_SEGMENT);	// set segment level

	int curr_size = 0;
	int ii = 0;
	for (string chr: chrs) {
		ii++;

		this->current_chr = chr;	// save currently processed chromosome (used in MC)

		if (debug_chromosomes_limit > 0 && ii > debug_chromosomes_limit) break;  // check chr number limit

		setLevel(LVL_SEGMENT);	// set segment level


		if (Settings::outputLevel >= 6) printv(current_level[chr], true, true, "curr seg");

		output(5, "position interaction blocks...\n");
		positionInteractionBlocks(current_level[chr]);

		setLevel(LVL_INTERACTION_BLOCK);	// set interaction blocks level
		if (Settings::outputLevel >= 6) printv(current_level[chr], true, true, "curr ib");

		// now go through all ibs
		curr_size = current_level[chr].size();
		for (int i = 0; i < curr_size; ++i) {

			output(1, " %s %d/%d\n", chr.c_str(), i+1, curr_size);

			// update active region
			active_region.assign(clusters[current_level[chr][i]].children.begin(), clusters[current_level[chr][i]].children.end());
			if (Settings::outputLevel >= 6)	printv(active_region, true, true, "active");

			// create anchor and subanchor heatmaps
			if (Settings::useSubanchorHeatmap || Settings::useAnchorHeatmap) {
				vector<int> anchors_gap_len;
				createSingletonSubanchorHeatmap(chr, anchors_gap_len);
			}

			//interpolatePosition(active_region, current_level[chr], i, true);
			for (size_t j = 0; j < active_region.size(); ++j) clusters[active_region[j]].pos = clusters[current_level[chr][i]].pos;

			// calc expected distances between anchors
			output(5, "calc exp distances for anchors\n");
			calcAnchorExpectedDistancesHeatmap();
			//heatmap_exp_dist_anchor.toFile(ftext("%sanchor_%s_%d_exp_dist.heat", output_dir.c_str(), chr.c_str(), i));

			// reconstruct anchor level model
			output(5, "reconstruct anchors\n");
			reconstructClusterArcsDistances(current_level[chr][i], i, false);

			//Heatmap hd = calcTrueDistancesHeatmapForRegion(active_region);
			//hd.toFile(ftext("%sanchor_%s_%d_true_dist.heat", output_dir.c_str(), chr.c_str(), i));

			//printActiveRegion();

			// now add new points (and interpolate them along the anchors)
			output(5, "densify\n");
			densifyActiveRegion(current_level[chr][i], true);

			// update info about anchors' neighbors
			if (Settings::useCTCFMotifOrientation) {
				output(5, "calcActiveAnchorsNeighbors\n");
				calcActiveAnchorsNeighbors();
			}

			// find average distance between sub-anchor beads
			if (Settings::useSubanchorHeatmap) {
				output(5, "get avg distance for subanchor beads\n");
				Heatmap subanchor_avg_dist;
				subanchor_avg_dist.init(active_region.size());

				// save current structure
				vector<vector3> str;
				for (size_t j = 0; j < active_region.size(); ++j) str.push_back(clusters[active_region[j]].pos);

				for (int j = 0; j < Settings::subanchorEstimateDistancesReplicates; ++j) {
					reconstructClusterArcsDistances(current_level[chr][i], i, true, false);

					Heatmap hd = calcTrueDistancesHeatmapForRegion(active_region);
					subanchor_avg_dist.add(hd);
					// restore previous structure
					for (size_t j = 0; j < active_region.size(); ++j) clusters[active_region[j]].pos = str[j];
				}
				subanchor_avg_dist.scale(1.0f / Settings::subanchorEstimateDistancesReplicates);
				//subanchor_avg_dist.toFile(ftext("%sdist/subanchor_dist_%s_%d_avg.heat", output_dir.c_str(), chr.c_str(), i));

				// create expected dist heatmap
				output(5, "create expected dist heatmap\n");
				createExpectedDistSubanchorHeatmap(subanchor_avg_dist);
				//heatmap_dist_subanchor.toFile(ftext("%sdist/subanchor_dist_%s_%d_exp.heat", output_dir.c_str(), chr.c_str(), i));
			}

			// reconstruct sub-anchor level model
			output(4, "reconstruct sub-anchors\n");
			reconstructClusterArcsDistances(current_level[chr][i], i, true, Settings::useSubanchorHeatmap);

			//			if (Settings::useSubanchorHeatmap) {
			//				Heatmap hd = calcTrueDistancesHeatmapForRegion(active_region);
			//				hd.toFile(ftext("%ssubanchor_dist_%s_%d_true.heat", output_dir.c_str(), chr.c_str(), i));
			//			}

			//exit(0);
			//getSnapshot(ftext("cluster_%d,_smoothed", i));
		}
	}
}

// position interaction blocks along the segments
// be default use spline interpolation to find correct position
// if there is only one segment we can use random walk, as we have no other reference
void LooperSolver::positionInteractionBlocks(std::vector<int> &segments) {

	if (segments.size() > 1) interpolateChildrenPositionSpline(segments, true);
	else {
		// TODO: create heatmap for ib?
		printf("only one segment, use random walk to position ib (step size=%f)\n", Settings::ibRandomWalkJumps);
		int n = clusters[segments[0]].children.size();
		rw_pos.set(0.0f, 0.0f, 0.0f);

		for (int j = 0; j < n; ++j) {
			rw_pos = displace(rw_pos, Settings::ibRandomWalkJumps, Settings::use2D);
			clusters[clusters[segments[0]].children[j]].pos = rw_pos;
		}
	}
}

// reconstruct the model on (sub)anchor level
// cluster - id of parent
// cluster_ind - parent index of current_level array (needed to get previous/next clusters for interpolation)
// smooth - if false then use anchor level parameters (ie. use proper energy function), otherwise sub-anchor level
// active_region should contain indices of all
// Note: function assumes that 'current_chr' stores the id of chromosome that 'cluster' belongs to
void LooperSolver::reconstructClusterArcsDistances(int cluster, int cluster_ind, bool smooth, bool use_subanchor_heatmap) {
	//printf("reconstruct cluster %d (current_region ind: %d, smooth=%d)\n", cluster, cluster_ind, smooth);
	//clusters[cluster].print();

	if (active_region.size() <= 1) return;  // nothing to position

	size_t active_size = active_region.size();

	if (active_size == 0) {
		error(ftext("active size == 0!\ncluster=%d, cluster ind=%d, smooth=%d", cluster, cluster_ind, smooth));
	}

	// check if everything is set up
	if (smooth) {
		if (use_subanchor_heatmap && heatmap_dist_subanchor.size != active_size)
			error(ftext("heatmap_dist_subanchor size mismatch (%d, should be %d)", heatmap_dist_subanchor.size, active_size));
	}
	else {
		if (heatmap_exp_dist_anchor.size != active_size)
			error(ftext("heatmap_exp_dist_anchor size mismatch (%d, should be %d)", heatmap_exp_dist_anchor.size, active_size));
	}

	// calc distances (for every point set proper distance to the next one)
	float dist;
	float noise_size = 0.0f;	// noise for subanchors
	float noise_size_small = 0.05f;	// noise for anchors

	for (size_t i = 0; i+1 < active_size; ++i) {
		int d = abs(clusters[active_region[i+1]].genomic_pos - clusters[active_region[i]].genomic_pos);
		dist = genomicLengthToDistance(d);
		clusters[active_region[i]].dist_to_next = dist;
		noise_size += dist;
		//printf("%d d=%d ; %d %d ; %f\n", i, d, clusters[active_region[i]].genomic_pos, clusters[active_region[i+1]].genomic_pos, dist);
	}

	// we will use distance between neighbors as a measure of noise - it will be used eg. to modify
	// magnitude of points displacement during MC simulation
	noise_size /= active_region.size();
	noise_size *= smooth ? Settings::noiseCoefficientLevelSubanchor : Settings::noiseCoefficientLevelAnchor;
	output(5, "noise size = %f   %d %d\n", noise_size, active_region.size(), active_size);

	// * initial structure
	// structure for a segment will typically consist of several tight loop clusters. Number of these clusters
	// depend on how anchors are connected by loops. We use different initial structures for anchor and subanchor level,
	// based on a fact that loops clusters anchors are close to each other (r), but distances between them may be large (R).
	// Loop size is large as well (R), so that loops span whole TD
	// here R is approximate size of TD (corresponding to current segment) and r distance between interacting anchors
	// - on anchor level:
	//	  - identify loop clusters (maximal subsets of anchors that are connected)
	//    - randomly position loop clusters 'origins' in distance R around the parent (segment) cluster position
	//    - for each loop cluster randomly position corresponding anchors in distance r from the loop cluster origin
	// - on subanchor level
	//    - (previously we interpolated subanchors between anchors beads, but this may not be the best idea,
	//      as now distances between anchors should be much smaller)
	//    - randomly position beads in distance R from ... something (eg. midpoint between neighboring anchor beads?)

	// TODO: how to best get values for R and r?
	// TODO: check if there are actually multiple loops clusters in segments (there are many arcs, so it may always be just one)

	output(5, "initial\n");
	std::vector<vector3> initial_structure;	// store initial position for all runs
	for (size_t i = 0; i < active_size; ++i) initial_structure.push_back(clusters[active_region[i]].pos);
	//getSnapshot(ftext("base_initial_structure_arcs%s", smooth ? "_smooth" : ""));

	double score, best_score = -1.0;
	std::vector<vector3> best_structure;

	// number of simulations for every segment
	int steps = Settings::simulationStepsLevelAnchor;
	if (smooth) {
		if (Settings::useSubanchorHeatmap && !use_subanchor_heatmap) steps = Settings::subanchorEstimateDistancesSteps;
		else steps = Settings::simulationStepsLevelSubanchor;
	}

	//smooth ? Settings::simulationStepsLevelArcsSmooth : Settings::simulationStepsLevelArcs;
	//if (!smooth && Settings::useSubanchorHeatmap && !use_subanchor_heatmap) steps /= 5;

	for (int k = 0; k < steps; ++k) {

		output(3, "MC segment id=%d, step %d/%d (smooth=%d)\n", cluster, k+1, steps, smooth);

		// TODO: generate it for all threads
		// set random initial position (use initial_structure and add some noise, but only to subanchors (ie. not-anchors))
		for (size_t i = 0; i < active_size; ++i) {
			clusters[active_region[i]].pos = initial_structure[i];
			if (!smooth || !clusters[active_region[i]].is_fixed) clusters[active_region[i]].pos +=
					random_vector(smooth ? noise_size : noise_size_small, Settings::use2D);
		}

		// run the simulation (and keep track of the current score)
		//if (smooth) {
		//	score = 1.0f;
		//}
		//else

		score = smooth ? MonteCarloArcsSmooth(noise_size, use_subanchor_heatmap) : MonteCarloArcs(noise_size);

		output(3, "score = %lf, best = %lf\n", score, best_score);

		// if current structure is better than best so far than save current structure as best
		if (score < best_score || best_score < 0.0) {
			best_structure.clear();
			for (size_t i = 0; i < active_size; ++i) best_structure.push_back(clusters[active_region[i]].pos);
			best_score = score;
		}
	}

	// restore best structure
	for (size_t i = 0; i < active_size; ++i) clusters[active_region[i]].pos = best_structure[i];
}

void LooperSolver::interpolateChildrenPosition(std::vector<int> &regions) {
	for (size_t i = 0; i < regions.size(); ++i) {
		interpolatePosition(clusters[regions[i]].children, regions, i);
	}
}

// interpolate points corresponding to 'regions'.
// we use this to calculate initial position for beads of subanchor level (ie. we interpolate them between anchor-level beads)
void LooperSolver::interpolatePosition(std::vector<int> &regions, std::vector<int> &clusters_template, int region_ind, bool noise) {
	int p1 = region_ind-1;
	int p2 = region_ind;
	size_t n = regions.size();
	float dst = 1.0f / n;
	float st = 0.5f + dst / 2.0;

	if (region_ind == 0) {
		p1 = region_ind;
		p2 = region_ind + 1;
		st = -0.5f + dst / 2.0f;
	}

	// calc noise to add - use fraction of avg distance between consecutive points.
	float avg_disp = 0.0f;
	if (noise) {
		if (clusters_template.size() > 1) {
			for (size_t i = 1; i < clusters_template.size(); ++i) avg_disp += (clusters[clusters_template[i]].pos-clusters[clusters_template[i-1]].pos).length();
			avg_disp /= (clusters_template.size() - 1);	// take average...
			avg_disp *= 0.25f;							// ... and the fraction of it
		}
		else avg_disp = 0.0f;	// if there is only one point do not add any noise (TODO: a global noise magnitude)
		output(6, "noise disp = %f\n", avg_disp);
	}

	for (size_t i = 0; i < n; ++i) {
		//printf("%d (%d %d) %f, set %d\n", i, p1, p2, st, active_region[i]);
		if (!clusters[regions[i]].is_fixed) {
			clusters[regions[i]].pos = interpolate(clusters[clusters_template[p1]].pos, clusters[clusters_template[p2]].pos, st);
			if (noise) clusters[regions[i]].pos += random_vector(avg_disp, Settings::use2D);
		}

		st += dst;

		// switch to the second interval
		if (st >= 1.0f && region_ind+1<(int)clusters_template.size()) {
			p1 = region_ind;
			p2 = region_ind + 1;
			st -= 1.0f;
		}
	}
}

void LooperSolver::interpolateChildrenPositionSpline(std::vector<int> &regions, bool use_genomic_dist) {
	size_t reg_cnt = regions.size();
	vector3 pos;

	// init points (for first region, we take P[-2]..P[1])
	vector3 pts[4];
	pts[1] = mirrorPoint(clusters[regions[0]].pos, clusters[regions[1]].pos);
	pts[0] = mirrorPoint(pts[1], clusters[regions[0]].pos);
	for (size_t i = 2; i < 4; ++i) pts[i] = clusters[regions[i-2]].pos;

	// calc two last points
	vector3 end_pt = mirrorPoint(clusters[regions[reg_cnt-1]].pos, clusters[regions[reg_cnt-2]].pos);
	vector3 end_pt2 = mirrorPoint(end_pt, clusters[regions[reg_cnt-1]].pos);
	vector<float> knots;
	float st, dst;

	// for every parent cluster
	for (size_t i = 0; i < reg_cnt; ++i) {

		size_t n = clusters[regions[i]].children.size();

		int switch_controls = -1;	// index of child after which we should switch control points by one

		// calculate spline parameters (knots)
		// by default use equidistant nodes
		knots.clear();
		if (use_genomic_dist) {

			// * first (old) version:
			// find start and end location of region (segment + half of flanking between-segments regions)
			// for every children get position relative to the region and transform it to the knot
			// the problem is that if right flanking region is much bigger than it will bias all the regions to the left (beacuse the
			// relative positions will get smaller)

			// * second version:
			// as before, but find positions relative to the corresponding halves - if the children is on the left of the segment center
			// than find position relative to the left flanking region and left half of the segment. This way if the right flanking region
			// is geting bigger and bigger it won't influence the anchors on the left side

			// get center, left and right boundary genomic location
			int center_loc = clusters[regions[i]].genomic_pos;
			int start_loc = clusters[regions[i]].start;
			int end_loc = clusters[regions[i]].end;

			// we need also to consider flanking between-segments space (we use center of this space)
			int left_flank = i==0 ? start_loc : clusters[regions[i-1]].end;

			size_t right_flank = i==regions.size()-1 ? end_loc : clusters[regions[i+1]].start;

			// calculate left and right position including flanking regions
			start_loc = (left_flank + start_loc) / 2;
			end_loc = (right_flank + end_loc) / 2;

			// calculate lengths of left and right subregion
			float left_length = (float)(center_loc - start_loc);
			float right_length = (float)(end_loc - center_loc);

			// for every child calculate its location relative to (start_loc, end_loc) and transform it to a knot [0.5->0.9->1.0->0.1->0.5]
			for (size_t j = 0; j < n; ++j) {
				int pos = clusters[clusters[regions[i]].children[j]].genomic_pos;
				float p = 1.0f;
				if (pos < center_loc) {
					p = (pos - start_loc) / left_length;  // location proportion [0, 1]
					p = 0.5f + p * 0.5;	// transform to [0.5, 1]
				}
				else {
					p = (pos - center_loc) / right_length;
					p = p * 0.5;	// transform to [0, 0.5]
					if (switch_controls == -1) switch_controls = j;	// this was the first child to cross 1.0
				}
				knots.push_back(p);
			}
		}
		else {
			dst = 1.0f / n;
			st = 0.5f + dst / 2.0;
			for (size_t j = 0; j < n; ++j) {
				knots.push_back(st);
				st += dst;
				if (st > 1.0f) {
					if (switch_controls == -1) switch_controls = (int)j;
					st -= 1.0f;
				}
			}
		}

		//printv(knots, true, true, "knots");

		if (knots.size() != n) error("knots number != n");

		// iterate over childs
		for (size_t j = 0; j < n; ++j) {
			if ((int)j == switch_controls || (switch_controls==-1 && j==n-1)) {
				for (int k = 0; k < 3; ++k) pts[k] = pts[k+1];
				if (i+2 == regions.size()) pts[3] = end_pt;
				else if (i+2 == regions.size()+1) pts[3] = end_pt2;
				else pts[3] = clusters[regions[i+2]].pos;
			}
			pos = interpolateSpline(knots[j], pts[0], pts[1], pts[2], pts[3]);
			clusters[clusters[regions[i]].children[j]].pos = pos;
		}
	}
}

double LooperSolver::MonteCarloArcs(float step_size) {

	std::cout << "Arcs started\n";

	double maxT = Settings::maxTemp;	  		// max temperature
	double dt = Settings::dtTemp;				// change of temperature
	double T = maxT;							// temperature
	double tp = 0.0;							// transition probability

	int i = 1, p, ind;
	vector3 tmp;

	bool ok;
	double score_prev, score_curr;		// global score (used to check if we can stop)
	double local_score_prev, local_score_curr;	// local score
	double milestone_score;		// score measured every N steps, used to see whether the solution significantly improves
	int success = 0;			// overall number of successes
	int milestone_success = 0;	// calc how many successes there were since last milestone
	int size = active_region.size();

	// the idea to keep the true score without recalculating it every time is to calculate total score once, on the beginning.
	// during the simulation we calculate local score before and after a move
	// new score is then equal to (old total score) - (local score before move) + (local score after move)

	// calculate total score
	score_curr = calcScoreDistancesActiveRegion();

	score_prev = score_curr;
	milestone_score = score_curr;

	output(4, "size = %d; initial score = %lf, step size = %f\n", size, score_curr, step_size);

	while (true) {

		p = random(size);		// select point to mutate
		ind = active_region[p];	// translate the index to 'clusters'

		if (clusters[ind].is_fixed) error("cluster fixed during arcs!\n");

		local_score_prev = calcScoreDistancesActiveRegion(p);

		tmp = random_vector(step_size, Settings::use2D);
		clusters[ind].pos += tmp;

		local_score_curr = calcScoreDistancesActiveRegion(p);

		score_curr = score_curr - local_score_prev + local_score_curr;

		ok = score_curr <= score_prev;

		if (!ok) {
			tp = Settings::tempJumpScale * exp(-Settings::tempJumpCoef * (score_curr/score_prev) / T);
			ok = withChance(tp);
			//printf("sc = %f vs. %f (%f), T=%f, tp=%f, ok=%d\n", score_curr, score_prev, score_curr/score_prev, T, tp, ok);
		}


		if (ok) {
			success++;
			milestone_success++;
		}
		else {
			// if we reject move then move back and restore old score
			clusters[ind].pos -= tmp;
			score_curr = score_prev;
		}

		T *= dt;  			// cooling

		// check if we should stop
		if (i % Settings::MCstopConditionSteps == 0) {
			output(7, "milestone: score = %lf, last = %lf (%lf), T=%lf successes = %d, last = %d\n", score_curr, milestone_score,
					score_curr/milestone_score, T, success, milestone_success);

			// stop if the improvement since the last milestone is less than threshold,
			// or if the score is too small (may happen when heatmap is small and points are perfectly arranged)
			if ((score_curr > Settings::MCstopConditionImprovement * milestone_score &&
					milestone_success < Settings::MCstopConditionMinSuccesses) || score_curr < 1e-5 || score_curr/milestone_score>0.9999) break;

			milestone_score = score_curr;
			milestone_success = 0;
		}

		score_prev = score_curr;

		i++;
	}

	return score_curr;
}

double LooperSolver::MonteCarloArcsSmooth(float step_size, bool use_subanchor_heatmap) {

	std::cout << "Arcs smooth started\n";


	double maxT = Settings::maxTempSmooth;	  	// max temperature
	double dt = Settings::dtTempSmooth;			// change of temperature
	double T = maxT;							// temperature
	double tp = 0.0;							// transition probability

	int i, p, ind;
	vector3 tmp;
	bool silent = false;// !use_subanchor_heatmap;
	bool ok;
	double score_prev, score_curr;		// global score (used to check if we can stop)
	double curr_score_structure = 0.0, prev_score_structure = 0.0, local_score_prev_structure = 0.0, local_score_curr_structure = 0.0;
	double curr_score_orientation = 0.0, prev_score_orientation = 0.0, local_score_prev_orientation = 0.0, local_score_curr_orientation = 0.0;
	double curr_score_heat = 0.0, prev_score_heat = 0.0, local_score_prev_heat = 0.0, local_score_curr_heat = 0.0;
	double milestone_score;		// score measured every N steps, used to see whether the solution significantly improves
	int success = 0;			// overall number of successes
	int milestone_success = 0;	// calc how many successes there were since last milestone
	int size = active_region.size();

	// when we move an anchor we may need to recalculate the orientation. then, if we reject the move, we need to restore it
	// prev_orientation will store old orientation, and orn_index denotes the index of anchor for which the orientation has changed
	int orn_index = -1;  // -1 means that there was no change of orientation. >=0 values denotes an 'anchor_orientation' index
	vector3 prev_orientation;
	int anchor_index; // tmp variable

	// mark anchors and their neighbors (used to see whether orientation score changed)
	// for every bead in active_region we assign a coded value
	// coding: 2 - right anchor neighbor, 1 - left neighbor
	// for anchors we use 3+index, ie. 3 for first anchor, 4 for second etc (used to quickly get orientation, this index refers to 'anchor_orientation')
	// negative values: if cluster 'i' is anchor, then cell i-2 keeps negative index of anchor neighboring from left, and i+2 - from right
	// this allows us to immediately see whether particular cluster is anchor, and what are neighbors (required to calculate score)
	// note: for this to work we need at least 4 subanchor clusters between anchors
	// sample result: [3 2 -6 0 0 1 4 2 -12 0 -6 1 5 2 -18 0 -12 1 6 2 -24 0 -18 1 7]
	// so active_region[5] is a left neighbor of the second anchor (with cluster-index: 4-3=1), neighboring anchors are in 0 and 12
	// anchors are denoted by 3,4,5,6,7.
	vector<vector3> anchor_orientation;  // store anchor orientation (we need to update this after anchors or neighbors moving)
	vector<int> cluster_type(size, 0);
	if (Settings::useCTCFMotifOrientation) {
		int last_anchor = -1;
		int anchor_ind = 0;
		for (int i=0; i<size; i++) {
			if (clusters[active_region[i]].is_fixed) {
				cluster_type[i] = 3 + anchor_ind;
				if (i>0) cluster_type[i-1] = 1;
				if (i+1<size) cluster_type[i+1] = 2;

				// save indices of left and right neighbors
				if (last_anchor >= 0) {
					cluster_type[i-2] = -last_anchor;
					cluster_type[last_anchor+2] = -i;
				}

				anchor_orientation.push_back(calcOrientation(i));  // calculate and save orientation

				last_anchor = i;
				anchor_ind++;
			}
		}
		//printv(cluster_type, true, true, "cluster type");
	}

	// init score
	curr_score_structure = calcScoreStructureSmooth(true, true);
	if (Settings::useCTCFMotifOrientation) curr_score_orientation = calcScoreOrientation(anchor_orientation);
	if (use_subanchor_heatmap && Settings::useSubanchorHeatmap) curr_score_heat = calcScoreSubanchorHeatmap();

	score_curr = curr_score_structure + curr_score_orientation + curr_score_heat;

	prev_score_orientation = curr_score_orientation;
	prev_score_structure = curr_score_structure;
	prev_score_heat = curr_score_heat;

	score_prev = score_curr;
	milestone_score = score_curr;

	if (!silent) output(4, "size = %d; initial score = %lf (%lf %lf %lf), step size=%f\n", size, score_curr, curr_score_structure, curr_score_orientation, curr_score_heat, step_size);

	i = 1;
	while (true) {

		p = random(size); // select point to mutate
		ind = active_region[p];

		if (clusters[ind].is_fixed) continue;

		// check if we moved an anchor neighbor, and update orientation if needed
		orn_index = -1;  // -1 means that there was no change of orientation. >=0 values denotes an 'anchor_orientation' index
		if (Settings::useCTCFMotifOrientation) {
			if (cluster_type[p] > 0) {
				// first find the index of anchor in active_region
				anchor_index = p;
				if (cluster_type[p] == 1) anchor_index = p+1;
				else if (cluster_type[p] == 2) anchor_index = p-1;

				// now extract the corresponding index in anchor_orientation
				orn_index = cluster_type[anchor_index] - 3;	// index is coded as (index+3)

				prev_orientation = anchor_orientation[orn_index];	// store current orientation in case we reject the move
				//anchor_orientation[orn_index] = calcOrientation(anchor_index);

				local_score_prev_orientation = calcScoreOrientation(anchor_orientation, orn_index);
			}
		}

		local_score_prev_structure = calcScoreStructureSmooth(p, true, true);
		if (use_subanchor_heatmap && Settings::useSubanchorHeatmap) local_score_prev_heat = calcScoreSubanchorHeatmap(p);

		tmp = random_vector(step_size, Settings::use2D);
		clusters[ind].pos += tmp;

		local_score_curr_structure = calcScoreStructureSmooth(p, true, true);
		if (orn_index != -1) {
			anchor_orientation[orn_index] = calcOrientation(anchor_index);
			local_score_curr_orientation = calcScoreOrientation(anchor_orientation, orn_index);
			assert(local_score_curr_orientation >= 0.0f);
			curr_score_orientation = curr_score_orientation + 2.0 * (local_score_curr_orientation - local_score_prev_orientation);
		}

		if (use_subanchor_heatmap && Settings::useSubanchorHeatmap) {
			local_score_curr_heat = calcScoreSubanchorHeatmap(p);
			curr_score_heat = curr_score_heat + 2.0 * (local_score_curr_heat - local_score_prev_heat);
		}

		curr_score_structure = curr_score_structure - local_score_prev_structure + local_score_curr_structure;
		//curr_score_structure = calcScoreStructureSmooth(true, true);

		score_curr = curr_score_structure + curr_score_orientation + curr_score_heat;

		ok = score_curr < score_prev;

		if (!ok && T > 0.0) {
			tp = Settings::tempJumpScaleSmooth * exp(-Settings::tempJumpCoefSmooth * (score_curr/score_prev) / T);
			ok = withChance(tp);
			//printf("%f %f %f %f   %d\n", score_curr, score_prev, T, tp, ok);
		}

		if (ok) {
			success++;
			milestone_success++;

			score_prev = score_curr;
			prev_score_structure = curr_score_structure;
			prev_score_orientation = curr_score_orientation;
			prev_score_heat = curr_score_heat;
		}
		else {
			// the move was rejected
			clusters[ind].pos -= tmp;		// restore previous position, ...
			score_curr = score_prev;  		// ... score ...
			curr_score_structure = prev_score_structure;
			curr_score_heat = prev_score_heat;

			if (orn_index >= 0) {
				curr_score_orientation = prev_score_orientation;
				anchor_orientation[orn_index] = prev_orientation; // ... and orientation (if needed)
			}
		}

		// check if we should stop
		if (i % Settings::MCstopConditionStepsSmooth == 0) {

			if (!silent) {
				output(7, "milestone: score = %lf (%lf %lf %lf), last = %lf (%lf), T=%lf successes = %d, last = %d\n", score_curr,
						curr_score_structure, curr_score_orientation, curr_score_heat, milestone_score, score_curr/milestone_score,
						T, success, milestone_success);			}

			if ((score_curr > Settings::MCstopConditionImprovementSmooth * milestone_score && milestone_success < Settings::MCstopConditionMinSuccessesSmooth) || score_curr < 1e-6)
				break;

			milestone_score = score_curr;
			milestone_success = 0;
		}

		T *= dt;  			// cooling

		i++;
	}

	//if (!silent) for (i=0; i<anchor_orientation.size(); i++) print_vector(anchor_orientation[i], "orient");
	return score_curr;
}

Heatmap LooperSolver::calculateRegionsDistances(const std::vector<int> clusters_ind) {
	int n = clusters_ind.size();
	Heatmap h(n);
	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			h.v[i][j] = (clusters[clusters_ind[i]].pos - clusters[clusters_ind[j]].pos).length();
			h.v[j][i] = h.v[i][j];
		}
	}
	return h;
}


// works only for chromosome and segment level
string LooperSolver::activeRegionIndexToChr(int ind) {
	for (string chr: chrs) {
		if (ind < active_region_to_chr[chr]) return chr;
	}
	return "";
}

std::set<int> LooperSolver::getCurrentHeatmapBreaks() {
	std::set<int> r;
	int curr_ind = 0;
	string chr;
	for (size_t j=0; j<chrs.size(); j++) {
		chr = chrs[j];
		curr_ind += current_level[chr].size();
		r.insert(curr_ind);
	}
	return r;
}

void LooperSolver::setDebuggingOptions(int chr_number_limit, int length_limit) {
	debug_chromosomes_limit = chr_number_limit;
	debug_length_limit = length_limit;
}

// calculate and save orientation
// orientation depends on the neighbors position (tangent) and motif orientation
// first find the tangent assuming orientation is R, and then, if it is L, just take the negative
vector3 LooperSolver::calcOrientation(int cind) {
	vector3 orn;
	if (cind==0) orn = clusters[active_region[cind+1]].pos - clusters[active_region[cind]].pos;
	else if (cind==(int)active_region.size()-1) orn = clusters[active_region[cind]].pos - clusters[active_region[cind-1]].pos;
	else orn = clusters[active_region[cind+1]].pos - clusters[active_region[cind-1]].pos;

	if (clusters[active_region[cind]].orientation == 'L') orn *= -1.0;

	orn.normalize();
	return orn;
}


void LooperSolver::showCoordinatesRange(int level) {

	vector3 min(1e6, 1e6, 1e6);
	vector3 max(-1e6, -1e6, -1e6);
	size_t dim = Settings::use2D ? 2 : 3;


	HierarchicalChromosome hc = getModel();
	hc.setLevel(level);	// go to the lowest level

	for (string chr: chrs) {
		for (size_t i = 0; i < hc.chr[chr].points.size(); ++i) {
			for (size_t j = 0; j < dim; ++j) {
				if (hc.chr[chr].points[i][j] < min[j]) min[j] = hc.chr[chr].points[i][j];
				if (hc.chr[chr].points[i][j] > max[j]) max[j] = hc.chr[chr].points[i][j];
			}
		}
	}


	if (dim == 2) {
		min[2] = 0.0f;
		max[2] = 0.0f;
	}

	print_vector(min, "min");
	print_vector(max, "max");
	print_vector(max-min, "diameter");

}

// if active_region then show coordinates range for the current region
// otherwise do it for all chromosomes for the subanchor level
void LooperSolver::showCoordinatesRange(bool use_active_region) {


	if (use_active_region) {

		vector3 min(1e6, 1e6, 1e6);
		vector3 max(-1e6, -1e6, -1e6);
		size_t dim = Settings::use2D ? 2 : 3;

		for (size_t i = 0; i < active_region.size(); ++i) {
			for (size_t j = 0; j < dim; ++j) {
				if (clusters[active_region[i]].pos[j] < min[j]) min[j] = clusters[active_region[i]].pos[j];
				if (clusters[active_region[i]].pos[j] > max[j]) max[j] = clusters[active_region[i]].pos[j];
			}
		}

		if (dim == 2) {
			min[2] = 0.0f;
			max[2] = 0.0f;
		}

		print_vector(min, "min");
		print_vector(max, "max");
		print_vector(max-min, "diameter");
	}
	else {
		showCoordinatesRange(5);
	}
}

void LooperSolver::calcAnchorsStats() {
	setLevel(LVL_SEGMENT);
	std::map<std::string, std::vector<int> > breaks;
	std::map<std::string, std::vector<int> > seg_length;
	std::map<std::string, std::vector<int> > inter_seg_length;
	vector<int> ib_in_segments; // number of ibs in a segment
	vector<int> ib_in_segments_chr;

	std::vector<int> seg_length_all;
	std::vector<int> seg_length_chr;

	int curr, next;

	for (string chr: chrs) {
		seg_length.clear();
		seg_length_chr.clear();

		ib_in_segments_chr.clear();

		breaks[chr].push_back(0);
		for (size_t i = 0; i < current_level[chr].size(); i++) {
			curr = current_level[chr][i];
			next = current_level[chr][i+1];
			seg_length[chr].push_back(clusters[curr].end - clusters[curr].start);
			seg_length_chr.push_back(clusters[curr].end - clusters[curr].start);
			seg_length_all.push_back(clusters[curr].end - clusters[curr].start);

			ib_in_segments_chr.push_back(clusters[curr].children.size());
			ib_in_segments.push_back(clusters[curr].children.size());

			if (i+1 < current_level[chr].size()) {
				breaks[chr].push_back((clusters[curr].end + clusters[next].start) / 2);
				inter_seg_length[chr].push_back(clusters[next].start - clusters[curr].end);
			}
		}
		breaks[chr].push_back(1e9);

		//saveToCSV(ftext("/d/tmp/motifs_gb/segment_len_%s.txt", chr.c_str()), seg_length_chr);

		//saveToCSV(ftext("/d/tmp/motifs_gb/ib_in_segments_%s.txt", chr.c_str()), ib_in_segments_chr);
	}

	printm(breaks, false, "split positions");
	printm(seg_length, true, "segments len");
	printf("\n");
	printm(inter_seg_length, true, "inter-segments len");
	printf("\n");

	//saveToCSV("/d/tmp/motifs_gb/segment_len.txt", seg_length_all);

	//saveToCSV("/d/tmp/motifs_gb/ib_in_segments.txt", ib_in_segments);


	setLevel(LVL_INTERACTION_BLOCK);
	vector<int> ib_length;   // ib length
	vector<int> ib_length_chr;
	vector<int> ib_dist;	// distance between neighboring anchors

	for (string chr: chrs) {

		setLevel(LVL_INTERACTION_BLOCK);
		int size = current_level[chr].size();
		printf("%s, size=%d\n", chr.c_str(), size);

		ib_length_chr.clear();
		for (int i = 0; i < size; ++i) {
			int ind = current_level[chr][i];
			ib_length.push_back(clusters[ind].end-clusters[ind].start);
			ib_length_chr.push_back(clusters[ind].end-clusters[ind].start);
			if (i+1<size) ib_dist.push_back(clusters[ind+1].start-clusters[ind].end);
			//			printf("%s %d %d %d %d\n", chr.c_str(), clusters[ind].start, clusters[ind].end, clusters[ind].base_start, clusters[ind].base_end);
		}

		//saveToCSV(ftext("/d/tmp/motifs_gb/ib_len_%s.txt", chr.c_str()), ib_length_chr);

	}


	//saveToCSV("/d/tmp/motifs_gb/ib_len.txt", ib_length);
	//saveToCSV("/d/motifs_gb/anchor_dist.txt", anchor_dist);

	exit(0);


	// check singletons - anchors relation
	// points of interest:
	// - how many falls completely inside one segment
	// - how many are inside one loop

	int between_segments = 0;   // starting between segments
	int inter_segmental = 0;
	int intra_segmental = 0;
	int intra_loops = 0;
	int inter_loops = 0;
	int starting_in_anchor = 0;
	int starting_in_loop = 0;

	setLevel(LVL_SEGMENT);

	// create a set of chromosomes to speed up checking if chromosome is to be processed
	std::set<std::string> chrs_set;
	for (std::string chr: chrs) chrs_set.insert(chr);

	char chr1[8], chr2[8];
	int sta, stb, enda, endb, sc;

	// read data from files
	for (size_t i = 0; i < arcs_singletons.size(); ++i) {
		printf("\nread [%s]\n", arcs_singletons[i].c_str());

		FILE *f = open(arcs_singletons[i].c_str(), "r");
		if (!f) exit(0);

		int k = 0;
		while (!feof(f)) {
			if (k++ % 1000000 == 0) printf(".");

			if (fscanf(f, "%s %d %d %s %d %d %d", chr1, &sta, &stb, chr2, &enda, &endb, &sc) < 7) continue;

			// check whether we should process current chromosomes
			if (chrs_set.find(chr1) == chrs_set.end()) continue;
			if (chrs_set.find(chr2) == chrs_set.end()) continue;

			// find indices of start and end of the arc
			int segm_ind_a = -1, segm_ind_b = -1;	// cluster-index of segment
			for (size_t j = 0; j < current_level[chr1].size(); ++j) {
				if (clusters[current_level[chr1][j]].contains(sta)) segm_ind_a = current_level[chr1][j];
				if (clusters[current_level[chr1][j]].contains(enda)) segm_ind_b = current_level[chr1][j];
			}

			if (segm_ind_a==-1 || segm_ind_b==-1) {
				between_segments++;
				continue;
			}

			if (segm_ind_a != segm_ind_b) {
				inter_segmental++;
				continue;
			}

			intra_segmental++;

			// singleton is inside one segment. check anchors
			int anch_ind_a = -1, anch_ind_b = -1;   // cluster-index of anchor within current segment
			bool inside_anchor_a, inside_anchor_b;   // is sigleton starting within anchor, or between-anchor region?
			for (size_t j = 0; j < clusters[segm_ind_a].children.size(); ++j) {
				int ch_ind = clusters[segm_ind_a].children[j];

				if (anch_ind_a == -1) {
					if (sta < clusters[ch_ind].start) {
						anch_ind_a = ch_ind;
						inside_anchor_a = false;
					}
					else if (clusters[ch_ind].contains(sta)) {
						anch_ind_a = ch_ind;
						inside_anchor_a = true;
					}
				}

				if (anch_ind_b == -1) {
					if (enda < clusters[ch_ind].start) {
						anch_ind_b = ch_ind;
						inside_anchor_b = false;
					}
					else if (clusters[ch_ind].contains(enda)) {
						anch_ind_b = ch_ind;
						inside_anchor_b = true;
					}
				}

				if (anch_ind_a != -1 && anch_ind_b != -1) break;
			}

			//			printf("anch ind = %d %d (%d %d)\n", anch_ind_a, anch_ind_b, inside_anchor_a, inside_anchor_b);
			//			if (anch_ind_a != -1) {
			//				if (anch_ind_a > 0) clusters[anch_ind_a-1].print();
			//				clusters[anch_ind_a].print();
			//			}
			//			if (anch_ind_b != -1) {
			//				if (anch_ind_b > 0) clusters[anch_ind_b-1].print();
			//				clusters[anch_ind_b].print();
			//			}

			if (anch_ind_a==-1 || anch_ind_b==-1) error("anchor id == -1");

			if (inside_anchor_a) starting_in_anchor++;
			else starting_in_loop++;
			if (inside_anchor_b) starting_in_anchor++;
			else starting_in_loop++;

			if (anch_ind_a==anch_ind_b || (abs(anch_ind_a-anch_ind_b)==1 && inside_anchor_a!=inside_anchor_b)) intra_loops++;
			else inter_loops++;
		}

		fclose(f);
	}

	printf("between segments: %d\n", between_segments);
	printf("intra/inter segments: %d %d\n", intra_segmental, inter_segmental);
	printf("intra/inter loops: %d %d\n", intra_loops, inter_loops);
	printf("start in anchor=%d loops=%d\n", starting_in_anchor, starting_in_loop);
}

Heatmap LooperSolver::calcTrueDistancesHeatmapForRegion(std::vector<int> &regions) {
	Heatmap h;
	int i,j;
	int n = regions.size();
	h.init(n);

	float d;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			d = (clusters[regions[i]].pos-clusters[regions[j]].pos).length();
			h.v[i][j] = h.v[j][i] = d;
		}
	}

	return h;
}


// create heatmap with expected distances between beads on subanchor level
// use heatmap with avg distance between beads and subanchor level singleton heatmap
// distances are based on avg distances, but we also consider subanchor heatmap
// assumes that 'heatmap_subanchor' is created
// result is saved to 'heatmap_dist_subanchor'
void LooperSolver::createExpectedDistSubanchorHeatmap(const Heatmap &subanchor_avg_dist) {
	if (subanchor_avg_dist.size != heatmap_subanchor.size) error("subanchor heatmaps size mismatch!");

	int n = subanchor_avg_dist.size;
	heatmap_dist_subanchor.init(n);

	//float min_val, max_val;
	//heatmap_subanchor.getRange(min_val, max_val);
	float avg_val = heatmap_subanchor.getAvg();
	printf("avg = %f\n", avg_val);
	if (avg_val < 1e-6) return;	// heatmap is empty, we can do nothing

	float s;
	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			s = (heatmap_subanchor.v[i][j] / avg_val) * Settings::subanchorHeatmapInfluence;
			if (s > 1.0f) s = 1.0f;
			heatmap_dist_subanchor.v[i][j] = subanchor_avg_dist.v[i][j] * (1.0f - s);
			heatmap_dist_subanchor.v[j][i] = heatmap_dist_subanchor.v[i][j];
		}
	}
}

// find all anchors in the active_region, and then their neighbors
// indices are relative to active_region
// save result to active_anchors_neighbors to speed up MC simulation
void LooperSolver::calcActiveAnchorsNeighbors() {
	int ai, arc;
	int other_end;

	active_anchors_neighbors.clear();
	active_anchors_neighbors_weight.clear();

	// create list of anchors
	vector<int> anchors;
	std::map<int, int> cluster_to_active_anchor_index;
	for (size_t i = 0; i < active_region.size(); ++i) {
		ai = active_region[i];
		if (clusters[ai].arcs.size() || clusters[ai].is_fixed) {
			cluster_to_active_anchor_index[ai] = anchors.size();
			anchors.push_back(i);
		}
	}

	// for every anchor iterate through all outcoming arcs
	// if they end in another anchor (and they should) then add them (ie. their active_region-anchor index) to result
	float w;
	for (size_t i = 0; i < anchors.size(); ++i) {
		ai = active_region[anchors[i]];
		for (size_t j = 0; j < clusters[ai].arcs.size(); ++j) {
			arc = clusters[ai].arcs[j];
			other_end = otherEnd(arc, ai);	// cluster-index of anchor

			active_anchors_neighbors[i].push_back(cluster_to_active_anchor_index[other_end]);
			w = sqrt(arcs.arcs[current_chr][arc].score);
			active_anchors_neighbors_weight[i].push_back(w);
		}
	}
}

// assume active_size contains anchors only
// save result to heatmap_exp_dist_anchor
void LooperSolver::calcAnchorExpectedDistancesHeatmap() {
	int ai = 0, other_end = 0, arc = 0, freq;
	double exp_dist = 0.0;

	// init with -1 values
	heatmap_exp_dist_anchor.init(active_region.size());
	heatmap_exp_dist_anchor.add(-1.0f);
	heatmap_exp_dist_anchor.clearDiagonal(1);
	//printf("size = %d\n", active_region.size());

	active_to_cluster_index.clear();

	std::map<int, int> cluster_to_active_index;
	for (size_t i = 0; i < active_region.size(); ++i) {
		cluster_to_active_index[active_region[i]] = i;
		active_to_cluster_index[i] = active_region[i];
	}

	for (size_t i = 0; i < active_region.size(); ++i) {
		ai = active_region[i];
		for (size_t j = 0; j < clusters[ai].arcs.size(); ++j) {
			arc = clusters[ai].arcs[j];

			other_end = otherEnd(arc, ai);
			if (other_end == -1) {
				printf("no other end (%d %d, %d %d)!\n", arc, ai, arcs.arcs[current_chr][arc].start, arcs.arcs[current_chr][arc].end);
				exit(0);
				continue;
			}

			if (other_end < ai) continue;	// only consider forward loops (to count each arc only once)

			// note: 'arc' is chromosome specific, but we know that we work on 'current_chr' right now
			freq = arcs.arcs[current_chr][arc].score;
			exp_dist = freqToDistance(freq, true);

			int a = cluster_to_active_index[ai];
			int b = cluster_to_active_index[other_end];

			heatmap_exp_dist_anchor.v[a][b] = exp_dist;
			heatmap_exp_dist_anchor.v[b][a] = exp_dist;
		}
	}

	// consider anchor heatmap
	if (Settings::useAnchorHeatmap) {

		if (heatmap_exp_dist_anchor.size != heatmap_anchor.size) error("anchor heatmaps size mismatch!");

		float min_val, max_val;
		//float avg_val = heatmap_anchor.getAvg();
		//printf("avg anchor heatmap val = %f\n", avg_val);
		heatmap_anchor.getRange(min_val, max_val);
		printf("max anchor heatmap val = %f\n", max_val);

		// only consider heatmap if it is not empty (we can't divide by max_val = 0)
		if (max_val > 1e-6) {
			int n = heatmap_anchor.size;
			float s;
			for (int i = 0; i < n; ++i) {
				for (int j = i+1; j < n; ++j) {
					if (heatmap_exp_dist_anchor.v[i][j] < 0.0f) continue;

					s = (heatmap_anchor.v[i][j] / max_val) * Settings::anchorHeatmapInfluence;
					if (s > 1.0f) s = 1.0f;
					heatmap_exp_dist_anchor.v[i][j] *= (1.0f - s);
					heatmap_exp_dist_anchor.v[j][i] = heatmap_exp_dist_anchor.v[i][j];
				}
			}
		}
	}
}

int LooperSolver::getGenomicLengthExcludingCentromere(string chr, int st, int end) {
	if (centromeres.regions.size() == 0) return end - st + 1;	// no info about centromeres
	for (size_t i = 0; i < centromeres.regions.size(); ++i) {
		if (centromeres.regions[i].chr == chr) {
			int cs = centromeres.regions[i].start;
			int ce = centromeres.regions[i].end;
			if (cs > end || ce < st) return end - st + 1;    //  p q  A  B
			if (cs < st) {
				if (ce < end) return end - ce + 1;    // p A  q  B
				return 0;    // p A    B  q
			}
			if (ce < end) return end - st - (ce-cs+1) + 1;    //   A   p  q  B
			return cs - st;    // A   p    B  q
		}
	}
	return end - st + 1;
}

void LooperSolver::printStructuresForGenomicPosition(string chr, int pos) {
	int curr_lev, clust;
	printf("look for position: %s %d\n", chr.c_str(), pos);
	// start with segments
	setLevel(LVL_SEGMENT);
	curr_lev = findCurrentLevelForGenomicPosition(chr, pos);
	if (curr_lev == -1) printf("no segment found!");
	else {
		clust = current_level[chr][curr_lev];
		printf("cluster id = %d, curr_level id = %d\n", clust, curr_lev);
		clusters[clust].print();

		setLevel(LVL_INTERACTION_BLOCK);
		curr_lev = findCurrentLevelForGenomicPosition(chr, pos);
		if (curr_lev == -1) printf("no ib found!");
		else {
			clust = current_level[chr][curr_lev];
			printf("cluster id = %d, curr_level id = %d\n", clust, curr_lev);
			clusters[clust].print();
		}
	}
}

int LooperSolver::findClusterForGenomicPosition(string chr, int pos) {
	for (size_t i = 0; i < current_level[chr].size(); ++i) {
		int c = current_level[chr][i];
		if (clusters[c].contains(pos)) return c;
	}
	return -1;
}

int LooperSolver::findCurrentLevelForGenomicPosition(string chr, int pos) {
	for (size_t i = 0; i < current_level[chr].size(); ++i) {
		int c = current_level[chr][i];
		if (clusters[c].contains(pos)) return i;
	}
	return -1;
}

void LooperSolver::output(int level, const char* format, ...) {
	if (level == 0 || level <= Settings::outputLevel) {
		va_list vl;
		va_start(vl, format);
		vprintf(format, vl);
		va_end(vl);
	}
}

void LooperSolver::inputArbitraryHeatmap(Heatmap h) {
	clusters.clear();
	active_region.clear();

	int n = h.size;
	string chr_name = "chrN";

	Cluster root_chr(1000, 1000*n+100);
	for (int i=1; i<=n; i++) root_chr.children.push_back(i);
	root_chr.parent = -1;
	clusters.push_back(root_chr);

	chr_root[chr_name] = 0;

	for (int i = 0; i < n; ++i) {
		Cluster c(1000*(i+1), 1000*(i+1)+100);
		c.parent = 0;
		clusters.push_back(c);
	}

	chrs.clear();
	chrs.push_back(chr_name);

	current_level.clear();
	current_level[chr_name].push_back(0);

	heatmap.clear();

	Heatmap tmp;
	heatmap.push_back(tmp);
	heatmap.push_back(h);
}

Heatmap LooperSolver::createIFHeatmapFromDistances(const Heatmap &heat) {
	Heatmap h(heat.size);
	for (size_t i = 0; i < heat.size; ++i) {
		for (size_t j = 0; j < heat.size; ++j) {
			h.v[i][j] = (i==j) ? 0.0f : distToFreqHeatmap(heat.v[i][j]);
		}
	}
	return h;
}

void LooperSolver::getChromosomeHeatmapBoundary(int p, int &start, int &end) {
	if (heatmap_chromosome_boundaries.size() == 0) return;
	if (p < 0 || p > heatmap_chromosome_boundaries[heatmap_chromosome_boundaries.size()-1]) return;

	for (size_t i = 0; i+1 < heatmap_chromosome_boundaries.size(); ++i) {
		if (heatmap_chromosome_boundaries[i] <= p && p < heatmap_chromosome_boundaries[i+1]) {
			start = heatmap_chromosome_boundaries[i];
			end = heatmap_chromosome_boundaries[i+1]-1;
			return;
		}
	}
}

void LooperSolver::reset() {
	for (size_t i=0; i<clusters.size(); i++) {
		clusters[i].pos.set(0.0, 0.0, 0.0);
		clusters[i].is_fixed = false;
		clusters[i].dist_to_next = 0.0;
	}
}

// remove the subanchor beads (remove the corresponding clusters, clean cluster::children etc)
void LooperSolver::removeSubanchorBeads() {
	// we use the property that the subanchor beads are placed on the very end of the 'clusters', and that for a given chromosome
	// the root is always the last non-subanchor entry. Thus, to locate all subanchors we only need to find the index of the last
	// root. All the following entries are subanchor beads.
	// To remove the subanchor beads we need to:
	// - remove the corresponding entries from 'clusters'
	// - remove all the links to them, ie. in the cluster::children list

	for (string chr: chrs) printf("%s %d\n", chr.c_str(), chr_root[chr]);

	// first, find the subanchor beads location
	uint i;
	int first_subanchor_bead = -1;
	for (i = 0; i < chrs.size(); ++i) {
		if (chr_root[chrs[i]] > first_subanchor_bead) first_subanchor_bead = chr_root[chrs[i]];
	}
	if (first_subanchor_bead == -1) error("could not find the subanchor beads starting index");
	//printf("root = %d\n", first_subanchor_bead);
	first_subanchor_bead++; // it was pointing to the last root, now to the first subanchor

	//while (first_subanchor_bead<clusters.size() &&
	//		clusters[first_subanchor_bead].start != clusters[first_subanchor_bead].end) first_subanchor_bead++;

	//printf("last non-subanchor:\n");
	//clusters[first_subanchor_bead-1].print();
	//printf("first subanchor:\n");
	//clusters[first_subanchor_bead].print();

	clusters.erase(clusters.begin()+first_subanchor_bead, clusters.end());	// remove the subanchor beads

	// remove their indicies from 'children' list
	for (i = 0; i < clusters.size(); ++i) {
		if (clusters[i].children.size() == 0) continue;

		for (int j = clusters[i].children.size()-1; j >= 0; --j) {
			if (clusters[i].children[j] >= first_subanchor_bead) clusters[i].children.erase(clusters[i].children.begin()+j);
		}
	}
}

vector3 LooperSolver::densityCoordToStructure(vector3 coord) {

	vector3 shift = density.origin - density.center;
	vector3 pos_map = vector3(coord.x+0.5f, coord.y+0.5f, coord.z+0.5f);
	vector3 pos_3d = (pos_map+shift) / Settings::densityScale;
	pos_3d.x *= 3.0f;
	pos_3d += density.center;
	return pos_3d;
}

vector3 LooperSolver::structureCoordToDensity(vector3 pos) {
	pos = (pos - density.center) * Settings::densityScale;
	pos.x /= 3.0f;
	pos += density.center - density.origin;
	return vector3((int)pos.x, (int)pos.y, (int)pos.z);
}

