/*
 * InteractionArcs.h
 *
 *  Created on: Apr 18, 2014
 *      Author: psz
 */

#ifndef INTERACTIONARCS_H_
#define INTERACTIONARCS_H_

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>

#include "../src/lib/common.h"
#include "Anchor.h"
#include "InteractionArc.h"
#include "Settings.h"
#include "BedRegion.h"
#include "BedRegions.h"

class InteractionArcs {
public:
	InteractionArcs();

	bool fromFile(string filename);
	void toFile(string filename);

	void clear();

	void selectRegion(BedRegion region); // restrict region to reconstruct (only anchors and clusters from this region will be loaded)

	void loadAnchorsData(string anchors_path);
	void loadPetClustersData(string pet_clusters_path, string factor_name, BedRegions predefined_segments);	// read arcs data, store result in 'raw_arcs'

	void print(int display_limit = 10);

	// update arcs to use anchor indicies rather than genomic position. input: 'anchors' and 'raw_arcs', output: 'arcs'
	void markArcs(bool ignore_missing);
	void parallelMarkArcs(bool ignore_missing);

	// because we filter long arcs there may be anchors without any arcs. this function removes such anchors
	// (should be called after marcArcs(), and then marcArcs() should be called again to update arcs information)
	void removeEmptyAnchors();

	void rewire();

	std::vector<string> chrs;		// list of chromosomes

	std::vector<string> factors;	// list of factors used

	// for each chromosome keep number of anchors and arcs
	std::unordered_map<string, int> anchors_cnt;
	std::unordered_map<string, int> arcs_cnt;

	std::unordered_map<string, std::vector<Anchor> > anchors;
	std::unordered_map<string, std::vector <InteractionArc> > arcs;		// list of arcs representing contacts (as indices for 'points')
	std::unordered_map<string, std::vector <InteractionArc> > raw_arcs;	// list of raw arcs representing contacts (as genomic positions)

	std::unordered_map<string, std::vector <InteractionArc> > long_arcs;	// list of long arcs (>pet_max_length), used to refine the segment level heatmap

	std::vector <float> expected_distances;

	BedRegion selected_region;
};

#endif /* INTERACTIONARCS_H_ */
