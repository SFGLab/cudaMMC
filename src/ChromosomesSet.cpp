/*
 * ChromosomesSet.cpp
 *
 *  Created on: May 25, 2014
 *      Author: psz
 */

#include "../include/ChromosomesSet.h"

ChromosomesSet::ChromosomesSet() {
	// TODO Auto-generated constructor stub

}


void ChromosomesSet::print() {
	printf("Set size = %lu\n", chromosome.size());

	for (size_t i = 0; i < chromosome.size(); ++i) {
	//	printf("%d %s\n", chromosome[i].points.size(), desc[i].c_str());
	}
}

void ChromosomesSet::add(std::map<std::string, Chromosome> chr) {
	add(chr, "<no_desc>");
}

void ChromosomesSet::add(std::map<std::string, Chromosome> chr, string desc) {
	chromosome.push_back(chr);
	if (desc == "") desc = "none";
	this->desc.push_back(desc);
}


void ChromosomesSet::toFile(string filename) {
	FILE *f;
	f = fopen(filename.c_str(), "w");
	if (f == NULL)
	{
		printf("Error opening file [%s]!\n", filename.c_str());
		return;
	}

	toFile(f);

	fclose(f);
}

void ChromosomesSet::toFile(FILE* file) {
	//printf("set = %d\n", chromosome.size());
	fprintf(file, "%lu\n", chromosome.size());
	for (size_t i=0; i<chromosome.size(); i++) {
		//fprintf(file, "%d\n%s\n", chromosome[i].points.size(), desc[i].c_str());
		fprintf(file, "%lu %s\n", chromosome[i].size(), desc[i].c_str());
		for (auto el: chromosome[i]) {
			fprintf(file, "%s %lu\n", el.first.c_str(), el.second.points.size());
			el.second.toFile(file);
			//chromosome[i].toFile(file);
		}
	}
}

void ChromosomesSet::fromFile(string filename) {
	FILE *f = open(filename, "r");
	if (f == NULL) return;
	fromFile(f);
	fclose(f);
}

void ChromosomesSet::fromFile(FILE* file) {
	chromosome.clear();

	char de[100], str_chr[10];
	int n, n_chr;
	fscanf(file, "%d", &n);

	int pts;
	for (int i=0; i<n; i++) {
		fscanf(file, "%d %100s", &n_chr, de);
		string s(de);
		desc.push_back(s);

		std::map<std::string, Chromosome> map_chr;
		for (int j = 0; j < n_chr; ++j) {
			fscanf(file, "%s %d", str_chr, &pts);
			Chromosome chr;
			chr.fromFile(file, pts);
			//chromosome.push_back(chr);
			map_chr[str_chr] = chr;
		}

		chromosome.push_back(map_chr);
	}
}
