/*
 * common.h
 *
 *  Created on: Aug 4, 2013
 *      Author: psz
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <map>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "mtxlib.h"

using namespace std;

#define epsilon 1e-6

#define LVL_CHROMOSOME 0
#define LVL_SEGMENT 1
#define LVL_INTERACTION_BLOCK 2
#define LVL_SUBANCHOR 3

string ftext(const string fmt_str, ...);

int random(int range);

float random_uniform();

float random(float range, bool negative = false);

vector3 random_vector(float max_size, bool in2D = false);
matrix44 random_rot_matrix(float max_angle = 180.0f);

float dist(const vector3 &v1, const vector3 &v2);
float angle(vector3 v1, vector3 v2);
float angle_norm(vector3 v1,
                 vector3 v2); // calc angle for already normalized vectors

vector3 displace(vector3, float displacement, bool in2D = false);
void displace_ref(vector3 &v, float displacement);

float distanceToInteractionFrequency(float distance);

float interpolate(float a, float b, float p);
vector3 interpolate(vector3 a, vector3 b, float p);

void random_shuffle(float *arr, int size);

void printv(const vector<int> &v, bool count = false, bool newline = false,
            const char *label = "");
void printv(const vector<float> &v, bool count = false, bool newline = false,
            const char *label = "");
void printv(const vector<string> &v, bool count = false, bool newline = false,
            const char *label = "");
void printv(vector<bool> v);
void printm(const std::map<int, int> &map, const char *label = "");
void printm(const std::map<int, vector<int>> &map, const char *label = "");
void printm(const std::map<int, vector<float>> &map, const char *label = "");
void printm(const std::map<std::string, int> &map, bool flat = false,
            bool newline = false, const char *label = "");
void printm(std::map<std::string, int> &map,
            const std::vector<std::string> &keys, bool flat = false,
            bool newline = false, const char *label = "");
void printm(const std::map<std::string, std::vector<int>> &map,
            bool flat = false, const char *label = "");

void print_vector(const vector3 &v, const char *label = "");

bool withChance(float chance);

void saveToCSV(string path, vector<int> v, bool add_index_column = false);
void saveToCSV(string path, vector<float> v, bool add_index_column = false);
void saveToCSV(string path, vector<bool> v, bool add_index_column = false);
void saveToCSV(string path, vector<vector<float>> v,
               bool add_index_column = false);
void saveToCSV(string path, vector<vector<int>> v,
               bool add_index_column = false);

bool file_exists(const std::string &name);
FILE *open(string path, string mode = "r");
// char* concat(const char *a, const char *b);
// char* ftext(const char *fmt, ...);

// string fstr(string s, ...);

int countWords(char *line);

// vector<float> interpolateInterval(float posl, float posr, float dl, float dr,
// float mn = 1.5f); vector<float> interpolatePoints(vector<float> points, float
// mn = 1.5f);

vector<int> interpolateInterval(int posl, int posr, int dl, int dr,
                                float mn = 1.5f);
vector<int> interpolatePoints(vector<int> points, float mn = 1.5f);

bool vector_contains(const std::vector<int> &v, int val);
void vector_insert_unique(std::vector<int> &v, int val);

void error(string msg);

vector3 interpolateSpline(float t, const vector3 &p1, const vector3 &p2,
                          const vector3 &p3, const vector3 &p4);
std::vector<vector3> interpolateSpline(const std::vector<vector3> &points,
                                       int splits);

std::vector<vector3> interpolateSplineCentripetal(std::vector<vector3> &points,
                                                  int splits);
vector3 interpolate(vector3 &p0, vector3 &p1, vector3 &p2, vector3 &p3,
                    double time[], double t);

vector3 mirrorPoint(vector3 fixed, vector3 pt);

void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim = ',');
std::vector<std::string> split_unknown(const std::string &s);

void split_file_path(const std::string &s, std::string &dir,
                     std::string &filename);

#endif /* COMMON_H_ */
