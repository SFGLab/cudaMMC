/*
 * Settings.cpp
 *
 *  Created on: Aug 4, 2013
 *      Author: psz
 */

#include "../include/Settings.h"

bool Settings::initialized;
bool Settings::use2D;

int Settings::debug;
int Settings::outputLevel;
bool Settings::randomWalk;

int Settings::milestoneFailsThreshold;
int Settings::cudaBlocksMultiplier;
int Settings::cudaThreadsPerBlock;

bool Settings::useCTCFMotifOrientation;
bool Settings::motifsSymmetric;
float Settings::motifOrientationWeigth;

bool Settings::useSubanchorHeatmap;
int Settings::subanchorEstimateDistancesReplicates;
int Settings::subanchorEstimateDistancesSteps;
float Settings::subanchorHeatmapInfluence;
float Settings::subanchorHeatmapDistWeight;

bool Settings::useAnchorHeatmap;
float Settings::anchorHeatmapInfluence;
float Settings::anchorHeatmapDistWeight;

int Settings::loopDensity;

int Settings::segmentSize;
int Settings::maxPETClusterLength;
float Settings::longPETClustersEffectPower;
float Settings::longPETClustersEffectScale;

bool Settings::rewiringEnabled;
int Settings::rewiringCertaintyThreshold;

float Settings::ibRandomWalkJumps;

std::string Settings::dataDirectory;
std::string Settings::dataAnchors;
std::string Settings::dataPetClusters;
std::string Settings::dataSingletons;
std::string Settings::dataSingletonsInter;
bool Settings::dataSplitSingletonFilesByChr;
std::string Settings::dataFactors;

std::string Settings::dataCentromeres;
std::string Settings::dataSegmentsSplit;

std::string Settings::dataSegmentHeatmap;

std::string Settings::templateSegment;
float Settings::templateScale;

std::string Settings::distHeatmap;
float Settings::distHeatmapScale;

bool Settings::useDensity;
std::string Settings::densityMapFile;
float Settings::densityScale;
float Settings::densityInfluence;
float Settings::densityWeight;

bool Settings::useTelomerePositions;
vector3 Settings::telomere_1_position;
vector3 Settings::telomere_2_position;

int Settings::simulationStepsLevelChr;
int Settings::simulationStepsLevelSegment;
int Settings::simulationStepsLevelAnchor;
int Settings::simulationStepsLevelSubanchor;

float Settings::noiseCoefficientLevelChr;
float Settings::noiseCoefficientLevelSegment;
float Settings::noiseCoefficientLevelAnchor;
float Settings::noiseCoefficientLevelSubanchor;

float Settings::heatmapInterScaling;
float Settings::heatmapDistanceHeatmapStretching;

float Settings::freqToDistHeatmapScale;
float Settings::freqToDistHeatmapPower;
float Settings::freqToDistHeatmapScaleInter;
float Settings::freqToDistHeatmapPowerInter;
float Settings::countToDistA;
float Settings::countToDistScale;
float Settings::countToDistShift;
float Settings::countToDistBaseLevel;
float Settings::genomicLengthToDistPower;
float Settings::genomicLengthToDistScale;
float Settings::genomicLengthToDistBase;

float Settings::springConstantStretch;
float Settings::springConstantSqueeze;
float Settings::springAngularConstant;

float Settings::springConstantStretchArcs;
float Settings::springConstantSqueezeArcs;

float Settings::maxTemp;
float Settings::dtTemp;
float Settings::tempJumpScale;
float Settings::tempJumpCoef;
float Settings::MCstopConditionImprovement;
int Settings::MCstopConditionMinSuccesses;
int Settings::MCstopConditionSteps;

float Settings::weightDistSmooth;
float Settings::weightAngleSmooth;
float Settings::maxTempSmooth;
float Settings::dtTempSmooth;
float Settings::tempJumpScaleSmooth;
float Settings::tempJumpCoefSmooth;
float Settings::MCstopConditionImprovementSmooth;
int Settings::MCstopConditionMinSuccessesSmooth;
int Settings::MCstopConditionStepsSmooth;

float Settings::maxTempHeatmap;
float Settings::dtTempHeatmap;
float Settings::tempJumpScaleHeatmap;
float Settings::tempJumpCoefHeatmap;
float Settings::MCstopConditionImprovementHeatmap;
int Settings::MCstopConditionMinSuccessesHeatmap;
int Settings::MCstopConditionStepsHeatmap;

float Settings::maxTempHeatmapDensity;
float Settings::dtTempHeatmapDensity;
float Settings::tempJumpScaleHeatmapDensity;
float Settings::tempJumpCoefHeatmapDensity;
float Settings::MCstopConditionImprovementHeatmapDensity;
int Settings::MCstopConditionMinSuccessesHeatmapDensity;
int Settings::MCstopConditionStepsHeatmapDensity;

bool Settings::useInputCache;

Settings::Settings() {
  if (!initialized)
    init();
}

void Settings::init() {
  initialized = true;
  debug = 0;
  outputLevel = 0;
  randomWalk = false;

  use2D = false;

  useInputCache = true;

  loopDensity = 5;

  milestoneFailsThreshold = 3;
  cudaBlocksMultiplier = 4;
  cudaThreadsPerBlock = 256;

  useCTCFMotifOrientation = false;
  motifsSymmetric = true;
  motifOrientationWeigth = 1.0f;

  useSubanchorHeatmap = false;
  subanchorEstimateDistancesReplicates = 5;
  subanchorEstimateDistancesSteps = 2;
  subanchorHeatmapInfluence = 0.5f;
  subanchorHeatmapDistWeight = 1.0f;

  useAnchorHeatmap = false;
  anchorHeatmapInfluence = 0.5f;
  anchorHeatmapDistWeight = 1.0f;

  segmentSize = 2000000;
  maxPETClusterLength = 1000000;
  longPETClustersEffectPower = 2.0f;
  longPETClustersEffectScale = 10.0f;

  rewiringCertaintyThreshold = 16;

  simulationStepsLevelChr = 2;
  simulationStepsLevelSegment = 2;
  simulationStepsLevelAnchor = 5;
  simulationStepsLevelSubanchor = 5;

  templateScale = 1.0;
  distHeatmapScale = 1.0;

  useDensity = false;
  densityScale = 1.0f;
  densityInfluence = 0.95f;
  densityWeight = 1.0f;

  ibRandomWalkJumps = 10.0f;

  noiseCoefficientLevelChr = 1.0f;
  noiseCoefficientLevelSegment = 0.1f;
  noiseCoefficientLevelAnchor = 0.5f;
  noiseCoefficientLevelSubanchor = 0.5f;

  heatmapInterScaling = 1.0f;
  heatmapDistanceHeatmapStretching = 2.0f;

  freqToDistHeatmapScale = 100.0f;
  freqToDistHeatmapPower = -0.333;
  freqToDistHeatmapScaleInter = 100.0f;
  freqToDistHeatmapPowerInter = -1.0;

  countToDistA = 0.5f;
  countToDistScale = 20.0f;
  countToDistShift = 1.0f;
  countToDistBaseLevel = 0.01f;

  genomicLengthToDistPower = 0.5f;
  genomicLengthToDistScale = 1.0f;
  genomicLengthToDistBase = 0.0f;

  maxTempHeatmap = 20.0f;
  dtTempHeatmap = 0.99995;
  tempJumpCoefHeatmap = 20.0f;
  tempJumpScaleHeatmap = 50.0f;
  MCstopConditionImprovementHeatmap = 0.995f;
  MCstopConditionMinSuccessesHeatmap = 5;
  MCstopConditionStepsHeatmap = 10000;

  maxTempHeatmapDensity = 20.0f;
  dtTempHeatmapDensity = 0.99995;
  tempJumpCoefHeatmapDensity = 20.0f;
  tempJumpScaleHeatmapDensity = 50.0f;
  MCstopConditionImprovementHeatmapDensity = 0.995f;
  MCstopConditionMinSuccessesHeatmapDensity = 5;
  MCstopConditionStepsHeatmapDensity = 10000;

  maxTemp = 20.0f;
  dtTemp = 0.99995;
  tempJumpCoef = 20.0f;
  tempJumpScale = 50.0f;
  MCstopConditionImprovement = 0.995f;
  MCstopConditionMinSuccesses = 5;
  MCstopConditionSteps = 10000;

  weightAngleSmooth = 1.0f;
  weightDistSmooth = 1.0f;
  maxTempSmooth = 20.0f;
  dtTempSmooth = 0.99995;
  tempJumpCoefSmooth = 20.0f;
  tempJumpScaleSmooth = 50.0f;
  MCstopConditionImprovementSmooth = 0.995f;
  MCstopConditionMinSuccessesSmooth = 5;
  MCstopConditionStepsSmooth = 10000;

  springConstantSqueeze = 0.1f;
  springConstantStretch = 0.1f;
  springAngularConstant = 0.1f;

  springConstantSqueezeArcs = 1.0f;
  springConstantStretchArcs = 1.0f;
}

void Settings::print(int level) {

  // printf("debug: %d\n", debug);
  printf("cache input: %s\n", useInputCache ? "yes" : "no");
  printf("random walk: %s\n", randomWalk ? "yes" : "no");
  printf("2D: %s\n", use2D ? "yes" : "no");

  printf("use motif orientation: %s (symmetric: %s)\n",
         useCTCFMotifOrientation ? "yes" : "no",
         motifsSymmetric ? "yes" : "no");
  printf("use singletons on subanchor level: %s\n",
         useSubanchorHeatmap ? "yes" : "no");

  printf("loop density: %d\n", loopDensity);

  printf("use density: %s\n", useDensity ? "yes" : "no");
  if (useDensity) {
    printf("   density map file: %s\n", densityMapFile.c_str());
    printf("   density scale: %f\n", densityScale);
    printf("   density influence: %f\n", densityInfluence);
    printf("   density weight: %f\n", densityWeight);
  }

  printf("use rewiring: %s\n", rewiringEnabled ? "yes" : "no");
  if (rewiringEnabled) {
    printf("   certainty threshold: %d\n", rewiringCertaintyThreshold);
  }

  printf("use telomere positions: %s\n", useTelomerePositions ? "yes" : "no");
  if (useTelomerePositions) {
    print_vector(telomere_1_position, "   telomere 1 position");
    print_vector(telomere_2_position, "   telomere 2 position");
  }

  printf("data\n");
  printf("   data directory: %s\n", dataDirectory.c_str());
  printf("   anchors: %s\n", dataAnchors.c_str());
  printf("   PET clusters: %s\n", dataPetClusters.c_str());
  printf("   singletons: %s\n", dataSingletons.c_str());
  printf("   split singleton files by chr: %d\n", dataSplitSingletonFilesByChr);
  printf("   singletons (inter): %s\n", dataSingletonsInter.c_str());
  printf("   factors: %s\n", dataFactors.c_str());
  printf("   centromeres: %s\n", dataCentromeres.c_str());
  printf("   segment split: %s\n", dataSegmentsSplit.c_str());
  printf("   segment heatmap: %s\n", dataSegmentHeatmap.c_str());

  printf("template structure: [%s]\n", templateSegment.c_str());
  printf("template scale: %f\n", templateScale);

  // printf("segment size: %d\n", segmentSize);

  printf("simulation steps, lvl1: %d\n", simulationStepsLevelChr);
  printf("simulation steps, lvl2: %d\n", simulationStepsLevelSegment);
  printf("simulation steps, arcs: %d\n", simulationStepsLevelAnchor);
  printf("simulation steps, smooth: %d\n", simulationStepsLevelSubanchor);

  //	printf("noise coef, lvl1: %f\n", noiseCoefficientLevel1);
  //	printf("noise coef, lvl2: %f\n", noiseCoefficientLevel2);
  //	printf("noise coef, arcs: %f\n", noiseCoefficientLevelArcs);
  //	printf("noise coef, smooth: %f\n", noiseCoefficientLevelArcsSmooth);

  printf("heatmap distance stretching: %f\n", heatmapDistanceHeatmapStretching);

  printf("ib random walk jump size: %f\n", ibRandomWalkJumps);
  printf("frequency-distance scale: %f\n", freqToDistHeatmapScale);
  printf("frequency-distance power: %f\n", freqToDistHeatmapPower);
  printf("frequency-distance power (inter): %f\n", freqToDistHeatmapPowerInter);
  printf("count-distance A: %f\n", countToDistA);
  printf("count-distance scale: %f\n", countToDistScale);
  printf("count-distance shift: %f\n", countToDistShift);
  printf("count-distance min dist: %f\n", countToDistBaseLevel);
  printf("genomic length (kb) to dist power: %f\n", genomicLengthToDistPower);
  printf("genomic length (kb) to dist scale: %f\n", genomicLengthToDistScale);
  printf("genomic length (kb) to dist base: %f\n", genomicLengthToDistBase);

  printf("spring stretch: %f\n", springConstantStretch);
  printf("spring squeeze: %f\n", springConstantSqueeze);
  printf("spring angular: %f\n", springAngularConstant);
  printf("spring stretch (arcs): %f\n", springConstantStretchArcs);
  printf("spring squeeze (arcs): %f\n", springConstantSqueezeArcs);

  printf("maxTemp: %f\n", maxTemp);
  printf("dTemp: %f\n", dtTemp);
  printf("jump coef: %f\n", tempJumpCoef);
  printf("jump scale: %f\n", tempJumpScale);
  printf("stop condition steps: %d\n", MCstopConditionSteps);
  printf("stop condition improvement threshold: %f\n",
         MCstopConditionImprovement);
  printf("stop condition successes threshold: %d\n",
         MCstopConditionMinSuccesses);

  printf("maxTemp (smooth): %f\n", maxTempSmooth);
  printf("dTemp: %f\n", dtTempSmooth);
  printf("jump coef: %f\n", tempJumpCoefSmooth);
  printf("jump scale: %f\n", tempJumpScaleSmooth);
  printf("stop condition steps: %d\n", MCstopConditionStepsSmooth);
  printf("stop condition improvement threshold: %f\n",
         MCstopConditionImprovementSmooth);
  printf("stop condition successes threshold: %d\n",
         MCstopConditionMinSuccessesSmooth);

  printf("maxTemp (heatmap): %f\n", maxTempHeatmap);
  printf("dTemp (heatmap): %f\n", dtTempHeatmap);
  printf("jump coef: %f\n", tempJumpCoefHeatmap);
  printf("jump scale: %f\n", tempJumpScaleHeatmap);
  printf("stop condition steps (heatmap): %d\n", MCstopConditionStepsHeatmap);
  printf("stop condition improvement threshold (heatmap): %f\n",
         MCstopConditionImprovementHeatmap);
  printf("stop condition successes threshold (heatmap): %d\n",
         MCstopConditionMinSuccessesHeatmap);
}

bool Settings::loadFromINI(std::string ini_path) {
  INIReader reader(ini_path);

  if (reader.ParseError() < 0) {
    printf("Can't load [%s]'\n", ini_path.c_str());
    return false;
  }

  cudaBlocksMultiplier =
      reader.GetInteger("cuda", "blocks_multiplier", cudaBlocksMultiplier);
  cudaThreadsPerBlock =
      reader.GetInteger("cuda", "num_threads", cudaThreadsPerBlock);
  milestoneFailsThreshold =
      reader.GetInteger("cuda", "milestone_fails", milestoneFailsThreshold);

  // debug = reader.GetInteger("main", "debug", debug);
  outputLevel = reader.GetInteger("main", "output_level", outputLevel);
  randomWalk = reader.GetBoolean("main", "random_walk", randomWalk);

  useInputCache = reader.GetBoolean("main", "cache_input", useInputCache);
  use2D = reader.GetBoolean("main", "use_2D", use2D);

  loopDensity = reader.GetInteger("main", "loop_density", loopDensity);

  useCTCFMotifOrientation = reader.GetBoolean(
      "motif_orientation", "use_motif_orientation", useCTCFMotifOrientation);
  motifsSymmetric = reader.GetBoolean("motif_orientation", "symmetric_motifs",
                                      motifsSymmetric);
  motifOrientationWeigth = (float)reader.GetReal("motif_orientation", "weight",
                                                 motifOrientationWeigth);

  useSubanchorHeatmap = reader.GetBoolean(
      "subanchor_heatmap", "use_subanchor_heatmap", useSubanchorHeatmap);
  subanchorEstimateDistancesReplicates =
      reader.GetInteger("subanchor_heatmap", "estimate_distances_replicates",
                        subanchorEstimateDistancesReplicates);
  subanchorEstimateDistancesSteps =
      reader.GetInteger("subanchor_heatmap", "estimate_distances_steps",
                        subanchorEstimateDistancesSteps);
  subanchorHeatmapInfluence = (float)reader.GetReal(
      "subanchor_heatmap", "heatmap_influence", subanchorHeatmapInfluence);
  subanchorHeatmapDistWeight = (float)reader.GetReal(
      "subanchor_heatmap", "heatmap_dist_weight", subanchorHeatmapDistWeight);

  useAnchorHeatmap = reader.GetBoolean("anchor_heatmap", "use_anchor_heatmap",
                                       useAnchorHeatmap);
  anchorHeatmapInfluence = (float)reader.GetReal(
      "anchor_heatmap", "heatmap_influence", anchorHeatmapInfluence);
  // anchorHeatmapDistWeight = (float)reader.GetReal("anchor_heatmap",
  // "heatmap_dist_weight", anchorHeatmapDistWeight);

  // segmentSize = reader.GetInteger("main", "segment_size", segmentSize);

  maxPETClusterLength =
      reader.GetInteger("main", "max_pet_length", maxPETClusterLength);
  longPETClustersEffectPower = (float)reader.GetReal(
      "main", "long_pet_power", longPETClustersEffectPower);
  longPETClustersEffectScale = (float)reader.GetReal(
      "main", "long_pet_scale", longPETClustersEffectScale);

  rewiringEnabled =
      reader.GetBoolean("rewiring", "rewiring_enabled", rewiringEnabled);
  rewiringCertaintyThreshold = reader.GetInteger(
      "rewiring", "rewiring_certainty_threshold", rewiringCertaintyThreshold);

  dataDirectory = reader.Get("data", "data_dir", dataDirectory);
  dataAnchors = reader.Get("data", "anchors", dataAnchors);
  dataPetClusters = reader.Get("data", "clusters", dataPetClusters);
  dataSingletons = reader.Get("data", "singletons", dataSingletons);
  dataSingletonsInter =
      reader.Get("data", "singletons_inter", dataSingletonsInter);
  dataSplitSingletonFilesByChr = reader.GetBoolean(
      "data", "split_singleton_files_by_chr", dataSplitSingletonFilesByChr);
  dataFactors = reader.Get("data", "factors", dataFactors);

  dataCentromeres = reader.Get("data", "centromeres", dataCentromeres);
  dataSegmentsSplit = reader.Get("data", "segment_split", dataSegmentsSplit);

  dataSegmentHeatmap =
      reader.Get("template", "segment_heatmap", dataSegmentHeatmap);

  templateSegment = reader.Get("template", "template_segment", templateSegment);
  templateScale =
      (float)reader.GetReal("template", "template_scale", templateScale);

  distHeatmap = reader.Get("template", "dist_heatmap", distHeatmap);
  distHeatmapScale =
      (float)reader.GetReal("template", "dist_heatmap_scale", distHeatmapScale);

  useDensity = reader.GetBoolean("density", "use_density", useDensity);
  densityMapFile = reader.Get("density", "density_data", densityMapFile);
  densityScale =
      (float)reader.GetReal("density", "density_scale", densityScale);
  densityInfluence =
      (float)reader.GetReal("density", "density_influence", densityInfluence);
  densityWeight =
      (float)reader.GetReal("density", "density_weight", densityWeight);

  useTelomerePositions =
      reader.GetBoolean("telomeres", "use_telomeres", useTelomerePositions);
  if (useTelomerePositions) {
    string p = reader.Get("telomeres", "telomere_1_position", "");
    sscanf(p.c_str(), "%f,%f,%f", &telomere_1_position.x,
           &telomere_1_position.y, &telomere_1_position.z);
    p = reader.Get("telomeres", "telomere_2_position", "");
    sscanf(p.c_str(), "%f,%f,%f", &telomere_2_position.x,
           &telomere_2_position.y, &telomere_2_position.z);
  }

  simulationStepsLevelChr =
      reader.GetInteger("main", "steps_lvl1", simulationStepsLevelChr);
  simulationStepsLevelSegment =
      reader.GetInteger("main", "steps_lvl2", simulationStepsLevelSegment);
  simulationStepsLevelAnchor =
      reader.GetInteger("main", "steps_arcs", simulationStepsLevelAnchor);
  simulationStepsLevelSubanchor =
      reader.GetInteger("main", "steps_smooth", simulationStepsLevelSubanchor);

  noiseCoefficientLevelChr =
      (float)reader.GetReal("main", "noise_lvl1", noiseCoefficientLevelChr);
  noiseCoefficientLevelSegment =
      (float)reader.GetReal("main", "noise_lvl2", noiseCoefficientLevelSegment);
  noiseCoefficientLevelAnchor =
      (float)reader.GetReal("main", "noise_arcs", noiseCoefficientLevelAnchor);
  noiseCoefficientLevelSubanchor = (float)reader.GetReal(
      "main", "noise_smooth", noiseCoefficientLevelSubanchor);

  heatmapInterScaling =
      (float)reader.GetReal("heatmaps", "inter_scaling", heatmapInterScaling);
  heatmapDistanceHeatmapStretching =
      (float)reader.GetReal("heatmaps", "distance_heatmap_stretching",
                            heatmapDistanceHeatmapStretching);

  ibRandomWalkJumps = (float)reader.GetReal("distance", "ib_random_walk_jumps",
                                            ibRandomWalkJumps);
  freqToDistHeatmapScale = (float)reader.GetReal("distance", "freq_dist_scale",
                                                 freqToDistHeatmapScale);
  freqToDistHeatmapPower = (float)reader.GetReal("distance", "freq_dist_power",
                                                 freqToDistHeatmapPower);
  freqToDistHeatmapScaleInter = (float)reader.GetReal(
      "distance", "freq_dist_scale_inter", freqToDistHeatmapScaleInter);
  freqToDistHeatmapPowerInter = (float)reader.GetReal(
      "distance", "freq_dist_power_inter", freqToDistHeatmapPowerInter);
  countToDistA =
      (float)reader.GetReal("distance", "count_dist_a", countToDistA);
  countToDistScale =
      (float)reader.GetReal("distance", "count_dist_scale", countToDistScale);
  countToDistShift =
      (float)reader.GetReal("distance", "count_dist_shift", countToDistShift);
  countToDistBaseLevel = (float)reader.GetReal(
      "distance", "count_dist_base_level", countToDistBaseLevel);
  genomicLengthToDistPower = (float)reader.GetReal(
      "distance", "genomic_dist_power", genomicLengthToDistPower);
  genomicLengthToDistScale = (float)reader.GetReal(
      "distance", "genomic_dist_scale", genomicLengthToDistScale);
  genomicLengthToDistBase = (float)reader.GetReal(
      "distance", "genomic_dist_base", genomicLengthToDistBase);

  springConstantStretch = (float)reader.GetReal("springs", "stretch_constant",
                                                springConstantStretch);
  springConstantSqueeze = (float)reader.GetReal("springs", "squeeze_constant",
                                                springConstantSqueeze);
  springAngularConstant = (float)reader.GetReal("springs", "angular_constant",
                                                springAngularConstant);
  springConstantStretchArcs = (float)reader.GetReal(
      "springs", "stretch_constant_arcs", springConstantStretchArcs);
  springConstantSqueezeArcs = (float)reader.GetReal(
      "springs", "squeeze_constant_arcs", springConstantSqueezeArcs);

  maxTempHeatmap = (float)reader.GetReal("simulation_heatmap",
                                         "max_temp_heatmap", maxTempHeatmap);
  dtTempHeatmap = (float)reader.GetReal("simulation_heatmap",
                                        "delta_temp_heatmap", dtTempHeatmap);
  tempJumpCoefHeatmap = (float)reader.GetReal(
      "simulation_heatmap", "jump_temp_coef_heatmap", tempJumpCoefHeatmap);
  tempJumpScaleHeatmap = (float)reader.GetReal(
      "simulation_heatmap", "jump_temp_scale_heatmap", tempJumpScaleHeatmap);
  MCstopConditionStepsHeatmap =
      reader.GetInteger("simulation_heatmap", "stop_condition_steps_heatmap",
                        MCstopConditionStepsHeatmap);
  MCstopConditionImprovementHeatmap = (float)reader.GetReal(
      "simulation_heatmap", "stop_condition_improvement_threshold_heatmap",
      MCstopConditionImprovementHeatmap);
  MCstopConditionMinSuccessesHeatmap = reader.GetInteger(
      "simulation_heatmap", "stop_condition_successes_threshold_heatmap",
      MCstopConditionMinSuccessesHeatmap);

  maxTempHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "max_temp_heatmap", maxTempHeatmapDensity);
  dtTempHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "delta_temp_heatmap", dtTempHeatmapDensity);
  tempJumpCoefHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "jump_temp_coef_heatmap",
      tempJumpCoefHeatmapDensity);
  tempJumpScaleHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "jump_temp_scale_heatmap",
      tempJumpScaleHeatmapDensity);
  MCstopConditionStepsHeatmapDensity = reader.GetInteger(
      "simulation_heatmap_density", "stop_condition_steps_heatmap",
      MCstopConditionStepsHeatmapDensity);
  MCstopConditionImprovementHeatmapDensity =
      (float)reader.GetReal("simulation_heatmap_density",
                            "stop_condition_improvement_threshold_heatmap",
                            MCstopConditionImprovementHeatmapDensity);
  MCstopConditionMinSuccessesHeatmapDensity =
      reader.GetInteger("simulation_heatmap_density",
                        "stop_condition_successes_threshold_heatmap",
                        MCstopConditionMinSuccessesHeatmapDensity);

  maxTemp = (float)reader.GetReal("simulation_arcs", "max_temp", maxTemp);
  dtTemp = (float)reader.GetReal("simulation_arcs", "delta_temp", dtTemp);
  tempJumpCoef =
      (float)reader.GetReal("simulation_arcs", "jump_temp_coef", tempJumpCoef);
  tempJumpScale = (float)reader.GetReal("simulation_arcs", "jump_temp_scale",
                                        tempJumpScale);
  MCstopConditionSteps = reader.GetInteger(
      "simulation_arcs", "stop_condition_steps", MCstopConditionSteps);
  MCstopConditionImprovement = (float)reader.GetReal(
      "simulation_arcs", "stop_condition_improvement_threshold",
      MCstopConditionImprovement);
  MCstopConditionMinSuccesses =
      reader.GetInteger("simulation_arcs", "stop_condition_successes_threshold",
                        MCstopConditionMinSuccesses);

  weightAngleSmooth = (float)reader.GetReal("simulation_arcs_smooth",
                                            "dist_weight", weightAngleSmooth);
  weightDistSmooth = (float)reader.GetReal("simulation_arcs_smooth",
                                           "angle_weight", weightDistSmooth);
  maxTempSmooth = (float)reader.GetReal("simulation_arcs_smooth", "max_temp",
                                        maxTempSmooth);
  dtTempSmooth = (float)reader.GetReal("simulation_arcs_smooth", "delta_temp",
                                       dtTempSmooth);
  tempJumpCoefSmooth = (float)reader.GetReal(
      "simulation_arcs_smooth", "jump_temp_coef", tempJumpCoefSmooth);
  tempJumpScaleSmooth = (float)reader.GetReal(
      "simulation_arcs_smooth", "jump_temp_scale", tempJumpScaleSmooth);
  MCstopConditionStepsSmooth =
      reader.GetInteger("simulation_arcs_smooth", "stop_condition_steps",
                        MCstopConditionStepsSmooth);
  MCstopConditionImprovementSmooth = (float)reader.GetReal(
      "simulation_arcs_smooth", "stop_condition_improvement_threshold",
      MCstopConditionImprovementSmooth);
  MCstopConditionMinSuccessesSmooth = reader.GetInteger(
      "simulation_arcs_smooth", "stop_condition_successes_threshold",
      MCstopConditionMinSuccessesSmooth);
  return true;
}
