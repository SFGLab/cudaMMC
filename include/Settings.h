/*
 * Settings.h
 *
 *  Created on: Aug 4, 2013
 *      Author: psz
 */

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "../src/lib/INIReader.h"
#include "../src/lib/common.h"
#include <stdio.h>
#include <string>

class Settings {

public:
  Settings();

  void init();
  void print(int level = 0);
  bool loadFromINI(std::string ini_path);

  static bool initialized;

  static bool use2D;
  static int debug;
  static bool randomWalk;
  static bool useInputCache;

  static int outputLevel;
  // 5 - milestone info

  static int loopDensity;

  // cuda
  static int milestoneFailsThreshold;
  static int cudaBlocksMultiplier;
  static int cudaThreadsPerBlock;

  static std::string dataDirectory;
  static std::string dataAnchors;
  static std::string dataPetClusters;
  static std::string dataSingletons;
  static std::string dataSingletonsInter;
  static std::string dataFactors;
  static bool dataSplitSingletonFilesByChr;

  static std::string dataCentromeres;
  static std::string dataSegmentsSplit;

  static std::string dataSegmentHeatmap;

  // structural template
  static std::string templateSegment;
  static float templateScale;

  // density
  static bool useDensity;
  static std::string densityMapFile;
  static float densityScale;
  static float densityInfluence;
  static float densityWeight;

  static bool useTelomerePositions;
  static vector3 telomere_1_position;
  static vector3 telomere_2_position;

  static float ibRandomWalkJumps;

  // distances matrix (segment level)
  static std::string distHeatmap;
  static float distHeatmapScale;

  // rewiring
  static bool rewiringEnabled;
  static int rewiringCertaintyThreshold;

  static bool useCTCFMotifOrientation;
  static bool motifsSymmetric;
  static float motifOrientationWeigth;

  static bool useAnchorHeatmap;
  static float anchorHeatmapInfluence;
  static float anchorHeatmapDistWeight;

  static bool useSubanchorHeatmap;
  static float subanchorHeatmapInfluence;
  static float subanchorHeatmapDistWeight;
  static int subanchorEstimateDistancesReplicates;
  static int subanchorEstimateDistancesSteps;
  static float subanchorLoopLoopInteractionsScale;

  static int maxPETClusterLength;
  static float longPETClustersEffectPower;
  static float longPETClustersEffectScale;

  static int segmentSize;

  static int simulationStepsLevelChr;
  static int simulationStepsLevelSegment;
  static int simulationStepsLevelAnchor;
  static int simulationStepsLevelSubanchor;

  static float noiseCoefficientLevelChr;
  static float noiseCoefficientLevelSegment;
  static float noiseCoefficientLevelAnchor;
  static float noiseCoefficientLevelSubanchor;

  // heatmaps
  static float heatmapInterScaling;
  static float
      heatmapDistanceHeatmapStretching; // how many averages the distance can be

  // map frequency (from singleton heatmap) to spatial distance
  static float freqToDistHeatmapScale;
  static float freqToDistHeatmapPower;
  static float freqToDistHeatmapScaleInter;
  static float freqToDistHeatmapPowerInter;

  // map count to spatial distance
  static float countToDistA;
  static float countToDistScale;
  static float countToDistShift;
  static float countToDistBaseLevel;

  // map genomic distance to spatial distance
  static float genomicLengthToDistPower;
  static float genomicLengthToDistScale;
  static float genomicLengthToDistBase;

  // springs
  static float springConstantStretch;
  static float springConstantSqueeze;
  static float springAngularConstant;

  static float springConstantStretchArcs;
  static float springConstantSqueezeArcs;

  // MC simulation
  static float maxTempHeatmap;
  static float dtTempHeatmap;
  static float tempJumpScaleHeatmap;
  static float tempJumpCoefHeatmap;
  static float MCstopConditionImprovementHeatmap;
  static int MCstopConditionMinSuccessesHeatmap;
  static int MCstopConditionStepsHeatmap;

  static float maxTempHeatmapDensity;
  static float dtTempHeatmapDensity;
  static float tempJumpScaleHeatmapDensity;
  static float tempJumpCoefHeatmapDensity;
  static float MCstopConditionImprovementHeatmapDensity;
  static int MCstopConditionMinSuccessesHeatmapDensity;
  static int MCstopConditionStepsHeatmapDensity;

  static float maxTemp;
  static float dtTemp;
  static float tempJumpScale;
  static float tempJumpCoef;
  static float MCstopConditionImprovement;
  static int MCstopConditionMinSuccesses;
  static int MCstopConditionSteps;

  static float weightDistSmooth;
  static float weightAngleSmooth;
  static float maxTempSmooth;
  static float dtTempSmooth;
  static float tempJumpScaleSmooth;
  static float tempJumpCoefSmooth;
  static float MCstopConditionImprovementSmooth;
  static int MCstopConditionMinSuccessesSmooth;
  static int MCstopConditionStepsSmooth;
};

#endif /* SETTINGS_H_ */
