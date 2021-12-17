#include <cuda.h>
#include <cuda_fp16.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>

#include <time.h>

#include "../include/LooperSolver.h"

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

struct gpu_settings {
  float dtTempHeatmap;
  float MCstopConditionImprovementHeatmap;
  bool use2D;
  float tempJumpCoefHeatmap;
  float tempJumpScaleHeatmap;
  int milestoneFailsThreshold;
};

struct __align__(8) half3 {
  __half x;
  __half y;
  __half z;
};

// __device__ void mutex_lock(unsigned int *mutex) {
//     unsigned int ns = 10;
//     while (atomicCAS(mutex, 0, 1) == 1) {
//         __nanosleep(ns);
//         if (ns < 256) {
//         ns *= 2;
//         }
//     }
// }

// __device__ void mutex_unlock(unsigned int *mutex) {
//     atomicExch(mutex, 0);
// }

__device__ int random(int range, curandState *state) {
  return curand(state) % range;
}

__device__ __half random(const float &range, bool negative,
                         curandState *state) {
  if (negative)
    return __float2half((2.0f * curand_uniform(state) - 1.0f) * range);
  return __float2half(range * curand_uniform(state));
}

__device__ bool withChance(const float &chance, curandState *state) {
  return curand_uniform(state) < chance;
}

__device__ void randomVector(half3 &vector, const float &max_size, bool &in2D,
                             curandState *state) {
  vector.x = random(max_size, true, state);
  vector.y = random(max_size, true, state);
  vector.z = in2D ? __float2half(0.0f) : random(max_size, true, state);
}

__device__ void addToVector(half3 &destination, const half3 &value) {
  destination.x = __hadd(destination.x, value.x);
  destination.y = __hadd(destination.y, value.y);
  destination.z = __hadd(destination.z, value.z);
}

__device__ void subtractValueFromVector(half3 &destination,
                                        const half3 &value) {
  destination.x = __hsub(destination.x, value.x);
  destination.y = __hsub(destination.y, value.y);
  destination.z = __hsub(destination.z, value.z);
}

__device__ void subtractVectors(half3 &destination, const half3 &value1,
                                const half3 &value2) {
  destination.x = __hsub(value1.x, value2.x);
  destination.y = __hsub(value1.y, value2.y);
  destination.z = __hsub(value1.z, value2.z);
}

__device__ float magnitude(const half3 &vector) {
  float x = __half2float(vector.x), y = __half2float(vector.y),
        z = __half2float(vector.z);
  return sqrtf(x * x + y * y + z * z);
}

__device__ void
getChromosomeHeatmapBoundary(const int p, int &start, int &end,
                             int *gpu_heatmap_chromosome_boundaries,
                             const int &chromosomeBoundariesSize) {
  if (chromosomeBoundariesSize == 0)
    return;
  if (p < 0 ||
      p > gpu_heatmap_chromosome_boundaries[chromosomeBoundariesSize - 1])
    return;

  for (size_t i = 0; i + 1 < chromosomeBoundariesSize; ++i) {
    if (gpu_heatmap_chromosome_boundaries[i] <= p &&
        p < gpu_heatmap_chromosome_boundaries[i + 1]) {
      start = gpu_heatmap_chromosome_boundaries[i];
      end = gpu_heatmap_chromosome_boundaries[i + 1] - 1;
      return;
    }
  }
}

__device__ float calcScoreHeatmapSingleActiveRegion(
    const int moved, half3 *__restrict__ clusters_positions,
    int *__restrict__ gpu_heatmap_chromosome_boundaries,
    float *__restrict__ heatmap_dist, const int &heatmapSize,
    const int &heatmapDiagonalSize, const int &activeRegionSize,
    const int &chromosomeBoundariesSize, const half3 &curr_vector,
    int &warpIdx) {
  double err = 0.0;
  float helper;
  int st = 0, end = activeRegionSize - 1;
  if (heatmapSize != activeRegionSize) {
    printf("heatmap sizes mismatch, dist size=%d, active region=%d",
           heatmapSize, activeRegionSize);
    return 0.0;
  }
  half3 temp_one;

  if (chromosomeBoundariesSize > 0) {
    getChromosomeHeatmapBoundary(moved, st, end,
                                 gpu_heatmap_chromosome_boundaries,
                                 chromosomeBoundariesSize);
  }

  for (int i = st; i <= end; ++i) {
    if (abs(i - moved) >= heatmapDiagonalSize) {
      if (i == moved || (helper = heatmap_dist[i * heatmapSize + moved]) < 1e-3)
        continue; // ignore 0 values
      // no warp divergence as all threads are from the same warp so the warpIdx
      // will evaluate the same
      subtractVectors(temp_one, *(clusters_positions + i),
                      (moved == warpIdx) ? curr_vector
                                         : *(clusters_positions + moved));
      helper = magnitude(temp_one) / helper - 1;
      err += helper * helper;
    }
  }
  return err;
}

__device__ float calcScoreHeatmapActiveRegion(
    const int moved, half3 *__restrict__ clusters_positions,
    int *__restrict__ gpu_heatmap_chromosome_boundaries,
    float *__restrict__ heatmap_dist, const int &heatmapSize,
    const int &heatmapDiagonalSize, const int &activeRegionSize,
    const int &chromosomeBoundariesSize, half3 &curr_vector, int &warpIdx) {
  float err = 0.0f;
  if (moved == -1) {
    for (int i = 0; i < activeRegionSize; ++i) {
      err += calcScoreHeatmapSingleActiveRegion(
          i, clusters_positions, gpu_heatmap_chromosome_boundaries,
          heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize,
          chromosomeBoundariesSize, curr_vector, warpIdx);
    }
    return err;
  }
  return calcScoreHeatmapSingleActiveRegion(
      moved, clusters_positions, gpu_heatmap_chromosome_boundaries,
      heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize,
      chromosomeBoundariesSize, curr_vector, warpIdx);
}

__global__ void setupKernel(curandState *state, time_t seed) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  curand_init(seed, tid, 0, &state[tid]);
}

__global__ void MonteCarloHeatmapKernel(
    curandState *__restrict__ state, bool *__restrict__ clusters_fixed,
    half3 *__restrict__ clusters_positions,
    int *__restrict__ heatmap_chromosome_boundaries, float T, const float score,
    gpu_settings settings, float *__restrict__ heatmap_dist, float step_size,
    const int numberOfClusters, const int activeRegionSize,
    const int heatmapSize, const int heatmapDiagonalSize,
    const int chromosomeBoundariesSize, bool *__restrict__ isDone) {
  int threadIndex = blockDim.x * blockIdx.x + threadIdx.x;
  int warpIdx = (threadIndex / warpSize) % activeRegionSize;
  int iterations = 0;
  int improvementMisses = 0;
  half3 curr_vector;
  half3 displacement;
  float score_curr = score;
  float score_prev = score_curr;
  float milestoneScore = score_curr;
  // __shared__ unsigned int mutex;

#if __CUDA_ARCH__ >= 800
  int i_winner;
  int i_score;
#else
  float f_winner;
#endif

  curandState localState = state[threadIndex];

  while (true) {
    curr_vector = clusters_positions[warpIdx];

#define N 512
#pragma unroll
    for (int i = 0; i < N; ++i) {
      if (clusters_fixed[warpIdx])
        return; // check if bead is fixed (for example telomere)

      randomVector(displacement, step_size, settings.use2D, &localState);
      addToVector(curr_vector, displacement);

      score_curr = calcScoreHeatmapSingleActiveRegion(
          warpIdx, clusters_positions, heatmap_chromosome_boundaries,
          heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize,
          chromosomeBoundariesSize, curr_vector, warpIdx);

      if ((score_curr <= score_prev) ||
          (T > 0.0f &&
           withChance(settings.tempJumpScaleHeatmap *
                          expf(-settings.tempJumpCoefHeatmap *
                               (score_curr * (1 / score_prev)) * (1 / T)),
                      &localState))) {
        score_prev = score_curr;
        continue;
      }
      score_curr = score_prev;
      subtractValueFromVector(curr_vector, displacement);
    }

    T *= settings.dtTempHeatmap;
    step_size *= 0.95;
    iterations += N;

    // find the best move
    __syncwarp();
#define FULL_MASK 0xffffffff
#if __CUDA_ARCH__ >= 800
    i_score = (int)score_curr;
    i_winner = __reduce_min_sync(FULL_MASK, i_score);
    if (i_winner == i_score) {
#else
    f_winner = score_curr;
#pragma unroll 5
    for (int offset = 16; offset > 0; offset /= 2)
      f_winner = fminf(f_winner, __shfl_down_sync(FULL_MASK, f_winner, offset));
    // propagate and check who's the lucky one
    __syncwarp();
    __shfl_sync(FULL_MASK, f_winner, 0);
    if (f_winner == score_curr) {
#endif
      // TODO doubling of the same work here -> memoization or specialized
      // function for this point
      // mutex_lock(&mutex);

      float curr = calcScoreHeatmapSingleActiveRegion(
          warpIdx, clusters_positions, heatmap_chromosome_boundaries,
          heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize,
          chromosomeBoundariesSize, curr_vector, warpIdx);

      float global = calcScoreHeatmapSingleActiveRegion(
          warpIdx, clusters_positions, heatmap_chromosome_boundaries,
          heatmap_dist, heatmapSize, heatmapDiagonalSize, activeRegionSize,
          chromosomeBoundariesSize, *(clusters_positions + warpIdx), warpIdx);

      if (curr < global) {
        clusters_positions[warpIdx] = curr_vector;
      }
      // mutex_unlock(&mutex);
    }
    __threadfence();

#if __CUDA_ARCH__ >= 800
    score_curr = i_winner;
#else
    score_curr = f_winner;
#endif
    score_prev = score_curr;

    // check if we should stop
    if (threadIndex == 0) {
      score_curr = calcScoreHeatmapActiveRegion(
          -1, clusters_positions, heatmap_chromosome_boundaries, heatmap_dist,
          heatmapSize, heatmapDiagonalSize, activeRegionSize,
          chromosomeBoundariesSize, *(clusters_positions + warpIdx), warpIdx);

      printf("Milestone score: %f \n", score_curr);

      if (score_curr >
          settings.MCstopConditionImprovementHeatmap * milestoneScore)
        ++improvementMisses;

      if (improvementMisses >= settings.milestoneFailsThreshold ||
          score_curr < 1e-04)
        *isDone = true;

      milestoneScore = score_curr;
    }

    if (*isDone)
      break;
  }
  state[threadIndex] = localState;
}

float LooperSolver::ParallelMonteCarloHeatmap(float step_size) {
  int threads = Settings::cudaThreadsPerBlock;
  int blocks = Settings::cudaBlocksMultiplier * active_region.size();

  double T = Settings::maxTempHeatmap; // set current temperature to max
  double score_curr;
  double score_density;
  string chr;
  int size = active_region.size();
  int heatmapSize = heatmap_dist.size;

  if (size <= 1)
    return 0.0f; // there is nothing to do

  bool *d_isDone;

  cudaMalloc((void **)&d_isDone, sizeof(bool));
  cudaMemset(d_isDone, 0, sizeof(bool));

  struct gpu_settings settings;
  settings.tempJumpScaleHeatmap = Settings::tempJumpScaleHeatmap;
  settings.tempJumpCoefHeatmap = Settings::tempJumpCoefHeatmap;
  settings.use2D = Settings::use2D;
  settings.dtTempHeatmap = Settings::dtTempHeatmap;
  settings.MCstopConditionImprovementHeatmap =
      Settings::MCstopConditionImprovementHeatmap;
  settings.milestoneFailsThreshold = Settings::milestoneFailsThreshold;

  thrust::host_vector<bool> h_clusters_fixed(active_region.size());
  thrust::host_vector<half3> h_clusters_positions_half(active_region.size());

  for (int i = 0; i < active_region.size(); ++i) {
    h_clusters_fixed[i] = clusters[i].is_fixed;
  }

  thrust::device_vector<int> d_heatmap_chromosome_boundaries(
      heatmap_chromosome_boundaries);
  thrust::device_vector<bool> d_clusters_fixed = h_clusters_fixed;
  thrust::device_vector<half3> d_clusters_positions(active_region.size());
  thrust::device_vector<float> d_heatmap_dist(heatmapSize * heatmapSize);

  for (int i = 0; i < heatmapSize; ++i) {
    thrust::copy(heatmap_dist.v[i], heatmap_dist.v[i] + heatmapSize,
                 d_heatmap_dist.begin() + heatmapSize * i);
  }

  for (int i = 0; i < active_region.size(); ++i) {
    h_clusters_positions_half[i].x =
        __float2half(clusters[active_region[i]].pos.x);
    h_clusters_positions_half[i].y =
        __float2half(clusters[active_region[i]].pos.y);
    h_clusters_positions_half[i].z =
        __float2half(clusters[active_region[i]].pos.z);
  }

  d_clusters_positions = h_clusters_positions_half;
  curandState *devStates;
  gpuErrchk(
      cudaMalloc((void **)&devStates, threads * blocks * sizeof(curandState)));

  // cuRand initialization
  setupKernel<<<blocks, threads>>>(devStates, time(NULL));
  gpuErrchk(cudaDeviceSynchronize());

  // calc initial score
  score_curr = calcScoreHeatmapActiveRegion();
  score_density = calcScoreDensity();
  output(2, "initial score: %lf (density=%lf)\n", score_curr, score_density);

  MonteCarloHeatmapKernel<<<blocks, threads,
                            active_region.size() * sizeof(float)>>>(
      devStates, thrust::raw_pointer_cast(d_clusters_fixed.data()),
      thrust::raw_pointer_cast(d_clusters_positions.data()),
      thrust::raw_pointer_cast(d_heatmap_chromosome_boundaries.data()), T,
      score_curr, settings, thrust::raw_pointer_cast(d_heatmap_dist.data()),
      0.75f * step_size, clusters.size(), active_region.size(), heatmapSize,
      heatmap_dist.diagonal_size, heatmap_chromosome_boundaries.size(),
      d_isDone);

  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaDeviceSynchronize());

  h_clusters_positions_half = d_clusters_positions;

  for (int i = 0; i < active_region.size(); ++i) {
    clusters[active_region[i]].pos.x =
        __half2float(h_clusters_positions_half[i].x);
    clusters[active_region[i]].pos.y =
        __half2float(h_clusters_positions_half[i].y);
    clusters[active_region[i]].pos.z =
        __half2float(h_clusters_positions_half[i].z);
  }

  score_curr = calcScoreHeatmapActiveRegion();

  printf("===========================================================\n");
  printf("FINAL SCORE IS %f \n", score_curr);

  cudaFree(d_isDone);
  cudaFree(devStates);

  return score_curr;
}
