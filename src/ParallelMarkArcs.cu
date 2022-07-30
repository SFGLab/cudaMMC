#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>

#include <time.h>

#include <InteractionArcs.h>

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

__device__ int anchor_length(int start, int end) {
  if (start == 0 && end == 0)
    return 0;
  return end - start + 1;
}

__device__ bool anchor_contains(int start, int end, int pos) {
  return pos >= start && pos <= end;
}

__global__ void kernel(int *__restrict__ anchors_starts,
                       int *__restrict__ anchors_ends,
                       int *__restrict__ raw_arcs_starts,
                       int *__restrict__ raw_arcs_ends, int *output_starts,
                       int *output_ends, const int anchors_count,
                       const int arcs_count) {
  int threadIndex = blockDim.x * blockIdx.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  if (threadIndex >= arcs_count)
    return;

  for (int i = threadIndex; i < arcs_count; i += stride) {
    int start = -1, end = -1;

    for (int j = 0; j < anchors_count; ++j) {
      if (anchor_length(anchors_starts[j], anchors_ends[j]) > 1) {
        if (anchor_contains(anchors_starts[j], anchors_ends[j],
                            raw_arcs_starts[i]))
          output_starts[i] = j;
        if (anchor_contains(anchors_starts[j], anchors_ends[j],
                            raw_arcs_ends[i]))
          output_ends[i] = j;
        if (start != -1 && end != -1)
          return;
      }
    }
  }
}

// we have anchors[], a sorted list of anchors, and raw_arcs[], a list of arcs
// (using genomic positions) we want to fill arc[] so that it contains a list of
// arcs with values referring to indices of 'anchors', and not genomic positions
void InteractionArcs::parallelMarkArcs(bool ignore_missing) {
  int threads = 1024;
  int blocks = 256;

  // for every chromosome:
  // go through all raw arcs, find the corresponding anchor, create new arc
  // (anchor based) and add it to list

  int cnt = 0;
  int last_start = -1;
  std::unordered_map<int, std::vector<InteractionArc>> tmp_arcs;

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

    std::vector<Anchor> current_anchors = anchors[chr];
    std::vector<InteractionArc> current_raw_arcs = raw_arcs[chr];

    thrust::host_vector<int> h_anchors_start(current_anchors.size());
    thrust::host_vector<int> h_anchors_end(current_anchors.size());
    thrust::host_vector<int> h_raw_arcs_start(cnt);
    thrust::host_vector<int> h_raw_arcs_end(cnt);
    thrust::host_vector<int> h_outputs_start;
    thrust::host_vector<int> h_outputs_end;

    thrust::device_vector<int> d_anchors_start;
    thrust::device_vector<int> d_anchors_end;
    thrust::device_vector<int> d_raw_arcs_start;
    thrust::device_vector<int> d_raw_arcs_end;
    thrust::device_vector<int> d_outputs_start(cnt);
    thrust::device_vector<int> d_outputs_end(cnt);

    thrust::fill(d_outputs_start.begin(), d_outputs_start.end(), -1);
    thrust::fill(d_outputs_end.begin(), d_outputs_end.end(), -1);

    for (int i = 0; i < anchors_cnt[chr]; ++i) {
      h_anchors_start[i] = current_anchors[i].start;
      h_anchors_end[i] = current_anchors[i].end;
    }

    for (int i = 0; i < cnt; ++i) {
      h_raw_arcs_start[i] = current_raw_arcs[i].start;
      h_raw_arcs_end[i] = current_raw_arcs[i].end;
    }

    d_anchors_start = h_anchors_start;
    d_anchors_end = h_anchors_end;
    d_raw_arcs_start = h_raw_arcs_start;
    d_raw_arcs_end = h_raw_arcs_end;

    kernel<<<blocks, threads>>>(
        thrust::raw_pointer_cast(d_anchors_start.data()),
        thrust::raw_pointer_cast(d_anchors_end.data()),
        thrust::raw_pointer_cast(d_raw_arcs_start.data()),
        thrust::raw_pointer_cast(d_raw_arcs_end.data()),
        thrust::raw_pointer_cast(d_outputs_start.data()),
        thrust::raw_pointer_cast(d_outputs_end.data()), anchors_cnt[chr], cnt);

    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    h_outputs_start = d_outputs_start;
    h_outputs_end = d_outputs_end;

    // we need i==cnt to process the remaining arcs
    for (int i = 0; i <= cnt; ++i) { // for every arc

      int st = -1, end = -1;

      if (i < cnt) {
        st = h_outputs_start[i];
        end = h_outputs_end[i];

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
