# cudaMMC - Parallel Multiscale Monte Carlo Approach to 3D Structure Modelling of the Human Genome

cudaMMC is a GPU-enhanced version of the 3D-GNOME software, which utilizes the Multiscale Monte Carlo approach for 3D chromatin structure modeling.

For those interested in the previous version, it's hosted on the [Bitbucket repository](https://bitbucket.org/3dome/3dgnome).

The methodology behind this tool is described in detail in the following publication:
**Szałaj, P., Tang, Z., Michalski, P., Pietal, M. J., Luo, O. J., Sadowski, M., ... & Plewczynski, D. (2016).** *An integrated 3-dimensional genome modeling engine for data-driven simulation of spatial genome organization.* _Genome research_, 26(12), 1697-1709.

cudaMMC has been integrated into the 3D-GNOME 3.0 web server. We invite you to explore it on the [3D-GNOME 3.0 web server](https://3dgnome.mini.pw.edu.pl/).
cudaMMC has been integrated into the 3D-GNOME 3.0 web server. We invite you to explore it on the [3D-GNOME 3.0 web server](https://3dgnome.mini.pw.edu.pl/). For a concise overview of our methodology, please refer to the help section on our web server: [3D-GNOME 3.0 web server - help](https://3dgnome.mini.pw.edu.pl/help/). In our article, we detail the usage and introduce new features for analyzing 3D chromatin models. Specifically, we discuss the impact of Structural Variation-driven changes on the 3D structure, which result in altered distances between enhancers and gene promoters. For more details, please see: [3D-GNOME 3.0: a three-dimensional genome modelling engine for analysing changes of promoter-enhancer contacts in the human genome](https://academic.oup.com/nar/article/51/W1/W5/7157515).



## System and hardware requirements
* UNIX based Operating System
* C++ 17 compliant compiler (gcc/clang)
* NVIDIA CUDA Toolkit 11.0 or higher with compatible NVIDIA driver
* NVIDIA GPU with Pascal architecture or newer (compute capability >= 6.0)
* CMake 3.13 or higher

You can install cudaMMC using docker or compile it without a container:

## Docker Instructions

Start by building the image with the following command:

`sudo docker build -t cudammc --progress=plain . 2>&1 | tee docker_build.log`

and then start it with

`sudo docker run -v ~/Projects:/Projects --rm -it --gpus all cudammc`

## Compilation
First, generate the makefile with CMake by executing the command below from the root directory of the package with your desired `path_to_build_folder` and `cuda_arch` (>= 60) values:
```
cmake -B ./<path_to_build_folder> -S ./ -DCUDA_ARCH="<cuda_arch>"
```
Example for GPU with compute capability >= 8.0:
```
cmake -B ./build -S ./ -DCUDA_ARCH="80"
```
If you wish to compile a multi-architecture binary, separate desired architectures with a semicolon eg. `-DCUDA_ARCH="60;70;80"`

Then, compile the program by executing the following command:
```
make --directory ./build
```

Finally, you can verify if the program has compiled properly by executing the cudaMMC binary:
```
./build/cudaMMC
```
You should see the usage message printed in the terminal.

Example execution command:
```
./build/cudaMMC -s config.ini -c chr14 -o ~/chiapet/chr14/ -n chr14_test
```

Program usage:
```
cudaMMC -s <path> -n <label>
    -a    action to perform (available values: 'create')
    -s    path to settings file
    -n    label to be used in names of output files
    -c    chromosomes (or region) to be included in the simulation (eg 'chr14:1:2500000', default: 'chr1-chr22,chrX')
    -o    output directory (with trailing '/'; must exists; default: './')
Debugging options
    -u    maximal number of chromosomes to reconstruct
    -l    length of chromosome fragment that will be reconstructed (in bp)
Notes
    PET clusters should have filenames of the form 'GM12878_CTCF_Rep1.withLinker.clusters.intra.PET<N>+.merged_anchors.simple.txt', where <M> corresponds to the flag -m
```

### Input and output files
Input of the algorithm consists of several files. 
First of all, there are files corresponding to different types of ChIA-PET data: anchors, PET clusters and singletons. As singleton files tend to be large it is possible to provide inter- and intrachromosomal singletons files separately (using **singletons_inter** and **singletons** options, respectively). This allows the program to skip the interchromosomal files reading if they are not needed (e.g. when a single chromosome is reconstructed). 
When the subanchor heatmaps are to be generated it is beneficial to create the intrachromosomal singletons files for every chromosome separately (one should denote this by setting the flag **split_singleton_files_by_chr**), and to use these per-chromosome files rather than the bigger, aggregated files. 
You can use script [extract_singletons.sh](extract_singletons.sh) for extracting singletons from files (run extract_singletons.sh -h).
You can also extract it using awk procedures like in following commands:
```
#!bash
mkdir chr;
for i in `seq 1 22` X; 
   do 
      ( cat data.txt | awk '{if ($1 == "chr'$i'") print $0}' > chr/data.txt.chr$i ) 
   done
```
, where data.txt is the original data file. The resulting files are created in the *chr* subdirectory of the data directory, and they have the same name as the original file but with a chromosome id as a suffix.



Additionally to the ChIA-PET files mentioned one can provide a BED file with centromeres' locations and a file with definition of split of the chromosomes into segments. 

Paths of the files are provided via the setting file in the **[data]** section. There are following options available:

* data_dir - path to the data directory (anchors, clusters, singletons and singletons_inter are all relative to this directory)
* anchors - name of the anchors file
* clusters - names of the clusters files (comma separated). Each file should contain clusters   
* singletons - names of the singleton files (comma separated)
* singletons_inter - names of interchromosomal singleton files (comma separated)
* split_singleton_files_by_chr - if *yes* then the chromosome-splitted intrachromosmal files are used 
* factors - names of protein factors used in the experiment (e.g. "CTCF,RNAPII")
* segment_split - path to the BED file with segments split definition
* centromeres - path to the BED file with centromeres' locations

Following is the format description of the corresponding input files.

### Anchor file
This file contain info about the anchors. For every anchor the chromosome and genomic position of start and end is provided. Optionally the CTCF motif orientation may be provided as 'L' (leftward) or 'R' (rightward). 'N' may be used if the orientation is unknown. File format is as follows:

```
#!text
chromosome_id anchor_start anchor_end [orientation]
```
A sample file (with motif orientation provided) may look like this:

```
#!text
chr1	838908	841157	R
chr1	911279	912011	L
chr1	918286	922335	R
chr1	967152	969271	R
chr1	997732	1000167	L
```

### Cluster files
These files contain info about the PET clusters. Each file corresponds to a single transcription factor used in the experiment. Every line of a file describes a single PET cluster and consists of the chromosome, start and end position for the anchors overlapping with the PET cluster. Additionally, a PET count is provided. File format is as follows:

```
#!text
chromosome_id_1 anchor_1_start anchor_1_end chromosome_id_2 anchor_2_start anchor_2_end PET_count
```
A sample file may look like this:

```
#!text
chr1	838908	841157	chr1	911279	912011	11
chr1	838908	841157	chr1	997732	1000167	7
chr1	918286	922335	chr1	997732	1000167	75
chr1	967152	969271	chr1	997732	1000167	52
chr1	967152	969271	chr1	1306262	1308124	7

```

### Singleton files
Singleton files possess the same format as the clusters files, but the genomic coordinates correspond to the paired-end read mapping locations rather than the overlapping anchors:
```
#!text
chromosome_id_1 read_1_start read_1_end chromosome_id_2 read_2_start read_2_end PET_count
```

```
#!text
chr9	77361535	77361615	chr9	78624184	78624253	1
chr2	74756819	74756889	chr2	201675832	201675878	1
chr12	51818625	51818729	chr12	70877421	70877572	1
chr4	169237541	169237655	chr4	181760704	181760799	1
chr6	129055926	129055968	chr6	129170407	129170485	1
```

## User manual

# cudaMMC Settings #
There is a number of settings available for the cudaMMC simulation. We will now shortly describe the most important ones with an intuitive meaning, the full list of options will follow.

We have added three additional parameters under [cuda] section that can be
specified in the config.ini configuration file. These are:
- num threads - this parameter specifies the number of CUDA threads to be launched per block
- blocks multiplier - stands for the number of CUDA blocks multiplied by the number of beads in the system
- milestone fails - this parameter determines after how many non improved milestone scores, the algorithm will stop

Rest parameters are the same as in 3D-GNOME version:

- freq_dist_power - this exponent describes the relation between singletons interaction frequency on a segment level and the physical distances between beads, significantly different values will yield different shapes. Different values were used in the literature. A value of -1.0 correspond to a simple inverse relation.
- freq_dist_scale - scaling of the segment level structure 
- genomic_dist_scale - responsible for the size of the chromatin loops
- use_motif_orientation - whether or not to use the CTCF motif orientation
- use_subanchor_heatmap - whether or not to use the subanchor heatmap to refine chromatin loops modeling

Below a description of all the settings available is provided.

##### [main]
* output_level - set the level of output messages (range=0..10)
* random_walk - if true then create a random walk structure on the segment level
* loop_density - number of subanchor beads that will be placed between the consecutive anchor beads
* use_2D - if true then the simulation is restricted to 2 dimensions
* max_pet_length - maximal length of the PET clusters used on the subanchor level (in bp)
* long_pet_power, long_pet_scale - describe how long PET cluster impact the segment heatmaps (scale*C^power, where C is PET count)

* steps_lvl1 (lvl2, arcs, smooth) - number of simulations on the corresponding levels  (lvl1 - chromosome level, lvl2 - segments, arcs - anchor, smooth - subanchor)
* noise_lvl1 (lvl2, arcs, smooth) - amount of noise used to create the initial structures on the corresponding levels

##### [data]
* data_dir - path of the directory with data files
* anchors - name of the anchor file
* clusters - names of the cluster files (comma separated)
* factors - names of the factors used in
* singletons - names of files with intrachromosomal singletons
* split_singleton_files_by_chr - flag denoting whether the files in 'singletons' were splitted by chromosome  
* singletons_inter - names of the files with interchromosomal singletons
* segment_split - absolute path to a BED file with the segment split info
* centromeres - absolute path to a BED file with the centromere locations

##### [distance]
* genomic_dist_power, genomic_dist_scale, genomic_dist_base - describe relationship between genomic distance and the physical distance between subanchor beads (3D dist = base+scale*d^power, where d is genomic distance in kb)
* freq_dist_scale, freq_dist_power - describe relationship between interaction frequency and physical distance, used to generate segment level expected distances matrix (3D dist = scale*F^power, where F is interaction frequency)
* freq_dist_scale_inter, freq_dist_power_inter - the same as freq_dist_scale, but used for the chromosome level. Allows to use different relation for segment and chromosome level.
* count_dist_a, count_dist_scale, count_dist_shift, count_dist_base_level - describe relationship between PET count and the physical distance between subanchor beads (3D dist = base+scale/e^[a*(shift+C)], where C is PET count)

##### [template]
* template_segment - path to the file the structural template file (if any)
* template_scale - scale for the structural template
* dist_heatmap - path to the file the structural template file (if any)
* dist_heatmap_scale - scale for the structural template

##### [motif_orientation]
* use_motif_orientation - whether to consider CTCF morif orientation during the simulation or not
* weight - the weight assigned to the motif orientation energy term

##### [anchor_heatmap]
* use_anchor_heatmap - whether or not to construct the anchor heatmap to refine the anchor beads placement
* heatmap_influence - the influence of the pairwise anchor distances matrix

##### [subanchor_heatmap]
* use_subanchor_heatmap - whether or not to construct the subanchor heatmap to refine the subanchor beads placement (i.e. chromatin loops shape and relative positions)
* estimate_distances_steps - number of structures created to obtain the expected distance matrix
* estimate_distances_replicates - number of simulation steps for every structure 
* heatmap_influence - the influence of the pairwise anchor distances matrix
* heatmap_dist_weight - the weight assigned to the expected distance matrix energy term

##### [heatmaps]
* inter_scaling - scaling factor applied to the interchromosomal contacts (segment level). This can be used  
* distance_heatmap_stretching - used to calculate the cap value for large 3D distances (cap = average * stretching) 

##### [springs]
* stretch_constant_arcs - weight assigned to the flexibility energy term on the anchor level, when the distance is higher than expected
* squeeze_constant_arcs - as above, but for distances smaller than expected
* stretch_constant - as stretch_constant_arcs, but on the subanchor level
* squeeze_constant - as squeeze_constant_arcs, but on the subanchor level
* angular_constant - weight assigned to the bending energy term (subanchor level)

##### [simulation_heatmap]
* max_temp_heatmap - initial temperature for the simulated annealing
* delta_temp_heatmap - temperature reduction between iterations
* jump_temp_scale_heatmap, jump_temp_coef_heatmap - parameters to scale the probability of accepting move with higher energy
* stop_condition_improvement_threshold_heatmap - improvement ratio that is required for the algorithm to stop (must be higher)
* stop_condition_successes_threshold_heatmap - if the number of accepted moves during a milestone is higher than this value than the algortihm continues
stop_condition_steps_heatmap - number of steps for each milestone

##### [simulation_arcs]
(the same as for [simulation_heatmap])
##### [simulation_arcs_smooth]
(the same as for [simulation_heatmap])


### Usage
Example usage of the software is described in the [README.md](example_data%2FREADME.md).

### Benchmark

We provide scripts to compare time performance and modeling quality between the previous version of the software, 3D-GNOME, and cudaMMC. Check out the [benchmark](benchmark) directory for more details.