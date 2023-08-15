Benchmarking Instructions
=========================

This README provides steps to perform a benchmark using our provided scripts and tools.

Steps:
------

1. **Install Necessary Packages**:
   Ensure you have all required packages by running:

   pip install -r requirements.txt


2. **Run Benchmark Script**:
Before running the benchmark script `benchmark_modelling_methods.py`, ensure you have a configuration similar to the one found in the cudaMMC repository under `example_data/stg_gpu_example_config.ini`. 

Use the `set_paths_in_config.py` script to set the appropriate paths in your config file.

**Example runs of the benchmark script**:

- For single runs of modeling:
  ```
  python benchmark_modelling_methods.py -c config_file -o output_with_logs_and_models -m 3dnome -i 10 -d autosomal
  ```

- For ensemble modeling:
  ```
  python benchmark_modelling_methods.py -c config_file -o output_with_logs_and_models_ensemble -m cudaMMC -i 10 -d chr1,chr14,chr21
  ```

Note: Replace `config_file`, `output_with_logs_and_models`, and other placeholders with the appropriate paths/names as per your setup.

In the next step, go through the analysis stored in [cudaMMC_benchmark_analysis.ipynb](cudaMMC_benchmark_analysis.ipynb).

Thank you for using our benchmarking tools!
