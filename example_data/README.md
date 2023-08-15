Example Modelling with cudaMMC
==============================

Follow these steps to successfully run the example modeling using cudaMMC:

1) Set Paths in Config File
---------------------------
Firstly, you need to ensure the paths in the configuration file are correctly set. You can do this by running:

python set_paths_in_config.py -d /absolute/path/to/the/example_data -c /absolute/path/to/the/example_data/stg_gpu_example_config.ini

**Note**: This script will modify the `stg_gpu_example_config.ini` file. It is advised to keep a backup of the original configuration file.

2) Run Modelling using cudaMMC
------------------------------
With the paths correctly set, you can now initiate the modeling using cudaMMC:

cudaMMC -s /path/to/the/example_data/stg_gpu_example_config.ini -c chr14:35138000-36160000 -o ./sim_out/

This will model the segment from `chr14:35138000-36160000` and store the output in the `./sim_out/` directory.

3) Convert Output to mmCIF Format
---------------------------------
The native output from cudaMMC is in the `.hcm` format. For visualization purposes, it is beneficial to convert this into the `.mmCIF` format which can then be viewed using software like UCFC Chimera.

**Prerequisites**:
- Ensure you have the `numpy` package installed. If not, you can install it using pip:

pip install numpy

To perform the conversion, run:
python convert_hcm2cif.py -i /sim_out/ -s 5000 -c /path/to/cudaMMC

You can now use [UCFC Chimera](https://www.cgl.ucsf.edu/chimera/download.html) to visualize the `.mmCIF` files generated.

---

We hope this guide assists you in running the example modeling smoothly. If you face any issues, please refer to the documentation or contact support.
