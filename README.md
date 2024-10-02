# **Installation**

## **1. `conda`**

Make sure you have your conda installed, conda env created and activated

- ### **From pre-built package**
    The pre-built package is located at `conda-build` directory. To install from this source:
    - Verify that current directory is gapzilla root directory (where `conda-recipe` folder and `README.md` are located)
        ```bash
        >>> tree -L 1
        .
        ├── conda-build
        ├── conda-recipe
        ├── docs
        ├── documentation.html
        ├── gapzilla
        ├── pyproject.toml
        ├── README.md
        ├── requirements.txt
        └── tsv_template.tsv

        ```
    - Run the following command:
        ```bash
        conda install -c ./conda-build gapzilla
        ```


- ### **Using `conda-build`**
    It is possible to build conda package locally from source. In order to do that:
    - `cd` into the root gapzilla folder.
    - Run the following command:
        ```bash
        conda build conda-recipe
        ```
        or, if you would like pre-compiled package to be stored at a certain path:
        ```bash
        conda build conda-recipe --output-folder path/to/output/folder
        ```


    ---
    **NOTE**

    Building process might take a while.

    ---
    - After building process is complete, your output package will be saved locally somewhere in `~/anaconda/conda-bld/linux-64/` if you used first build method or at path you've provided
    - To actually install package if no path was provided run:
        ```bash
        conda install --use-local gapzilla
        ```
        or, if you provided custom path on previous step:
        ```bash
        conda install -c path/to/output/folder gapzilla
        ```
## **2. `pip`**

This type of installation can be very useful for those who are used to working with [python venv](https://docs.python.org/3/library/venv.html). This package is not available in PyPI, so one have to make manual installation. 

- First make sure you have your venv created and activated. 
- Then `cd` into the root gapzilla folder (where `pyproject.toml` and `README.md` files are located).
- If you plan on editing source files of this repo, run:

    ```bash
    pip3 install -e .
    ```
    Otherwise:
    ```bash
    pip3 install .
    ```
    The difference is that `-e` option [provides a symbolic link](https://stackoverflow.com/a/59667164/19559362) to this package in pip env, so that changes user makes here will immediately take effect on the package.

    ---
    **NOTE**

    It is possible to install package via `pip` into conda env. However, this might be [not a very good idea](https://www.anaconda.com/blog/using-pip-in-a-conda-environment) if you have a complicated conda env.

    ---

## **3. `sys.path`**
This method is the least preferred one and is usually used as a temporary solution. To make package visible in python one can append it to `sys.path` directly in python session by executing:
```python
import sys
sys.path.append("path/to/gapzilla")
```

# **Usage**

## **CLI**

Command-line arguments are available via help message:
```bash
>>> gapzilla -h
usage: gapzilla [-h] [--path_to_gbk [N ...]] [--output_file_name OUTPUT_FILE_NAME] [--suffix_file_name SUFFIX_FILE_NAME] [--path_to_output_folder PATH_TO_OUTPUT_FOLDER] [--min_gap_length MIN_GAP_LENGTH]
                [--max_gap_length MAX_GAP_LENGTH] [--min_flanks_length MIN_FLANKS_LENGTH] [--max_flanks_length MAX_FLANKS_LENGTH] [--mfe_threshold_hpa MFE_THRESHOLD_HPA] [--mfe_threshold_hpt MFE_THRESHOLD_HPT]
                [--hairpin_similarity_thres HAIRPIN_SIMILARITY_THRES] [--border_shift BORDER_SHIFT] [--threads THREADS] [--verbosity VERBOSITY] [--avoid_plotting [N ...]] [--custom_config N]
                [--table_input TABLE_INPUT]

Find RNA hairpin-flanked gaps in GBK file for potential insert sites

options:
-h, --help            show this help message and exit
--path_to_gbk [N ...], -p [N ...]
                        Path to the gbk input file (default: None)
--output_file_name OUTPUT_FILE_NAME, -f OUTPUT_FILE_NAME
                        Output file name. If not specified, will save input file name and add suffix (default: )
--suffix_file_name SUFFIX_FILE_NAME, -s SUFFIX_FILE_NAME
                        Specify filename suffix to add (default: output)
--path_to_output_folder PATH_TO_OUTPUT_FOLDER, -o PATH_TO_OUTPUT_FOLDER
                        Path to the gbk output folder (default: )
--min_gap_length MIN_GAP_LENGTH, -min_gap MIN_GAP_LENGTH
                        Minimum gap length (default: 0)
--max_gap_length MAX_GAP_LENGTH, -max_gap MAX_GAP_LENGTH
                        Maximum gap length (default: 150)
--min_flanks_length MIN_FLANKS_LENGTH, -min_flanks MIN_FLANKS_LENGTH
                        Minimum flanking CDS length (default: 700)
--max_flanks_length MAX_FLANKS_LENGTH, -max_flanks MAX_FLANKS_LENGTH
                        Maximum flanking CDS length (default: 10000000)
--mfe_threshold_hpa MFE_THRESHOLD_HPA, -hpa MFE_THRESHOLD_HPA
                        MFE (minimal free energy) threshold for RNA hairpins (default: -7)
--mfe_threshold_hpt MFE_THRESHOLD_HPT, -hpt MFE_THRESHOLD_HPT
                        MFE threshold to filter most stable RNA hairpins (default: -15)
--hairpin_similarity_thres HAIRPIN_SIMILARITY_THRES, -similarity HAIRPIN_SIMILARITY_THRES
                        Threshold for filtering similar hairpins. Lower thres -> more hairpins dropped (default: 0.9)
--border_shift BORDER_SHIFT, -bs BORDER_SHIFT
                        Distance to the right and to the left from gap to search for haiirpins (default: 75)
--threads THREADS, -t THREADS
                        Maximum number of threads to use. Takes all available cores, if not specified (default: <all available cores>)
--verbosity VERBOSITY, -v VERBOSITY
                        Specify verbosity of output (0 - silent, 1 - max) (default: 1)
--avoid_plotting [N ...], -ap [N ...]
                        Specify instances that will not be plotted. Example to drop hpa: -ap all_hairpins_f all_hairpins_r. Full list of plottable instances: gaps, top_hairpins_f, top_hairpins_r, all_hairpins_f,
                        all_hairpins_r, insertion_sites (default: [])
--custom_config N, -cc N
                        Path to custom yaml config file that will override defaults and any other commandline args (default: None)
--table_input TABLE_INPUT, -tsv TABLE_INPUT
                        Path to custom tsv table that contains parameters for multiple samples (default: None)
```

The only required argument is `--path_to_gbk` or `-p`. 

Arguments can also be passed via custom config `.yaml` file or `.tsv` table.

### **Creating custom config file**
It is possible to create custom config file. Say, we would like to pass another `suffix_name`, MFE threshold for hpt (`mfe_threshold_hpt`) and not to plot `hpa` hairpins at all. The provided custom config then should look like this:
```yaml
suffix_file_name: "res"
mfe_threshold_hpt: -20
avoid_plotting: "all_hairpins_f all_hairpins_r"
```
It's not necessary to provide full config, only parameters to be changed should be passed. Full list of default parameters for config can be found at `config.py` file.

### **Creating custom table input**
Also `.tsv` table can be used as input. It should have the following format:
- column names should be *exact* names of parameters to be set
- row has all the parameters for corresponding table
- not neccessarily all cells should be filled

**Example:**

| path_to_gbk       | suffix_file_name | min_gap_length |
|-------------------|------------------|----------------|
| path/to/file1.gbk | out              | 50             |
| path/to/file2.gbk |                  | 75             |
| path/to/file3.gbk | res              |                |

Full template with all possible columns is available at `tsv_template.tsv`

## **Python package**
This package can also be imported as a regular python package. Some functions from this package may be useful in other projects too.
```python
from gapzilla import process_gbk, find_top_hairpins
```

# **Output**

Script outputs a `.gbk` file with additional annotations that include:
- `gap`: intergenic gap, that meets required length and flankong regions lengths
- `hpt`: most stable RNA hairpins detected
- `hpa`: all other hairpins detected that do not match stability criteria
- `ins{n}`: possible insertion sites, where `n` is a score from 1 to 4. Score 4 can only be achieved when site is surrounded by 4 stable hairpins (`hpt`s): two on the left (forward and reverse) and two on the right.

Output can be opened by any tool that allows view of `.gbk` files like UGENE or SnapGene.

# **Requirements**
```
python>=3.10
biopython
pyyaml
pandas
tqdm
viennarna
```






