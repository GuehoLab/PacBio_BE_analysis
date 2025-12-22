# PacBio_BE_analysis

This pipeline provides tools for demultiplexing PacBio HiFi reads, analyzing mutations with CRISPRlungo, and visualizing substitution patterns.

## 1. Installation

### requirements

- python3
- regex
- seaborn
- regex
- CRISPRLungo

First, create the analysis environment using the provided `environment.yml` file. Then, install the package in editable mode.

```bash
# Create and activate the mamba environment
mamba env create -f environment.yml -n pacbio_be_analysis_env
conda activate pacbio_be_analysis_env

# Install the package
pip install -e .

```
After activate environment, CRISPRLungo also should be installed (https://github.com/pinellolab/CRISPRlungo).

## 2. Usage

The analysis consists of three main steps: demultiplexing, mutation analysis, and data organization.

### Step 1: Demultiplexing

Categorize PacBio reads based on your index table.

```bash
PacBe_demultiplex {hifi.fastq.gz} {paf file} {index_table file} {output_dir}
```

### Step 2: CRISPRlungo Analysis

Run the core mutation analysis. This command uses all mutation types and a wide window size for comprehensive detection.

```bash
CRISPRlungo {reference fasta} {demultiplexed fastq} \
    --control {control file} \
    {lungo_output dir} \
    {target} \
    --using_all_mutations \
    --window 10000
```

### Step 3: Organize & Visualize

Generate organized data summaries and final visualization plots.

```bash
PacBe_organize_data {lungo_output dir} {target sequence}
```
After the pipeline finishes, you will find the results in the {lungo_output} directory. The primary visualization file is:
- subsitution.png: A heatmap representing base substitution patterns.
  - C > T transitions are indicated in Blue.
  - G > A transitions are indicated in Red.

## 4. Example

To verify the installation and see the pipeline in action, you can run the provided example script. The example dataset is designed to complete within 1 minute.

```bash
cd example
./run_example.sh
```



