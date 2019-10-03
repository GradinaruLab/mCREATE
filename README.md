# mCREATE

Example code for generating analyses and plots for the multiplexed CREATE paper
by Kumar et al. (2019). The code is divided into three parts: alignment,
variant analysis, and clustering. The scripts given in the Juypter notebooks
are "demo" versions with example data or input files in the relative
"example_data" sub-folder.

## Alignment

### Prerequisites

Alignment uses Python 3 and Jupyter notebooks, and depends on the following
Python packages:
```
jupyter
pepars
protfarm
```

Pepars and Protfarm are available at their respective GitHub repositories:
- Pepars: [https://github.com/GradinaruLab/pepars](https://github.com/GradinaruLab/pepars)
- Protfarm: [https://github.com/GradinaruLab/protfarm](https://github.com/GradinaruLab/protfarm)

### Scripts

#### align_and_export.ipynb

This scripts initializes a Protfarm data workspace, downloads example FASTQ files, aligns them against a template, and calculates and exports the count and enrichment data. For more details on Protfarm workspaces and enrichment, see the [Protfarm](https://github.com/GradinaruLab/protfarm) repository.

## Variant Analysis

Variant analysis contains a set of scripts that operate on count and enrichment
Excel sheets, as generated by aligning and exporting data via Protfarm. The
time required to run these scripts are in the order of 2-5 mins depending on
the length of the dataset.

### Prerequisites

Variant analysis uses Python 3 and Jupyter notebooks, and requires the following
Python packages:
```
jupyter
pepars
```

### Scripts

#### mCREATE_supplementary_analysis.ipynb
This is a compilation of codes performing multiple functions. 
Each cell is a separate function and the relevant functions are mentioned on top of every cell. 
The users are required to choose the functions/cells they want and run them separately.

#### Primer_generator_for_11mer_clones_587-590AAV9.ipynb
This script generates 2 overlapping forward and reverse primers for AAV variants where 11 codons are modified from the parent 
(as described in the RavindraKumar et al 2019 manuscript, Methods section (unpublished)). These primers are used for cloning AAV variants.

#### Reverse_primer_generator_11mer-spikein_synthetic_lib.ipynb
This script generates 1 reverse primer and 1 reverse primer with alternate mammalian codon for synthesizing 11 codon modified 
Spike-in library oligopool. This is then sent to Twist Bioscience for oligopool synthesis.

#### Reverse_primer_generator_7mer_synthetic_lib.ipynb
This script generates 1 reverse primer and 1 reverse primer with alternate mammalian codon for synthesizing 7 codon modified 
588-89 library oligopool. This is then sent to Twist Bioscience for oligopool synthesis.

#### heatmap_588-89+library.ipynb (previous heatmap generation method - see heatmap_7mer.ipynb for updated version)
This script generates a heatmap of amino acid distribution in the diversified region. 
The raw data is subjected to codon frequency normalization and the standard score is plotted as a heatmap 
(with options to generate heatmaps before and after normalization).

#### heatmap_7mer.ipynb
This script generates a heatmap of differences in amino acid distribution between two random sets of sequences (generated by different simulated templates).


## Clustering

Clustering analyses use sequence and count data as output from Protfarm, or from manual analysis, to investigate families of variants.

### Prerequisites

Clustering analysis requires MATLAB (tested on MATLAB R2017b (9.3.0.7135779), 64-bit (maci64)).

### Scripts

#### hammDist.m
This is a Matlab function that generates a .txt list of paired sequences with their reverse hamming distance, 
i.e. number of positions with identical AA. Output may be visualized in Cytoscape software.

#### clustering_example.m
This is a script that runs hammDist (above) using the included example data
