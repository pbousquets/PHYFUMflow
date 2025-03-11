# Phyfum

[![PyPI - Version](https://img.shields.io/pypi/v/phyfum.svg)](https://pypi.org/project/phyfum)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/phyfum.svg)](https://pypi.org/project/phyfum)
[![PHYFUM CI/CD Pipeline](https://github.com/pbousquets/PHYFUMflow/actions/workflows/main.yml/badge.svg)](https://github.com/pbousquets/PHYFUMflow/actions/workflows/main.yml)
![GitBook Badge](https://img.shields.io/badge/Gitbook-docs-orange?style=flat&logo=gitbook&link=https%3A%2F%2Fphyfum.gitbook.io%2Ftutorial)

 
> [!TIP]
>
> Visit our [GitBook](https://phyfum-1.gitbook.io/tutorial/) for a detailed tutorial of Phyfum

--- 
## Quick start

Phyfum allows two different workflows. If you are working with raw data (IDAT files), you can run phyfum in __complete__ mode. In this mode, phyfum will preprocess the files with [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html). If needed and if both tumor and normal samples are available, it will also run a copy number analysis with [rascal](https://github.com/crukci-bioinformatics/rascal) to blacklist fCpGs located within copy-number-altered regions, which do not behave as the model expects.

An example run for this workflow would look like this:

```{bash}
input_dir="epic_array_dir" #path to input directory with idat files and proper folder structure
output_dir="experiment1" #path to output directory
patient_info="${input_dir}/sample_sheet.csv" #path to csv file with experiment design

phyfum run --input ${input_dir}\
 --output ${output_dir}\
 --workdir ${output_dir}\
 --patientinfo ${patient_info}\
 --patient-col patient\
 --age-col age\
 --patient-col patient\
 --sample-col sample\
 --sample-type-col group\
 --stemcells 3-10-3 
```

If you have already pre-processed the data and have the beta values, you can run phyfum in __trees__ mode. The pipeline will simply deploy the XMLcreator tool to format the input data as expected by [our modified version](https://github.com/pbousquets/PHYFUM) of [BEAST](https://beast.community/) and run the inference.

```{bash}
input_dir="beta_dir"
input=${input_dir}/"exampleBeta.csv" #path to input file with beta values
output_dir="onlybetas" #path to output directory
patient_info="${input_dir}/meta.csv" #path to csv file with metadata

phyfum run --input ${input}\
 --output ${output_dir}\
 --workdir ${output_dir}\
 --patientinfo ${patient_info}\
 --patient-col patient\
 --age-col age\
 --patient-col patient\
 --sample-type-col group\
 --stemcells 3-10-3 
```

Phyfum auto-detects what kind of input is provided and selects automatically the optimal workflow.

## Installation

A docker image of [Phyfum](https://hub.docker.com/repository/docker/pbousquets/phyfum/general) is available, and is our __recommended way to use the tool:__

```{bash}
docker pull pbousquets/phyfum
```

The commands above can be ran as:
```{bash}
input_dir="epic_array_dir" #path to input directory with idat files and proper folder structure
output_dir="experiment1" #path to output directory
patient_info="${input_dir}/sample_sheet.csv" #path to csv file with experiment design

docker run --rm -it -v ${input_dir}:${input_dir} -v ${output_dir}:${output_dir}\
pbousquets/phyfum --input ${input_dir}\
 --output ${output_dir}\
 --workdir ${output_dir}\
 --patientinfo ${patient_info}\
 --patient-col patient\
 --age-col age\
 --patient-col patient\
 --sample-col sample\
 --sample-type-col group\
 --stemcells 3-10-3 
```

```{bash}
input_dir="beta_dir"
input=${input_dir}/"exampleBeta.csv" #path to input file with beta values
output_dir="onlybetas" #path to output directory
patient_info="${input_dir}/meta.csv" #path to csv file with metadata

docker run --rm -it -v ${input_dir}:${input_dir} -v ${output_dir}:${output_dir}\
pbousquets/phyfum --input ${input}\
 --output ${output_dir}\
 --workdir ${output_dir}\
 --patientinfo ${patient_info}\
 --patient-col patient\
 --age-col age\
 --patient-col patient\
 --sample-type-col group\
 --stemcells 3-10-3 
```

### Manual installation

Prior to installing phyfum, you'll need to install [our modified version of BEAST](https://github.com/pbousquets/PHYFUM) to enable fCpGs to be analyzed under the framework. Then, make sure you have installed __python3__ and __R (>4.0.0)__ and simply run:

```console
pip install phyfum
```

In order to preprocess IDAT files, we use minfi, conumee and rascal, as well as some tidyverse packages. Missing dependencies will automatically be installed during the first run with Phyfum, so it may take longer than usual to run. You can also install them yourself with:

```{r}
if (!require("pacman")) install.packages("pacman")
p_load(optparse, cli, conumee, minfi, parallel, tibble, tidyr, dplyr, data.table, gtools)
p_load_gh("crukci-bioinformatics/rascal")
```

## Preparing the sample sheet / metadata

Phyfum relies on the Array Sample sheet for the __complete__ workflow and a custom metadata file for the __trees__  workflow. In any case, the file must be a comma-separated file (.csv). 

- __Sample sheet__. When running the `complete` workflow, we recommend passing the array sample sheet. Custom columns can be added to specify parameters that are required by the pipeline (sample age, age_at_diagnosis, etc.). Additionally, if the user wanted to remove any sample from the analysis, the corresponding row in the sample sheet can be filtered out to exclude it from the analysis. 

  The pipeline will try to find how many "normal" or "control" samples exist to use them as controls for the CNV pipeline. You can provide the column name with the argument `--sample-type-col`. If no normals are found, this part of the pipeline will be skipped.

- __Custom metadata__. Sample-wise file providing information about the sample age, patient, age_at_diagnosis, etc. It doesn't require anything special as long as it is in CSV format. In order to identify what the columns are, you can use the arguments `--patient-col`, `--sample-col` and `age-col`, if the column names in your file are different from the defaults.

> Both the custom metadata and the sample sheet are passed through `--patientinfo`.


## License

`phyfum` is distributed under the terms of the [CC-BY-NC-SA](https://spdx.org/licenses/MIT.html) license.

