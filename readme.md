# IsoNet

## Introduction
`IsoNet` is an R package for isotopologue similarity networking to discover unknown reactions or metabolites.

The docker image (https://hub.docker.com/r/zhulab/isonet-r) contains entire environment for running `IsoNet`. For convenience and taking fully use of `IsoNet`, users can pull it and run `IsoNet` just as following.

## What is IsoNet-r

`IsoNet-r` is a Docker environment to build isotopologue similarity network with `IsoNet` R package. It is based on the [`r-base`](https://hub.docker.com/_/r-base/) docker.

## Pulling image

Users can pull the IsoNet-r image with the following script

```bash
docker pull zhulab/IsoNet-r
```

## Data preparation

The data folder should contain isotopologue pattern table (.csv), the identification table of labeled metabolite (.csv) and MS2 sprectra (.msp). Demo files could be downloaded from [https://www.zhulab.cn/usr/uploads/misc/isonet/IsoNet_demo_data.zip].


## R script preparation
To run the data processing, an R script named [run.R](extra/run.R) should be placed in the data folder.
Here we provide an example. Users only need to change tracer_table_file,  ms2_file,  id_table, and sample_names . Other parameters are recommended parameters.
  
```R
library(IsoNet)
wd <- '.'
setwd(wd)

constructMIDNet(tracer_table_file = "./Result/isotopologue_pattern_table.csv",
                ms2_file = "ms2_data.msp",
                id_table = "./Result/ID_labeled_metabolites.csv",
                sample_names = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6"),
                dir_path = '.',
                mid_cutoff = 0.7,
                mid_fc = 20,
                mid_isoDegree = 0.1,
                mid_min_motifLen = 0.5,
                ms2_score_cutoff = 0.5,
                mid_max_motif = 1000L,
                ignore_max = TRUE,
                scoring_approach = 'gnps',
                mass_diff_freq_cut_off = 4L)
```

**Import parameters**

- `tracer_table_file`: the file name of isotopologue pattern table which contains the labeled fraction of each isotopologue.
- `ms2_file`: the file name of ms2 file.
- `id_table`: the identification table name of labeled metabolites.
- `sample_names`: The sample names in tracer_table_file that be used to calculate isotopologue pattern similarity.
- `mid_cutoff`: the isotopogue pattern similarity score cutoff
- `ms2_score_cutoff`：the MSMS similarity score cutoff
- `mid_fc`: the labeled fraction radio between the highest and the second highest isotopologues to class the type II and type III.
- `mid_isoDegree`: The cutoff for the minimum labeling fraction of highest isotopologue in the motif.
- `mid_min_motifLen`：In Type III, The minimum carbon number of generated motif. For example, setting it to 0.5 indicates that the carbon number of the motif should be at least 50% longer than that of the metabolite with more carbons.
- `mid_max_motif`：In Type III, the maximum number of generated motifs.
- `ignore_max`：In Type I, whether to consider the position of the highest isotopologue. Setting it to TRUE indicates that it is not considered, and the calculation is performed directly.
- `scoring_approach`：The algorithm used for scoring the MSMS simlairty.
- `mass_diff_freq_cut_off`: The threshold for delta masses used for edge annotation based on its frequency in the KEGG database.


## Run data processing work with IsoNet-r image

- go to your data folder (e.g., data)

```base
cd data
```

- run docker using following code (*User should be permitted to run docker service*)

```bash
# MUST keep the code exactly as it is!
docker run -it --rm -v "$PWD":/data -u $(id -u ${USER}):$(id -g ${USER}) zhulab/IsoNet-r Rscript run.R
```

- wait till data processing work done

- Explaining `docker run` arguments
  
  - `-v "$PWD":/home/${USER}`: mapping current directory as home directory in docker container
  
  - `-u $(id -u ${USER}):$(id -g ${USER})`: using current user to run the container
  
  - `Rscript ~/run.R`: run run.R in container home directory with `Rscript`  command

## The result 

The main results are listed following:
  - labeled metabolite pairs named "Labeled metabolite pairs.csv". It can be imported into Cytoscape to built a network.


# License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a>
  
  This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
