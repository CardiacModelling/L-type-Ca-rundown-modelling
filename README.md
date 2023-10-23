# Modelling the effect of experimental conditions that influence rundown of L-type calcium current

This repository contains all data, protocols, codes, and figure generating files used in this study.

## Table of Contents
- [run.py](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/blob/main/run.py) runs all figure generating files in the [figures](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/tree/main/figures) and [supporting](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/tree/main/supporting) directory
- [helpers.py](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/blob/main/helpers.py) provides supporting functions used by the figure generating files
- [parameters.py](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/blob/main/parameters.py) has the values of the parameters used 
- The [resources](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/tree/main/resources) directory has model and other files from previous experiments used as input here
- The [supporting](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/tree/main/supporting) directory has output generating files used for [ExtendedData.ipynb](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/tree/main/ExtendedData.ipynb) 
- The [figures](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/tree/main/figures) directory has the files that generate all figures in the paper and their output
- Supporting material for the paper is provided at [ExtendedData.ipynb](https://github.com/CardiacModelling/L-type-Ca-rundown-modelling/tree/main/ExtendedData.ipynb)


## Software
- **Python** The code requires Python (version 3.3+) and the Python dependencies `myokit` and `pints`
  - The Python packages used in this project can be installed by running `pip install -r requirements.txt` in an environment with `Python 3.10.11`

## Hardware
  - Some of the scripts require high processing, and the results here were generated on 20 CPU cores of a x86_64 Intel(R) Xeon(R) (3.00GHz) system
