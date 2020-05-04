# ASL Pipeline for the Human Connectome Project
This repository contains the ASL processing pipeline scripts for the Human Connectome Project.

## Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)

## Prerequisites
The HCP list some prerequisites for their pipelines: https://github.com/Washington-University/HCPpipelines/wiki/Installation-and-Usage-Instructions.

The prerequisites specific to this pipeline, along with links to their installation pages, are listed below:
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
- [Workbench](https://www.humanconnectome.org/software/get-connectome-workbench)

## Installation
To install hcpasl in the current Python environment:

```
git clone https://github.com/ibme-qubic/hcp-asl.git
cd hcp-asl
python -m pip install --user .
```

## Usage
Once installed, the pipeline may be used as a command-line script as follows:

```
hcp_asl ${SubjectDirectory} ${mt_scaling_factors}
```