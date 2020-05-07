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

If available, the user can supply a gradient coefficients file for use in gradient 
distortion correction as follows:

```
hcp_asl ${SubjectDirectory} ${mt_scaling_factors} -g ${grad_coeffs}
hcp_asl ${SubjectDirectory} ${mt_scaling_factors} --grads ${grad_coeffs}
```

The distortion correction script can also be called directly:

```
hcp_asl_distcorr ${StudyDirectory} ${SubjectNumber} (-g ${grad_coeffs})
```

If the gradient coefficients are not supplied, the script will perform the other 
motion correction and registration steps without including gradient distortion 
correction.