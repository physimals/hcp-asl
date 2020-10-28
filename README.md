# ASL Pipeline for the Human Connectome Project
This repository contains the ASL processing pipeline scripts for the Human Connectome Project.

## Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)

## Prerequisites
The HCP list some prerequisites for their pipelines: https://github.com/Washington-University/HCPpipelines/wiki/Installation-and-Usage-Instructions.

The prerequisites specific to this pipeline, along with links to their installation pages, are listed below:
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) (version >= 6.0.4)
- [Workbench](https://www.humanconnectome.org/software/get-connectome-workbench)

FSL version 6.0.4 is required to use new features in fabber, added for this pipeline. Alternatively, you can try to update these files separately and provide the fabberdir argument - more support about this will be provided in future.

## Installation
A dependency of the hcpasl pipeline requires igl which is installable via conda. It is advised that the pipeline is installed in a conda environment which has igl installed, for example following the steps below:

```
git clone https://github.com/ibme-qubic/hcp-asl.git
cd hcp-asl
conda install igl -c conda-forge
python -m pip install --user .
```

## Usage
Once installed, the pipeline may be run as a command-line script as follows:

```
hcp_asl --studydir ${StudyDir} --subid ${Subjectid} --mtname ${ScalingFactors} -g ${GradientCoeffs} -s ${T1} --sbrain ${T1Brain} --surfacedir ${fs32kDirectory} --mbpcasl ${mbpcasl} --fmap_ap ${SEFM_AP} --fmap_pa ${SEFM_PA}
```

If --surfacedir is provided, the script will look for the surfaces required for ribbon-constrained projection within this directory. Alternatively, the user can instead specify the filenames of each of these 6 surfaces individually.

The filepaths passed to the script may be relative or absolute. A more detailed explanation of the arguments can be found by running:

```
hcp_asl --help
```