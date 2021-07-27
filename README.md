# ASL Pipeline for the Human Connectome Project
This repository contains the ASL processing pipeline scripts for the Human Connectome Project.

## Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)

## Prerequisites
The HCP list some prerequisites for their pipelines: https://github.com/Washington-University/HCPpipelines/wiki/Installation-and-Usage-Instructions.

The prerequisites specific to this pipeline, along with links to their installation pages, are listed below:
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) (version >= 6.0.4)
- [Workbench](https://www.humanconnectome.org/software/get-connectome-workbench) (version >= 1.5.0)
- [HCP Pipelines](https://github.com/Washington-University/HCPpipelines)
- [HCP's gradunwarp](https://github.com/Washington-University/gradunwarp)

FSL version 6.0.4 is required to use new features in fabber, added for this pipeline.

Workbench >= v1.5.0 is required. The `-weighted` option in `wb_command`'s `-volume-to-surface-mapping` introduced in version 1.5.0 is used in projections to the surface.

The HCP Pipelines must be installed and the environment variable `HCPPIPEDIR` set in order for the (Sub-)Cortical LUTs to be used in SE-based bias correction.

## Installation
It is advised that the pipeline is installed in a conda environment, for example following the steps below:

```
conda create -n hcpasl pip cython numpy python=3.7
conda activate hcpasl
git clone https://github.com/physimals/hcp-asl.git
cd hcp-asl
python -m pip install --user .
```

## Usage
Once installed, the pipeline may be run as a command-line script as follows:

```
hcp_asl --studydir ${StudyDir} --subid ${Subjectid} --outdir ${IntermediateOutputDir}--mtname ${ScalingFactors} -g ${GradientCoeffs} -s ${T1} --sbrain ${T1Brain} --mbpcasl ${mbpcasl} --fmap_ap ${SEFM_AP} --fmap_pa ${SEFM_PA} --cores ${n_cores} --interpolation ${InterpType} --pvcorr --wmparc ${wmparc.nii.gz} --ribbon ${ribbon.nii.gz} --territories_atlas ${Atlas} --territories_labels ${AtlasLabels} --wbdir ${wbDir} -v
```

The filepaths passed to the script may be relative or absolute. A more detailed explanation of the arguments can be found by running:

```
hcp_asl --help
```