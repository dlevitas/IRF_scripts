# IRF_scripts

This repository contains code for Indiana University's [Imaging Research Facilty (IRF)](https://irf.indiana.edu/), including automatic DICOM to BIDS conversion, and setup for BIDS apps (MRIQC, fMRIPrep, QSIprep, and xcpEngine).  

This repository is idiosyncratic to the IRF, which houses a Siemens Prisma scanner and stores data on IU's [Carbonate HPC](https://kb.iu.edu/d/aolp). Therefore, this repository is primarily meant for IRF users, as the scripts only work for Siemens Prisma imaging data on the Carbonate HPC architecture. Furthermore, Carbonate uses slurm as it's management & job scheduler system. 

These scripts are meant to be run on IU's [Research Desktop (RED)](https://kb.iu.edu/d/apum) or Carbonate. Running them elsewhere will require significant modifications, which is not recommend. The repository can be downloaded to users' RED or Carbonate accounts with the following command: *git clone https://github.com/dlevitas/IRF_scripts* 

A brief description of each script:

1). *IRF_BIDS_converter_ReproIn.py* - converts DICOMS from Siemens scanner that adhere to the [ReproIn setup](https://github.com/ReproNim/reproin) for easier, more robust BIDS conversion. As of spring 2021, the IRF **highly** recommends setting up new studies that follow ReproIn. 

2). *IRF_BIDS_converter_non-ReproIn.py* - converts DICOMS from Siemens scanner that do not adhere to ReproIn, to BIDS format. This script is not as robust, as it attempts to account for a wide range of protocol setups.

3). *fmriprep.sh* - provides setup for [fMRIPrep](https://readthedocs.org/projects/fmriprep/) anatomical & functional pre-processing and submits slurm job(s) on IU's HPC.

4). *mriqc.sh* - provides setup for [MRIQC](https://mriqc.readthedocs.io/en/stable/) quality assurance and submits slurm job(s) on IU's HPC.

5). *qsiprep.sh* - provides setup for [QSIprep](https://qsiprep.readthedocs.io/en/latest/) DWI and submits slurm job(s) on IU's HPC.

6). *xcpEngine.sh* - provides setup for [xcpEngine](https://xcpengine.readthedocs.io/) post-processing of neuroimaging data.

