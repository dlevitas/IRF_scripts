#!/bin/bash

#Notes:
# 1). This script runs the BIDS app "MRIQC" on specified subject(s)
# 2). Users don't have write access to this file, so copy it somewhere else 
# 		to make the changes to the mriqc command itself if desired.
# 3). The mriqc command itself is found towards the bottom of this script; 
# 		however, you may wish to adjust or add/remove certain options that 
# 		better fit the purposes of your study
# 4). Refer to https://mriqc.readthedocs.io/en/stable/running.html for 
# 		guidance on which options you wish to use
# 5). If you are having issues:
# 		feel free to reach out to me (dlevitas@iu.edu), 
# 		or post your questions/issues to https://neurostars.org, 
# 		or https://github.com/poldracklab/fmriprep/issues,
# 		or isainnis@iu.edu


# User inputs:
bids_root_dir=/N/slate/dlevitas/BIDS_tutorial
subjects=01
nodes=1 # Number of nodes to use on HPC job submission(s)
cores=2 # Number of threads to use on HPC job submission(s)
mem_mb=15000 # Memory to allocate (in mb)
walltime=2:00:00 # Walltime for job
email=dlevitas@iu.edu # email for sending updates on the job submission(s)
name=mriqc # Name for job submissions
mriqc_version= # Can be left blank. If so, will default to version 0.15.1
analysis_level=participant # participant or group. Specifies which level of analysis you wish to run
overwrite=no # yes or no. If yes, will remove files from specific analysis level and restart


## Begin:
echo ""
date +"Starting mriqc_slate_v2.sh"
echo ""

current_dir=$PWD

#Check user inputs
echo ""
if [[ -z "$bids_root_dir" || -z "$subjects" || -z "$nodes" || -z "$cores" || -z "$mem_mb" || -z "$walltime" || -z "$email" || -z "$name" || -z "$analysis_level" || -z "$overwrite" ]]; then
	echo "Make sure user inputs variables are defined. Exiting now"
	return
fi

# Ensure derivatives directory exists
if [ ! -d $bids_root_dir/derivatives ]; then
	mkdir $bids_root_dir/derivatives
fi


# Convert mem_mb to gb
mem_gb=`echo $((mem_mb / 1000))`

# Check mriqc version
if [ -z "$mriqc_version" ]; then
	mriqc_version=0.15.1
fi

if [ $analysis_level == participant ]; then

	# Loop through each subject
	for s in ${subjects[*]}
	do

		#overwrite
		if [ $overwrite == yes ]; then
			rm -rf $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/*
		fi

		# Add Singularity to PATH
		if [ -z `command -v singularity` ]; then
			module load singularity/3.6.4
		fi

		# Add Freesurfer to PATH
		if [ -z `command -v freesurfer` ]; then
			module load freesurfer/6.0.0 #not latest version, but don't use freesurfer/7.1.0
		fi

		# Generate MRIQC singularity image if not present
		if [ ! -f /N/dcwan/projects/irf/containers/mriqc-${mriqc_version}.simg ]; then
			echo ""
			echo "Building singularity image for MRIQC-${mriqc_version}"
			echo ""

			export SINGULARITY_TMPDIR="/N/slate/$(whoami)"
			export SINGULARITY_CACHEDIR="/N/slate/$(whoami)"
			singularity build /N/dcwan/projects/irf/containers/mriqc-${mriqc_version}.simg docker://poldracklab/mriqc:${mriqc_version}
		fi

		# Make mriqc directory and participant directory (if not already present) in derivatives folder
		if [ ! -d $bids_root_dir/derivatives/mriqc-${mriqc_version} ]; then
			mkdir $bids_root_dir/derivatives/mriqc-${mriqc_version}
		fi

		if [ ! -d $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s} ]; then
			mkdir $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}
		fi

		# Create plugin.yml file, which helps mriqc with memory allocation
		# If plugin.yml file already exists, delete it and create new
		if [ -f $bids_root_dir/mriqc_plugin_${s}.yml ]; then
			rm -rf $bids_root_dir/mriqc_plugin_${s}.yml
		fi
	    # Create new .yml
		touch $bids_root_dir/mriqc_plugin_${s}.yml
		echo "plugin: LegacyMultiProc" >> $bids_root_dir/mriqc_plugin_${s}.yml
		echo "plugin_args: {maxtasksperchild: $nodes, memory_gb: $mem_gb, n_procs: $cores, raise_insufficient: false}" >> $bids_root_dir/mriqc_plugin_${s}.yml

		# Make folder in user's slate directory for storing slurm error and text reports
		if [ ! -d /N/slate/$(whoami)/slurm_output ]; then
			mkdir /N/slate/$(whoami)/slurm_output
			chmod -R 777 /N/slate/$(whoami)/slurm_output
		fi

		# Run MRIQC
		echo ""
		echo "Running mriqc-${mriqc_version} on sub-${s}"
		echo ""

		mriqc_command="unset PYTHONPATH; singularity run -B ${bids_root_dir}:${bids_root_dir} /N/dcwan/projects/irf/containers/mriqc-${mriqc_version}.simg \
		$bids_root_dir $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s} \
		participant 													\
		--participant_label $s 											\
		--n_proc $cores													\
		--hmc-fsl														\
		--correct-slice-timing											\
		--mem_gb $mem_gb 												\
		--float32 														\
		--ants-nthreads 2                                               \
		--no-sub												        \
		--use-plugin $bids_root_dir/mriqc_plugin_${s}.yml 				\
		-w $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}"

		touch $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/mriqc_${s}.sh

		echo "#!/bin/bash" > $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/mriqc_${s}.sh
		echo "" >> $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/mriqc_${s}.sh
		echo $mriqc_command >> $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/mriqc_${s}.sh

		cd $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}

		# Submit slurm job(s)
		sbatch -p general \
		-J ${name}-${s} \
		-o /N/slate/$(whoami)/slurm_output/slurm_${name}-${s}.txt \
		-e /N/slate/$(whoami)/slurm_output/slurm_${name}-${s}.txt \
		--mail-type=ALL \
		--mail-user=$email \
		--nodes=$nodes \
		--cpus-per-task=$cores \
		--ntasks-per-node=1 \
		--time=$walltime \
		--mem=$mem_mb \
		$bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/mriqc_${s}.sh

		cd $current_dir
	done

else
	# overwite group-level files
	if [ $overwrite == yes ]; then
		rm -rf $bids_root_dir/derivatives/mriqc-${mriqc_version}/group*
	fi

	unset PYTHONPATH; singularity run -B ${bids_root_dir}:${bids_root_dir} /N/dcwan/projects/irf/containers/mriqc-${mriqc_version}.simg \
	$bids_root_dir $bids_root_dir/derivatives/mriqc-${mriqc_version} \
	group

	# # Remove temporary files (only do this once all data has gone through participant level; needed for group level reports)
	# rm -rf $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/workflow_enumerator
	# rm -rf $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/logs
	# rm -rf $bids_root_dir/derivatives/mriqc-${mriqc_version}/sub-${s}/sub-${s}*
fi
