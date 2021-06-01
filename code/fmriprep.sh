#!/bin/bash

#Notes:
# 1). This script runs the BIDS app "fMRIPrep" on specified subject(s)
# 2). Users don't have write access to this file, so copy it somewhere else 
# 		to make the changes to the fmriprep command itself if desired.
# 3). The fmriprep command itself is found towards the bottom of this script; 
# 		however, you will want to adjust or add/remove certain options that 
# 		better fit the purposes of your study
# 4). Refer to https://fmriprep.readthedocs.io/en/stable/running.html for 
# 		guidance on which options you wish to use
# 5). Be aware that fmriprep only performs preprocessing steps that the neuroimaging community has deemed absolutely necessary
#		As an example, spatial smoothing is not performed with fmriprep (unless --use-aroma is specified), 
# 		so spatial smoothing would need to be done using a separate program afterwards, if you so choose
# 6). fmriprep generates a ton of temporary files in the work directory (-w) /N/slate/$(whoami)/fmriprep-${fmriprep_version}_temp 
#		If fmriprep successfully finishes, these temp files should be removed. Otherwise, re-run fmriprep and it will pick up from where the crash occurred 
# 7). If you are having issues:
# 		feel free to reach out to me (dlevitas@iu.edu), 
# 		or post your questions/issues to https://neurostars.org, 
# 		or https://github.com/nipreps/fmriprep/issues,
# 		or isainnis@iu.edu


## User inputs:
bids_root_dir=/N/dcwan/projects/irf/BIDS/BIDS_tutorial
subjects=01
nodes=1 # Number of nodes to use on HPC job submission(s)
cores=4 # Number of processors to use on HPC job submission(s). 16 is the max
mem_mb=25000 # Memory to allocate
walltime=10:00:00 # Walltime for job
email=dlevitas@iu.edu # email for sending updates on the job submission(s)
name=fmriprep # Name for job submissions
fmriprep_version= # Can be left blank. If so, will default to version 20.2.0, which is the long term support (LTS)
output_dir= # can be left blank, in which case the fmriprep output will be stored in $bids_root_dir/derivatives/fmriprep
overwrite=yes # yes or no. If yes, will remove files and restart


# Begin:
echo ""
date +"Starting fmriprep_slate_v2.sh"
echo ""

current_dir=$PWD

# Ensure derivatives directory exists
if [ ! -d $bids_root_dir/derivatives ]; then
	mkdir $bids_root_dir/derivatives
fi

#Check user inputs
echo ""
if [[ -z "$bids_root_dir" || -z "$subjects" || -z "$nodes" || -z "$cores" || -z "$mem_mb" || -z "$walltime" || -z "$email" || -z "$name" || -z "$overwrite" ]]; then
	echo "Make sure user inputs variables are defined. Exiting now"
	return
fi

if [ -z "$fmriprep_version" ]; then
	fmriprep_version=20.2.0
fi

if [ -z "$output_dir" ]; then
	output_dir=$bids_root_dir/derivatives
fi

if [ ! -d $output_dir/fmriprep ]; then
	mkdir $output_dir/fmriprep
fi

# Convert mem_mb to gb
mem_gb=`echo $((mem_mb / 1000))`

# Loop through each subject
for s in ${subjects[*]}
do

	# Remove old folders if overwrite == yes
	if [ $overwrite == yes ]; then
		rm -rf $bids_root_dir/derivatives/fmriprep_${s}*.o*
		rm -rf $bids_root_dir/derivatives/fmriprep/sub-${s}*
		rm -rf /$bids_root_dir/derivatives/fmriprep/dataset_description.json
	fi

	# Add .bidsignore file to tell fmriprep which directories/files to skip during validation
	if [ ! -f $bids_root_dir/.bidsignore ]; then
		touch $bids_root_dir/.bidsignore
	fi 

	# Add Freesurfer to PATH
	if [ -z `command -v freesurfer` ]; then
		module unload perl
		module load freesurfer/6.0.0 # Don't use version 7.1.0 b/c seems to throw error
	fi

	export FS_LICENSE=$FREESURFER_HOME/license.txt

	# Add text to README (required for BIDS validator)
	if [ $(wc -l <$bids_root_dir/README) -eq 0 ]; then
		echo "Running fmriprep-${fmriprep_version}" >> $bids_root_dir/README
	fi
	
	# Reduce memory consumption for fMRIPrep
	reduced_mem=`echo $((mem_mb-5000))` #remove a little less than what was required in job (buffer space)

	# Create plugin.yml file, which helps fmriprep with memory allocation
	# If plugin.yml file already exists, delete it and make anew
	if [ -f $bids_root_dir/fmriprep_plugin_${s}.yml ]; then
		rm -rf $bids_root_dir/fmriprep_plugin_${s}.yml
	fi
	
	# Create new yaml file
	touch $bids_root_dir/fmriprep_plugin_${s}.yml
	echo "plugin: LegacyMultiProc" >> $bids_root_dir/fmriprep_plugin_${s}.yml
	echo "plugin_args: {maxtasksperchild: $nodes, memory_gb: $mem_gb, n_procs: $cores, raise_insufficient: false}" >> $bids_root_dir/fmriprep_plugin_${s}.yml

	# Add singularity to PATH
	if [ -z `command -v singularity` ]; then
		module load singularity/3.6.4
	fi

	# Generate fmriprep singularity image if one doesn't exist in /N/dcwan/projects/irf/containers
	if [ ! -f /N/dcwan/projects/irf/containers/fmriprep-${fmriprep_version}.simg ]; then
		echo ""
		echo "Building singularity image for fmriprep-${fmriprep_version}"
		echo ""

		export SINGULARITY_TMPDIR="$bids_root_dir"
		export SINGULARITY_CACHEDIR="$bids_root_dir"
		singularity build /N/dcwan/projects/irf/containers/fmriprep-${fmriprep_version}.simg docker://nipreps/fmriprep:${fmriprep_version}
	fi

	# Make temp directory in user's slate directory for files
	if [ ! -d /N/slate/$(whoami)/fmriprep-${fmriprep_version}_temp ]; then
		mkdir /N/slate/$(whoami)/fmriprep-${fmriprep_version}_temp
		chmod -R 777 /N/slate/$(whoami)/fmriprep-${fmriprep_version}_temp
	fi

	# Make folder in user's slate directory for storing slurm error and text reports
	if [ ! -d /N/slate/$(whoami)/slurm_output ]; then
		mkdir /N/slate/$(whoami)/slurm_output
		chmod -R 777 /N/slate/$(whoami)/slurm_output
	fi

	# Run fmriprep
	if [ ! -d $bids_root_dir/derivatives/fmriprep/sub-${s}/figures ]; then
		echo ""
		echo "Running fmriprep-${fmriprep_version} on sub-${s}"
		echo ""

		fmriprep_command="unset PYTHONPATH; singularity run \
		-B ${bids_root_dir}:${bids_root_dir} \
		-B /N/slate/$(whoami)/fmriprep-${fmriprep_version}_temp:/N/slate/$(whoami)/fmriprep-${fmriprep_version}_temp \
		-B ${output_dir}:${output_dir} \
		/N/dcwan/projects/irf/containers/fmriprep-${fmriprep_version}.simg \
		$bids_root_dir $output_dir												\
			participant															\
			--participant-label $s 												\
			--skip-bids-validation 												\
			--md-only-boilerplate												\
			--fs-license-file $FREESURFER_HOME/license.txt						\
			--fs-no-reconall 													\
			--output-spaces MNI152NLin2009cAsym:res-2							\
			--nthreads $cores													\
			--stop-on-first-crash 												\
			--resource-monitor 													\
			--low-mem 															\
			--mem_mb $reduced_mem 												\
			--use-plugin $bids_root_dir/fmriprep_plugin_${s}.yml 				\
			--verbose 															\
			-w /N/slate/$(whoami)/fmriprep-${fmriprep_version}_temp"  


		touch $bids_root_dir/derivatives/fmriprep_${s}.sh
		touch /N/slate/$(whoami)/slurm_output/slurm_${name}-${s}.txt
		touch /N/slate/$(whoami)/slurm_output/slurm_${name}-${s}.err
		chmod -R 777 $bids_root_dir/derivatives

		echo "#!/bin/bash" > $bids_root_dir/derivatives/fmriprep_${s}.sh
		echo "" >> $bids_root_dir/derivatives/fmriprep_${s}.sh
		echo $fmriprep_command >> $bids_root_dir/derivatives/fmriprep_${s}.sh

		cd $bids_root_dir/derivatives

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
		$bids_root_dir/derivatives/fmriprep_${s}.sh


		cd $current_dir
	fi
done


