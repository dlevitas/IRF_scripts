#!/bin/bash


# User inputs:
bids_root_dir=/N/dcwan/projects/irf/BIDS/BIDS_tutorial
subjects=01
nodes=1 # Number of nodes to use on HPC job submission(s)
cores=4 # Number of processors to use on HPC job submission(s). 16 is the max
mem_mb=25000 # Memory to allocate
walltime=10:00:00 # Walltime for job
email=dlevitas@iu.edu # email for sending updates on the job submission(s)
name=qsiprep # Name for job submissions
qsiprep_version= # Can be left blank. If so, will default to version 20.2.0, which is the long term support (LTS)
output_dir= # can be left blank, in which case the qsiprep output will be stored in $bids_root_dir/derivatives/qsiprep
overwrite=yes # yes or no. If yes, will remove files and restart

# Begin:

#Check user inputs
echo ""
if [[ -z "$bids_root_dir" || -z "$subjects" || -z "$nodes" || -z "$cores" || -z "$mem_mb" || -z "$walltime" || -z "$email" || -z "$name" || -z "$overwrite" ]]; then
	echo "Make sure user inputs variables are defined. Exiting now"
	return
fi

# Convert mem_mb to gb
mem_gb=`echo $((mem_mb / 1000))`

# Check mriqc version
if [ -z "$qsiprep_version" ]; then
	qsiprep_version=0.13.0RC2
fi

if [ -z "$output_dir" ]; then
	output_dir=$bids_root_dir/derivatives
fi
chmod -R 777 $output_dir

if [ ! -d $output_dir/qsiprep ]; then
	mkdir $output_dir/qsiprep
fi

for s in ${subjects[*]}
do

	# Add Singularity to PATH
	if [ -z `command -v singularity` ]; then
		module load singularity/3.6.4
	fi

	# Generate MRIQC singularity image if need be
	if [ ! -f /N/dcwan/projects/irf/containers/qsiprep-${qsiprep_version}.sif ]; then
		echo ""
		echo "Building singularity image for qsiprep-${qsiprep_version}"
		echo ""

		export SINGULARITY_TMPDIR="/N/slate/$(whoami)"
		export SINGULARITY_CACHEDIR="/N/slate/$(whoami)"
		singularity build /N/dcwan/projects/irf/containers/qsiprep-${qsiprep_version}.sif docker://pennbbl/qsiprep:${qsiprep_version}
	fi

	# Create plugin.yml file, which helps mriqc with memory allocation
	if [ -f $bids_root_dir/qsiprep_plugin_${s}.yml ]; then
		rm -rf $bids_root_dir/qsiprep_plugin_${s}.yml
	fi

	touch $bids_root_dir/qsiprep_plugin_${s}.yml
	echo "plugin: LegacyMultiProc" >> $bids_root_dir/qsiprep_plugin_${s}.yml
	echo "plugin_args: {maxtasksperchild: $nodes, memory_gb: $mem_gb, n_procs: $cores, raise_insufficient: false}" >> $bids_root_dir/qsiprep_plugin_${s}.yml

	# Make temp directory in user's slate directory for files
	if [ ! -d /N/slate/$(whoami)/qsiprep-${qsiprep_version}_temp ]; then
		mkdir /N/slate/$(whoami)/qsiprep-${qsiprep_version}_temp
		chmod -R 777 /N/slate/$(whoami)/qsiprep-${qsiprep_version}_temp
	fi

	# Make folder in user's slate directory for storing slurm error and text reports
	if [ ! -d /N/slate/$(whoami)/slurm_output ]; then
		mkdir /N/slate/$(whoami)/slurm_output
		chmod -R 777 /N/slate/$(whoami)/slurm_output
	fi

	# Run MRIQC
	if [ ! -d $bids_root_dir/derivatives/qsiprep/sub-${s} ]; then
		echo ""
		echo "Running qsiprep on participant $s"
		echo ""
		qsiprep_command="unset PYTHONPATH; singularity run \
		-B ${bids_root_dir}:${bids_root_dir} \
		-B /N/slate/$(whoami)/qsiprep-${qsiprep_version}_temp:/N/slate/$(whoami)/qsiprep-${qsiprep_version}_temp \
		-B ${output_dir}:${output_dir} \
		 /N/dcwan/projects/irf/containers/qsiprep-${qsiprep_version}.sif 	\
		$bids_root_dir $bids_root_dir/derivatives				 		\
		participant 													\
		--participant_label $s 											\
		--skip_bids_validation 											\
		--output-space T1w			 									\
		--combine_all_dwis 												\
		--output-resolution 1											\
		--nthreads $cores												\
		--fs-license-file $FREESURFER_HOME/license.txt					\
		--stop-on-first-crash 											\
		--resource-monitor 												\
		--use-plugin $bids_root_dir/qsiprep_plugin_${s}.yml 			\
		-w /N/slate/$(whoami)/qsiprep-${qsiprep_version}_temp"

		touch $bids_root_dir/derivatives/qsiprep_${s}.sh
		touch /N/slate/$(whoami)/slurm_output/slurm_${name}-${s}.txt
		touch /N/slate/$(whoami)/slurm_output/slurm_${name}-${s}.err
		# chmod -R 777 $bids_root_dir/derivatives

		echo "#!/bin/bash" > $bids_root_dir/derivatives/qsiprep_${s}.sh
		echo "" >> $bids_root_dir/derivatives/qsiprep_${s}.sh
		echo $qsiprep_command >> $bids_root_dir/derivatives/qsiprep_${s}.sh

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
		$bids_root_dir/derivatives/qsiprep_${s}.sh

		cd $current_dir
	fi
done