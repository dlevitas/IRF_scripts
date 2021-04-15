#!/bin/bash

#Notes:
#1) This script automatically sets up things to run xcpEngine (e.g. creates cohort file(s))
#2) Users will modify/edit "user inputs" before running
#3) Script will submit slurm job(s) at the end
#4) This script only works for functional connectivity pipeline, using resting-state data
#5) Only one standard output space accepted. Either specify it as a user input, or script will find it


#User inputs (required):
bids_root_dir=/N/dcwan/projects/irf/BIDS/BIDS_tutorial
subjects=(01) #Can leave blank, in which case, script will find all fmriprep preprocessed subjects
fmriprep_data_dir=${bids_root_dir}/derivatives/fmriprep # directory where the friprep output is stored
space=(MNI152NLin2009cAsym) # specifywhich fmriprep standard output space do you wish to use for xcpEngine
design_file=fc-24p_gsr.dsn # Just the design file name, not the entire path. For choices, see https://github.com/PennBBL/xcpEngine/tree/master/designs
nodes=1 #keep at one
cores=1 # number of cores for xcpEngine (1 should be fine)
mem_mb=15000 # how much job memory you want
walltime=1:30:00 # how long you wish for the xcpEngine job to run for
email=dlevitas@iu.edu
name=xcpEngine # Name for job submissions

#User inputs (recommended, but not required):
xcpEngine_version=1.2.3
output_dir=${bids_root_dir}/derivatives/xcpEngine-${xcpEngine_version}
overwrite=yes #yes or no. If set to yes, will overwrite all xcpEngine work and restart

#Begin:

#Check user inputs
echo ""
if [[ -z "$bids_root_dir" || -z "$fmriprep_data_dir" || -z "$design_file" || -z "$nodes" || -z "$cores" || -z "$mem_mb" || -z "$walltime" || -z "$email" ]]; then
	echo "Make sure user inputs variables are defined. Exiting now"
	return
fi

if [ -z "$xcpEngine_version" ]; then
	xcpEngine_version=1.2.3
	echo "'xcpEngine_version' variable was not specified. Setting to $xcpEngine_version"
	echo ""
fi

if [ -z "$output_dir" ]; then
	output_dir=${bids_root_dir}/derivatives/xcpEngine-${xcpEngine_version}
	echo "'output_dir' variable was not specified. Setting to $output_dir"
	echo ""
fi

if [[ *"$design_file"* == *"anat"* ]]; then
  echo "Anatomical pipeline not configured for this script. Exiting now"
	return
fi

if [ -z "$overwrite" ]; then
	overwrite=no
	echo "'overwrite' variable was not specified. Setting to $overwrite"
	echo ""
fi

if [ ! -d $output_dir ]; then
	mkdir $output_dir
fi

if [ $overwrite == 'yes' ]; then
	rm -rf $output_dir/*
fi

#Add singularity to PATH
if [ -z `command -v singularity` ]; then
	module load singularity/3.6.4
fi

#Generate fmriprep singularity image if need be
if [ ! -f /N/dcwan/projects/irf/containers/xcpEngine-${xcpEngine_version}.simg ]; then
	echo ""
	echo "Building singularity image for xcpEngine version ${xcpEngine_version}"
	echo ""

	export SINGULARITY_TMPDIR="/N/slate/$(whoami)"
	export SINGULARITY_CACHEDIR="/N/slate/$(whoami)"
	singularity build /N/dcwan/projects/irf/containers/xcpEngine-${xcpEngine_version}.simg docker://pennbbl/xcpengine:${xcpEngine_version}
fi

#Copy xcpEngine stuff into user's slate directory
if [ ! -d /N/slate/$(whoami)/xcpEngine-github ]; then
	echo "Copying xcpEngine-github folder to: /N/slate/$(whoami)"
	cp -r /N/dcwan/projects/irf/xcpEngine-github /N/slate/$(whoami)
fi
design_file="/N/slate/$(whoami)/xcpEngine-github/designs/${design_file}"

#Get subject IDs, and check to see if data is multi-session
if [ -z "$subjects" ]; then
	subjects=(`ls -d $fmriprep_data_dir/sub-*/ | xargs -n 1 basename`)
fi
sesdircount=$(find $fmriprep_data_dir/sub-*/ -name ses-* -type d | wc -l)

rm -rf $output_dir/temp*/*

if [ "$sesdircount" -gt 0 ]; then
	path=$fmriprep_data_dir/sub-*/ses-*/func
else
	path=$fmriprep_data_dir/sub-*/func
fi

if [[ -z "$space" ]]; then
	space=($(find $path -name *desc-preproc_bold.nii.gz -type f | sed 's|.*\(space-\)||g' | cut -d _ -f 1 | sort -u))
	if [ $space -gt 1 ]; then
		echo "More than one standard output space identified. Please select only one for xcpEngine. Exiting now"
		return
	fi
fi

# task=($(find $path -name *desc-preproc_bold.nii.gz -type f | sed 's|.*\(task-\)||g' | cut -d _ -f 1 | sort -u))
task=`find $path -name *desc-preproc_bold.nii.gz -type f | sed 's|.*\(task-\)||g' | cut -d _ -f 1 | sort -u`
if [[ ! "${task[@]}" =~ "rest" ]]; then
	echo "no resting state data found, cannot perform xcpEngine. Exiting now"
	return
fi

#Iterate through subject ID list (and sessions if applicable) to complete cohort file(s)
if [ ! -f $output_dir/cohort.csv ]; then
	touch $output_dir/cohort.csv
	mkdir $output_dir/temp
else
	printf '' | tee $output_dir/cohort.csv
fi
cohort_file=$output_dir/cohort.csv
temp_dir=$output_dir/temp

if [ "$sesdircount" -eq 0 ]; then
	echo "id0,id1,img" >> $cohort_file
else
	echo "id0,id1,id2,img" >> $cohort_file
fi

for sub in ${subjects[*]}
do
	if [ "$sesdircount" -gt 0 ]; then
		sessions=(`ls -d $fmriprep_data_dir/sub-${sub}/ses-*/ | xargs -n 1 basename`)

		for ses in ${sessions[*]}
		do

			reg_preproc=($(find $fmriprep_data_dir/sub-${sub}/ses-${ses}/func -name *task-rest*space-${space}*desc-preproc_bold.nii.gz -type f | sort -n | sed "s|$fmriprep_data_dir/||g"))

			if [ ${#reg_preproc[@]} -gt 0 ]; then
				for i in ${!reg_preproc[*]}
				do
					run=$((i+1))
					echo "sub-$sub,ses-$ses,run-${run},${reg_preproc[$i]}" >> $cohort_file
				done
			else
				echo ""
				echo "There do not appear to be any functional preprocessed data for subject $sub, session $ses"
				echo ""
			fi
		done

	else
		reg_preproc=($(find $fmriprep_data_dir/sub-${sub}/func -name *task-rest*space-${space}*desc-preproc_bold.nii.gz -type f | sort -n | sed "s|$fmriprep_data_dir/||g"))

		if [ ${#reg_preproc[@]} -gt 0 ]; then
			for i in ${!reg_preproc[*]}
			do
				run=$((i+1))
				echo "sub-$sub,run-${run},${reg_preproc[$i]}" >> $cohort_file
			done

		else
			echo ""
			echo "There do not appear to be any functional preprocessed data for subject $sub"
			echo ""
		fi
	fi
done

#Set up Singularity command for running xcpEngine on HPC
xcpEngine_command="unset PYTHONPATH; singularity run -B /N/slate/$(whoami):/N/slate/$(whoami) \
				/N/dcwan/projects/irf/containers/xcpEngine-${xcpEngine_version}.simg \
				-d $design_file \
				-c $cohort_file \
				-i $temp_dir \
				-r $fmriprep_data_dir \
				-t 2 \
				-o $output_dir"

if [ -f $output_dir/xcpEngine_slurm.sh ]; then
	rm -rf $output_dir/xcpEngine_slurm.sh
fi

touch $output_dir/xcpEngine_slurm.sh
touch /N/slate/$(whoami)/slurm_output/slurm_${name}.txt
touch /N/slate/$(whoami)/slurm_output/slurm_${name}.err
chmod -R 777 $bids_root_dir/derivatives

echo "#!/bin/bash" > $output_dir/xcpEngine_slurm.sh
echo "" >> $output_dir/xcpEngine_slurm.sh
echo $xcpEngine_command >> $output_dir/xcpEngine_slurm.sh

# Submit slurm job
sbatch -p general \
-J ${name} \
-o /N/slate/$(whoami)/slurm_output/slurm_${name}.txt \
-e /N/slate/$(whoami)/slurm_output/slurm_${name}.txt \
--mail-type=ALL \
--mail-user=$email \
--nodes=$nodes \
--cpus-per-task=$cores \
--ntasks-per-node=1 \
--time=$walltime \
--mem=$mem_mb \
$output_dir/xcpEngine_slurm.sh