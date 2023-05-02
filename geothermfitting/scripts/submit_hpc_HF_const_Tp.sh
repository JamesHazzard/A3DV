#!/bin/bash

# Define home directory
HOMEDIR=/rds/general/user/jah220/home/a3dv/geothermfitting
cd $HOMEDIR/scripts
cp $HOMEDIR/config/jameshpc.ini $HOMEDIR/scripts/config.ini
date=$(date '+%Y-%m-%d_%H-%M-%S')
echo "date = ${date}" >> $HOMEDIR/scripts/config.ini
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')

# Define chunk size (i.e. number of individual LAB runs in a chunk)
chunksize=20

# Define total number of runs
minnumber=0
maxnumber=999

#Specify tomography model
model=ANT-20

# Make output folders
SRCDIR=/rds/general/user/jah220/ephemeral/a3dv
mkdir -p $SRCDIR
SRCDIR=${SRCDIR}/geothermfitting
mkdir -p $SRCDIR
SRCDIR=${SRCDIR}/data_output
mkdir -p $SRCDIR

#Make job submission script
cat << \EOF > concept
#!/bin/bash
#PBS -lselect=1:ncpus=1:mem=12gb
#PBS -lwalltime=72:0:0

# Load personal conda environment, make sure necessary depencies and libraries are installed in this environment
module load anaconda3/personal
source activate py39

# Get input/output directories
HOMEDIR=/rds/general/user/jah220/home/a3dv/geothermfitting
SRCDIR=/rds/general/user/jah220/ephemeral/a3dv/geothermfitting

# Change into directory pertaining to chunk
cd $SRCDIR
rm -rf batch_start_foo_end_bar
mkdir batch_start_foo_end_bar
cp $HOMEDIR/scripts/config.ini batch_start_foo_end_bar
cp $HOMEDIR/scripts/make_HF_const_Tp.py batch_start_foo_end_bar
cp $HOMEDIR/scripts/fit_geotherm.f batch_start_foo_end_bar
cp $HOMEDIR/scripts/Makefile batch_start_foo_end_bar
cp $HOMEDIR/scripts/libseis_min.py batch_start_foo_end_bar
cp $HOMEDIR/scripts/lib_anelasticity.py batch_start_foo_end_bar
cd batch_start_foo_end_bar
make clean
make all

for i in $(seq foo 1 bar); do
  python3 make_HF_const_Tp.py -v distribution -i ${i}
done

EOF

# Loop for creating chunks
for startrun in $(seq $minnumber $chunksize $maxnumber); do

    step=$(echo ${chunksize}-1 | bc -l)
    endrun=$(echo ${startrun} + ${step} | bc -l)
	echo "First run in chunk is $startrun; final run is $endrun, model is $model"

    # Write variables into job submission script
    cp concept junk
	sed "s/foo/$startrun/g" junk | sed "s/bar/$endrun/g" | sed "s/zoo/$model/g" > send_HF_const_Tp_${startrun}_${endrun}_${model}.sh

    # Make script executable
	chmod +x send_HF_const_Tp_${startrun}_${endrun}_${model}.sh

    # Submit chunk of LAB runs
	qsub ./send_HF_const_Tp_${startrun}_${endrun}_${model}.sh

    # Wait 5 seconds before submitting next chunk of runs
    wait

    # Clean up directory
    rm junk

done
