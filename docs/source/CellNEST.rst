CellNEST Analysis of CCI
=========================================== 

The following script demonstrates how we use CellNEST to conduct the CCI analysis.

- The CellNEST documentation, tutorial and package could be found here: https://github.com/schwartzlab-methods/CellNEST

- The paper could be found here: https://doi.org/10.1038/s41592-025-02721-3

CellNEST could only be run on GPU node. Here we demonstrate how to run on High Computing Cluster (HPC) using singularity by submitting LSF jobs.

1. Preprocessing Step
-------------------------

.. code-block:: r

   #BSUB -J CellNEST_IPMN_preprocess
   #BSUB -W 12:00
   #BSUB -q gpu
   #BSUB -gpu "num=1:gmem=16"
   #BSUB -n 4
   #BSUB -R "rusage[mem=16GB]"
   #BSUB -M 16
   #BSUB -o logs/cellnest_preprocess.%J.out
   #BSUB -e logs/cellnest_preprocess.%J.err
   #BSUB -cwd /home/lku/CellNEST
   #BSUB -B 
   #BSUB -N 
   #BSUB -u lku@mdanderson.org

   module load singularity

   singularity exec --nv /home/lku/cellnest_image.sif \
     bash cellnest preprocess \
     --data_name='IPMN_LG2' \ 
     --data_from='data/IPMN_LG2/'


2. Train the CellNEST Model 
-------------------------

.. code-block:: r

   #!/bin/bash
   #BSUB -J LG2_train
   #BSUB -W 48:00
   #BSUB -q gpu
   #BSUB -gpu "num=1:gmem=16"
   #BSUB -n 4
   #BSUB -R "rusage[mem=32GB]"
   #BSUB -M 32
   #BSUB -o logs/IPMN/LG2_train.%J.out
   #BSUB -e logs/IPMN/LG2_train.%J.err
   #BSUB -cwd /home/lku/CellNEST
   #BSUB -B
   #BSUB -N
   #BSUB -u lku@mdanderson.org

   module load singularity

   for RUNID in 1 2 3 4 5
   do
       echo "Starting run_id=${RUNID} at $(date)"
       singularity exec --nv \
           --bind /home/lku:/home/lku \
           /home/lku/cellnest_image.sif \
           bash cellnest run \
               --data_name='IPMN_LG2' \
               --num_epoch 80000 \
               --model_name='IPMN_LG2' \
               --run_id=${RUNID} \
               --manual_seed='yes' \
               --seed=${RUNID} \
               > logs/IPMN_train/LG2_output_run_${RUNID}.log 2>&1
       echo "Finished run_id=${RUNID} at $(date)"
   done



3. Postprocessing Step   
-------------------------

.. code-block:: bash

   #!/bin/bash
   #BSUB -J IPMN_postprocess_LG2
   #BSUB -W 24:00
   #BSUB -q gpu
   #BSUB -gpu "num=1:gmem=16"
   #BSUB -n 4
   #BSUB -R "rusage[mem=32GB]"
   #BSUB -M 32
   #BSUB -o logs/cellnest_postprocess.%J.out
   #BSUB -e logs/cellnest_postprocess.%J.err
   #BSUB -cwd /home/lku/CellNEST
   #BSUB -B
   #BSUB -N
   #BSUB -u lku@mdanderson.org

   module load singularity

   singularity exec --nv \
       --bind /home/lku:/home/lku \
       /home/lku/cellnest_image.sif \
       bash cellnest postprocess \
           --data_name='IPMN_LG2' \
           --model_name='IPMN_LG2' \
           --total_runs=5





