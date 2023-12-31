This is a simple script / Jupyter notebook which runs BLASTn on extendors from raw SPLASH output. 

The script takes two command-line arguments: 
    1. Path to a working directory. This can be a directory that exists, or a subdirectory to be created within an existing directory. 
    2. Path to SPLASH output. This should have columns of the form of 'most_freq_target_N' and 'cnt_most_freq_target_N'. 
  
The script uses Python's multiprocessing library to run BLAST on chunks of extendors in the SPLASH output. 
    1. In the working directory, FASTAs and BLAST output files are created as intermediates. 
    2. In the working directory, once BLAST has completed, a file called SPLASH_output_BLAST_merged.tsv is created. 
    3. This file reports all BLAST hits having e-value < 0.05 for all anchors' extendors. The table is presented with one row for each 
pair of (anchor, target, BLAST hit). 

The script can be run using a command of this form: 
sbatch /oak/stanford/groups/horence/george/blast_internal_parallelization/blast_internal_parallelization.sh /oak/stanford/groups/horence/george/blast_internal_parallelization/test_Nov1_2023 /oak/stanford/groups/horence/george/blast_internal_parallelization/test/test_data.tsv

To change the BLAST script used, the user may modify the path to the BLAST bash script on line 12 in blast_internal_parallelization.py.

A virtual environment specific to this procedure is provided.
