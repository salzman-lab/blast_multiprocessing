This is a simple script / Jupyter notebook which runs BLASTn on extendors from raw SPLASH output. 

The script takes two command-line arguments: 
    1. Path to a working directory. This can be a directory that exists, or a subdirectory to be created within an existing directory. 
    2. Path to SPLASH output. This should have columns entitled 'most_freq_target_N' and 'cnt_most_freq_target_N'. 
  
The script uses Python's multiprocessing library to run BLAST on chunks of extendors in the SPLASH output. 
In the working directory, FASTAs and BLAST output files are created as intermediates. 
In the working directory, once BLAST has completed, a file called SPLASH_output_BLAST_merged.tsv is created. 
This file reports all BLAST hits having e-value < 0.05 for all anchors' extendors. The table is presented with one row for each 
pair of (anchor, target, BLAST hit). 

The script can be run using a command of this form: 

sbatch /oak/stanford/groups/horence/george/blast_internal_parallelization/blast_internal_parallelization.sh /oak/stanford/groups/horence/george/blast_internal_parallelization/test_Nov1_2023 /oak/stanford/groups/horence/george/blast_internal_parallelization/test/test_data.tsv
