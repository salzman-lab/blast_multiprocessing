import numpy as np
import pandas as pd
import os
import multiprocessing as mp
import glob
import sys

def blast(fasta, dummy):
    """
    Simple helper to run BLAST. 
    """
    os.system('bash /oak/stanford/groups/horence/george/splash_utils/blast.sh ' + fasta)
    return 

def main():
    
    ## The user sets a path to a SPLASH input and a working directory.
    splash_input = sys.argv[1]
    working_dir = sys.argv[2]

    ## Create and/or move to the working directory.  
    try:
        os.mkdir(working_dir)
        os.chdir(working_dir)
    except FileExistsError:
        os.chdir(working_dir)

    ## Load SPLASH output, available at the user-specified path. 
    splash_input = pd.read_csv(splash_input,sep='\t')

    ## Determine how many target columns we have. 
    most_freq_columns = [i for i in splash_input.columns if 'most_freq_target' in i and 'cnt' not in i]

    ## Instantiate a table having columns for anchor, target 1, and target 1 count. 
    ## Rename columns to 'target' and 'count'.
    df = splash_input[['anchor','most_freq_target_1','cnt_most_freq_target_1']].rename(columns={'most_freq_target_1':'target','cnt_most_freq_target_1':'target_count'})
    df['target_rank'] = 1

    ## Load and rename the remaining target and target count columns. 
    ## Append these columns to produce a long-form dataframe having rows for each anchor-target pair. 
    for i in range(2, len(most_freq_columns) + 1): 
        fd = splash_input[['anchor','most_freq_target_'+str(i),'cnt_most_freq_target_'+str(i)]].rename(columns={'most_freq_target_'+str(i):'target','cnt_most_freq_target_'+str(i):'target_count'})
        fd['target_rank'] = i
        df = pd.concat([df,fd])
    df = df[df['target'] != '-'].reset_index(drop=True)

    ## Introduce the SPLASH statistics back to the long-form dataframe. 
    df = df.merge(splash_input[[i for i in splash_input.columns if 'most_freq_target' not in i]],how='left')

    ## Get a Pandas series representing extendors. 
    cs = df['anchor'] + df['target']

    ## Get the number of CPUs available on this Jupyter notebook. 
    workers = int(os.environ['SLURM_JOB_CPUS_PER_NODE']) 

    ## Split the compactors into chunks of number == num CPUs. 
    cs_inds = np.array_split(cs.index,workers)

    ## Write FASTAs, one FASTA per available CPU.
    for i in range(len(cs_inds)):
        fasta = open('11blast_'+str(i)+'.fasta','a')
        sel = cs.iloc[cs_inds[i]]
        for ind in sel.index:
            fasta.write('>'+str(i)+'_'+str(ind)+'\n')
            fasta.write(cs[ind]+'\n')
        fasta.close()

    ## Use Python multiprocessing to submit a subprocess for each FASTA. 
    if __name__ == "__main__":
        with mp.Pool(workers) as p:
            outs = p.starmap(blast, [(i,0) for i in glob.glob('11*fasta')])       

    ## Load the FASTA input and BLAST output tables. 
    ## Except cases where the BLAST output is empty. 
    outs = glob.glob('*BLAST_out.txt')
    ins = glob.glob('*.fasta') 

    ## Read the first FASTA. 
    fasta = pd.read_csv(ins[0],header=None,sep='\t',engine='python')

    ## Read the remaining FASTAs, appending to the first. 
    for i in range(1,len(ins)):
        fasta1 = pd.read_csv(ins[i],header=None,sep='\t',engine='python')
        fasta = pd.concat([fasta,fasta1])

    ## Read the first BLAST output. 
    blast = pd.read_csv(outs[0],header=None,sep='\t',engine='python')

    ## Read the remaining BLAST outputs, appending to the first. 
    for i in range(1,len(outs)):
        try:
            blast1 = pd.read_csv(outs[i],header=None,sep='\t',engine='python')
            blast = pd.concat([blast,blast1])
        except pd.errors.EmptyDataError:
            pass

    ## Define the BLAST output columns concordantly with the specified 'fmt' in the BLAST command.  
    blast.columns = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseqid sgi sacc slen staxids'.split(' ')
    headers = [i[1:] for i in fasta[0] if i[0] == '>']
    seqs = [i for i in fasta[0] if i[0] != '>']
    fast = pd.DataFrame({'qseqid':headers,'sequence':seqs})
    blast['qseqid'] = blast['qseqid'].astype(str)

    ## Select BLAST hits having e-values < 0.05 and merge hits left onto a table generated from the input FASTA. 
    blast = blast[blast['evalue']<0.05]

    ## Merge the BLAST results onto the FASTA, thus joining BLAST information with query sequence information. 
    blast = fast.merge(blast,how='left')

    ## Extract the anchor and target; assume anchor is length 27.
    blast['anchor'] = blast['sequence'].str[:27]
    blast['target'] = blast['sequence'].str[27:]

    ## Merge SPLASH output with all (e-value < 0.05) BLAST results and write the file. 
    df = df.merge(blast,how='left')
    df.to_csv('SPLASH_output_BLAST_merged.tsv',sep='\t')
    
    return 
## every man is a piece of the continent, a part of the 
main()
