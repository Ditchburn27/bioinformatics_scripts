import os
import pandas as pd
import sys

# this script will print CellRanger count command line for any RL provided
batches = ['RL2103_ga22_v3', 'RL2107_ga24_v3', 'RL2121_ga34_v3', 'RL1777_2d_v3',   'RL1612_34d_v2',  'RL2100_86d_v3',
           'RL2104_118d_v3', 'RL2108_179d_v3', 'RL2122_301d_v3', 'RL2125_422d_v3', 'RL2105_627d_v3', 'RL1613_2yr_v2', 
           'RL1786_2yr_v3',  'RL2129_3yr_v3',  'RL2109_4yr_v3',  'RL2106_6yr_v3',  'RL1614_8yr_v2',  'RL2110_10yr_v3',
           'RL2126_10yr_v3', 'RL2127_12yr_v3', 'RL2130_14yr_v3', 'RL2102_16yr_v3', 'RL2131_17yr_v3', 'RL2123_20yr_v3', 
           'RL2128_20yr_v3', 'RL2132_25yr_v3', 'RL2124_40yr_v3']

fastq_df = pd.read_csv( "/dd_userdata/usrdat03/userdata/cherring/ref_data/fastq_locations.csv")
rls = sys.argv[1:]

fastq_dict = {}
sample_dict = {}
id_dict = {}
for r_itr in rls:
    # get full paths to each sequencing run, usually just 1 run but some have 2+
    exp_paths = fastq_df.loc[fastq_df['Library_ID']==r_itr,'Path'].values
    full_fastq_paths = []
    sample_ids = []
    for exp_itr in exp_paths:
        fastq_path = exp_itr + '/fastqs/Neuronal_Maturation'
        dirs_itr = os.listdir(fastq_path)
        dir_nm = [ii for ii in dirs_itr if r_itr in ii]
        full_fastq_paths.append(fastq_path+'/'+str(dir_nm[0]))
        sample_ids.append(str(dir_nm[0]))
    # check it full paths matach number of exp_paths
    if len(exp_paths)!=len(full_fastq_paths):
        print('Check PATHS!!!!!!')
    fastq_dict[r_itr] = full_fastq_paths
    sample_dict[r_itr] = sample_ids
    id_dict[r_itr] = [ii for ii in batches if r_itr in ii]


for rl in rls:
    fastqs = ','.join(fastq_dict[rl])
    samples = ','.join(sample_dict[rl])
    print(f"cellranger count --fastqs={fastqs} --sample={samples} --id={id_dict[rl][0]} --localmem=400 --localcores=15 --transcriptome=/home/lditchburn/working_data_02/refdata-gex-GRCh38-2020-A")


#if __name__ == "__main__":
#    main()