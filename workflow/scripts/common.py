from os.path import join as opj, exists
import sys

def get_samples(df, data_dir):
    """
    Parses the sample list and generates a sample dictionary and a pool
    dictionary with information on how the workflow should handle input.

    :param df: pandas DataFrame with samples
    :param data_dir: Input data directory from config
    :return: sample_dict, pool_dict dictionaries
    """
    sample_dict = {}
    pool_dict = {}
    for i in df.index:
        r = df.loc[i]
        pool_id = r["sample"]
        try:
            run_id = r.run_id
            sample_id = "{}.{}".format(pool_id, run_id)
        except AttributeError:
            sample_id = r["sample"]
        # Set start location for sample
        R1 = opj(data_dir, "{}_R1.fastq.gz".format(sample_id))
        R2 = opj(data_dir, "{}_R2.fastq.gz".format(sample_id))
        # See if the sample has an SRA id that can be used to download it
        if "sra_id" in df.columns:
            sample_dict[sample_id] = r.sra_id
        # If not, check that the starting file exists under data dir
        else:
            if not exists(R1) or not exists(R2):
                sys.exit("""
                
###### WORKFLOW ERROR ######
Sample {} has no SRA id and missing at least one of
R1: {}
R2: {}
###### WORKFLOW ERROR ######
                
                """.format(sample_id, R1, R2))
        # Add files to pool
        if pool_id not in pool_dict.keys():
            pool_dict[pool_id] = {"R1": [], "R2": []}
        pool_dict[pool_id]["R1"].append(R1)
        pool_dict[pool_id]["R2"].append(R2)
    return sample_dict, pool_dict
