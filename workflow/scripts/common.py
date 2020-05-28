from os.path import join as opj, exists
import sys

def exit(sample_id, R1, R2):
    sys.exit("""

###### WORKFLOW ERROR ######
Sample {} has no SRA id and missing at least one of
R1: {}
R2: {}
###### WORKFLOW ERROR ######

            """.format(sample_id, R1, R2))


def get_samples(df, data_dir, results_dir):
    """
    Parses the sample list and generates a sample dictionary and a pool
    dictionary with information on how the workflow should handle input.

    :param df: pandas DataFrame with samples
    :param data_dir: Input data directory from config
    :return: sample_dict, pool_dict dictionaries
    """
    sample_dict = {}
    dirnames = []
    for i in df.index:
        r = df.loc[i]
        sample_id = r["sample"]
        # Check if there's a dirname for sample
        try:
            dirname = r.dirname
        except AttributeError:
            dirname = ""
        if dirname not in sample_dict.keys():
            sample_dict[dirname] = {sample_id: {}}
        else:
            sample_dict[dirname][sample_id] = {}
        dirnames.append(dirname)
        # Set start location for sample
        R1 = opj(data_dir, dirname, "{}_R1.fastq.gz".format(sample_id))
        R2 = opj(data_dir, dirname, "{}_R2.fastq.gz".format(sample_id))
        sample_dict[dirname][sample_id]["R1"] = R1
        sample_dict[dirname][sample_id]["R2"] = R2
        # See if the sample has an SRA id that can be used to download it
        print(R1)
        if "sra_id" in df.columns:
            sample_dict[dirname][sample_id]["acc"] = r.sra_id
        # If not, check that the starting file exists under data dir
        else:
            if not exists(R1) or not exists(R2):
                exit(sample_id, R1, R2)
    return sample_dict, list(set(dirnames))
