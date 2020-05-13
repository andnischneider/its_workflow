from os.path import join as opj


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
            sample_id = r.sample
        if "sra_id" in df.columns:
            sample_dict[sample_id] = r.sra_id
        if pool_id not in pool_dict.keys():
            pool_dict[pool_id] = {"R1": [], "R2": []}
        pool_dict[pool_id]["R1"].append(
            opj(data_dir, "{}_R1.fastq.gz".format(sample_id)))
        pool_dict[pool_id]["R2"].append(
            opj(data_dir, "{}_R2.fastq.gz".format(sample_id)))
    return sample_dict, pool_dict
