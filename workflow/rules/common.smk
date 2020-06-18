from snakemake.utils import validate
from scripts.common import get_samples
from os.path import join as opj
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://ubuntu:18.04"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
if config["test"]:
    config["spots"] = "-X 1000"

sample_df = pd.read_csv(config["samples"], sep="\t")
validate(sample_df, schema="../schemas/samples.schema.yaml")
samples, dirnames, mocks = get_samples(sample_df,
                                config["data_dir"],
                                config["mock"])