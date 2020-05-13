configfile: "config.yml"
from os.path import join as opj
from src.parse_samples import get_samples

wildcard_constraints:
    run = "\d+"

# Use only a small set of reads if testing the workflow
if config["test"]:
    config["maxSpotId"] = "-X 1000"
else:
    config["maxSpotId"] = ""
samples, pools, config = get_samples(config["sample_file"], config)