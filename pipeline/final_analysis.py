import pandas as pd
import os
import glob

base_dir = '/home/egarcia/workspace/github/sv-code/pipeline'


def intersect_stats():
    for intersect_file in glob.glob(os.path.join(base_dir, "results", "analysis", "raw", "*_vs_*")):
        data = pd.read_csv(intersect_file, sep='\t', low_memory=False)
        data =  data.drop_duplicates(subset=["IDX_SRC"], keep="last")
        print(intersect_file + "--> " + str(len(data.index)))


intersect_stats()
