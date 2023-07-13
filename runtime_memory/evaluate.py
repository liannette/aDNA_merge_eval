import pandas as pd

filenames = ["AdapterRemoval", 
             "bbmerge", 
             "ClipAndMerge", 
             "fastp", 
             "leeHom", 
             "SeqPrep", 
             "seqtk_adna_trim"]

# calculate summary statistics for time (s) and memory (USS) usage
for filename in filenames:
    path = "output/benchmarks/" + filename + ".tsv"
    df = pd.read_csv(path, sep="\t")
    df = df[["s", "max_uss"]].describe()
    df.to_csv("output/evaluation/{}.csv".format(filename))