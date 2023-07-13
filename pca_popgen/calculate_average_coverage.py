import pandas as pd

def weighted_average_m1(distribution, weights):
  
    numerator = sum([distribution[i]*weights[i] for i in range(len(distribution))])
    denominator = sum(weights)
    
    return round(numerator/denominator,2)


df = pd.read_csv("output/alignment/AdapterRemoval/HG002_markdup_coverage.tsv", sep="\t")
distribution = df["meandepth"]
weights = df["endpos"]
print(weighted_average_m1(distribution, weights))

df = pd.read_csv("output/alignment/AdapterRemoval/HG005_markdup_coverage.tsv", sep="\t")
distribution = df["meandepth"]
weights = df["endpos"]
print(weighted_average_m1(distribution, weights))