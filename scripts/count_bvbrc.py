import pandas as pd

df = pd.read_csv("data/raw/bvbrc_amr.txt", sep="\t", low_memory=False,
                 quoting=3, dtype=str)
df["genome_name"] = df["genome_name"].str.replace('"', "", regex=False).str.strip()
df["antibiotic"] = df["antibiotic"].str.strip().str.lower()

kp = df[df["genome_name"].str.contains("Klebsiella pneumoniae", na=False)]
print("Klebsiella rows:", len(kp), " genomes:", kp["genome_id"].nunique())
print("Lab methods:\n", kp["laboratory_typing_method"].value_counts(dropna=False).head(8))

kp = kp[kp["resistant_phenotype"].isin(["Resistant", "Susceptible"])]
print("Rows with R/S call:", len(kp), " genomes:", kp["genome_id"].nunique())

tab = (kp.drop_duplicates(["genome_id", "antibiotic"])
         .groupby(["antibiotic", "resistant_phenotype"]).size()
         .unstack(fill_value=0))
tab["min_class"] = tab.min(axis=1)
print(tab.sort_values("min_class", ascending=False).head(20))
