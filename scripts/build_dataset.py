import pandas as pd

df = pd.read_csv("data/raw/bvbrc_amr.txt", sep="\t", low_memory=False,
                 quoting=3, dtype=str)
df["genome_name"] = df["genome_name"].str.replace('"', "", regex=False).str.strip()
df["antibiotic"] = df["antibiotic"].str.strip().str.lower()

kp = df[df["genome_name"].str.contains("Klebsiella pneumoniae", na=False)]
kp = kp[kp["laboratory_typing_method"].notna()]
kp = kp[kp["resistant_phenotype"].isin(["Resistant", "Susceptible"])]

drugs = ["ciprofloxacin", "ceftriaxone", "gentamicin"]
kp = kp[kp["antibiotic"].isin(drugs)]

# drop genome/drug pairs with conflicting labels
n_lab = kp.groupby(["genome_id", "antibiotic"])["resistant_phenotype"].nunique()
bad = n_lab[n_lab > 1].index
print("Conflicting genome/drug pairs removed:", len(bad))
kp = kp[~kp.set_index(["genome_id", "antibiotic"]).index.isin(bad)]

kp = kp.drop_duplicates(["genome_id", "antibiotic"])
kp["label"] = (kp["resistant_phenotype"] == "Resistant").astype(int)
wide = kp.pivot(index="genome_id", columns="antibiotic", values="label")

print("\nPer drug (1=Resistant, 0=Susceptible):")
for d in drugs:
    s = wide[d].dropna()
    print(f"  {d}: n={len(s)}  R={int(s.sum())}  S={int((1-s).sum())}")
print("Genomes with a label for all 3 drugs:", int(wide.notna().all(axis=1).sum()))
print("Genomes with at least one label:", len(wide))

wide.to_csv("data/processed/labels_wide.csv")
wide.index.to_series().to_csv("data/processed/genome_list.txt", index=False, header=False)
print("\nSaved data/processed/labels_wide.csv and genome_list.txt")
