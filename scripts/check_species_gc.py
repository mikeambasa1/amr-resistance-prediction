import glob, os
import pandas as pd

rows = []
for p in sorted(glob.glob("data/genomes/*.fna")):
    gc = at = 0
    with open(p) as f:
        for line in f:
            if line[0] == ">":
                continue
            s = line.strip().upper()
            gc += s.count("G") + s.count("C")
            at += s.count("A") + s.count("T")
    rows.append((os.path.basename(p)[:-4], gc + at, round(100 * gc / (gc + at), 2)))

df = pd.DataFrame(rows, columns=["genome_id", "bases", "gc_percent"])
df.to_csv("data/processed/gc_check.csv", index=False)
print(df[["bases", "gc_percent"]].describe().round(2))
odd = df[(df.gc_percent < 56.0) | (df.gc_percent > 58.5) | (df.bases > 6_500_000) | (df.bases < 4_800_000)]
print("\nOutliers (GC outside 56-58.5% or size outside 4.8-6.5 Mb):", len(odd), "of", len(df))
print(odd.sort_values("gc_percent").head(25).to_string(index=False))
