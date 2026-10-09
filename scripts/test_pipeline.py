import os, glob
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import balanced_accuracy_score, f1_score

os.chdir(os.path.expanduser("~/amr-project"))

def fasta_stats(path):
    lens, cur = [], 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur: lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur: lens.append(cur)
    lens.sort(reverse=True)
    total = sum(lens); half = 0; n50 = 0
    for l in lens:
        half += l
        if half >= total / 2:
            n50 = l
            break
    return total, len(lens), n50

ids = sorted(os.path.basename(p)[:-4] for p in glob.glob("data/amr/*.tsv"))
print("Annotated genomes found:", len(ids))

# 1. Quality check
rows = []
for i in ids:
    t, n, n50 = fasta_stats(f"data/genomes/{i}.fna")
    rows.append((i, t, n, n50))
qc = pd.DataFrame(rows, columns=["genome_id", "total_bp", "contigs", "n50"]).set_index("genome_id")
qc["pass"] = qc.total_bp.between(4_800_000, 6_500_000) & (qc.contigs <= 500)
qc.to_csv("data/processed/qc_report_test.csv")
print("Passed quality check:", int(qc["pass"].sum()), "of", len(qc))
print(qc[["total_bp", "contigs", "n50"]].describe().round(0).loc[["min", "50%", "max"]])
good = list(qc.index[qc["pass"]])

# 2. Gene/mutation table (AMR rows only; genomes with no hits become all zeros)
hits = {}
for i in good:
    df = pd.read_csv(f"data/amr/{i}.tsv", sep="\t")
    hits[i] = set(df.loc[df["Type"] == "AMR", "Element symbol"])
X = pd.DataFrame(0, index=good, columns=sorted(set().union(*hits.values())))
for i, s in hits.items():
    X.loc[i, list(s)] = 1
X = X.loc[:, X.sum() >= 3]
X.to_csv("data/processed/features_test.csv")
print("Feature table:", X.shape[0], "genomes x", X.shape[1], "features")

# 3. Labels
labels = pd.read_csv("data/processed/labels_wide.csv", dtype={"genome_id": str}).set_index("genome_id")
labels = labels.loc[X.index]
print("Label counts (1 = resistant, 0 = susceptible):")
for d in labels.columns:
    print(" ", d, labels[d].value_counts().to_dict())

# 4. Sanity check: do known mutations match the label?
for g in ["gyrA_S83Y", "gyrA_S83I"]:
    if g in X.columns:
        print("\nCiprofloxacin label among genomes with", g, ":",
              labels.loc[X[g] == 1, "ciprofloxacin"].value_counts().to_dict())

# 5. Quick practice model (a test only, not a real result)
for d in ["ciprofloxacin", "ceftriaxone", "gentamicin"]:
    y = labels[d].dropna().astype(int)
    Xd = X.loc[y.index]
    minority = y.value_counts().min() if y.nunique() == 2 else 0
    if minority < 10:
        print(f"\n{d}: skipped (smaller class has only {minority})")
        continue
    cv = StratifiedKFold(n_splits=min(5, minority), shuffle=True, random_state=0)
    rf = RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=0)
    pred = cross_val_predict(rf, Xd, y, cv=cv)
    print(f"\n{d}: n={len(y)}  balanced accuracy={balanced_accuracy_score(y, pred):.2f}  F1={f1_score(y, pred):.2f}")
    rf.fit(Xd, y)
    top = pd.Series(rf.feature_importances_, index=Xd.columns).sort_values(ascending=False).head(5)
    print("  top features:", ", ".join(top.index))
