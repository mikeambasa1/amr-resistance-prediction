from pathlib import Path
import joblib
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "dashboard" / "data"
MODELS = ROOT / "models"
DRUGS = ["ceftriaxone", "ciprofloxacin", "gentamicin"]

st.set_page_config(page_title="K. pneumoniae AMR", layout="wide")


@st.cache_data
def load_tables():
    t = {n: pd.read_csv(DATA / f"{n}.csv") for n in
         ["prevalence", "gene_frequency", "gene_weights",
          "model_results", "model_results_grouped"]}
    X = pd.read_csv(DATA / "gene_matrix_filtered.csv", index_col=0, dtype={"genome_id": str})
    Y = pd.read_csv(DATA / "labels_clean.csv", index_col=0, dtype={"genome_id": str})
    X.index = X.index.astype(str)
    Y.index = Y.index.astype(str)
    return t, X, Y


@st.cache_resource
def load_model(drug):
    return joblib.load(MODELS / f"{drug}_logreg.joblib")


t, X, Y = load_tables()

st.title("Klebsiella pneumoniae antibiotic resistance")
st.caption(f"{len(X):,} genomes from BV-BRC | AMRFinderPlus gene calls | "
           "three antibiotics | logistic regression models")

tab1, tab2, tab3 = st.tabs(["Overview", "Gene explorer", "Predict"])

with tab1:
    st.subheader("Resistance prevalence")
    p = t["prevalence"]
    st.bar_chart(p.set_index("drug")[["resistant", "susceptible"]])
    st.dataframe(p, hide_index=True)
    st.subheader("Model performance (5-fold cross-validation)")
    res = t["model_results"][["drug", "model", "roc_auc", "bal_acc", "avg_prec", "f1"]]
    st.dataframe(res.round(3), hide_index=True)
    st.subheader("Lineage-aware cross-validation")
    st.dataframe(t["model_results_grouped"], hide_index=True)
    st.caption("The dataset is enriched for clinical isolates, so prevalence here "
               "does not represent K. pneumoniae in the general population.")

with tab2:
    drug = st.selectbox("Antibiotic", DRUGS, key="explore_drug")
    n = st.slider("Number of genes", 5, 25, 12)
    w = (t["gene_weights"].query("drug == @drug")
         .sort_values("weight", ascending=False).head(n))
    st.subheader("Genes pushing the prediction towards resistant")
    st.bar_chart(w.set_index("gene")["weight"])
    f = t["gene_frequency"].query("drug == @drug").merge(
        w[["gene", "weight"]], on="gene")
    st.dataframe(f.sort_values("weight", ascending=False), hide_index=True)
    st.caption("Weights show association with resistance in the model, "
               "not proof that a gene causes it.")

with tab3:
    st.subheader("Predict resistance from a gene profile")
    start = st.radio("Start from", ["Empty profile", "A real genome from the dataset"],
                     horizontal=True)
    default, gid = [], "none"
    if start != "Empty profile":
        gid = st.selectbox("Genome ID", list(X.index))
        default = list(X.columns[X.loc[gid] == 1])
        st.write("Observed phenotype (1 = resistant):",
                 Y.loc[gid].astype(int).to_dict())
    genes = st.multiselect("Genes present (everything else is treated as absent)",
                           list(X.columns), default=default, key=f"genes_{gid}")
    v = pd.DataFrame(0, index=[0], columns=X.columns)
    for g in genes:
        v.loc[0, g] = 1
    cols = st.columns(3)
    for c, d in zip(cols, DRUGS):
        b = load_model(d)
        score = float(b["model"].predict_proba(v[b["genes"]])[0, 1])
        c.metric(d.capitalize(), f"{score:.2f}",
                 "predicted resistant" if score > 0.5 else "predicted susceptible",
                 delta_color="off")
    st.caption("Scores are class-weighted model scores, not calibrated probabilities "
               "(above 0.5 = predicted resistant). Real genomes carry many "
               "co-occurring genes, so a few ticked genes can score lower than "
               "the same genes in a full profile. For research and teaching only, "
               "not for clinical decisions.")
