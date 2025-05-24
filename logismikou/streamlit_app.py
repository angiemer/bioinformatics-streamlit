import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import streamlit as st

# --- Page config ---
st.set_page_config(page_title="🔬 scRNA-seq Explorer", layout="wide")
st.title("🔬 Single-cell RNA-seq Explorer")
st.markdown("Upload and explore single-cell RNA-seq `.h5ad` files interactively.")

# --- Tabs Layout ---
tab1, tab2, tab3 = st.tabs(["Data Analysis", "Differential Expression", "Team Info"])

with tab1:
    uploaded_file = st.file_uploader("📂 Upload `.h5ad` file (max 300MB)", type="h5ad")

    if uploaded_file:
        file_size_mb = len(uploaded_file.getvalue()) / (1024 * 1024)
        if file_size_mb > 300:
            st.error(f"❌ File too large ({file_size_mb:.1f} MB > 300 MB). Please upload a smaller file.")
        else:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
                tmp.write(uploaded_file.read())
                tmp.flush()
                adata = sc.read(tmp.name)

            st.success("File loaded successfully!")

            st.markdown("### 🧬 Metadata Preview")
            st.dataframe(adata.obs.head(), height=200)

            if "X_umap" in adata.obsm:
                st.markdown("### 🗺️ UMAP Plot")
                color_by = st.selectbox("Select metadata column to color UMAP:", options=adata.obs.columns)
                fig_umap = sc.pl.umap(adata, color=color_by, return_fig=True)
                fig_umap.set_size_inches(6, 4)
                st.pyplot(fig_umap)
            else:
                st.warning("⚠️ UMAP coordinates not found.")

            if "X_pca" in adata.obsm:
                st.markdown("### 📊 PCA Plot")
                fig_pca = sc.pl.pca(adata, color=color_by, return_fig=True)
                fig_pca.set_size_inches(6, 4)
                st.pyplot(fig_pca)

with tab2:
    if uploaded_file and "rank_genes_groups" in adata.uns:
        st.markdown("### Differential Gene Expression")

        deg_result = adata.uns["rank_genes_groups"]
        degs_df = pd.DataFrame({
            "Gene": deg_result["names"]["case"],
            "p-value": deg_result["pvals"]["case"],
            "Adjusted p-value": deg_result["pvals_adj"]["case"],
            "Log2 Fold Change": deg_result["logfoldchanges"]["case"],
        })

        with st.sidebar:
            st.markdown("### ⚙️ Volcano Plot Settings")
            pval_thresh = st.slider("p-value threshold", 0.001, 0.1, 0.05)
            logfc_thresh = st.slider("Log2 Fold Change threshold", 0.0, 2.0, 1.0)

        degs_df["-log10 p-value"] = -np.log10(degs_df["p-value"])
        degs_df["Status"] = "Not Significant"
        degs_df.loc[(degs_df["Log2 Fold Change"] > logfc_thresh) & (degs_df["p-value"] < pval_thresh), "Status"] = "Upregulated"
        degs_df.loc[(degs_df["Log2 Fold Change"] < -logfc_thresh) & (degs_df["p-value"] < pval_thresh), "Status"] = "Downregulated"

        st.markdown("#### Top 10 DE Genes")
        st.dataframe(degs_df.head(10), height=220)

        st.markdown("#### Volcano Plot")
        plt.figure(figsize=(7, 5))
        sns.scatterplot(
            data=degs_df,
            x="Log2 Fold Change",
            y="-log10 p-value",
            hue="Status",
            palette={"Upregulated": "#bb0c00", "Downregulated": "#00AFBB", "Not Significant": "grey"},
            alpha=0.7
        )
        plt.axhline(y=-np.log10(pval_thresh), color="gray", linestyle="--")
        plt.axvline(x=-logfc_thresh, color="gray", linestyle="--")
        plt.axvline(x=logfc_thresh, color="gray", linestyle="--")
        plt.title("Volcano Plot of DE Genes", fontsize=14)
        plt.tight_layout()
        st.pyplot(plt)

    elif uploaded_file:
        st.warning(" No DE results found in the file.")
    else:
        st.info(" Upload a `.h5ad` file first under the 'Analysis' tab.")

with tab3:
    st.markdown("### Project Team")
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("#### Members:")
        st.write("- Μερκούρη Αγγελική: Streamlit διεπαφή, usecase diagram, Latex")
        st.write("- Χρυσούλα Δήμα: streamlit διεπαφή, class diagram, Latex")

    with col2:
        st.markdown("#### Project Info:")
        st.write("**Course:** Τεχνολογία Λογισμικού")
        st.write("**Year:** 2025")
        st.write("**Institution:** Ιόνιο Πανεπιστήμιο")