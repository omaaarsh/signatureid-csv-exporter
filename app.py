import streamlit as st
import pandas as pd
import numpy as np
import os
import tempfile
import zipfile
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
import gseapy as gp
from gseapy.plot import barplot
import io

# --- STREAMLIT CONFIG ---
st.set_page_config(layout="wide")
st.title("🧬 Batch Gene Signature Analysis from Zipped Drug Treatments")
st.markdown("""
Upload a **ZIP file** that contains multiple **sub-ZIPs**, each representing a **drug**.  
Each sub-ZIP should contain CSVs for multiple **cell lines** treated with the same drug.

✅ Each CSV must include: `ID_geneid`, `Name_GeneSymbol`, `Value_LogDiffExp`, `Significance_pvalue`  
✅ Naming format: `DrugName cell_line_name.xls - DrugName cell_line_name.xls.csv`
""")

# --- PARAMETERS ---
log2fc_threshold = st.sidebar.slider("log₂ Fold Change Threshold", min_value=0.5, max_value=3.0, value=1.0, step=0.1)
consistency_threshold = st.sidebar.slider("Cell Line Consistency Threshold (%)", min_value=30, max_value=100, value=50, step=10)
pvalue_threshold = st.sidebar.slider("Significance P-Value Threshold", min_value=0.0, max_value=0.2, value=0.2, step=0.01)


# --- ANALYSIS FUNCTION ---
def run_analysis(cell_line_data: dict, drug_name: str):
    report_lines = []
    report_lines.append(f"📊 Summary for {drug_name}")

    # --- Merge ---
    gene_dfs = []
    for cell_line, df in cell_line_data.items():
        temp = df[['ID_geneid', 'Name_GeneSymbol', 'Value_LogDiffExp']].copy()
        temp = temp.rename(columns={'Value_LogDiffExp': cell_line})
        gene_dfs.append(temp)

    merged_df = reduce(lambda left, right: pd.merge(left, right, on=['ID_geneid', 'Name_GeneSymbol'], how='outer'), gene_dfs)
    merged_df.set_index(['ID_geneid', 'Name_GeneSymbol'], inplace=True)

    # --- Filter ---
    cell_line_count = len(cell_line_data)
    min_consistent = int(cell_line_count * (consistency_threshold / 100))

    upregulated_mask = (merged_df > log2fc_threshold).sum(axis=1)
    downregulated_mask = (merged_df < -log2fc_threshold).sum(axis=1)

    consistently_up_genes = upregulated_mask[upregulated_mask >= min_consistent].index
    consistently_down_genes = downregulated_mask[downregulated_mask >= min_consistent].index

    report_lines.append(f"Total genes in matrix: {merged_df.shape[0]}")
    report_lines.append(f"Genes upregulated in ≥{consistency_threshold}% cell lines: {len(consistently_up_genes)}")
    report_lines.append(f"Genes downregulated in ≥{consistency_threshold}% cell lines: {len(consistently_down_genes)}\n")

    st.subheader("📊 Summary")
    st.write(f"**Total genes in matrix:** {merged_df.shape[0]}")
    st.write(f"**Genes upregulated in ≥{consistency_threshold}% cell lines:** {len(consistently_up_genes)}")
    st.write(f"**Genes downregulated in ≥{consistency_threshold}% cell lines:** {len(consistently_down_genes)}")

    # --- Heatmap ---
    consistent_genes = consistently_up_genes.union(consistently_down_genes)
    heatmap_data = merged_df.loc[consistent_genes].fillna(0)

    if not heatmap_data.empty:
        row_linkage = linkage(pdist(heatmap_data, metric='correlation'), method='average')
        col_linkage = linkage(pdist(heatmap_data.T, metric='correlation'), method='average')

        st.subheader("🧯 Hierarchical Clustering Heatmap")
        fig = sns.clustermap(
            heatmap_data,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            cmap='RdBu_r',
            center=0,
            linewidths=0.5,
            figsize=(10, 10)
        )
        st.pyplot(fig.fig)

        # Save heatmap
        heatmap_file_path = os.path.join(tempfile.gettempdir(), f"{drug_name}_heatmap.png")
        fig.savefig(heatmap_file_path)

        with open(heatmap_file_path, "rb") as img_file:
            st.download_button(
                label=f"📸 Download Heatmap for {drug_name}",
                data=img_file,
                file_name=f"{drug_name}_heatmap.png",
                mime="image/png"
            )
    else:
        st.warning("⚠️ No consistent genes for heatmap.")
        report_lines.append("No consistent genes for heatmap.\n")

    # --- Enrichment ---
    st.subheader("🔬 Pathway Enrichment")

    def run_enrichment(gene_set, label):
        report_lines.append(f"Top enriched pathways in {label} genes:\n")
        if len(gene_set) == 0:
            report_lines.append(f"- No genes available for enrichment in {label} set.\n")
            return
        symbols = [gene[1] for gene in gene_set]
        try:
            enr = gp.enrichr(
                gene_list=symbols,
                gene_sets=["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"],
                organism='Human',
                outdir=None,
                cutoff=0.05
            )
            if enr.results.empty:
                st.info(f"No significant enrichment found for {label}.")
                report_lines.append(f"- No significant enrichment found for {label}.\n")
            else:
                st.write(f"**Top enriched pathways in {label} genes:**")
                st.dataframe(enr.results.head(10))
                fig, ax = plt.subplots(figsize=(10, 6))
                barplot(enr.results.sort_values('Adjusted P-value').head(10), title=f"{label} Enrichment", ax=ax)
                st.pyplot(fig)

                top = enr.results.sort_values('Adjusted P-value').head(10)
                report_lines.append(top[['Term', 'Adjusted P-value']].to_string(index=False))
                report_lines.append("\n")
        except Exception as e:
            report_lines.append(f"⚠️ Error during enrichment for {label}: {str(e)}\n")
            st.error(f"Error during enrichment for {label}: {e}")

    run_enrichment(consistently_up_genes, "Upregulated")
    run_enrichment(consistently_down_genes, "Downregulated")

    # --- Download Report ---
    report_text = "\n".join(report_lines)
    report_bytes = io.BytesIO(report_text.encode("utf-8"))
    st.download_button(
        label=f"📝 Download Report for {drug_name}",
        data=report_bytes,
        file_name=f"{drug_name}_report.txt",
        mime="text/plain"
    )

    # --- Download Merged CSV ---
    csv_download = merged_df.reset_index().to_csv(index=False).encode('utf-8')
    st.download_button(
        label=f"📥 Download Merged Gene Expression Data for {drug_name}",
        data=csv_download,
        file_name=f"{drug_name}_expression_matrix.csv",
        mime='text/csv'
    )


# --- ZIP UPLOAD ---
st.subheader("📁 Upload ZIP of Drug ZIPs")
uploaded_main_zip = st.file_uploader("Upload main ZIP file (contains drug ZIPs)", type="zip")

if uploaded_main_zip:
    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = os.path.join(tmpdir, "main.zip")
        with open(zip_path, "wb") as f:
            f.write(uploaded_main_zip.read())

        with zipfile.ZipFile(zip_path, 'r') as zipf:
            zipf.extractall(tmpdir)

        subzips = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.endswith(".zip")]

        if not subzips:
            st.error("❌ No sub-zip files found in uploaded archive.")
        else:
            for subzip_path in subzips:
                drug_name = os.path.basename(subzip_path).replace(".zip", "")
                st.markdown(f"---\n### 💊 Processing Drug: `{drug_name}`")

                extract_path = os.path.join(tmpdir, drug_name)
                with zipfile.ZipFile(subzip_path, 'r') as subzip:
                    subzip.extractall(extract_path)

                csv_files = [os.path.join(extract_path, f) for f in os.listdir(extract_path) if f.endswith(".csv")]
                cell_line_data = {}

                for csv_path in csv_files:
                    try:
                        file_name = os.path.basename(csv_path)
                        if " - " in file_name:
                            cell_line = file_name.split(" - ")[0].replace(drug_name, "").replace(".xls", "").strip()
                            df = pd.read_csv(csv_path)
                            df = df[df['Significance_pvalue'] <= pvalue_threshold]
                            cell_line_data[cell_line] = df
                    except Exception as e:
                        st.warning(f"⚠️ Failed to load {file_name}: {e}")

                if cell_line_data:
                    run_analysis(cell_line_data, drug_name)
                else:
                    st.warning(f"No valid cell line data found in `{drug_name}` ZIP.")
