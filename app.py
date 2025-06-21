import streamlit as st
import pandas as pd
import os
import io

st.title("🔬 SignatureId Matcher & CSV Exporter🧑‍💻")

st.write("""
Upload your metadata file and data file.  
The app will match rows by `SignatureId` and generate separate downloadable `.csv` files.
""")

# === Upload files ===
meta_file = st.file_uploader("📁 Upload Metadata File (.csv, .tsv, .xls)", type=['csv', 'tsv', 'xls'])
data_file = st.file_uploader("📁 Upload Data File (.csv, .tsv, .xls)", type=['csv', 'tsv', 'xls'])

if meta_file and data_file:
    # === Load metadata ===
    try:
         meta_df = pd.read_csv(meta_file, sep='\t')
         filtered_metaDate=meta_df.drop_duplicates(subset='Tissue', keep='first')
    except Exception as e:
        st.error(f"❌ Failed to read metadata: {e}")
        st.stop()

    # === Load data file ===
    try:
        data_df = pd.read_csv(data_file, sep='\t')
    except Exception as e:
        st.error(f"❌ Failed to read data file: {e}")
        st.stop()
        
    # === Validate columns ===
    required_meta_cols = {"SignatureId", "Perturbagen", "Tissue", "CellLine"}
    result = filtered_metaDate[["SignatureId", "Perturbagen","Tissue", "CellLine"]]
    perturbagen = result['Perturbagen'].iloc[0]
    if not required_meta_cols.issubset(filtered_metaDate.columns):
        st.error(f"❌ Metadata file must contain columns: {required_meta_cols}")
        st.stop()

    if 'signatureID' not in data_df.columns:
        st.error("❌ Data file must contain column: 'signatureID'")
        st.stop()

    st.success("✅ Files loaded successfully!")

    # Output zip buffer
    from zipfile import ZipFile
    from io import BytesIO
    zip_buffer = BytesIO()

    with ZipFile(zip_buffer, "a") as zipf:
        for _, row in filtered_metaDate.iterrows():
            sig_id = str(row['SignatureId']).strip()
            perturbagen = str(row['Perturbagen']).strip().replace(" ", "_")
            tissue = str(row['Tissue']).strip().replace(" ", "_")
            cell_line = str(row['CellLine']).strip().replace(" ", "_")

            matched = data_df[data_df['signatureID'] == sig_id]

            if matched.empty:
                continue

            filename = f"{perturbagen} {tissue} {cell_line} - {perturbagen} {tissue} {cell_line}.csv"
            csv_bytes = matched.to_csv(index=False).encode('utf-8')
            zipf.writestr(filename, csv_bytes)

    st.success("✅ All matching rows processed!")

    # Download the ZIP
    st.download_button(
        label="⬇️ Download All Results (ZIP)",
        data=zip_buffer.getvalue(),
        file_name=f"\{perturbagen}_Results.zip",
        mime="application/zip"
    )
