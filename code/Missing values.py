import pandas as pd
path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/Exports'#'C:/Users/ariad/OneDrive/Desktop/Proyecto/Exports'
import glob
import os

# Load your metadata file (example: CSV)
# Replace 'metadata.csv' with your file path
metadata_files = sorted(glob.glob(path + "/datMeta.unionexon.*"+".csv"))
expression_files = sorted(glob.glob(path + "/datExpr.HTSC.unionexon.*"+"filtered.csv"))

def get_dataset_name(file_path):
    #namee = os.path.basename(file_path).replace("datExpr.HTSC.unionexon.", "").replace(".filtered.csv", "")
    namee = os.path.basename(file_path).replace("datMeta.unionexon.", "").replace(".csv", "").replace("datExpr.HTSC.unionexon.", "").replace(".filtered", "")
    if namee == "C" or namee == "CBL":
        return "Vermis" 
    elif namee == "F" or namee == "FRT":
        return "Frontal"    
    elif namee == "T" or namee == "TEM":
        return "Temporal"

for meta_file in metadata_files:
    dataset_name = get_dataset_name(meta_file)
    print("\n\n====== Processing Dataset: {} ======\n".format(dataset_name))
    
    metadata = pd.read_csv(meta_file, index_col=0)
    selected_columns = ['ASD.CTL', 'Age', 'BrainBank', 'Seizures', 'Sex', 'SeqBatch']
    new_meta = metadata[selected_columns].copy()
    # Count missing values per column
    missing_counts = new_meta.isnull().sum()

    # Display both counts and percentages
    missing_summary = pd.DataFrame({
        'Missing Values': missing_counts,
        'Percentage': (missing_counts / len(new_meta)) * 100
    })
    print("Total number of samples:", len(new_meta))
    print(missing_summary)

for expr_file in expression_files:
    dataset_name = get_dataset_name(expr_file)
    print("\n\n====== Processing Dataset: {} ======\n".format(dataset_name))
    
    expression_df = pd.read_csv(expr_file, index_col=0)
    # Count missing values per column
    missing_counts = expression_df.isnull().sum()

    # Display both counts and percentages
    missing_summary = pd.DataFrame({
        'Missing Values': missing_counts,
        'Percentage': (missing_counts / len(expression_df)) * 100
    })
    print("Total number of samples:", len(expression_df.T))
    print(missing_summary.sort_values(by='Missing Values', ascending=False))