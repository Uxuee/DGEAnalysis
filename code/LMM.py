import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import os
#path = 'D:/Proyect/Exports'
#path = os.path.join(os.path.dirname(__file__), 'Exports')  
#path = 'C:/Users/ariad/OneDrive/Desktop/Proyect WCQN/Exports2'#'C:/Users/ariad/OneDrive/Desktop/Proyecto/Exports'
path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/Exports'  # Adjust this path as needed
import glob
import pickle
np.random.seed(123)  # for reproducibility

# Get all relevant expression/metadata pairs from folder
expression_files = sorted(glob.glob(path + "/datExpr.*"+".csv"))
metadata_files = sorted(glob.glob(path + "/datMeta.*"+".csv"))

def get_dataset_name(file_path):
    namee = os.path.basename(file_path).replace("datExpr.HTSC.unionexon.", "").replace(".filtered.csv", "")
    #namee = os.path.basename(file_path).replace("datExpr_", "").replace("_adjusted.csv", "")
    if namee == "C" or namee == "CBL":
        return "Vermis" 
    elif namee == "F" or namee == "FRT":
        return "Frontal"    
    elif namee == "T" or namee == "TEM":
        return "Temporal"

def encode_categorical_metadata(meta):
    """
    Encodes categorical metadata for machine learning:
    - Applies binary label encoding to selected columns.
    - Performs one-hot encoding for multiclass columns.
    - Simplifies comorbidity and cause of death.
    
    Parameters:
        meta (pd.DataFrame): The original metadata DataFrame.
    
    Returns:
        pd.DataFrame: Encoded metadata with categorical variables prepared for ML.
    """
    # Copy only the relevant columns to avoid modifying the original metadata
    encoded_meta = meta.copy()

    # 1. Binary label encoding
    binary_map = {
        'Sex': {'M': 0, 'F': 1},
        'BrainBank': {'ATP': 0, 'NICHD': 1},
        'ASD.CTL': {'CTL': 0, 'ASD': 1},
        'Seizures': {'No': 0, 'Yes': 1},
        'Pyschiatric.Medications': {'No': 0, 'Yes': 1},
        'DeathCategory': {'Sudden': 0,'Prolonged': 1}
    }
    for col, mapping in binary_map.items():
        if col in encoded_meta.columns:
            encoded_meta[col] = encoded_meta[col].map(mapping)

    return encoded_meta

print("Script has started...")

output_path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/LMMResults' 


for expr_file, meta_file in zip(expression_files, metadata_files):
    dataset_name = get_dataset_name(expr_file)
    print(dataset_name)

for expr_file, meta_file in zip(expression_files, metadata_files):
    dataset_name = get_dataset_name(expr_file)
    print("\n\n====== Processing Dataset: {} ======\n".format(dataset_name))
    
    expression_df = pd.read_csv(expr_file, index_col=0)
    metadata = pd.read_csv(meta_file, index_col=0)
    metadata['Age2'] = metadata['Age'] ** 2

    # Transpose expression: samples x genes
    expression_df = expression_df.T

    # Ensure alignment: keep only shared samples
    common_samples = metadata.index.intersection(expression_df.index)
    expression_df = expression_df.loc[common_samples]
    metadata = metadata.loc[common_samples]
    metadata= encode_categorical_metadata(metadata)
    metadata['Diagnosis'] = metadata['ASD.CTL']

    print("Metadata shape: {}".format(metadata.shape))
    print("Expression shape: {}".format(expression_df.shape))
    print("Number of samples: {}".format(len(common_samples)))

    print("Processing region: {}".format(expr_file))

    results = []

    for gene in expression_df.columns:
        gene_expr = expression_df[gene]
        
        # Combine gene expression and metadata
        df = metadata[['Diagnosis', 'Age', 'Age2', 'BrainBank', 'Seizures', 'Sex', 'SeqBatch']].copy()
        df['Y'] = gene_expr
        
        # Impute missing values for numeric columns with mean, categorical with mode
        for col in ['Y', 'Diagnosis', 'Age', 'Age2', 'BrainBank', 'Seizures', 'Sex', 'SeqBatch']:
            if df[col].dtype.kind in 'biufc':  # Numeric types
                df[col] = df[col].fillna(df[col].mean())
            else:  # Categorical types
                df[col] = df[col].fillna(df[col].mode()[0])

        # Now you don't need to drop NA rows
        # df = df.dropna(subset=['Y', 'Diagnosis', 'Age', 'Age2', 'BrainBank', 'Seizures', 'Sex', 'SeqBatch'])        

        if len(df) < 20:
            continue  # Skip if not enough samples to fit model
        
        ## HERE IS THE LMM
        try:
            model = smf.mixedlm(
                "Y ~ 1 + Diagnosis + Age + Age2 + BrainBank + Seizures + Sex",
                df,
                groups=df["SeqBatch"]
            ).fit()

            
            # Get the coefficient and p-value for Diagnosis (ASD vs CTL)
            coef = model.params['Diagnosis']
            pval = model.pvalues['Diagnosis']
                        
            # Standardize the effect size: divide by the standard deviation of the gene's expression
            std_effect = coef / gene_expr.std()

            results.append({'Gene': gene, 'EffectSize': coef, 'StdEffect': std_effect, 'P-Value': pval})
            
        except Exception as e:
            print("Skipped {} due to: {}".format(gene, e))

    results_df = pd.DataFrame(results)
    results_df = results_df.dropna()

    if results_df.empty:
        print("No valid models for dataset: {}".format(dataset_name))
        continue 

    # Apply Benjamini/Yekutieli
    results_df['FDR'], results_df['pBH'], _, _ = multipletests(results_df['P-Value'], method='fdr_by')

    # Select significant genes
    sig_genes = results_df[(results_df['pBH'] < 0.05) & (abs(results_df['StdEffect']) > 0.8)]

    # Save results
    results_df.to_csv(os.path.join(output_path,'{}_LMM.csv'.format(dataset_name+"_" )), index=True)
    sig_genes.to_csv(os.path.join(output_path,'{}_LMM_SG.csv'.format(dataset_name +"_")), index=True)
    print("Results saved for {}".format(dataset_name +"_" ))
    print("Processing completed for dataset: {}".format(dataset_name))

    # Ensure -log10(p-values) and volcano annotations
    results_df['-log10(pBH)'] = -np.log10(results_df['pBH'])
    results_df['Significant'] = (results_df['pBH'] < 0.05) & (abs(results_df['StdEffect']) > 0.8)

    results_df = results_df.replace([np.inf, -np.inf], np.nan).dropna()
    # Optional: Set plotting limits for clarity
    xlim = (-3, 3)
    ylim = (0, results_df['-log10(pBH)'].max() + 2)

    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=results_df,
        x='EffectSize', y='-log10(pBH)',
        hue='Significant',
        palette={True: 'red', False: 'grey'},
        edgecolor=None,
        alpha=0.7
    )

    plt.title("Volcano Plot: {}".format(dataset_name))
    plt.xlabel("Effect Size (β for Diagnosis)")
    plt.ylabel("-log10 Adjusted p-value (FDR)")
    plt.axvline(x=-0.6, color='blue', linestyle='--', linewidth=1)
    plt.axvline(x=0.6, color='blue', linestyle='--', linewidth=1)
    plt.axhline(y=-np.log10(0.05), color='green', linestyle='--', linewidth=1)

    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend(title='Significant', loc='upper right')
    plt.tight_layout()

    # Save or show
    plt.savefig(os.path.join(output_path, "{}_VolcanoPlot.png".format(dataset_name)), dpi=300)
    plt.close()