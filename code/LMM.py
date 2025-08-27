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
path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/Exports'  # Adjust this path as needed
import glob
import pickle
np.random.seed(123)  # for reproducibility
import scipy.stats as stats
import statsmodels.api as sm

# Get all relevant expression/metadata pairs from folder
expression_files = sorted(glob.glob(path + "/datExpr.*"+".csv"))
metadata_files = sorted(glob.glob(path + "/datMeta.*"+".imputed.csv"))

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
    # Centering and squaring Age
    metadata['Age_c'] = (metadata['Age'] - metadata['Age'].mean())/(2*metadata['Age'].std())  # Centering Age
    metadata['Age2']  = (metadata['Age_c']) ** 2   # quadratic term

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
        df = metadata[['Diagnosis', 'Age_c', 'Age2', 'BrainBank', 'Sex', 'SeqBatch','Seizures']].copy()
        df['Y'] = gene_expr

        if len(df) < 20:
            continue  # Skip if not enough samples to fit model
        
        ## HERE IS THE LMM
        try:
            model = smf.mixedlm(
                "Y ~ 1 + Diagnosis + Age_c + Age2 + BrainBank + Seizures + Sex",
                df,
                groups=df["SeqBatch"]
            ).fit()

            # Only plot and export for the first 5 genes
            if len(results) < 5:
                # KDE plot of residuals
                fig_kde = plt.figure(figsize=(10, 6))
                ax_kde = sns.kdeplot(model.resid, fill=True, linewidth=1, color='blue', label='Residuals KDE')
                x = np.linspace(model.resid.min(), model.resid.max(), 100)
                ax_kde.plot(x, stats.norm.pdf(x, model.resid.mean(), model.resid.std()), color='black', linestyle='--', label='Normal PDF')
                ax_kde.set_title(f"KDE Plot of Model Residuals: {gene}")
                ax_kde.set_xlabel("Residuals")
                ax_kde.legend()
                plt.tight_layout()
                kde_path = os.path.join(output_path, f"{dataset_name}_{gene}_resid_KDE.png")
                plt.savefig(kde_path, dpi=200)
                plt.close(fig_kde)

                # Q-Q plot
                fig_qq = plt.figure(figsize=(8, 8))
                ax_qq = fig_qq.add_subplot(111)
                sm.qqplot(model.resid, dist=stats.norm, line='s', ax=ax_qq)
                ax_qq.set_title(f"Q-Q Plot of Model Residuals: {gene}")
                plt.tight_layout()
                qq_path = os.path.join(output_path, f"{dataset_name}_{gene}_resid_QQ.png")
                plt.savefig(qq_path, dpi=200)
                plt.close(fig_qq)

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
        x='StdEffect', y='-log10(pBH)',
        hue='Significant',
        palette={True: 'red', False: 'grey'},
        edgecolor=None,
        alpha=0.7
    )

    plt.title(f"Volcano Plot: {dataset_name}")
    plt.xlabel("Standardized Effect Size (β for Diagnosis)")
    plt.ylabel("-log10 Adjusted p-value (FDR)")
    plt.axvline(x=-0.8, color='blue', linestyle='--', linewidth=1)  # threshold for effect size
    plt.axvline(x=0.8, color='blue', linestyle='--', linewidth=1)
    plt.axhline(y=-np.log10(0.05), color='green', linestyle='--', linewidth=1)

    plt.xlim(-3, 3)
    plt.ylim(0, results_df['-log10(pBH)'].max() + 2)
    plt.yscale('symlog', linthresh=0.01)  # Use symmetric log scale 
    plt.legend(title='Significant', loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f"{dataset_name}_VolcanoPlot.png"), dpi=300)
    plt.close()