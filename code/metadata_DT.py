print("Loading packages...")
import os
path = 'C:/Users/ariad/OneDrive/Desktop/Proyect/Exports'
#path= 'C:/Users/palom/Desktop/Proyect/Exports'
#path = os.path.join(os.path.dirname(__file__), 'Exports')   
import glob
import pandas as pd
import numpy as np
from itertools import repeat
np.random.seed(123)  # for reproducibility

print("Script has started...")


# Get all relevant expression/metadata pairs from folder
expression_files = sorted(glob.glob(os.path.join(path, "*DT*.csv")))
metadata_files = sorted(glob.glob(path + "/*_Metadata.csv"))
metadata_files=sorted([x for item in metadata_files for x in repeat(item, 3)])

def get_dataset_name(file_path):
    return os.path.basename(file_path).replace(".csv", "")

for expr_file, meta_file in zip(expression_files, metadata_files):
    dataset_name = get_dataset_name(expr_file)
    print(dataset_name)

categorical_vars = ["Sex", "SeqBatch", "BrainBank", "ASD.CTL", 
  "Seizures", "Pyschiatric.Medications", "Comorbidity.notes..other.than.seizures.", 
  "Primary.Cause.of.Death"]

continuous_vars = ["Age", "PMI", "Brain.Weight", "pH", "IQ",
  "TotalReadsInBam", "After.rmdup.samtools", "Num.dup.samtools", 
  "Unique.Reads.samtools", "PropExonicReads.HTSC",
  "TotalReads.picard", "Aligned.Reads.picard", "HQ.Aligned.Reads.picard",
  "PF.All.Bases.picard", "Coding.Bases.picard", "UTR.Bases.picard", 
  "Intronic.Bases.picard", "Intergenic.bases.picard",
  "Median.CV.Coverage.picard", "Median.5prime.Bias.picard", 
  "Median.3prime.Bias.picard", "Median.5to3prime.Bias.picard",
  "AT.Dropout.picard", "GC.Dropout.picard"
]


for expr_file, meta_file in zip(expression_files, metadata_files):
    dataset_name = get_dataset_name(expr_file)
    print("\n\n====== Processing Dataset: {} ======\n".format(dataset_name))

    dataset = pd.read_csv(expr_file, index_col=0)
    meta = pd.read_csv(meta_file, index_col=0)

    sample_ids = dataset.columns

    cleaned_meta = meta[categorical_vars + continuous_vars+["BrainID"]+['Detailed.Diagnosis']].copy()
    cleaned_meta.drop_duplicates(subset=['BrainID'], keep='first', inplace=True)
    cleaned_meta.set_index('BrainID', inplace=True)
    cleaned_meta = cleaned_meta.loc[sample_ids]

    
    output_path = os.path.join(path,'{}_Metadata.csv'.format(dataset_name)) 
    cleaned_meta.to_csv(output_path, index=True)
    
    print("Saved Continuous Results for {}".format(dataset_name))
    print("Results saved to: {}".format(output_path))
