# Importing libraries
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import pandas as pd
import os
#path = 'D:/Proyect/Exports'
#path = os.path.join(os.path.dirname(__file__), 'Exports')  
#path = 'C:/Users/ariad/OneDrive/Desktop/Proyect WCQN/Exports2'#'C:/Users/ariad/OneDrive/Desktop/Proyecto/Exports'
path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/Exports'  # Adjust this path as needed
import glob


# Load your dataset
metadata_files = sorted(glob.glob(path + "/datMeta.*"+".csv"))


def get_dataset_name(file_path):
    namee = os.path.basename(file_path).replace("datMeta.unionexon.", "").replace(".csv", "")
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


# ...existing code...

for meta_file in metadata_files:
    dataset_name = get_dataset_name(meta_file)
    print(dataset_name)
    metadata = pd.read_csv(meta_file, index_col=0)
    metadata = encode_categorical_metadata(metadata)

    # Use only relevant columns
    features = ['ASD.CTL', 'Age', 'BrainBank', 'Sex']
    # Ensure no missing values in features for prediction
    metadata = metadata.dropna(subset=features, how='any')

    # Split into known and unknown Seizures
    known = metadata[metadata['Seizures'].notna()]
    unknown = metadata[metadata['Seizures'].isna()]

    X_train = known[features]
    y_train = known['Seizures'].astype(int)  # Make sure it's int for classification

    # Train the model
    rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_model.fit(X_train, y_train)

    # Predict for missing Seizures
    if not unknown.empty:
        X_pred = unknown[features]
        y_pred = rf_model.predict(X_pred)
        # Save predictions
        imputed = unknown.copy()
        imputed['Seizures_predicted'] = y_pred
        # Convert numeric predictions back to string labels if needed
        label_map = {0: "No", 1: "Yes"}
        imputed_str = [label_map[val] for val in imputed['Seizures_predicted']]
        with open(os.path.join(output_path, f"{dataset_name}_rf_imputed_seizures.csv"), "w") as f:
            f.write(",".join(imputed_str))
        print(f"Imputed {len(y_pred)} missing Seizures for {dataset_name}")
    else:
        print(f"No missing Seizures to impute for {dataset_name}")

    # Optionally, evaluate on known data (train set)
    y_train_pred = rf_model.predict(X_train)
    acc = accuracy_score(y_train, y_train_pred)
    report_text = classification_report(y_train, y_train_pred)
    print("Training Accuracy:", acc)
    print("Training Classification Report:\n", report_text)

    # Export accuracy and report
    with open(os.path.join(output_path, f"{dataset_name}_rf_accuracy.txt"), "w") as f:
        f.write(f"Training Accuracy: {acc}\n")
    with open(os.path.join(output_path, f"{dataset_name}_rf_classification_report.txt"), "w") as f:
        f.write(report_text)
# ...existing code...
    # Now, update the original metadata with imputed values
    if not unknown.empty:
        metadata.loc[imputed.index, 'Seizures'] = imputed['Seizures_predicted']
    # Save the updated metadata
    metadata.to_csv(os.path.join(output_path, f"{dataset_name}_updated_metadata.csv"))
    print(f"Updated metadata saved for {dataset_name}")