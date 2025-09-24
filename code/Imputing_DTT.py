# Importing libraries
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import os
#path = 'D:/Proyect/Exports'
#path = os.path.join(os.path.dirname(__file__), 'Exports')  
#path = 'C:/Users/ariad/OneDrive/Desktop/Proyect WCQN/Exports2'#'C:/Users/ariad/OneDrive/Desktop/Proyecto/Exports'
path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/Exports'  # Adjust this path as needed
import glob
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier



# Load your dataset
metadata_files = sorted(glob.glob(path + "/DT_datMeta_*"+".csv"))
expression_files = sorted(glob.glob(path + "/DT_datExpr_*"+".csv"))


def get_dataset_name(file_path):
    namee = os.path.basename(file_path).replace("DT_datMeta_", "").replace("_unionexon.csv", "")
    #namee = os.path.basename(file_path).replace("datExpr_", "").replace("_adjusted.csv", "")
    if namee == "C_T" or namee == "CBL":
        return "VermisTemporal" 
    elif namee == "C_F" or namee == "FRT":
        return "VermisFrontal"    
    elif namee == "F_T" or namee == "TEM":
        return "FrontalTemporal"

print("Script has started...")

output_path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/imputed' 

# ...existing code...

for meta_file,expr_file in zip(metadata_files,expression_files):
    dataset_name = get_dataset_name(meta_file)
    print(dataset_name)
    metadata = pd.read_csv(meta_file, index_col=0)
    expression_df = pd.read_csv(expr_file, index_col=0)  # genes x samples

    # Transpose expression so rows = samples, columns = genes
    expression_df = expression_df.T

    # Encode Seizures column ("Yes"/"No" → 1/0)
    metadata['Seizures'] = metadata['Seizures'].map({'Yes': 1, 'No': 0})

    # Align samples
    common_samples = metadata.index.intersection(expression_df.index)
    expression_df = expression_df.loc[common_samples]
    metadata = metadata.loc[common_samples]

    # Split into known and unknown Seizures
    known = metadata[metadata['Seizures'].notna()]
    unknown = metadata[metadata['Seizures'].isna()]

    X_train =  expression_df.loc[known.index]
    y_train = known['Seizures'].astype(int)  # Make sure it's int for classification

    # Train the model
    classifiers = {
        "RandomForest": RandomForestClassifier(n_estimators=250, random_state=42),
        "LogisticRegression": LogisticRegression(max_iter=1000, random_state=42),
        "SVC": SVC(kernel='rbf', probability=True, random_state=42),
        "KNeighbors": KNeighborsClassifier(n_neighbors=5)
    }

    best_score = 0
    best_name = None
    best_model = None

    for name, clf in classifiers.items():
        scores = cross_val_score(clf, X_train, y_train, cv=5, scoring='accuracy')
        mean_score = scores.mean()
        print(f"{name} cross-validated accuracy: {mean_score:.4f}")
        if mean_score > best_score:
            best_score = mean_score
            best_name = name
            best_model = clf

    print(f"Best classifier: {best_name} with accuracy {best_score:.4f}")
    rf_model = best_model
    rf_model.fit(X_train, y_train)

    # Predict for missing Seizures
    if not unknown.empty:
        X_pred = expression_df.loc[unknown.index]
        y_pred = rf_model.predict(X_pred)
        # Save predictions
        imputed = unknown.copy()
        imputed['Seizures'] = y_pred

    # Optionally, evaluate on known data (train set)
    cv_scores = cross_val_score(rf_model, X_train, y_train, cv=5, scoring='accuracy')
    print("Cross-validated accuracy:", cv_scores.mean())
    acc = cv_scores.mean()

    # Export accuracy and report
    with open(os.path.join(output_path, f"DT_{dataset_name}_rf_accuracy.txt"), "w") as f:
        f.write(f"Training Accuracy: {acc}\n")

    if not unknown.empty:
        metadata.loc[imputed.index, 'Seizures'] = imputed['Seizures']
    # Save the updated metadata
    metadata.to_csv(os.path.join(output_path, f"DT_datMeta_{dataset_name}_imp.csv"), index=True)
    print(f"Updated metadata saved for {dataset_name}")


