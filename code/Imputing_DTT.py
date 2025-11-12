# Importing libraries
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import os
#path = 'D:/Proyect/Exports'
#path = os.path.join(os.path.dirname(__file__), 'Exports')  
#path = 'C:/Users/ariad/OneDrive/Desktop/Proyect WCQN/Exports2'#'C:/Users/ariad/OneDrive/Desktop/Proyecto/Exports'
path = 'C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/Exports'  # Adjust this path as needed
import glob
from sklearn.model_selection import cross_val_score, GridSearchCV
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

    # Define hyperparameter grids for each classifier
    param_grids = {
        "RandomForest": {
            'n_estimators': [100, 250],
            'max_depth': [None, 10],
            'min_samples_split': [2, 5]
        },
        "LogisticRegression": {
            'C': [0.1, 1.0, 10.0],
            'max_iter': [1000]
        },
        "SVC": {
            'C': [0.1, 1.0, 10.0],
            'kernel': ['rbf'],
            'probability': [True]
        },
        "KNeighbors": {
            'n_neighbors': [3, 5, 7],
            'weights': ['uniform', 'distance']
        }
    }

    # Define classifiers
    classifiers = {
        "RandomForest": RandomForestClassifier(random_state=42),
        "LogisticRegression": LogisticRegression(random_state=42),
        "SVC": SVC(random_state=42),
        "KNeighbors": KNeighborsClassifier()
    }

    # Define multiple scoring metrics
    scoring = {
        'accuracy': 'accuracy',
        'precision': 'precision',
        'recall': 'recall',
        'f1': 'f1'
    }

    best_score = 0
    best_name = None
    best_model = None

    for name, clf in classifiers.items():
        # Use GridSearchCV with multi-metric scoring
        # Set refit='accuracy' to specify which metric to use for selecting best model
        grid_search = GridSearchCV(
            clf, 
            param_grids[name], 
            cv=5, 
            scoring=scoring,
            refit='accuracy',  # Explicitly set refit parameter for multi-metric scoring
            n_jobs=-1
        )
        grid_search.fit(X_train, y_train)
        
        mean_score = grid_search.cv_results_['mean_test_accuracy'][grid_search.best_index_]
        print(f"{name} cross-validated accuracy: {mean_score:.4f}")
        print(f"  Best parameters: {grid_search.best_params_}")
        
        if mean_score > best_score:
            best_score = mean_score
            best_name = name
            best_model = grid_search.best_estimator_

    print(f"Best classifier: {best_name} with accuracy {best_score:.4f}")
    rf_model = best_model

    # Predict for missing Seizures
    if not unknown.empty:
        X_pred = expression_df.loc[unknown.index]
        y_pred = rf_model.predict(X_pred)
        # Save predictions
        imputed = unknown.copy()
        imputed['Seizures'] = y_pred

    # Evaluate on known data using cross-validation
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


