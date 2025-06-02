import os
from tqdm import tqdm
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from sklearn.ensemble import RandomForestClassifier
from pulearn import BaggingPuClassifier

import warnings
warnings.simplefilter('ignore')

def Parameter_tuning(X, y):
    """
    Parameter tuning for PU learning using BaggingPuClassifier with RandomForestClassifier.
    This function evaluates the performance of the model with different numbers of estimators and max sample times.
    Args:
        X (DataFrame-like): Input features for the model.
        y (array-like): Labels for the input data, where 1 indicates positive samples, 0 indicates unlabeled samples, and -1 indicates negative samples.
    Returns:
        tuple: A tuple containing the following
            df_TPR (pandas.DataFrame): DataFrame summarizing True Positive Rate (TPR) for each validation.
            df_TNR (pandas.DataFrame): DataFrame summarizing True Negative Rate (TNR) for each validation.
            df_MCC (pandas.DataFrame): DataFrame summarizing Matthews correlation coefficient (MCC) for each validation.
    """

    if type(X) != np.ndarray:
        X = X.to_numpy()
    
    if type(y) != np.ndarray:
        y = y.to_numpy()

    p_index = np.where(y == 1)[0]
    n_index = np.where(y == -1)[0]

    y[y == -1] = 0
    
    clf = RandomForestClassifier(random_state=0, class_weight='balanced')
    n_estimators_options = [100, 200, 300, 400, 500]
    max_samples_times = [t for t in range(1, 21)]

    TPR_li1 = []
    TNR_li1 = []
    MCC_li1 = []

    for n_estimators_ in tqdm(n_estimators_options):
        
        TPR_li2 = []
        TNR_li2 = []
        MCC_li2 = []
        
        for time in max_samples_times:
            # Validation_P
            TP = 0
            total_P = 0

            for i in p_index:
                y_train = y.copy()
                y_train[i] = 0
                X_test = X[i].reshape(1, -1)
                y_test = y[i]
                
                pu_estimator = BaggingPuClassifier(estimator=clf, n_estimators=n_estimators_,
                                                   max_samples=sum(y_train)*time, random_state=0, n_jobs=-1)
                
                pu_estimator.fit(X, y_train)
                y_pred = pu_estimator.predict(X_test)
                del pu_estimator
                
                if y_pred[0] == y_test:
                    TP += 1
                
                total_P += 1
            
            FN = total_P - TP
            TPR = TP / total_P
            TPR_li2.append(TPR)
            
            # Validation_N
            X_test = X[n_index]
            y_test = y[n_index]
            pu_estimator = BaggingPuClassifier(estimator=clf, n_estimators=n_estimators_,
                                               max_samples=sum(y)*time, random_state=0, n_jobs=-1)
            
            pu_estimator.fit(X, y)
            y_pred = pu_estimator.predict(X_test)
            del pu_estimator
            
            TN = len(y_test) - sum(y_pred)
            FP = sum(y_pred)
            total_N = len(y_test)
            TNR = TN / total_N
            TNR_li2.append(TNR)

            MCC = (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
            MCC_li2.append(MCC)

        TPR_li1.append(TPR_li2)
        TNR_li1.append(TNR_li2)
        MCC_li1.append(MCC_li2)
    
    df_TPR = pd.DataFrame(TPR_li1, columns=max_samples_times, index=n_estimators_options)
    df_TNR = pd.DataFrame(TNR_li1, columns=max_samples_times, index=n_estimators_options)
    df_MCC = pd.DataFrame(MCC_li1, columns=max_samples_times, index=n_estimators_options)

    return df_TPR, df_TNR, df_MCC

def Fixed_parameter_validation(X, y, ylds_li, smiles_li, condition, descriptor,
                               n_estimators, max_sample_times):
    """
    This function evaluates the performance of the model at fixed parameters. 
    In addition, it also outputs images of compounds that were wrongly predicted, etc.
    Args:
        X (DataFrame-like): Input features for the model.
        y (array-like): Labels for the input data, where 1 indicates positive samples, 0 indicates unlabeled samples, and -1 indicates negative samples.
        ylds_li (array-like): Dimer yield for each data (compound) as substrate under the target reaction condition.
        smiles_li (array-like): SMILES for each data (compound).
        condition (str): Name of the target reaction condition.
        descriptor (str): Descriptor Name.
        n_estimators (int): Parameter 'n_estimators' in BaggingPuClassifier
        max_sample_times (int): Set the parameter 'max_samples' of BaggingPuClassifier to the set number x the number of Positive samples.

    Returns:
        pandas.DataFrame: DataFrame summarizing the DPS and other data for each compound output in the validation.
    """
    
    result_path = f'Fixed_parameter_validation/{descriptor}'
    uncorrect_positive_path = f'{result_path}/uncorrect_compounds/positive'
    uncorrect_negative_path = f'{result_path}/uncorrect_compounds/negative'
    os.makedirs(result_path, exist_ok=True)
    os.makedirs(f'{result_path}/uncorrect_compounds', exist_ok=True)
    os.makedirs(uncorrect_positive_path, exist_ok=True)
    os.makedirs(uncorrect_negative_path, exist_ok=True)
    os.makedirs(f'{uncorrect_positive_path}/{condition}', exist_ok=True)
    os.makedirs(f'{uncorrect_negative_path}/{condition}', exist_ok=True)

    if type(X) != np.ndarray:
        X = X.to_numpy()

    if type(y) != np.ndarray:
        y = y.to_numpy()
    
    if type(ylds_li) != np.ndarray:
        ylds_li = ylds_li.to_numpy()
    
    if type(smiles_li) != np.ndarray:
        smiles_li = smiles_li.to_numpy()

    p_index = np.where(y == 1)[0]
    n_index = np.where(y == -1)[0]

    y[y == -1] = 0
    clf = RandomForestClassifier(random_state=0, class_weight='balanced')
    
    # Validation_p
    TP = 0
    total_p = 0
    dps_li_p = []

    for i in p_index:
        y_train = y.copy()
        y_train[i] = 0
        X_test = X[i].reshape(1, -1)
        y_test = y[i]
        
        pu_estimator = BaggingPuClassifier(estimator=clf, n_estimators=n_estimators,
                                           max_samples=sum(y_train)*max_sample_times, random_state=0, n_jobs=-1)
        pu_estimator.fit(X, y_train)
        y_pred = pu_estimator.predict(X_test)
        dps = pu_estimator.oob_decision_function_[i, 1]
        dps_li_p.append(dps)
        del pu_estimator

        if y_pred[0] == y_test:
            TP += 1
        
        else:
            uncorrect_mol = Chem.MolFromSmiles(smiles_li[i])
            Draw.MolToFile(uncorrect_mol, f'{uncorrect_positive_path}/{condition}/{smiles_li[i]}.png')
        
        total_p += 1
    
    TPR = TP / total_p
    print(f'TPR of {condition} is {TPR}')

    df_dy_p = pd.DataFrame({'smiles':smiles_li[p_index], 'ylds':ylds_li[p_index], 'DPS':dps_li_p})

    # Validation_n
    if len(n_index) != 0:
        X_test = X[n_index]
        y_test = y[n_index]
        
        pu_estimator = BaggingPuClassifier(estimator=clf, n_estimators=n_estimators,
                                           max_samples=sum(y)*max_sample_times, random_state=0, n_jobs=-1)
        pu_estimator.fit(X, y)
        y_pred = pu_estimator.predict(X_test)
        dps_li_n = pu_estimator.oob_decision_function_[:, 1][n_index]
        del pu_estimator

        for i, n_i in enumerate(n_index):
            if y_pred[i] != y_test[i]:
                uncorrect_mol = Chem.MolFromSmiles(smiles_li[n_i])
                Draw.MolToFile(uncorrect_mol, f'{uncorrect_negative_path}/{condition}/{smiles_li[n_i]}.png')
                
        TNR = (len(y_test) - sum(y_pred))/len(y_test)
        print(f'TNR of {condition} is {TNR}')
        
        df_dy_n = pd.DataFrame({'smiles':smiles_li[n_index], 'ylds':ylds_li[n_index], 'DPS':dps_li_n})
    
    try:
        df_dy = pd.concat([df_dy_p, df_dy_n])
    
    except:
        df_dy = df_dy_p
    
    return df_dy

def Predict_unlabeled(X, y, n_estimators, max_sample_times):
    """
    Args:
        X (DataFrame-like): Input features for the model.
        y (array-like): Labels for the input data, where 1 indicates positive samples, 0 indicates unlabeled samples, and -1 indicates negative samples.
        n_estimators (int): Parameter 'n_estimators' in BaggingPuClassifier
        max_sample_times (int): Set the parameter 'max_samples' of BaggingPuClassifier to the set number x the number of Positive samples.

    Returns:
        tuple: A tuple containing the following
            y_pred (numpy.ndarray): An array containing the predicted value for each compound (1 for Positive, 0 for Negative)
            dps_li (numpy.ndarray): Array containing the DPS for each compound
    """
    
    if type(X) != np.ndarray:
        X = X.to_numpy()

    if type(y) != np.ndarray:
        y = y.to_numpy()
    
    y[y == -1] = 0

    clf = RandomForestClassifier(random_state=0, class_weight='balanced')
    pu_estimator = BaggingPuClassifier(estimator=clf, n_estimators=n_estimators,
                                       max_samples=sum(y)*max_sample_times, random_state=0, n_jobs=-1)
    pu_estimator.fit(X, y)
    y_pred = pu_estimator.predict(X)
    dps_li = pu_estimator.oob_decision_function_[:, 1]
    del pu_estimator

    return y_pred, dps_li