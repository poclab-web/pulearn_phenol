import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap

def pca(df_data, n, whiten=True, fn1= 'Dimension_reduction/PCA_results.csv.', fn2 = 'Dimension_reduction/PCA_correlation.csv'):
    """
    The function to perform Principal Component Analysis (PCA) on given data.

    Args:
        df_data (pandas.DataFrame): Data to be used.
        n (int): Number of new axes.
        whiten (bool): Parameter of PCA. Whether to whiten or not. Defaults to True.
        fn1 (str): Path of the CSV file summarizing the principal components obtained. Defaults to 'Dimension_reduction/PCA_results.csv'.
        fn2 (str): Path of the CSV file summarizing "the correlation coefficients between the principal components obtained and each feature. Defaults to 'Dimension_reduction/PCA_correlation.csv'.

    Returns:
        tuple: A tuple containing the following
            df (pandas.DataFrame): DataFrame summarizing the principal components obtained.
            pc_corr (pandas.DataFrame): DataFrame summarizing "the correlation coefficients between the principal components obtained and each feature.
            ratio (float): Cumulative contribution ratio
    """
    
    model_pca = PCA(n_components=n, whiten=whiten, random_state=0)
    model_pca.fit(df_data)
    new = model_pca.transform(df_data)
    new_df = pd.DataFrame(new)
    new_df.columns = [f'PC{i+1}' for i in new_df.columns]
    new_df.index = df_data.index

    df = pd.concat([df_data, new_df], axis=1)
    pc_corr = df.corr().iloc[:df_data.shape[1], df_data.shape[1]:]

    df.to_csv(f'{fn1}.csv')
    pc_corr.to_csv(f'{fn2}.csv')

    ratio_ar = model_pca.explained_variance_ratio_
    ratio = sum(ratio_ar[0:n])*100

    return df, pc_corr, ratio

def UMAP(df_data, n, n_neighbors=15, min_dist=0.1, metric='euclidean', fn = 'Dimension_reduction/UMAP_results.csv'):
    """
    The function to perform Uniform Manifold Approximation and Projection (UMAP) on given data.

    Args:
        df_data (pandas.DataFrame): Data to be used.
        n (int): Number of new axes.
        n_neighbors (int): Parameter of UMAP. Defaults to 15.
        min_dist (float): Parameter of UMAP. Defaults to 0.1.
        metric (str): Parameter of UMAP. Defaults to 'euclidean'.
        fn (str): Path of the CSV file summarizing the results. Defaults to 'Dimension_reduction/UMAP_results.csv'.

    Returns:
        pandas.DataFrame: DataFrame summarizing the results.
    """
    
    reducer = umap.UMAP(n_components=n, n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=0)
    embedding = reducer.fit_transform(df_data)
    new_df = pd.DataFrame(embedding)
    new_df.columns = [f'UMAP{i+1}' for i in new_df.columns]
    new_df.index = df_data.index

    df = pd.concat([df_data, new_df], axis=1)
    df.to_csv(f'{fn}.csv')

    return df