"""
Package Name: CAT Bridge (Compounds And Trancrips Bridge)
Author: Bowen Yang
email: by8@ualberta
Version: 0.1.0
Description: A detailed description of your package, what it does, and how to use it. 

For more detailed information on specific functions or classes, use the help() function on them. For example:
help(your_package_name.your_function_name)


#  ██████╗ █████╗ ████████╗    ██████╗ ██████╗ ██╗██████╗  ██████╗ ███████╗
# ██╔════╝██╔══██╗╚══██╔══╝    ██╔══██╗██╔══██╗██║██╔══██╗██╔════╝ ██╔════╝
# ██║     ███████║   ██║       ██████╔╝██████╔╝██║██║  ██║██║  ███╗█████╗  
# ██║     ██╔══██║   ██║       ██╔══██╗██╔══██╗██║██║  ██║██║   ██║██╔══╝  
# ╚██████╗██║  ██║   ██║       ██████╔╝██║  ██║██║██████╔╝╚██████╔╝███████╗
#  ╚═════╝╚═╝  ╚═╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝╚═════╝  ╚═════╝ ╚══════╝
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from bioinfokit import analys, visuz
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from PIL import Image
import networkx as nx
from sklearn.metrics.pairwise import cosine_similarity


from scipy.stats import spearmanr
from scipy.stats import pearsonr
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
from sklearn.preprocessing import MinMaxScaler
import subprocess
import math
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import pandas as pd
from tslearn.clustering import TimeSeriesKMeans
from tslearn.preprocessing import TimeSeriesScalerMeanVariance

import openai
import getpass



"""
***********  1. Pre-Processing ***********
#  ██████╗ █████╗ ████████╗    ██████╗ ██████╗ ██╗██████╗  ██████╗ ███████╗
# ██╔════╝██╔══██╗╚══██╔══╝    ██╔══██╗██╔══██╗██║██╔══██╗██╔════╝ ██╔════╝
# ██║     ███████║   ██║       ██████╔╝██████╔╝██║██║  ██║██║  ███╗█████╗  
# ██║     ██╔══██║   ██║       ██╔══██╗██╔══██╗██║██║  ██║██║   ██║██╔══╝  
# ╚██████╗██║  ██║   ██║       ██████╔╝██║  ██║██║██████╔╝╚██████╔╝███████╗
#  ╚═════╝╚═╝  ╚═╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝╚═════╝  ╚═════╝ ╚══════╝
""" 

# *********** Read File ***********
def read_upload(file_name):
    """
    Read the uploaded file and return a dataframe
    supported file format: csv, txt, tsv, xls, xlsx
    """
    #if the last 3 letters of the file name is csv, we use csv module to read it
    if file_name[-3:] == 'csv':
        df = pd.read_csv(file_name, index_col=0)
    #if the last 3 letters of the file name is txt, we use csv module to read it
    elif file_name[-3:] == 'txt':
        df = pd.read_csv(file_name, sep='\t', index_col=0)
    #if the last 3 letters of the file name is tsv, we use csv module to read it
    elif file_name[-3:] == 'tsv':
        df = pd.read_csv(file_name, sep='\t', index_col=0)
    elif file_name[-3:] == 'xls':
        df = pd.read_excel(file_name, index_col=0)
    elif file_name[-4:] == 'xlsx':
        df = pd.read_excel(file_name, index_col=0)
    else:
        print('File format is not supported')
    return df


#get target compounds
def get_target(name, df):
    """
    find the target compound in the dataframe and return the row of the target
    """
    if name in df.index:
        name = df.loc[name].tolist()
    else:
        print('Target is not in the index of the dataframe')
    return name



#def normalize function, that allow log2, log10, and z-score
def normalize(df, method):
    """
    Normalize the data using the specified method.
    """
    if method == 'log2':
        df = np.log2(df + 1)
    elif method == 'log10':
        df = np.log10(df + 1)
    elif method == 'z-score':
        df = (df - df.mean()) / df.std()
    return df

  
    
#  ************* Scaling ***************
def scaling(df, method):
    """
    Scale the data using the specified method.
    """
    if method == 'min-max':
        scaler = MinMaxScaler()
        df = scaler.fit_transform(df)
    elif method == 'pareto':
        df = (df - df.mean()) / np.sqrt(df.std())


# ************* One Column Scaling ***************
def scale_column(df, column_name):
    """
    Scale the values of a column using the MinMaxScaler.
    """
    # Create a scaler object
    scaler = MinMaxScaler()

    # Create a copy of the DataFrame to avoid modifying the original one
    df_scaled = df.copy()

    # Reshape the data to fit the scaler
    data = df[column_name].values.reshape(-1, 1)

    # Fit and transform the data
    df_scaled[column_name] = scaler.fit_transform(data)

    return df_scaled



# ************* Biological Replicates ***************
def repeat_aggregation_max(df, design):
    """
    Aggregate biological replicates by taking the maximum value run for each row (gene/compound).
    
    Parameters:
        df (pandas DataFrame): The DataFrame to be aggregated.
        design (pandas DataFrame): The experimental design DataFrame.
    
    Returns:
        new_df (pandas DataFrame): The aggregated DataFrame.
    """
    # Calculate the number of suffixes
    num_suffixes = int(len(design.index) / len(design['group'].unique()))

    # Extract base column names and maintain their order
    base_cols = [col.rsplit('_', 1)[0] for col in df.columns if '_' in col]
    base_cols = sorted(set(base_cols), key=base_cols.index)

    # Create sub DataFrames and store their means
    mean_dfs = {}
    sub_dfs = {}
    for suffix in np.arange(1, num_suffixes + 1).astype(str):
        sub_df = df.filter(regex=f'_{suffix}$')
        sub_dfs[suffix] = sub_df
        mean_dfs[suffix] = sub_df.mean(axis=1)

    # Create a DataFrame to store which sub DataFrame has the highest mean for each row
    max_mean_df = pd.DataFrame(mean_dfs).idxmax(axis=1)

    # Create a new DataFrame, preserving the original index
    new_df = pd.DataFrame(index=df.index)

    # For each base column name
    for base_col in base_cols:
        for idx in new_df.index:
            # Determine which sub DataFrame to pull from for this row
            sub_df_idx = max_mean_df.loc[idx]

            # Get value from the appropriate sub DataFrame
            new_df.loc[idx, base_col] = sub_dfs[sub_df_idx].loc[idx, f"{base_col}_{sub_df_idx}"]


    new_df.columns = [str(col) + '_1' for col in new_df]
    mapping_dict = design['group'].to_dict()
    new_df.columns = new_df.columns.map(lambda x: mapping_dict[x] if x in mapping_dict else x)
    
    return new_df



def repeat_aggregation_mean(df: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate biological replicates by taking the mean value run for each row (gene/compound).
    
    Parameters:
        df (pandas DataFrame): The DataFrame to be aggregated.
        design (pandas DataFrame): The experimental design DataFrame.
        
    Returns:
        new_df (pandas DataFrame): The aggregated DataFrame.
    """
    # Generate new column names from design DataFrame
    new_column_names = {i: design.loc[i]['group'] for i in df.columns}

    # Rename columns and average by new column names
    df.rename(columns=new_column_names, inplace=True)
    df = df.T.groupby(level=0).mean().T
    
    return df



# ************* Merge ***************
def merge_dataframes(dataframes):
    """
    Merge multiple dataframes based on the 'Name' column.

    Parameters:
        dataframes (list): A list of pandas DataFrames to be merged.

    Returns:
        merged_dataframe (pandas DataFrame): The merged DataFrame.
    """
    # Check if the input list is not empty
    if not dataframes:
        raise ValueError("The input list of dataframes is empty.")

    # Merge the dataframes one by one
    merged_dataframe = dataframes[0]
    for df in dataframes[1:]:
        merged_dataframe = merged_dataframe.merge(df, on='Name', how='outer')

    return merged_dataframe 







   
   
"""
************** 2 COMPUTE ******************
#  ██████╗ █████╗ ████████╗    ██████╗ ██████╗ ██╗██████╗  ██████╗ ███████╗
# ██╔════╝██╔══██╗╚══██╔══╝    ██╔══██╗██╔══██╗██║██╔══██╗██╔════╝ ██╔════╝
# ██║     ███████║   ██║       ██████╔╝██████╔╝██║██║  ██║██║  ███╗█████╗  
# ██║     ██╔══██║   ██║       ██╔══██╗██╔══██╗██║██║  ██║██║   ██║██╔══╝  
# ╚██████╗██║  ██║   ██║       ██████╔╝██║  ██║██║██████╔╝╚██████╔╝███████╗
#  ╚═════╝╚═╝  ╚═╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝╚═════╝  ╚═════╝ ╚══════╝
""" 

# ************* 2.1 Correlation Score ***************

def compute_spearman(df, target):
    # Initialize lists to store results
    index_list = []
    corr_list = []
    p_value_list = []

    # Loop through each row in the dataframe
    for idx, row in df.iterrows():
        # Compute Spearman correlation and p-value
        corr, p_value = spearmanr(row, target)

        # Append results to the lists
        index_list.append(idx)
        corr_list.append(corr)
        p_value_list.append(p_value)

    # Create a new dataframe to store results
    results_df = pd.DataFrame({
        'Name': index_list,
        'Spearman': corr_list
        #'P-value': p_value_list
    })

    return results_df



"""
Pearson 
df is a dataframe, target is a list
Pearson's correlation coefficient is the covariance of the two variables divided by the product of their standard deviations. The form of the definition involves a "product moment", that is, the mean (the first moment about the origin) of the product of the mean-adjusted random variables; hence the modifier product-moment in the name.
"""
def compute_pearson(df, target):
    # Initialize lists to store results
    index_list = []
    corr_list = []
    p_value_list = []

    # Loop through each row in the dataframe
    for idx, row in df.iterrows():
        # Compute Pearson correlation and p-value
        corr, p_value = pearsonr(row, target)

        # Append results to the lists
        index_list.append(idx)
        corr_list.append(corr)
        p_value_list.append(p_value)

    # Create a new dataframe to store results
    results_df = pd.DataFrame({
        'Name': index_list,
        'Pearson': corr_list
        #'P-value': p_value_list
    })

    return results_df





# ***************** Granger *********************
def compute_granger(df, target, maxlag=1):
    # Initialize lists to store results
    index_list = []
    p_value_list = []

    # Loop through each row in the dataframe
    for idx, row in df.iterrows():
        # Skip if the row or the target list has constant values
        if np.std(row.values) == 0 or np.std(target) == 0:
            continue

        # Combine the row and target into a DataFrame
        data = pd.concat([pd.Series(target), pd.Series(row.values)], axis=1)

        # Perform the Granger causality test
        try:
            result = grangercausalitytests(data, maxlag=maxlag, verbose=False)

            # Extract the p-value of the F test for the maximum lag
            p_value = result[maxlag][0]['ssr_ftest'][1]
            p_value = 1-p_value

            # Append results to the lists
            index_list.append(idx)
            p_value_list.append(p_value)
        except Exception as e:
            #add NA if there is an error
            index_list.append(idx)
            p_value_list.append(np.nan)
            
            # print(f"Error at index {idx}: {str(e)}")

    # Create a new dataframe to store results
    results_df = pd.DataFrame({
        'Name': index_list,
        'Granger': p_value_list
    })

    return results_df


# *************** 2.2 FC Score ***************
# ********* Noontide ************
def find_noontide(df, row_name):
    """
    Find the column with the highest value in the row, and return the column, and the column after it (name).
    If the column with the highest value is the last column, then return the column with the second highest value, and the column after it (name).
    
    Parameters:
        df (pandas DataFrame): The DataFrame has been aggregated.
        row_name (str): The name of the row to be detected.
        
    Returns:
        column_n (str): The name of the column with the highest value.
    """
    # Get the row with the specified name
    row = df.loc[row_name]

    # Identify the column with the highest value
    column_n = row.idxmax()
    col_idx = list(df.columns).index(column_n)

    # If column with max value is last column, then find the second highest column
    if col_idx == len(df.columns) - 1:
        row[column_n] = row.min() # set value in column_n to minimum value of row
        column_n = row.idxmax()   # get column with max value now

    # Find column after column_n
    col_idx = list(df.columns).index(column_n)
    if col_idx < len(df.columns) - 1: # Make sure it's not the last column
        column_n_plus_1 = df.columns[col_idx + 1]
    else:
        raise Exception("There is no column after the column with the maximum value.")

    # Keep only column n and column n+1
    df_filtered = df[[column_n, column_n_plus_1]]

    return df_filtered.columns


# ********* df for fc ************
def df_for_fc(df1, target, df2, design):
    """
    Gnerate the design matrix and matrix for computing the FC score.
    
    Parameters:
        df1 (pandas DataFrame): The DataFrame has been aggregated (processed_metabo).
        target (str): The name of the row to be detected (Capsaicin).
        df2 (pandas DataFrame): The DataFrame for fc computing (gene).
        design (pandas DataFrame): study design, the samle and group information.
    """
    noontide = find_noontide(df1, target)
    design_fc = design[design['group'].isin(noontide)]
    matrix_fc = df2[design_fc.index]

    # Saving to CSV files instead of returning
    design_fc.to_csv('result/design_fc.csv')
    matrix_fc.to_csv('result/matrix_fc.csv')
    
    

def no_repeat_fc(df, noontide):
    df['log2FoldChange'] = df[noontide[0]]/df[noontide[1]]
    df['log2FoldChange'] = df['log2FoldChange'].apply(lambda x: math.log2(x))
    #make the range of log2FoldChange from 0-1
    scaler = MinMaxScaler()
    # Apply the scaler to the 'log2FoldChange' column
    df['log2FoldChange'] = scaler.fit_transform(df[['log2FoldChange']])
    df = df[['log2FoldChange']]
    return df



# ********* FC Compute ********
def fc_comp():
    result = subprocess.run(['Rscript', 'FC.R'], stdout=subprocess.PIPE)
    fc = read_upload('result/fc.csv')
    fc.index.name = 'Name'
    #for value in fc['log2FoldChange'], do scaling to makeit range from 0-1
    scaler = MinMaxScaler()
    # Apply the scaler to the 'log2FoldChange' column 
    fc['log2FoldChange'] = -1 * fc['log2FoldChange']
    fc['log2FoldChange'] = scaler.fit_transform(fc[['log2FoldChange']])
    return fc



# ***************  compute score ******************
def score(df):
    """
    Compute a score based on specified columns of a dataframe, and rank the rows based on the score.
    """
    # Exclude the 'Name' column and compute row sums
    df['Score'] = df.drop(columns=['Name']).sum(axis=1)
    
    # Rank rows from high to low based on the sum, '1' being the highest rank
    df['Rank'] = df['Score'].rank(method='min', ascending=False)
    
    # Sort dataframe by rank
    df_sorted = df.sort_values('Rank')
    # df_sorted.set_index('Rank', inplace=True)
    # df_sorted = df_sorted[['Name', 'Score', 'Rank']]
    
    return df_sorted


keywords_scores = {'ase': 0.2, 'enzyme': 0.2, 'synthase': 0.2}
def annotation_score(df, keywords=keywords_scores):
    """
    This function replaces the value in the 'Description' column of the input dataframe with the highest score
    from the keywords dictionary if a word in the description ends with any of the keywords. If a description is NaN,
    it is replaced with 0.1. If no keyword is found in a description, the description is replaced with 0.

    Parameters:
    df (pandas.DataFrame): The input dataframe.
    keywords (dict): A dictionary where the keys are the keywords and the values are the scores. Defaults to keywords_scores.

    Returns:
    df (pandas.DataFrame): The modified dataframe.
    """
    def replace_in_description(description):
        try:
            if pd.isna(description):  # Check if the value is NaN
                return 0.1
            words = description.split()  # Split the description into words
            max_score = None  # Initialize max score
            for word in words: 
                for keyword, score in keywords.items(): 
                    if word.endswith(keyword):  # If the word ends with a keyword
                        if max_score is None or score > max_score:  # If score is higher than current max score
                            max_score = score  # Update max score
            if max_score is not None:  # If a keyword was found
                return max_score  # Return max score
            return 0  # Return 0 if no keywords are found
        except Exception as e:
            print(f"An error occurred: {e}")
            return description  # Return the original description in case of an error

    df['Description'] = df['Description'].apply(replace_in_description)
    return df


# Clustering
def ts_clustering(df, n_clusters):
    """
    Cluster the time series data in the input dataframe using the TS-KMeans algorithm.
    
    Parameters:
        df (pandas.DataFrame): The input dataframe.
        n_clusters (int): The number of clusters to create.
    """
    # Convert DataFrame to NumPy array for compatibility with tslearn
    data = df.values

    # Rescale the time series data so that their mean is 0 and their standard deviation is 1
    scaler = TimeSeriesScalerMeanVariance(mu=0., std=1.)
    data_scaled = scaler.fit_transform(data)

    # Create the KMeans model
    km = TimeSeriesKMeans(n_clusters=n_clusters, metric="euclidean", max_iter=5, random_state=0)

    # Fit the model to the data
    km.fit(data_scaled)

    # Get the cluster labels for each time series
    labels = km.labels_

    # Add the labels as a new column in the original DataFrame
    df['Cluster'] = labels
    df = df['Cluster']

    return df


















#********* 3 PLOT FUNCTION ********   
#  ██████╗ █████╗ ████████╗    ██████╗ ██████╗ ██╗██████╗  ██████╗ ███████╗
# ██╔════╝██╔══██╗╚══██╔══╝    ██╔══██╗██╔══██╗██║██╔══██╗██╔════╝ ██╔════╝
# ██║     ███████║   ██║       ██████╔╝██████╔╝██║██║  ██║██║  ███╗█████╗  
# ██║     ██╔══██║   ██║       ██╔══██╗██╔══██╗██║██║  ██║██║   ██║██╔══╝  
# ╚██████╗██║  ██║   ██║       ██████╔╝██║  ██║██║██████╔╝╚██████╔╝███████╗
#  ╚═════╝╚═╝  ╚═╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝╚═════╝  ╚═════╝ ╚══════╝
    
    
    
# Line
def plot_line(target, title):
    plt.figure(figsize=(10, 3))
    plt.plot(target)
    plt.title(title)
    plt.show()
    return

# Heatmap
def plot_heatmap(dataframe, palette='viridis', figsize=(10, 8), row_threshold=50, save_path=None):
    #remov the rows that have all 0 values
    dataframe = dataframe.loc[~(dataframe==0).all(axis=1)]
    if dataframe.shape[0] > row_threshold:
        row_labels = False
    else:
        row_labels = True
    
    sns.clustermap(dataframe, cmap=palette, yticklabels=row_labels, xticklabels=True, figsize=figsize)
    #set x axis label
    plt.xlabel(dataframe.columns.name)
    
    if save_path:
        plt.savefig(save_path)
    
    plt.show()
    return 


# Hexbin
def plot_hexbin(data, x_axis, y_axis, gridsize=20):
    sns.jointplot(x=x_axis, y=y_axis, data=data, kind='hex', gridsize=20)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.show()


# Line Heatmap
def plot_line_heatmap(df, name_value, cmap='vlag'):
    """
    Plot a line plot and a heatmap of the specified row in the input dataframe.
    
    Parameters:
        df (pandas.DataFrame): The input dataframe.
        name_value (str): The value in the 'Name' column of the row to plot.
        cmap (str): The name of the colormap to use for the heatmap. Defaults to 'vlag'.
    """
    # Find the row where 'Name' is equal to name_value
    row = df.loc[df['Name'] == name_value]

    # Extract the desired columns and convert them to a list
    desired_columns = ['Granger', 'log2FoldChange', 'Pearson', 'Spearman']
    list1 = row[desired_columns].values.tolist()[0]

    # Get the other columns by dropping the ones already in list1
    remaining_columns = df.drop(columns=['Name'] + desired_columns)
    list2 = remaining_columns.loc[remaining_columns.index == row.index[0]].values.tolist()[0]

    # Reshape list1 into a 2D array for the heatmap
    list1_2d = np.array(list1).reshape(-1, 1)

    # Create subplots with adjusted sizes
    fig = plt.figure(figsize=(10, 2))
    grid = plt.GridSpec(1, 10, hspace=0.2, wspace=0.2)  # We'll use 8 columns in total

    # Plot the line graph on the left (using 7 out of 8 columns)
    ax1 = plt.subplot(grid[:9])  # equivalent to grid[0, :7]
    ax1.plot(list2, color='navy')
    ax1.set_title(name_value)
    ax1.set_xticks(range(len(remaining_columns.columns)))
    ax1.set_xticklabels(remaining_columns.columns, rotation=90)

    # Plot the heatmap on the right (using 1 out of 8 columns) without color bar
    ax2 = plt.subplot(grid[9:])  # equivalent to grid[0, 7:]
    sns.heatmap(list1_2d, ax=ax2, cmap=cmap, cbar=False, yticklabels=desired_columns)
    ax2.yaxis.tick_right()
    ax2.xaxis.set_visible(False)  # hide the x-axis
    ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0)  # set rotation to 0

    # Show the plotc
    plt.tight_layout()
    plt.show()


# Line Heatmap for all data
def plot_data(df, name_value, ax1, ax2, cmap='vlag'):
    """
    plot the line graph and heatmap for the given gene name
    
    df: dataframe
    name_value: the name of the gene
    ax1: the axis for the line graph, ax1 = plt.subplot(grid[:9])
    ax2: the axis for the heatmap, ax2 = plt.subplot(grid[9:])
    cmap: the color map for the heatmap
    """
    # Find the row where 'Name' is equal to name_value
    row = df.loc[df['Name'] == name_value]

    # Extract the desired columns and convert them to a list
    desired_columns = ['Granger', 'log2FoldChange', 'Pearson', 'Spearman']
    list1 = row[desired_columns].values.tolist()[0]

    # Get the other columns by dropping the ones already in list1
    remaining_columns = df.drop(columns=['Name'] + desired_columns)
    list2 = remaining_columns.loc[remaining_columns.index == row.index[0]].values.tolist()[0]

    # Reshape list1 into a 2D array for the heatmap
    list1_2d = np.array(list1).reshape(-1, 1)

    # Plot the line graph on the left 
    ax1.plot(list2, color='navy')
    ax1.set_title(name_value)
    ax1.set_xticks(range(len(remaining_columns.columns)))
    ax1.set_xticklabels(remaining_columns.columns, rotation=90)

    # Plot the heatmap on the right without color bar
    sns.heatmap(list1_2d, ax=ax2, cmap=cmap, cbar=False, yticklabels=desired_columns)
    ax2.yaxis.tick_right()
    ax2.xaxis.set_visible(False)  # hide the x-axis
    ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0)  # set rotation to 0


def plot_all_data(df):
    """
    plot the line graph and heatmap for the top 10 genes in the dataframe
    """
    df = df.head(10)
    num_rows = df.shape[0]
    fig = plt.figure(figsize=(10, num_rows * 2))  # adjust the figure height based on the number of rows
    grid = plt.GridSpec(num_rows, 10, hspace=0.5, wspace=0.2)

    for i, row in df.iterrows():
        ax1 = plt.subplot(grid[i, :9])  # line graph
        ax2 = plt.subplot(grid[i, 9:])  # heatmap
        plot_data(df, row['Name'], ax1, ax2)
        
        # Hide the x labels except for the last line plot
        if i < num_rows - 1:
            ax1.set_xticklabels([]) 
        else:
            ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

    plt.tight_layout()
    plt.show()




# Volcano Plot
def plot_volcano(path, lfc_threshold, padj_threshold):
    # Import and preprocess data
    fc_for_volcano = pd.read_csv(path)
    fc_for_volcano.reset_index(inplace=True)
    fc_for_volcano.rename(columns={'index':'Name'}, inplace=True)
    gene_exp = fc_for_volcano
    gene_exp = gene_exp.dropna(subset=['padj'])

    # Check if at least 10 gene names exist
    genenames = gene_exp['Name'].head(10) if len(gene_exp['Name']) >= 10 else None

    # Create plot
    visuz.GeneExpression.volcano(df=gene_exp, 
                                lfc='log2FoldChange', pv='padj', sign_line=True,
                                lfc_thr=(lfc_threshold, lfc_threshold), pv_thr=(padj_threshold, padj_threshold),
                                plotlegend=True, legendpos='upper right', legendanchor=(1.46,1),
                                color=('maroon','gainsboro','steelblue'), theme='whitesmoke',
                                valpha=1, dotsize=5,
                                geneid = 'Name',
                                genenames = tuple(genenames) if genenames is not None else None,
                                )

    # Save and display image
    # plt.savefig('result/volcano_plot.png')
    img = Image.open('volcano.png')
    
    img = Image.open('volcano.png')  # replace with your image file path if not in the same directory

    # Create a figure and a set of subplots
    fig, ax = plt.subplots()

    # Display the image
    ax.imshow(img)

    # Remove the axis
    ax.axis('off')

    # Show the figure
    plt.show()

    # return gene_exp



def plot_ts_clusters(result, processed_gene, palette_name='Paired'):
    result = result[['Name', 'Cluster']]
    data = merge_dataframes([result, processed_gene])

    # Drop the 'Name' column
    data = data.drop(columns=['Name'])

    # Reset index, group by 'Cluster', and set the index back
    data = data.reset_index()
    grouped = data.groupby('Cluster')

    # Set color palette
    palette = sns.color_palette(palette_name, 12)  # The palette_name palette has maximum 12 distinct colors

    # Iterate over groups (clusters) and plot each one
    for name, group in grouped:
        group = group.set_index('index')  # set the index back to 'index'
        group = group.drop(columns='Cluster')  # drop the 'Cluster' column

        plt.figure(figsize=(10, 6))
        for i, feature in enumerate(group.index):
            plt.plot(group.columns, group.loc[feature], color=palette[i % len(palette)])

        plt.xlabel('Time')
        plt.ylabel('Value')
        plt.title(f'Cluster: {name}', weight='bold')
        plt.show()



# PCA
def plot_pca(gene, design, n_clusters):
    # transpose your DataFrame, as PCA works on the features (columns), not on the samples (rows)
    gene_transposed = gene.T

    # Perform PCA on your data
    pca = PCA(n_components=2)  # here we ask for the first two principal components
    pca_result = pca.fit_transform(gene_transposed)

    # convert the PCA result to a DataFrame
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])

    # add a 'sample' column
    pca_df['sample'] = gene_transposed.index

    # get group information from the design dataframe
    pca_df = pca_df.merge(design[['group']], left_on='sample', right_index=True)

    if n_clusters:
        # Perform K-means clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(pca_df[['PC1', 'PC2']])
        pca_df['Cluster'] = kmeans.labels_

        # create the figure and axis objects
        fig, ax = plt.subplots(figsize=(8, 6))

        # plot the points on the scatterplot
        scatter = sns.scatterplot(x="PC1", y="PC2", hue="group", data=pca_df, palette="Paired", s=100, alpha=0.7, ax=ax)

        # for each cluster, add a circle at the mean coordinates with radius proportional to the standard deviation
        for cluster in set(kmeans.labels_):
            cluster_points = pca_df[pca_df['Cluster'] == cluster][['PC1', 'PC2']]
            # calculate mean and standard deviation for the cluster
            cluster_mean = cluster_points.mean().values
            cluster_std = cluster_points.std().values
            # add a circle at the mean coordinates with radius=stddev
            circle = Circle(cluster_mean, np.linalg.norm(cluster_std), alpha=0.1)
            ax.add_artist(circle)

        # hide the legend
        ax.get_legend().remove()
    else:
        # create the figure and axis objects
        fig, ax = plt.subplots(figsize=(8, 6))

        # plot the points on the scatterplot without clustering
        scatter = sns.scatterplot(x="PC1", y="PC2", hue="group", data=pca_df, palette="Paired", s=100, alpha=0.7, ax=ax)

        # hide the legend
        ax.get_legend().remove()

    # annotate points on the graph with the sample names
    for i, sample in enumerate(pca_df['sample']):
        plt.annotate(sample, (pca_df.iloc[i].PC1, pca_df.iloc[i].PC2), color='gray')

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA with K-means clustering', fontweight='bold')

    plt.show()




def plot_network(data, target_index, num_nodes):
    # Compute the similarity
    similarity = cosine_similarity(data)
    similarity_df = pd.DataFrame(similarity, index=data.index, columns=data.index)
    
    # Get the top num_nodes similar nodes for the target
    target_similarities = similarity_df.loc[target_index].sort_values(ascending=False)[1:num_nodes+1].to_dict()

    # Create a network graph
    G = nx.Graph()

    # Add nodes and edges
    for node, similarity in target_similarities.items():
        G.add_edge(target_index, node, weight=similarity)
    
    # Draw the network
    plt.figure(figsize=(8,8))
    pos = nx.spring_layout(G)
    colors = ['red' if node == target_index else 'moccasin' for node in G.nodes()]
    
    nx.draw_networkx_nodes(G, pos, node_color=colors)
    nx.draw_networkx_labels(G, pos, font_size=8)
    nx.draw_networkx_edges(G, pos, edge_color='gray', width=[G[u][v]['weight'] for u,v in G.edges()], alpha=0.7)

    plt.title(f'Top {num_nodes} similar nodes to the target', fontsize=10)
    plt.axis('off')  # to turn off the frame
    plt.show()








# AI
def Yuanfang(df, target):
    # annotation = read_upload(annotation_file)
    # df = merge_dataframes([df, annotation])
    df = df.head(20)
    
    hits = df['Description'].to_list()
    hits = [str(item) for item in hits]
    hits = ', '.join(hits)
    
    q = hits + '\n\n\nWhich one may be involved in the synthesis of ' + target + '?'
    
    openai_api_key = getpass.getpass("Please enter your OpenAI API Key: ")
    openai.api_key = openai_api_key

    messages = [
        {"role": "system", "content": "You are a biological chemist and can explain biological mechanisms"},
        {"role": "user", "content": q}
    ]

    completion = openai.ChatCompletion.create(
        model = "gpt-3.5-turbo",
        temperature = 0.8,
        max_tokens = 2000,
        messages = messages
    )
    
    print(' ')
    print(completion.choices[0].message.content)
    print(' ')
    print(' ')
    print('NOTICE: The output was produced by the large language model GPT 3.5 turbo, so it should only be regarded as a source of inspiration.')









# ************* Pipeline *************
def pipeline(gene_file, metabo_file, design_file, annotation_file, target, cluster_count, aggregation_func=None):
    """
    This function processes gene expression data, performs computations, and returns results.

    Parameters:
    - gene_file (str): The filename of the gene count data.
    - metabo_file (str): The filename of the metabolome data.
    - design_file (str): The filename of the experimental design data (can be None).
    - annotation_file (str): The filename of the gene annotation data (can be None).
    - target (str): The target metabolite for the analysis.
    - cluster_count (int): The number of clusters for time series clustering.
    - aggregation_func (function): The function to be used for data aggregation.

    Returns:
    - result (pd.DataFrame): A pandas DataFrame containing the processed results, with annotations and clustering information if provided.
    """
    # Read data
    gene = read_upload(gene_file)
    metabo = read_upload(metabo_file)

    if design_file is not None:
        # If there is a design file
        design = read_upload(design_file)

        # Process data
        processed_gene = aggregation_func(gene, design)
        processed_metabo = aggregation_func(metabo, design)

        # Get target data
        t = get_target(target, processed_metabo)

        # Compute Granger causality
        granger = compute_granger(processed_gene, t, 1)

        # Prepare dataframe for fold change calculation
        df_for_fc(processed_metabo, target, gene, design)

        # Compute fold change
        fc = fc_comp()

        # Merge data
        data = merge_dataframes([granger, fc])

    else:
        # If there is no design file
        t = get_target(target, metabo)
        noontide = find_noontide(metabo, target)
        granger = compute_granger(gene, t)

        # Compute fold change
        fc = no_repeat_fc(gene, noontide)
        data = merge_dataframes([granger, fc])
        
    # Compute score
    result = score(data)

    # Perform clustering
    cluster = ts_clustering(gene, cluster_count)

    if annotation_file is not None:
        # If there is an annotation file
        annotation = read_upload(annotation_file)
        data = merge_dataframes([data, annotation])
        data = annotation_score(data)
        result = merge_dataframes([result, annotation, cluster])
    else:
        result = merge_dataframes([result, cluster])

    result.set_index('Rank', inplace=True)

    return result





# print('''
      
      
#  ██████╗ █████╗ ████████╗    ██████╗ ██████╗ ██╗██████╗  ██████╗ ███████╗
# ██╔════╝██╔══██╗╚══██╔══╝    ██╔══██╗██╔══██╗██║██╔══██╗██╔════╝ ██╔════╝
# ██║     ███████║   ██║       ██████╔╝██████╔╝██║██║  ██║██║  ███╗█████╗  
# ██║     ██╔══██║   ██║       ██╔══██╗██╔══██╗██║██║  ██║██║   ██║██╔══╝  
# ╚██████╗██║  ██║   ██║       ██████╔╝██║  ██║██║██████╔╝╚██████╔╝███████╗
#  ╚═════╝╚═╝  ╚═╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝╚═════╝  ╚═════╝ ╚══════╝
# ''')
# print(r'''
#                                    ___
#                       |\---/|  / )|Guo|
#           ------------;     |-/ / |Lab|
#                       )     (' /  `---'
#           ===========(       ,'==========
#           ||  _      |      |      
#           || ( (    /       ;
#           ||  \ `._/       /
#           ||   `._        /|
#           ||      |\    _/||
#         __||_____.' )  |__||____________
#          ________\  |  |_________________
#                   \ \  `-.
#                    `-`---'  
                   



# CAT Bridge is a Python-based, cross-platform tool designed for integrated analysis of transcriptomes and metabolites. It is adept at processing time series data, aiding in the discovery of genes contributing to specific metabolite synthesis. 

# For more information on our other research outputs, welcome to visit the Guo Lab's website at: http://www.guo-lab.site.

# For any queries about this software, feel free to contact us at by8@ualberta.ca.

# ''')