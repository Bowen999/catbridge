"""
Package Name: your_package_name
Version: your_package_version
Description: A detailed description of your package, what it does, and how to use it. 

For more detailed information on specific functions or classes, use the help() function on them. For example:
help(your_package_name.your_function_name)
"""


import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec


from scipy.stats import spearmanr
from scipy.stats import pearsonr
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
from sklearn.preprocessing import MinMaxScaler
import subprocess

import openai



"""
***********  1. Pre-Processing ***********
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



def compute_score(data, col1, desc_col=None, col2=None, keywords=None):
    """
    Compute a score based on specified columns of a dataframe. Can also convert a description
    column to a score based on specified keywords. 

    Parameters:
    data (pandas.DataFrame): The dataframe to process.
    col1 (str): The name of the first column to include in the score.
    desc_col (str, optional): The name of the description column to convert to a score. Defaults to None.
    col2 (str, optional): The name of a second column to include in the score. Defaults to None.
    keywords (list, optional): List of keywords to use in converting the description to a score. Defaults to None.

    Returns:
    pandas.DataFrame: The processed dataframe, sorted by score in descending order.
    """
    data[col2] = data[col2].astype(float)
    data[col1] = data[col1].astype(str)
    try:
        if col2 is not None:
            data['score'] = data[col1] + data[col2]
        else:
            data['score'] = data[col1]
    except KeyError as e:
        print(f"Error: {e} not found in dataframe.")
        return None
    
    if desc_col is not None:
        try:
            data[desc_col] = data[desc_col].astype(str)
            if keywords is not None:
                data = data.apply(convert_description_to_score, args=(desc_col, keywords), axis=1)
        except KeyError as e:
            print(f"Error: {e} not found in dataframe.")
            return None
        
    data = data.sort_values(by='score', ascending=False)
    try:
        data = data.drop(columns=['index', 'level_0'])
    except:
        pass
    return data


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








#********* 3 PLOT FUNCTION ********
#plot a line plot for the target
def plot_line(target, title):
    plt.figure(figsize=(10, 3))
    plt.plot(target)
    plt.title(title)
    plt.show()
    return


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



def plot_hexbin(data, x_axis, y_axis, gridsize=20):
    sns.jointplot(x=x_axis, y=y_axis, data=data, kind='hex', gridsize=20)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.show()



def plot_line_heatmap(df, name_value, cmap='vlag'):
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




# AI
def smart_aleck(question):
    openai.api_key = "my_api_key"

    messages = [
        {"role": "system", "content": "You are a biological chemist and can explain biological mechanisms, and also translate your answer to Chinese"},
        {"role": "user", "content": question}
    ]

    completion = openai.ChatCompletion.create(
        model = "gpt-3.5-turbo",
        temperature = 0.8,
        max_tokens = 2000,
        messages = messages
    )

    return completion.choices[0].message.content




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