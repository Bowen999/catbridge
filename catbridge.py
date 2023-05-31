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

import openai


# *********** Read File ***********
def read_upload(file_name):
    #if the last 3 letters of the file name is csv, we use csv module to read it
    if file_name[-3:] == 'csv':
        df = pd.read_csv(file_name, index_col=0)
    #if the last 3 letters of the file name is txt, we use csv module to read it
    elif file_name[-3:] == 'txt':
        df = pd.read_csv(file_name, sep='\t', index_col=0)
    #if the last 3 letters of the file name is tsv, we use csv module to read it
    elif file_name[-3:] == 'tsv':
        df = pd.read_csv(file_name, sep='\t', index_col=0)
    return df


#target compounds
def read_target(name, df):
    if name in df.index:
        name = df.loc[name].tolist()
    else:
        print('Target is not in the index of the dataframe')
    return name








# ************* Normalize ************
#def normalize function, that allow log2, log10, and z-score
def normalize(df, method):
    if method == 'log2':
        df = np.log2(df + 1)
    elif method == 'log10':
        df = np.log10(df + 1)
    elif method == 'z-score':
        df = (df - df.mean()) / df.std()
    return df

  
    
#  ************* Scaling ***************
def scaling(df, method):
    if method == 'min-max':
        scaler = MinMaxScaler()
        df = scaler.fit_transform(df)
    elif method == 'pareto':
        df = (df - df.mean()) / np.sqrt(df.std())
        
        
        

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


def scale_column(df, column_name):
    # Create a scaler object
    scaler = MinMaxScaler()

    # Create a copy of the DataFrame to avoid modifying the original one
    df_scaled = df.copy()

    # Reshape the data to fit the scaler
    data = df[column_name].values.reshape(-1, 1)

    # Fit and transform the data
    df_scaled[column_name] = scaler.fit_transform(data)

    return df_scaled












#*********PLOT FUNCTION********
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

    # Show the plot
    plt.tight_layout()
    plt.show()




def plot_data(df, name_value, ax1, ax2, cmap='vlag'):
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

   
   
   
   
   
   
   
#*********** COMPUTE ******************
"""
Spearman
"""
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




"""
***************** Granger *********************
"""
def compute_granger(df, target, maxlag=2):
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


def compute_score(data, col1, col2, column_name, keywords):
    data['score'] = data[col1] + data[col2]
    data[column_name] = data[column_name].astype(str)
    data = data.apply(convert_description_to_score, args=(column_name, keywords), axis=1)
    data = data.sort_values(by='score', ascending=False)
    try:
        data = data.drop(columns=['index', 'level_0'])
    except:
        pass
    return data


def convert_description_to_score(row, column_name, keywords):
    total_score = 0
    for word in row[column_name].split():
        for keyword in keywords:
            if word.endswith(keyword):
                total_score += keywords[keyword]
                total_score = min(total_score, 0.2)  # Cap the total score at 0.2
    row['score'] += total_score
    return row






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




print('''
      
      
 ██████╗ █████╗ ████████╗    ██████╗ ██████╗ ██╗██████╗  ██████╗ ███████╗
██╔════╝██╔══██╗╚══██╔══╝    ██╔══██╗██╔══██╗██║██╔══██╗██╔════╝ ██╔════╝
██║     ███████║   ██║       ██████╔╝██████╔╝██║██║  ██║██║  ███╗█████╗  
██║     ██╔══██║   ██║       ██╔══██╗██╔══██╗██║██║  ██║██║   ██║██╔══╝  
╚██████╗██║  ██║   ██║       ██████╔╝██║  ██║██║██████╔╝╚██████╔╝███████╗
 ╚═════╝╚═╝  ╚═╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝╚═════╝  ╚═════╝ ╚══════╝
''')
print(r'''
                                   ___
                      |\---/|  / )|Guo|
          ------------;     |-/ / |Lab|
                      )     (' /  `---'
          ===========(       ,'==========
          ||  _      |      |      
          || ( (    /       ;
          ||  \ `._/       /
          ||   `._        /|
          ||      |\    _/||
        __||_____.' )  |__||____________
         ________\  |  |_________________
                  \ \  `-.
                   `-`---'  
                   



CAT Bridge is a Python-based, cross-platform tool designed for integrated analysis of transcriptomes and metabolites. It is adept at processing time series data, aiding in the discovery of genes contributing to specific metabolite synthesis. 

For more information on our other research outputs, welcome to visit the Guo Lab's website at: http://www.guo-lab.site.

For any queries about this software, feel free to contact us at by8@ualberta.ca.

''')