import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import numpy as np
import pandas as pd


def plot_grouped_bar_chart(motifs, group_proportions, group_labels, chart_title, colors=['#4daf4a', '#377eb8', 'darkturquoise']):
    """
    Plot a grouped bar chart with multiple attributes and their proportions for each motif.
    
    Parameters:
        motifs (list): List of motifs to be displayed on the x-axis.
        group_proportions (dict): Dictionary containing attribute names as keys and their corresponding proportions as values.
        group_labels (dict): Dictionary containing attribute names as keys and their corresponding labels for bar annotations.
        chart_title (str): Title for the chart.
        colors (list, optional): List of colors for the bars. Defaults to a colorblind-friendly palette.

    Returns:
        None
    """
    x = np.arange(len(motifs))  # the label locations
    width = 0.32  # Adjust width for better spacing
    multiplier = 1 - (len(group_proportions) - 1) / 2

    fig, ax = plt.subplots(figsize=(10, 6))

    for i, (attribute, measurement) in enumerate(group_proportions.items()):
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute, color=colors[i])
        labels = group_labels[attribute]
        ax.bar_label(rects, labels, padding=3, fontsize=7.5)  # Smaller font for labels
        multiplier += 1

    ax.set_ylabel('Proportion of Proteins (%)', fontsize=12)
    ax.set_title(chart_title, fontsize=16)
    ax.set_xticks(x + width * (len(group_proportions) - 1) / 2)  # Adjust the ticks to center them properly
    ax.set_xticklabels(motifs, fontsize=12)
    ax.legend(loc='upper right', bbox_to_anchor=(0.95, 0.99), ncol=1, fontsize=12, frameon=False)  # Clean legend
    ax.set_ylim(0, 100)

    # Adjust grid, remove excessive lines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, color='gray', linestyle='--', linewidth=0.5)  # Subtle gridlines

    # Adjust x-axis limits to reduce space on sides
    ax.set_xlim(-0.3, len(motifs) -0.1)
    
    plt.tight_layout()
    
    return fig, ax


def plot_histogram_with_hue(df1, df2=None, hue_col=None, title=None, color=None):
    """
    Plot a bar plot with or without hue based on input DataFrames.

    Parameters:
        df1 (pd.DataFrame): First DataFrame for plotting.
        df2 (pd.DataFrame, optional): Second DataFrame for plotting with hue. Defaults to None.
        hue_col (str, optional): Column name for hue distinction. Defaults to None.
        title (str, optional): Title for the plot. Defaults to None.
        color (str, optional): Color for the bars without hue. Defaults to None.

    Returns:
        None
    """
        
    if df2 is not None:
        df1[hue_col] = 'farnesylated'
        df2[hue_col] = 'geranylgeranylated'
        res = pd.concat([df1, df2])
        ax = sns.barplot(x='x', y='y', data=res, hue=hue_col)
    else:
        ax = sns.barplot(x='x', y='y', data=df1, color=color)
    
    plt.xticks(rotation=90)

    ax.set_ylabel('N of proteins')
    ax.set_xlabel('position of most C-terminal peptide')
    
    if title:
        ax.set_title(title)

    plt.grid(axis='y')


def plot_histogram_4dfs(df1, df2, df3, df4, hue_col, title):
    """
    Plot a bar plot with hue for four input DataFrames.

    Parameters:
        df1 (pd.DataFrame): First DataFrame for plotting.
        df2 (pd.DataFrame): Second DataFrame for plotting.
        df3 (pd.DataFrame): Third DataFrame for plotting.
        df4 (pd.DataFrame): Fourth DataFrame for plotting.
        hue_col (str): Column name for hue distinction.
        title (str): Title for the plot.

    Returns:
        None
    """
    
    df1[hue_col] = 'CI farnesylated'
    df2[hue_col] = 'CI geranylgeranylated'
    df3[hue_col] = 'CI control'
    df4[hue_col] = 'TE'
    res = pd.concat([df1, df2, df3, df4])
    
    ax = sns.barplot(x='x', y='y', data=res, hue=hue_col)

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    
    plt.xticks(rotation=90)

    ax.set_ylabel('% of proteins')
    ax.set_xlabel('position of most C-terminal peptide')
    ax.set_title(title)

    
def plot_peptide_positions(df1, df2, moiety, motif, ax):
    """
    Plot the positions of peptides based on input DataFrames.

    Parameters:
        df1 (pd.DataFrame): First DataFrame for plotting.
        df2 (pd.DataFrame): Second DataFrame for plotting.
        moiety (str): Type of protein moiety (e.g., 'farnesylated').
        motif (str): Motif name (e.g., 'CAAX').
        ax (matplotlib.axes._subplots.AxesSubplot): Axes object for plotting.

    Returns:
        None
    """
    
    df1[f'peptides of {moiety} proteins with {motif}'] = 'total extract'
    df2[f'peptides of {moiety} proteins with {motif}'] = 'click-it'

    res = pd.concat([df1, df2])
    sns.barplot(x='x', y='y', data=res, hue=f'peptides of {moiety} proteins with {motif}', ax=ax)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)  # Rotate x-axis labels
    ax.set_ylabel('N of proteins')
    ax.set_xlabel('position of most C-terminal peptide')
    ax.set_title('Positions of peptides')

    ax.set_ylim(0, 35)  # Set y-axis limit

    ax.grid(axis='y')
