import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="Plots the pca")
    
    # required arguments
    parser.add_argument(
        "-evec", action="store", type=str, required=True, dest="evec_file", 
        help='path to the evec file from the smartPCA output (position of '
        'each individual along eigenvectors)')
    parser.add_argument(
        "-eval", action="store", type=str, required=True, dest="eval_file", 
        help='path to the eval file from the smartPCA output (the ordered '
        'eigenvalues corresponding to the eigenvectors)')
    parser.add_argument(
        "-pop", action="store", type=str, required=True, dest="poplist_file", 
        help='path to the grouped population list tsv file')
    parser.add_argument(
        "-ind", action="store", type=str, required=True, dest="HG_ind_file", 
        help="eigenstrat ind file with the HG002 and HG005 samples")
    parser.add_argument(
        "-o", action="store", type=str, required=False, dest="outdir",
        help="path of the output directory")  

    args = parser.parse_args()
    return args
    
    
def load_eigenvectors(evec_file):
    """ Loads the evec file from the smartPCA output, which contains the 
    coordinates in PC space (eigenvectors). """
    col_names = ["Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", 
                 "PC8", "PC9", "PC10", "Population"]
    df_eigenvectors = pd.read_csv(evec_file, 
                         delim_whitespace=True, 
                         header=0, 
                         names=col_names)
    return df_eigenvectors


def load_eigenvalues(eval_file):
    """ Loads the eval file from the smartPCA output, which contains the 
    ordered eigenvalues corresponding to the eigenvectors. """
    df_eigenvalues = pd.read_csv(eval_file, 
                                 header=0, 
                                 names=["eigenvalue"])
    df_eigenvalues.insert(0, 'PC', range(1, 1 + len(df_eigenvalues)))
    return df_eigenvalues


def load_grouped_eurasian_populations(poplist_file):
    """ Reads the tsv file of the grouped populations from the human origins 
    dataset """
    col_names = ["Population", "color_index", "symbol_index"]
    df_eurasian_pop = pd.read_csv(poplist_file, 
                                  names=col_names, 
                                  sep="\t")
    return df_eurasian_pop


def load_HG_populations(HG_ind_file):
    """Reads the eigenstrat ind file to get the all populations of the HG002 
    and HG005 samples. Based on the population name, columns for coverage 
    depth, color_index, and symbol_index are added for plotting. """
    col_names = ["Name", "Sex", "Population"]
    df_HG_samples = pd.read_csv(HG_ind_file, 
                                sep="\s+", 
                                names=col_names)
    # remove the reference/giab samples
    df_HG_samples = df_HG_samples.loc[~df_HG_samples['Name'].str.contains("giab")]
    # add columns to specify sample, trim/merging tool and subsample fraction
    cols = ["Sample", "Tool", "Fraction"]
    df_HG_samples[cols] = df_HG_samples["Population"] \
                              .str.split("-", expand=True)
    # translate the subsample fractions into coverage depths
    fraction2depth = {"06":"0.25X",
                      "12":"0.5X",
                      "24":"1X",
                      "48":"2X",
                      "96":"4X",}
    df_HG_samples["CovDepth"] = df_HG_samples["Fraction"] \
                                    .apply(lambda x: fraction2depth[x])
    
    # Add color and symbol indeces for plotting
    df_HG_samples["color_index"] = df_HG_samples["Sample"].factorize()[0]
    df_HG_samples["symbol_index"] = df_HG_samples["Tool"].factorize()[0]
    
    # Only take the columns and rows that are specific for each population,
    # with all subsample seeds being in one population
    final_colums = ["Population", "Sample", "Tool", "CovDepth", "color_index", 
                    "symbol_index"]
    df_HG_populations = df_HG_samples[final_colums] \
                            .drop_duplicates() \
                            .sort_values(by=["Sample", "Tool"])
    return df_HG_populations


def add_explained_variance(df_eigenvalues):
    """ Adds a column that specifies explained variance per principal component 
    and a column that specifies the cumulative explained variance """
    sum_eval = df_eigenvalues["eigenvalue"].sum()
    df_eigenvalues = df_eigenvalues.assign(
        explained_var=lambda x: x.eigenvalue/sum_eval*100)
    df_eigenvalues = df_eigenvalues.assign(
        cumulative_var=lambda x: np.cumsum(x.explained_var))
    return df_eigenvalues


def plot_individual_explained_variance(df_eigenvalues, num_pc, outdir):
    """ Plot the individual explained variance of the first num_pc PCs""" 
    plt.figure(figsize=(10, 5))
    plt.bar(df_eigenvalues[:10].PC, 
            df_eigenvalues[:10].explained_var, 
            alpha=0.5, 
            label='used to calculate within-cluster variance')
    plt.bar(df_eigenvalues[10:num_pc].PC, 
            df_eigenvalues[10:num_pc].explained_var, 
            alpha=0.7, 
            label='not used to calculate within-cluster variance')
    # Set x axis limits
    ax = plt.gca()
    ax.set_xlim(0, num_pc+1)
    # add labels
    plt.ylabel('Individual explained variance [%]')
    plt.xlabel('Principal component')
    #plt.title(f"Explained variance of the first {num_pc} principal components "
    #          f"(Total number of principal components: {len(df_eigenvalues)})")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(f"{outdir}/pca_explained_variance_individual.png", 
                dpi='figure', 
                format="png")


def plot_cumulative_explained_variance(df_eigenvalues, outdir):
    """ Plot the cumulative explained variance of all PCs""" 
    plt.figure(figsize=(10, 5))
    plt.step(df_eigenvalues.PC, 
             df_eigenvalues.cumulative_var, 
             where='mid',
             label='Cumulative explained variance')
    # Set axes limits
    ax = plt.gca()
    ax.set_xlim(0, len(df_eigenvalues))
    ax.set_ylim(0, 101)
    # add labels
    plt.ylabel('Explained variance, cumulative [%]')
    plt.xlabel('Principal component')
    plt.tight_layout()
    plt.savefig(f"{outdir}/pca_explained_variance_cumulative.png", 
                dpi='figure', 
                format="png")


def _plot_eurasian_populations(df_HO_populations, df_eigenvectors):
    """Plot the eurasian populations from the Affymetrix Human Origins Array"""
    symbols = [".", "o", "v", "^", "<", ">", "D", "s", "*", "+", "x", "X", "P",
               "d", "1", "2", "3", "4", "p", "h", "$∗$",]
    colors = [u'#8dd3c7', 
              u'#bebada', 
              u'#fb8072', 
              u'#80b1d3', 
              u'#fdb462',
              u'#b3de69', 
              u'#e5c494', 
              u'#bc80bd', 
              u'#fccde5']
    for index, row in df_HO_populations.iterrows():
        df = df_eigenvectors[df_eigenvectors.Population == row["Population"]]
        color = colors[row["color_index"]]
        marker = symbols[row["symbol_index"]]
        plt.plot(-df["PC1"], 
                 -df["PC2"], 
                 color=color,
                 marker=marker,
                 fillstyle="none",
                 linestyle="none",
                 label=row["Population"]
                 )
    

def _plot_HG_samples(df_HG_populations, df_eigenvectors):
    """ Plot the HG002 and HG005 samples"""
    colors = ("saddlebrown", "dimgrey")
    symbols = ("^", "8", "s", "p", "P", "*", "h")
    
    for index, row in df_HG_populations.iterrows():
        df = df_eigenvectors[df_eigenvectors.Population == row["Population"]]
        color = colors[row["color_index"]]
        marker = symbols[row["symbol_index"]]
        label = f'{row["Sample"]}, {row["Tool"]}'
        plt.plot(-df["PC1"], 
                 -df["PC2"], 
                 color=color,
                 marker=marker, 
                 linestyle="none",
                 label=label,
                 )
    


def plot_pca(df_eigenvectors, df_eurasian_populations, df_HG_populations, 
             df_eigenvalues, cov_depth, outdir):
    """ Plot PC1 agaist PC2, flipping the x axis to make the correlation to 
    Geography more obvious"""
    
    plt.figure(figsize=(11, 11))
    
    # Plot the datapoints
    _plot_eurasian_populations(df_eurasian_populations, df_eigenvectors)
    _plot_HG_samples(df_HG_populations, df_eigenvectors)
    
    # Get the explained variance for PC1 and PC2
    PC1_variance = float(df_eigenvalues[df_eigenvalues.PC==1]["explained_var"])
    PC2_variance = float(df_eigenvalues[df_eigenvalues.PC==2]["explained_var"])
    
    # add labels
    plt.xlabel(f"PC1 (variance explained: {round(PC1_variance,2)}%)")
    plt.ylabel(f"PC2 (variance explained: {round(PC2_variance,2)}%)")
    plt.legend(loc=(1.1, 0), ncol=3)
    plt.title(f"Coverage depth {cov_depth}")
    lgd = plt.legend(loc=(1.1, 0), ncol=3)
    plt.savefig(f"{outdir}/pca_all_{cov_depth}.png", 
                dpi='figure', 
                format="png", 
                bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


def plot_pca_all_depths(df_eigenvectors, df_eurasian_populations, 
                        df_HG_populations, df_eigenvalues, outdir):
    """ Creates a PCA plot for each coverage depth (0.25X, 0.5X, 1X, 2X, 4X)"""
    for cov_depth in ["0.25X", "0.5X", "1X", "2X", "4X"]:
        cov_depth_mask = df_HG_populations["CovDepth"] == cov_depth
        df_HG_populations_one_depth = df_HG_populations[cov_depth_mask]
        plot_pca(df_eigenvectors, 
                 df_eurasian_populations, 
                 df_HG_populations_one_depth, 
                 df_eigenvalues,
                 cov_depth, 
                 outdir)


def _plot_HG_sample_zoom(df_HG_populations, df_pca, ax, sample):
    """ Plot the samples (HG002 or HG005) from the fragmented, trimmed and   
    merged HG002 and HG005 samples"""
    df = df_HG_populations[df_HG_populations.Sample == sample]
    for index, row in df.iterrows():
        df_pop = df_pca[df_pca.Population == row["Population"]]
        label = row["Tool"]
        ax.plot(df_pop["PC1"], 
                df_pop["PC2"], 
                marker="o",
                linestyle="none",
                label=label,
                )
    ax.set_title(sample)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    plt.draw()
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        
        
def plot_pca_zoom(df_eigenvectors, df_HG_populations, cov_depth, outdir):
    """ flipping the axes to make the correlation to Geography more obvious
    """
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(11, 5))
    
    _plot_HG_sample_zoom(df_HG_populations, df_eigenvectors, ax0, "HG002")
    _plot_HG_sample_zoom(df_HG_populations, df_eigenvectors, ax1, "HG005")

    plt.legend(loc=(1.1, 0.3), ncol=1)
    fig.suptitle(f"Coverage depth {cov_depth}")
    fig.tight_layout()
    plt.savefig(f"{outdir}/pca_zoom_{cov_depth}.png", 
                dpi='figure', 
                format="png")
    
    
def plot_pca_all_depths_zoom(df_eigenvectors, df_HG_populations, outdir):
    """ Creates a PCA plot for each coverage depth (0.25X, 0.5X, 1X, 2X, 4X)"""
    for cov_depth in ["0.25X", "0.5X", "1X", "2X", "4X"]:
        cov_depth_mask = df_HG_populations["CovDepth"] == cov_depth
        df_HG_populations_one_depth = df_HG_populations[cov_depth_mask]
        plot_pca_zoom(df_eigenvectors, 
                      df_HG_populations_one_depth,
                      cov_depth,
                      outdir)
        

def within_cluster_variance(df):
    # Select only the columns starting with "PC"
    pc_cols = [col for col in df.columns if col.startswith('PC')]
    pc_df = df[pc_cols]
    # Calculate the centroid of the cluster
    centroid = pc_df.mean()
    # Calculate the variance of the distances from each point to the centroid
    variance = pc_df.sub(centroid).pow(2).sum().sum() / pc_df.size
    return variance


def get_df_of_within_cluster_variances(df_HG_populations, df_eigenvectors):
    """Calculates the within-cluster variance"""
    variance_rows = []
    for index, row in df_HG_populations.iterrows():
        # Get the population
        population = row["Population"]
        sample, tool, fraction = population.split("-")
        fraction2depth = {"06":"0.25X", 
                        "12":"0.5X", 
                        "24":"1X", 
                        "48":"2X", 
                        "96":"4X",}
        depth = fraction2depth[fraction]
        # calculate the variance
        df = df_eigenvectors[df_eigenvectors.Population == population]
        variance = within_cluster_variance(df)
        # store results
        variance_rows.append([sample, tool, depth, variance])
    
    # Create the pandas dataframe with all cluster variances
    columns = ["Sample", "Tool", "Depth", "Variance"]
    df_within_cluster_variance = pd.DataFrame(variance_rows, columns=columns)
    df_within_cluster_variance =  df_within_cluster_variance.sort_values(
        ['Sample', 'Depth'])
    
    return df_within_cluster_variance


def plot_within_cluster_variance(df_within_cluster_variance, outdir):
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))

    for sample in df_within_cluster_variance["Sample"].unique():
        mask_sample = df_within_cluster_variance["Sample"] == sample
        
        if sample == "HG002":
            ax = ax0
        elif sample == "HG005":
            ax = ax1
        else: 
            ax = None
            
        for tool in df_within_cluster_variance["Tool"].unique():
            mask_tool = df_within_cluster_variance["Tool"] == tool
            
            df = df_within_cluster_variance[mask_sample & mask_tool]
            ax.plot(df["Depth"], 
                    df["Variance"],
                    label=tool,
                    marker=".",
                    )
        ax.set_title(sample)
        ax.set_ylabel("Within-cluster Variance (first 10 PCs)")
        ax.set_xlabel("Coverage depth")
        ax.legend(loc='best')
        
    fig.tight_layout()
    plt.savefig(f"{outdir}/pca_within_cluster_variance.png", 
                dpi='figure', 
                format="png")


if __name__ == "__main__":
    
    args = parse_arguments()
    
    # load files
    df_eigenvectors = load_eigenvectors(args.evec_file)
    df_eigenvalues = load_eigenvalues(args.eval_file)
    df_eurasian_populations = load_grouped_eurasian_populations(
        args.poplist_file)
    df_HG_populations = load_HG_populations(args.HG_ind_file)
    
    # explained variance by each PC
    df_eigenvalues = add_explained_variance(df_eigenvalues)
    plot_individual_explained_variance(df_eigenvalues, 25, args.outdir)
    plot_cumulative_explained_variance(df_eigenvalues, args.outdir)
    
    # PC1 against PC2
    plot_pca_all_depths(df_eigenvectors, df_eurasian_populations, 
                        df_HG_populations, df_eigenvalues, args.outdir)
    plot_pca_all_depths_zoom(df_eigenvectors, df_HG_populations, args.outdir)
    
    # Within cluster variance
    df_within_cluster_variance = get_df_of_within_cluster_variances(
        df_HG_populations, df_eigenvectors)
    plot_within_cluster_variance(df_within_cluster_variance, args.outdir)
    
    