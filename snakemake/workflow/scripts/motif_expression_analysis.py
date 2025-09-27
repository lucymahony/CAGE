# Takes motif count per promoter 
# and expression df and looks for relationships

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, spearmanr, sem, ttest_ind
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.model_selection import cross_val_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
from scipy.stats import mannwhitneyu, sem
import re
import os
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


from scipy.stats import fisher_exact, sem

def map_gene_ids(gene_df, correspondence_df_file_path, from_version='v2.1', to_version='v1.1', gene_column='gene'):
    """
    Map gene IDs from one IWGSC RefSeq version to another using a correspondence table.
    """
    correspondence_df = pd.read_csv(correspondence_df_file_path, sep=' ')

    # Clean gene names 
    gene_df = gene_df.copy()
    correspondence_df = correspondence_df.copy()
    gene_df[gene_column] = gene_df[gene_column].str.strip()
    correspondence_df[from_version] = correspondence_df[from_version].astype(str).str.strip()
    correspondence_df[to_version] = correspondence_df[to_version].astype(str).str.strip()

    # Merge on source version
    mapped = pd.merge(
        gene_df,
        correspondence_df[[from_version, to_version]],
        left_on=gene_column,
        right_on=from_version,
        how='left'
    )
    num_total = gene_df.shape[0]
    num_unmapped = mapped[to_version].isna().sum()
    print(f"{num_unmapped} of {num_total} genes could not be mapped ({(num_unmapped / num_total) * 100:.2f}%)")
    mapped = mapped.dropna(subset=[to_version])

    # Replace gene column with mapped values
    mapped[gene_column] = mapped[to_version]
    mapped = mapped.drop(columns=[from_version, to_version])

    return mapped



def plot_motif_distribution(motif_matrix, output_file_path):
    """Plot histogram of total motif counts per gene."""
    motif_matrix['total_motifs'] = motif_matrix.drop('gene', axis=1).sum(axis=1)
    sns.histplot(motif_matrix['total_motifs'], bins=20)
    plt.title("Distribution of Total Motif Counts per Gene")
    plt.xlabel("Total Motifs per Gene")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_file_path)




def run_motif_continuous_enrichment(motif_matrix, tissue_specific, alpha=0.05):
    motif_matrix['gene'] = motif_matrix['gene'].str.strip()
    tissue_specific['gene'] = tissue_specific['gene'].str.strip()

    all_results = []
    tissues = tissue_specific['tissue'].unique()

    for tissue in tissues:
        tissue_specific_genes = tissue_specific[tissue_specific['tissue'] == tissue]['gene'].tolist()
        motif_long_df = motif_matrix.melt(id_vars='gene', var_name='motif', value_name='count')
        motif_long_df = motif_long_df.merge(tissue_specific[['gene', 'tissue']], on='gene', how='left')
        motif_long_df['group'] = np.where(motif_long_df['tissue'] == tissue, tissue, 'Other')

        results = []
        for motif in motif_matrix.columns[1:]:
            subset = motif_long_df[motif_long_df['motif'] == motif]
            ts_counts = subset[subset['group'] == tissue]['count']
            other_counts = subset[subset['group'] == 'Other']['count']

            if len(ts_counts) > 1 and len(other_counts) > 1:
                stat, pval = mannwhitneyu(ts_counts, other_counts, alternative='greater')
                results.append({
                    'tissue': tissue,
                    'motif': motif,
                    'u_stat': stat,
                    'raw_pval': pval
                })

        df_results = pd.DataFrame(results)
        if not df_results.empty:
            print(f'the number of motifs is {len(df_results)}')
            
            bonf_alpha = alpha / len(df_results)
            print(f'So the adjusted alpha is {bonf_alpha}')
            df_results['pval_adj'] = df_results['raw_pval']
            df_results['significant'] = df_results['pval_adj'] < bonf_alpha

            def star(p):
                if p <0.0014705882352941176:
                    return '*'
                else:
                    return ''

            df_results['star'] = df_results['raw_pval'].apply(star)
            all_results.append(df_results)

    return pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()

def plot_motif_enrichment_barplot(motif_matrix, tissue_specific, stats_df, outdir='.'):
    custom_palette = {
    'Other': '#878d96',  # grey
    'root': '#fa4d56',   
    'spike': '#ee5396',  
    'flag_leaf': '#a56eff',
    'grain': '#009d9a',
    'widespread': '#1192e8' }

    tissues = stats_df['tissue'].unique()

    for tissue in tissues:
        df = stats_df[stats_df['tissue'] == tissue]
        sig_motifs = df[df['significant']]['motif'].tolist()
        if not sig_motifs:
            continue

        # Prepare long-form data
        motif_long = motif_matrix.melt(id_vars='gene', value_vars=sig_motifs,
                                       var_name='motif', value_name='count')
        motif_long = motif_long.merge(tissue_specific[['gene', 'tissue']], on='gene', how='left')
        motif_long['group'] = np.where(motif_long['tissue'] == tissue, tissue, 'Other')

        # Compute summary stats
        summary = motif_long.groupby(['motif', 'group'])['count'].agg(mean='mean', sem=sem).reset_index()

        # Use seaborn barplot
        motif_order = sorted(sig_motifs)
        fig, ax = plt.subplots(figsize=(max(10, len(sig_motifs) * 0.6), 6))
        barplot = sns.barplot(
            data=summary,
            x='motif',
            y='mean',
            hue='group',
            ax=ax,
            order=motif_order,
            palette=custom_palette,
            err_kws={'linewidth': 0}  # disables internal error bars
        )

        # Compute group offsets
        group_levels = summary['group'].unique()
        n_groups = len(group_levels)
        bar_width = 0.8 / n_groups  # split full bar width

        motif_to_x = {motif: i for i, motif in enumerate(motif_order)}

        for i, motif in enumerate(motif_order):
            for j, group in enumerate(group_levels):
                row = summary[(summary['motif'] == motif) & (summary['group'] == group)]
                if row.empty:
                    continue
                y = row['mean'].values[0]
                yerr = row['sem'].values[0]
                x = i - 0.4 + j * bar_width + bar_width / 2
                ax.errorbar(x, y, yerr=yerr, fmt='none', ecolor='black', capsize=4, lw=1)

        # Add significance star
        for i, motif in enumerate(motif_order):
            row = df[df['motif'] == motif]
            if not row.empty and row['star'].notna().values[0]:
                stars = row['star'].values[0]
                if stars:
                    y_vals = summary[summary['motif'] == motif]['mean'] + summary[summary['motif'] == motif]['sem']
                    ymax = y_vals.max()
                    offset = ymax * 0.1
                    ax.plot([i - 0.2, i + 0.2], [ymax + offset] * 2, color='black', lw=1)
                    ax.text(i, ymax + offset * 1.1, stars, ha='center', va='bottom', fontsize=14)

        legend_labels = {
            'grain': 'Grain',
            'root': 'Root',
            'spike': 'Spike',
            'flag_leaf': 'Flag leaf',
            'Other': 'All other genes'
        }

        # Apply renamed labels
        handles, labels = ax.get_legend_handles_labels()
        new_labels = [legend_labels.get(l, l) for l in labels]
        ax.legend(handles, new_labels, frameon=False, title=None, fontsize=14)
        ax.set_ylabel('Motif Count', fontsize=20)
        ax.set_xlabel('')
        ax.set_xticks(range(len(motif_order)))
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.set_xticklabels(motif_order, rotation=45, ha='right')
        legend = ax.get_legend()
        legend.set_title(None)  
        legend.get_frame().set_linewidth(0) 

        plt.tight_layout()
        sns.despine()
        plt.savefig(f'{outdir}/barplot_{tissue}.png')
        plt.close()


# =====================
# CORRELATION ANALYSIS
# =====================

def correlation_with_expression(motif_matrix, expression_df, top_num=5, outdir='.'):
    motif_matrix['gene'] = motif_matrix['gene'].str.strip()
    expression_df['gene'] = expression_df['gene'].str.strip()

    # 2. Merge motif matrix and expression
    merged = motif_matrix.merge(expression_df, on='gene', how='inner')

    results = []
    for motif in motif_matrix.columns[1:]:  # skip 'gene'
        motif_vals = merged[motif]
        expr_vals = merged['gene']
        corr, pval = spearmanr(motif_vals, expr_vals)
        results.append({
            'motif': motif,
            'correlation': corr,
            'pvalue': pval
        })
    corr_df = pd.DataFrame(results)

    # --- Plot 1: Barplot of top motifs by abs(correlation) ---
    top_corrs = corr_df.reindex(corr_df['correlation'].abs().sort_values(ascending=False).index).head(top_num)
    print(f'TOP {top_num} correlations are:')
    print(top_corrs)
    plt.figure(figsize=(10, 6))
    barplot = sns.barplot(data=top_corrs, x='correlation', y='motif')
    for i, (corr, pval) in enumerate(zip(top_corrs['correlation'], top_corrs['pvalue'])):
        barplot.text(
            corr + 0.005 if corr > 0 else corr - 0.005,
            i,
            f"p = {pval:.3g}",
            color='black',
            ha='left' if corr > 0 else 'right',
            va='center',
            fontsize=10)
    plt.title('')
    plt.ylabel('')
    plt.xlabel('Spearman rank correlation')
    plt.xlim(-0.03, 0.03)
    plt.axvline(0, color='grey', linestyle='--')
    sns.despine()
    plt.savefig(os.path.join(outdir, f'motif_expression_top_{top_num}_correlations_barplot.png'))
    plt.close()

    return pd.DataFrame(results)


# =====================
# CLASSIFICATION ANALYSIS
# =====================

def motif_classification(motif_matrix, tissue_specific, target_tissue):
    """Use motifs to classify whether a gene is specific to a target tissue (logistic regression)."""
    tissue_specific['is_specific'] = (tissue_specific['tissue'] == target_tissue).astype(int)
    y_map = tissue_specific.set_index('gene')['is_specific']
    motif_sub = motif_matrix[motif_matrix['gene'].isin(y_map.index)].set_index('gene')
    y = y_map.loc[motif_sub.index]
    X = motif_sub.reset_index(drop=True)

    model = LogisticRegression(max_iter=1000)
    auc_scores = cross_val_score(model, X, y, cv=5, scoring='roc_auc')
    return auc_scores.mean()


# =====================
# PCA VISUALIZATION
# =====================

def pca_motif_clustering(motif_matrix, tissue_specific):
    """Perform PCA on motif matrix and plot 2D projection colored by tissue-specificity."""
    X = StandardScaler().fit_transform(motif_matrix.drop('gene', axis=1))
    pca = PCA(n_components=2)
    components = pca.fit_transform(X)
    df_pca = pd.DataFrame(components, columns=['PC1', 'PC2'])
    df_pca['gene'] = motif_matrix['gene']
    df_pca = pd.merge(df_pca, tissue_specific[['gene', 'tissue']], on='gene', how='left')

    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df_pca, x='PC1', y='PC2', hue='tissue', alpha=0.7)
    plt.title("PCA of Motif Profiles Colored by Tissue Specificity")
    plt.tight_layout()
    plt.savefig('pca')


def tsne_motif_clustering(motif_matrix, tissue_specific, sample_size=300, perplexity=30, random_state=42):
    """
    Perform t-SNE on motif matrix and plot 2D projection colored by tissue-specificity.
    
    Parameters:
    - motif_matrix: DataFrame with 'gene' column and motif counts
    - tissue_specific: DataFrame with 'gene' and 'tissue' assignments
    - sample_size: Number of widespread genes to sample for balance
    - perplexity: t-SNE perplexity parameter
    - random_state: Random seed for reproducibility
    """
    # Preprocessing
    motif_data = motif_matrix.drop(columns='gene')
    scaled_data = StandardScaler().fit_transform(motif_data)

    # Merge with tissue labels
    merged = motif_matrix[['gene']].copy()
    merged[['gene']] = motif_matrix[['gene']]
    merged = merged.merge(tissue_specific[['gene', 'tissue']], on='gene', how='left')
    merged['tissue'] = merged['tissue'].fillna('widespread')

    # Downsample widespread for visualization clarity
    widespread = merged[merged['tissue'] == 'widespread']
    others = merged[merged['tissue'] != 'widespread']
    if len(widespread) > sample_size:
        widespread = widespread.sample(n=sample_size, random_state=random_state)
    df_sampled = pd.concat([widespread, others])
    data_sampled = motif_matrix[motif_matrix['gene'].isin(df_sampled['gene'])]

    # t-SNE
    scaled_sampled = StandardScaler().fit_transform(data_sampled.drop(columns='gene'))
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=random_state)
    tsne_result = tsne.fit_transform(scaled_sampled)

    # Assemble plot DataFrame
    tsne_df = pd.DataFrame(tsne_result, columns=['TSNE1', 'TSNE2'])
    tsne_df['gene'] = data_sampled['gene'].values
    tsne_df = tsne_df.merge(df_sampled[['gene', 'tissue']], on='gene', how='left')

    # Plot
    plt.figure(figsize=(10, 6))
    palette = {
        'root': '#E69F00',
        'spike': '#56B4E9',
        'grain': '#009E73',
        'flag_leaf': '#CC79A7',
        'widespread': 'lightgray'
    }
    sns.scatterplot(data=tsne_df, x='TSNE1', y='TSNE2', hue='tissue', palette=palette, alpha=0.7, edgecolor=None)
    plt.title("t-SNE of Motif Profiles Colored by Tissue Specificity")
    plt.legend(title='', frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('tsne.pdf')

# =====================
# EXPRESSION PREDICTION
# =====================

def expression_prediction(motif_matrix, expression, tissue="flag_leaf"):
    """Use motif counts to predict expression levels in a specific tissue using Ridge regression."""
    expression_avg = (expression.groupby(['gene', 'tissue'])['expression_value'].mean().reset_index()
)

    expression_wide = expression_avg.pivot(index="gene", columns="tissue", values="expression_value").reset_index()
    merged = pd.merge(motif_matrix, expression_wide, on='gene')
    X = merged.drop(columns=['gene', tissue])
    y = merged[tissue]

    model = RidgeCV(alphas=[0.1, 1.0, 10.0])
    scores = cross_val_score(model, X, y, cv=5, scoring='r2')
    return scores.mean()


#
# For each motif print the mean and std counts per gene 

motif_matrix = pd.read_csv('motif_count_matrix.csv')
# Clean the motif names by removing 1- 2- e.c.t 
motif_cols = motif_matrix.columns[1:]
cleaned_cols = [re.sub(r'^\d+-', '', col) for col in motif_cols]
motif_matrix.columns = ['gene'] + cleaned_cols


list_motifs = [col for col in motif_matrix.columns if col != 'gene']
expression_matrix = pd.read_csv('../../../tmp/snakemake_intermediate_data/cadenza_mean_expression_per_tissue.csv')
tissue_specific_table = pd.read_csv('../../../tmp/snakemake_intermediate_data/cadenza_tissue_specific_list_z1.468448433935507.csv')
conversion_table = '../../../input_data/iwgsc_refseq_all_correspondances.csv'

# Expression matrix and tissue specific table, need to be converted from v1.1 to v2.1 

expression_matrix =  map_gene_ids(expression_matrix, conversion_table, from_version='v1.1', to_version='v2.1', gene_column='gene')
tissue_specific_table = map_gene_ids(tissue_specific_table, conversion_table, from_version='v1.1', to_version='v2.1', gene_column='gene')

#print(f'motif_matrix {motif_matrix}')
#print(f'expression{expression_matrix}')
#print(f'tissue_specific {tissue_specific_table}')
#print(f"There are the following tissues {tissue_specific_table['tissue'].unique()}")
#tissue_types = tissue_specific_table['tissue'].unique()
#for tissue in tissue_types:
#    print(f"There are {(tissue_specific_table[tissue_specific_table['tissue'] == tissue].shape[0])} genes assigned to {tissue}")
#
#for motif in list_motifs:
#    print(f'Motif {motif} has mean {motif_matrix[motif].mean()} and std {motif_matrix[motif].std()}')

# Print how many genes are assigned to 1 tissue ?

#plot_motif_distribution(motif_matrix, 'histograms.pdf')

#stats_df = run_motif_continuous_enrichment(motif_matrix, tissue_specific_table, alpha=0.05)
#plot_motif_enrichment_barplot(motif_matrix, tissue_specific_table, stats_df, '../../../results/')


print(expression_matrix.head())
average_expression_all_tissues= expression_matrix.groupby('gene', as_index=False)['expression_value'].mean()

corr = correlation_with_expression(motif_matrix, average_expression_all_tissues)




#auc = motif_classification(motif_matrix, tissue_specific, target_tissue="spike")
#tissue_specific = tissue_specific_table[tissue_specific_table['tissue']!='widespread']
#pca_motif_clustering(motif_matrix, tissue_specific)
#for tissue in ['root' ,'spike','grain','flag_leaf']:
#    r2 = expression_prediction(motif_matrix, expression_matrix, tissue=tissue)
#    print(tissue, r2)
#tsne_motif_clustering(motif_matrix, tissue_specific_table)

