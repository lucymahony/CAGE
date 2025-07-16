import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet
import pyranges as pr
import os
from tqdm import tqdm
import sys

def get_overlapping_regions(df, match_on, reduce_regions=False):
    """
    For one sample, return a set of regions defined by the tuple of values in 'match_on'.
    If reduce_regions is False, each row becomes one tuple.
    
    If reduce_regions is True, the function uses PyRanges to merge intervals that overlap
    (on a per-group basis) and then returns the corresponding set of region tuples.
    
    Note: The grouping is done on the columns that appear in the intersection of match_on 
          and {"genes", "strand", "Chromosome"}, so that merging is done within each gene/chromosome/strand.
    """
    if not reduce_regions:
        return set([tuple(row[col] for col in match_on) for _, row in df.iterrows()])
    
    # Rename columns to match PyRanges expectations
    df = df.rename(columns={"seqnames": "Chromosome", "start": "Start", "end": "End"})

    # Create a PyRanges object; include only the columns we need.
    # Note: we add the "genes" column only if it is requested.
    cols = ['Chromosome', 'Start', 'End', 'strand']
    if "genes" in match_on:
        cols.append("genes")
    pr_obj = pr.PyRanges(df[cols])

    # Group the PyRanges dataframe on the intersection of (genes, strand, Chromosome)
    # (order is set in the following list)
    grouping_cols = ["genes", "strand", "Chromosome"]
    present_groups = list(set(match_on) & set(grouping_cols))
    grouped = []
    for key, group in pr_obj.df.groupby(present_groups):
        reduced = pr.PyRanges(group).merge()
        reduced_df = reduced.df.copy()
        # Restore *all* fields from the group key
        if not isinstance(key, tuple):
            key = (key,)
        for col, val in zip(present_groups, key):
            reduced_df[col] = val

        grouped.append(reduced_df)
    
    reduced_df = pd.concat(grouped)
    print(f"Reduced DataFrame: {reduced_df}")
    
    # Ensure all columns in match_on are present in reduced_df.
    for col in match_on:
        if col not in reduced_df.columns:
            reduced_df[col] = None  # Add missing columns with default value

    # Return a set of region tuples in the order provided in match_on.
    return set([tuple(row[col] for col in match_on) for _, row in reduced_df.iterrows()])


def regions_overlap(region_sample, region_global, match_on):
    """
    Checks whether a region from a sample (region_sample) overlaps a global region (region_global).
    
    Both region_sample and region_global are tuples with fields in order given by match_on.
    We assume that the coordinate columns are "Start" and "End". Overlap is defined as:
      Two regions (a,b) and (c,d) overlap if a <= d and c <= b.
    
    Also, we require that the non-coordinate columns (e.g. 'genes', 'Chromosome', 'strand') match.
    
    IMPORTANT: Here we assume intervals are inclusive of the end.
    """
    # Check that the non-coordinate columns match
    for col in match_on:
        if col in ("Start", "End"):
            continue
        # Compare the corresponding values by index.
        if region_sample[match_on.index(col)] != region_global[match_on.index(col)]:
            return False

    start_sample = int(region_sample[match_on.index("Start")])
    end_sample   = int(region_sample[match_on.index("End")])
    start_global = int(region_global[match_on.index("Start")])
    end_global   = int(region_global[match_on.index("End")])
    # Overlap condition (inclusive): sample overlaps global if sample.start <= global.end and global.start <= sample.end
    return (start_sample <= end_global) and (start_global <= end_sample)


def global_merge_regions(sample_sets, match_on):
    """
    Given a dictionary of sample region sets (each a set of tuples in order given by match_on),
    merge overlapping regions globally (across samples) by grouping on the non-coordinate columns.
    
    For example, if match_on = ["genes", "Chromosome", "Start", "End", "strand"], then the grouping keys are
    all except "Start" and "End". For each group, intervals are merged if they overlap.
    
    Returns a list of merged region tuples (in the order of match_on).
    """
    # Determine the grouping keys (non-coordinate columns)
    grouping_keys = [col for col in match_on if col not in ("Start", "End")]
    print(f"Grouping keys: {grouping_keys}")
    # Determine their indices in match_on
    group_indices = [match_on.index(key) for key in grouping_keys]
    
    # Collect intervals from all samples
    groups = {}
    for sample, regions in sample_sets.items():
        for region in regions:
            # Group key: tuple of values from region for the grouping keys
            key = tuple(region[i] for i in group_indices)
            start = int(region[match_on.index("Start")])
            end = int(region[match_on.index("End")])
            groups.setdefault(key, []).append((start, end))
    
    merged_global_regions = []
    for group_key, intervals in groups.items():
        # Sort intervals by start coordinate
        intervals.sort(key=lambda x: x[0])
        merged = []
        current_start, current_end = intervals[0]
        for st, en in intervals[1:]:
            # Merge if the next interval starts before or at the current interval's end (inclusive)
            if st <= current_end:
                current_end = max(current_end, en)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = st, en
        merged.append((current_start, current_end))
        
        # For each merged interval, rebuild a tuple in the order specified by match_on.
        # We create a dict for non-coordinate fields from group_key.
        group_dict = {col: group_key[i] for i, col in enumerate(grouping_keys)}
        for st, en in merged:
            new_region = []
            for col in match_on:
                if col == "Start":
                    new_region.append(st)
                elif col == "End":
                    new_region.append(en)
                else:
                    new_region.append(group_dict.get(col, None))
            merged_global_regions.append(tuple(new_region))
    return merged_global_regions

def upset_plot_regions(
    sample_files_dict,
    output_path,
    promoter_only=True,
    match_on=["genes", "Chromosome", "Start", "End", "strand"],
    merge_replicates=True,
    file_sep=",",
    reduce_regions=False):
    """
    Generates an UpSet plot of overlapping regions (e.g., promoter-tagged regions) across samples,
    with accurate handling of strand and interval merging.

    Returns a summary dict with counts of shared and unique regions.
    """
    sample_sets = {}

    # Process each sample's replicate files
    for sample, paths in sample_files_dict.items():
        if merge_replicates:
            dfs = [pd.read_csv(f, sep=file_sep) for f in paths]
            merged_df = pd.concat(dfs, ignore_index=True)

            # Normalize strand field just in case
            merged_df["strand"] = merged_df["strand"].astype(str).str.strip()

            if promoter_only:
                merged_df = merged_df[merged_df['annotation'] == 'promoter']

            sample_sets[sample] = get_overlapping_regions(merged_df, match_on, reduce_regions=reduce_regions)
        else:
            for i, f in enumerate(paths):
                rep_df = pd.read_csv(f, sep=file_sep)
                rep_df["strand"] = rep_df["strand"].astype(str).str.strip()

                if promoter_only:
                    rep_df = rep_df[rep_df['annotation'] == 'promoter']

                sample_id = f"{sample}_rep{i+1}"
                sample_sets[sample_id] = get_overlapping_regions(rep_df, match_on, reduce_regions=reduce_regions)

    # Merge overlapping regions globally
    print("Global merging of regions across samples...")
    global_merged = global_merge_regions(sample_sets, match_on)
    print(f"Number of global merged regions: {len(global_merged)}")

    # Build a boolean matrix indicating whether each sample overlaps each global region
    overlap_data = {}
    for sample, regions in sample_sets.items():
        sample_bool = []
        for global_region in tqdm(global_merged, desc=f"Checking overlaps for {sample}"):
            overlaps = any(regions_overlap(sample_region, global_region, match_on) for sample_region in regions)
            sample_bool.append(overlaps)
        overlap_data[sample] = sample_bool

    # DataFrame indexed by global region, with boolean presence per sample
    overlap_df = pd.DataFrame(
        overlap_data,
        index=pd.MultiIndex.from_tuples(global_merged, names=match_on)
    )

    # For UpSet: count each unique sample-overlap pattern
    boolean_series = overlap_df.groupby(list(overlap_df.columns), as_index=True).size()

    # Count regions shared in all samples
    shared_all = overlap_df[overlap_df.all(axis=1)]
    print(f"Number of regions shared in all samples: {len(shared_all)}")

    # Count regions unique to each sample
    unique_counts = {}
    for sample in sample_sets:
        is_unique = (overlap_df[sample]) & (~overlap_df.drop(columns=sample).any(axis=1))
        count = is_unique.sum()
        unique_counts[sample] = count
        print(f"Number of unique regions in {sample}: {count}")

    # Plot using UpSet
    upset = UpSet(boolean_series)
    upset.plot()
    plt.title("")
    plt.savefig(output_path, dpi=900)
    plt.close()
    print(f"UpSet plot saved to {output_path}")

    # Return summary stats
    return {
        "shared_all": len(shared_all),
        "unique_per_sample": unique_counts,
        "global_region_count": len(global_merged),
    }


def test_datasets():
    """
    Create test CSV files with sample data for Sample1 and Sample2.
    
    In this test data, Sample2 has two regions (280-291 and 290-300) that overlap.
    With the new global merge, these will be merged with any overlapping interval from another sample.
    """
    base_dir = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/test_datasets/"
    # Ensure the directory exists.
    os.makedirs(base_dir, exist_ok=True)
    
    files = ['sample1_rep1.csv', 'sample1_rep2.csv']
    for file in files:
        file_path = os.path.join(base_dir, file)
        with open(file_path, 'w') as f:
            f.write("\"seqnames\",\"start\",\"end\",\"width\",\"strand\",\"score\",\"dominant_ctss.seqnames\",\"dominant_ctss.pos\",\"dominant_ctss.strand\",\"dominant_ctss.score\",\"nr_ctss\",\"min_density\",\"max_density\",\"annotation\",\"genes\",\"nearest_gene\",\"gene\"\n")
            f.write("\"Chr1A\",100,110,11,\"+\",1.0,\"Chr1A\",100,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
            f.write("\"Chr1A\",200,210,11,\"+\",1.0,\"Chr1A\",200,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
            f.write("\"Chr1A\",300,310,11,\"+\",1.0,\"Chr1A\",300,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
    
    file = 'sample2_rep1.csv'
    with open(os.path.join(base_dir, file), 'w') as f:
        f.write("\"seqnames\",\"start\",\"end\",\"width\",\"strand\",\"score\",\"dominant_ctss.seqnames\",\"dominant_ctss.pos\",\"dominant_ctss.strand\",\"dominant_ctss.score\",\"nr_ctss\",\"min_density\",\"max_density\",\"annotation\",\"genes\",\"nearest_gene\",\"gene\"\n")
        f.write("\"Chr1A\",50,60,11,\"+\",1.0,\"Chr1A\",50,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
        f.write("\"Chr1A\",200,210,11,\"+\",1.0,\"Chr1A\",200,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
        f.write("\"Chr1A\",300,310,11,\"+\",1.0,\"Chr1A\",300,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
    
    file = 'sample2_rep2.csv'
    with open(os.path.join(base_dir, file), 'w') as f:
        f.write("\"seqnames\",\"start\",\"end\",\"width\",\"strand\",\"score\",\"dominant_ctss.seqnames\",\"dominant_ctss.pos\",\"dominant_ctss.strand\",\"dominant_ctss.score\",\"nr_ctss\",\"min_density\",\"max_density\",\"annotation\",\"genes\",\"nearest_gene\",\"gene\"\n")
        f.write("\"Chr1A\",105,106,11,\"-\",1.0,\"Chr1A\",100,\"-\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
        f.write("\"Chr1A\",200,210,11,\"+\",1.0,\"Chr1A\",200,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
        f.write("\"Chr1A\",280,291,11,\"+\",1.0,\"Chr1A\",280,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
        f.write("\"Chr1A\",290,300,11,\"+\",1.0,\"Chr1A\",290,\"+\",1.0,1,0.001,1.0,\"promoter\",\"TraesCS1A03G0000400.1\",\"TraesCS1A03G0000400.1\",1\n")
    
    # Prepare the dictionary of sample file paths.
    sample_files_dict = {
        "Sample1": [os.path.join(base_dir, f"sample1_rep{i}.csv") for i in [1, 2]],
        "Sample2": [os.path.join(base_dir, f"sample2_rep{i}.csv") for i in [1, 2]],
    }
    output_path = os.path.join(base_dir, "upset_plot.png")
    
    summary = upset_plot_regions(
        sample_files_dict,
        output_path,
        promoter_only=True,
        match_on=["genes", "Chromosome", "Start", "End", "strand"],
        merge_replicates=True,
        file_sep=",",
        reduce_regions=True)

    print("\n--- Summary ---")
    print(f"Shared in all samples: {summary['shared_all']}")
    for sample, count in summary["unique_per_sample"].items():
        print(f"{sample} unique: {count}")
    

if __name__ == "__main__":
    #test_datasets()
    # no merge reps
    base_dir = str(sys.argv[1]) # "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps/"
    clustering_name = str(sys.argv[2]) #"paraclu_md100_single_0.3"  # There is no leaf only root.
    #
    #cadenza_vs_fielder = {
    #    "cadenza": [
    #        os.path.join(base_dir, f"TC_GR_CR{i}_{clustering_name}_myCAGEset_my_annot.csv") for i in [1, 2, 3, 4, 5]
    #    ],
    #    "fielder": [
    #        os.path.join(base_dir, f"TC_GR_SRO{i}_{clustering_name}_myCAGEset_my_annot.csv") for i in [1, 2, 3]
    #    ]
    #}
    #upset_plot_regions(
    #    sample_files_dict=cadenza_vs_fielder,
    #    output_path=os.path.join(base_dir,f"{clustering_name}_upset_plot.png"),
    #    promoter_only=True,
    #    match_on=["genes", "Chromosome", "Start", "End", "strand"],
    #    merge_replicates=True,
    #    file_sep=",",
    #    reduce_regions=True)
    #
    fielder_reps = {'fielder_rep1': [os.path.join(base_dir, f"TC_GR_SRO1_{clustering_name}_myCAGEset_my_annot.csv")],
                        'fielder_rep2': [os.path.join(base_dir, f"TC_GR_SRO2_{clustering_name}_myCAGEset_my_annot.csv")],
                        'fielder_rep3': [os.path.join(base_dir, f"TC_GR_SRO3_{clustering_name}_myCAGEset_my_annot.csv")]}
    upset_plot_regions(
        sample_files_dict=fielder_reps,
        output_path=os.path.join(base_dir,f"{clustering_name}_fielder_reps_upset_plot.png"),
        promoter_only=True,
        match_on=["genes", "Chromosome", "Start", "End", "strand"],
        merge_replicates=False,
        file_sep=",",
        reduce_regions=True)
    #base_dir = str(sys.argv[1]) 
    #clustering_name = str(sys.argv[2])     
    #cadenza_vs_fielder = {
    #    "cadenza": [
    #        os.path.join(base_dir, f"TC_GR_CR_{clustering_name}_myCAGEset_my_annot.csv")
    #    ],
    #    "fielder": [
    #        os.path.join(base_dir, f"TC_GR_SRO_{clustering_name}_myCAGEset_my_annot.csv") 
    #    ]
    #}
    #upset_plot_regions(
    #    sample_files_dict=cadenza_vs_fielder,
    #    output_path=os.path.join(base_dir,f"{clustering_name}_upset_plot.png"),
    #    promoter_only=True,
    #    match_on=["genes", "Chromosome", "Start", "End", "strand"],
    #    merge_replicates=True,
    #    file_sep=",",
    #    reduce_regions=True)


