##This code was written by Dr. Christine Sun. 

import pandas as pd
import sys

# --------------------------
# COMMAND-LINE ARGUMENTS
# --------------------------
if len(sys.argv) != 4:
    print("Usage: python script.py <hic_input_csv> <iphop_input_csv> <name_suffix>")
    sys.exit(1)

hic_path = sys.argv[1]        # Hi-C scores input file
iphop_path = sys.argv[2]      # iPHoP scores input file
suffix = sys.argv[3]          # Custom suffix for output files

# --------------------------
# LOAD INPUTS
# --------------------------
hic_raw = pd.read_csv(hic_path)
iphop_raw = pd.read_csv(iphop_path)

# --------------------------
# TAXONOMIC LEVELS
# --------------------------
tax_levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def ensure_tax_columns(df):
    """Ensure all taxonomic levels exist in the dataframe."""
    for level in tax_levels:
        if level not in df.columns:
            df[level] = pd.NA
    return df

# --------------------------
# COMPARISON FUNCTION
# --------------------------
def compare_taxonomy(df1, df2):
    """Compare taxonomy between two DataFrames by virus_name and return match summary."""
    df1 = ensure_tax_columns(df1)
    df2 = ensure_tax_columns(df2)

    shared_viruses = sorted(set(df1['virus_name']) & set(df2['virus_name']))
    match_rows = []

    for virus in shared_viruses:
        sub1 = df1[df1['virus_name'] == virus]
        sub2 = df2[df2['virus_name'] == virus]

        match_row = {'virus_name': virus}
        for level in tax_levels:
            set1 = set(sub1[level].dropna())
            set2 = set(sub2[level].dropna())

            if not set1 and not set2:
                match_row[f"{level}_match"] = pd.NA  # Both missing → NA
            else:
                match_row[f"{level}_match"] = len(set1 & set2) > 0  # True or False

        match_rows.append(match_row)

    return pd.DataFrame(match_rows)

# --------------------------
# 1) TOP SCORE ONLY
# --------------------------
def get_top_score(df, score_col='score'):
    # Keep all rows tied for the top score per virus_name
    return df[df[score_col] == df.groupby('virus_name')[score_col].transform('max')]

hic_top = get_top_score(hic_raw)
iphop_top = get_top_score(iphop_raw)
df_top = compare_taxonomy(hic_top, iphop_top)
df_top.to_csv(f"virus_matches_top_{suffix}.csv", index=False)

# --------------------------
# 2) FILTERED (original logic)
# --------------------------

# Hi-C: keep score >= 0.5, then top 20% per virus
hic_filt = hic_raw[hic_raw['score'] >= 0.5]
def filter_top_20_percent(group, pct=0.20):
    max_score = group['score'].max()
    threshold = max_score * (1 - pct)
    return group[group['score'] >= threshold]

hic_filt = hic_filt.groupby('virus_name', group_keys=False).apply(filter_top_20_percent)

# iPHoP: keep 90–100, then within 2 points of max
iphop_filt = iphop_raw[(iphop_raw['score'] >= 90) & (iphop_raw['score'] <= 100)]
def filter_within_2_points(group, max_diff=2):
    max_score = group['score'].max()
    threshold = max_score - max_diff
    return group[group['score'] >= threshold]

iphop_filt = iphop_filt.groupby('virus_name', group_keys=False).apply(filter_within_2_points)

df_filtered = compare_taxonomy(hic_filt, iphop_filt)
df_filtered.to_csv(f"virus_matches_filtered_{suffix}.csv", index=False)

# --------------------------
# 3) ALL DATA (no filtering)
# --------------------------
df_all = compare_taxonomy(hic_raw, iphop_raw)
df_all.to_csv(f"virus_matches_all_{suffix}.csv", index=False)

# --------------------------
# DONE
# --------------------------
print("✅ Done. Output files created:")
print(f" - virus_matches_top_{suffix}.csv")
print(f" - virus_matches_filtered_{suffix}.csv")
print(f" - virus_matches_all_{suffix}.csv")
