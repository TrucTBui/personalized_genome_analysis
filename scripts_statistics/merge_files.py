import os
import pandas as pd

def merge_tsv_files_in_directory(dir):
    all_files = [f for f in os.listdir(dir) if f.endswith('.tsv')]
    if not all_files:
        print(f"No .tsv files found in directory: {dir}")
        return None
    dataframes = []
    for file in all_files:
        # if the file name contains 'merged', skip it
        if 'merged' in file:
            continue
        file_path = os.path.join(dir, file)
        df = pd.read_csv(file_path, sep='\t')
        dataframes.append(df)
    merged_df = pd.concat(dataframes, ignore_index=True)
    # reset the index of the merged DataFrame
    merged_df.reset_index(drop=True, inplace=True)
    return merged_df

FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
]
for person in FAMILY:
    dir = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Statistics_no_dup/{person}/"
    merged_df = merge_tsv_files_in_directory(dir)
    if merged_df is not None:
        output_file = os.path.join(dir, f"{person}_merged.txt")
        merged_df.to_csv(output_file, sep='\t', index=False)
        print(f"Merged file for {person} saved to {output_file}")
    