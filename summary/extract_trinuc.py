import os
import argparse
import pandas as pd

def extract_values(base_folder, suffix, output):
    rows = []

    for root, dirs, files in os.walk(base_folder):
        for fname in files:
            if fname.endswith(suffix):
                sample = fname.replace(suffix, "")
                methylation_path = os.path.join(root, fname)

                try:
                    # Read as 3 columns: trinuc, count, percent
                    methy_df = pd.read_csv(
                        methylation_path, sep='\t', header=None,
                        names=['trinuc', 'count', 'percent']
                    )

                    # Extract metrics safely (handle missing trinucs)
                    def get_percent(trinuc):
                        vals = methy_df.loc[methy_df['trinuc'] == f"{trinuc}:", 'percent']
                        return float(vals.values[0].strip('%')) if len(vals) > 0 else None

                    noncpg = get_percent('ACT')
                    endo   = get_percent('ACG')
                    exo    = get_percent('GCT')

                    rows.append([sample, noncpg, endo, exo])

                except Exception as e:
                    print(f"⚠️ Skipping {fname}: {e}")

    df = pd.DataFrame(rows, columns=['sample', 'noncpg', 'endo', 'exo'])
    df.to_csv(output, sep='\t', index=False)
    print(f"✅ Done! Found {len(rows)} samples. Saved to {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Recursively extract trinuc methylation (ACT/ACG/GCT) from chrM QC files."
    )
    parser.add_argument("--folder", required=True, help="Top-level folder containing .trinuc_methy.chrM.txt files")
    parser.add_argument("--suffix", required=True, help="File suffix to remove for sample names, e.g. .hg38.calmd.trinuc_methy.chrM.txt")
    parser.add_argument("--output", required=True, help="Output TSV file name")
    args = parser.parse_args()

    extract_values(args.folder, args.suffix, args.output)


    