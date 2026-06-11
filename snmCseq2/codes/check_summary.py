import sys
import pandas as pd

f = sys.argv[1]
d = pd.read_csv(f)
print(f"== {f} ==")
print("cells:", d.shape[0], "metrics:", d.shape[1])

# pick a few representative QC columns if present
for col in ["FinalDNAReads", "R1UniqueMappedReads", "R2UniqueMappedReads",
            "R1MappingRate", "R2MappingRate", "DNAReadsYield",
            "CHNRate", "CGNRate", "CCCRate", "InputReads", "R1InputReads"]:
    if col in d.columns:
        s = pd.to_numeric(d[col], errors="coerce")
        print(f"{col:22s} median={s.median():.4g}  min={s.min():.4g}  max={s.max():.4g}")
