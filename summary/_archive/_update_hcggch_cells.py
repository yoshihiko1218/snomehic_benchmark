import json

path = "qc.ipynb"
nb = json.load(open(path))

load_src = '''# Per-cell detected HCG / GCH loci -- CONSISTENT method.
#   HCG (all 5 methylation methods): destrand (collapse +/- per CpG) + remove GCG
#     by reference genome. Calling = Bismark methylation extraction for
#     scnome/smallwood/snmCseq2, allcools for snmCseq3 (bhmem), BisSNP-1.0.1 NOMe
#     for scnomehic (the exception). No MAPQ-filter inconsistency.
#   GCH (GpC accessibility): NOMe methods ONLY -> scnome (Bismark NOMe.GpC),
#     scnomehic (BisSNP NOMe GCH.6plus2). N/A for the bisulfite-only methods.
# Built by each <tech>/codes/compute_hcg.py + summary/make_dataset_summaries.py.
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

gch_hcg = pd.read_csv("summary/gch_hcg_counts/all_methods.summary.txt", sep="\\t")
for c in ["GCH_n", "HCG_n"]:
    gch_hcg[c] = pd.to_numeric(gch_hcg[c], errors="coerce")

print("cells per dataset:\\n", gch_hcg["dataset"].value_counts())
print("\\nmedian loci per dataset:")
print(gch_hcg.groupby("dataset")[["HCG_n", "GCH_n"]].median())
gch_hcg.head()
'''

# cell 55 = load; replace its source, clear outputs
c55 = nb["cells"][55]
assert c55["cell_type"] == "code"
c55["source"] = load_src.splitlines(keepends=True)
c55["outputs"] = []
c55["execution_count"] = None

# cell 56 = plot; just refresh the two titles/ylabels to reflect GCG-removed HCG
c56 = nb["cells"][56]
s = "".join(c56["source"])
s = s.replace('title="Detected HCG (CpG) Loci per Cell"',
              'title="Detected HCG Loci per Cell (GCG removed)"')
s = s.replace('ylabel="HCG loci (count)"', 'ylabel="HCG loci (count)"')
s = s.replace('title="Detected GCH Loci per Cell"',
              'title="Detected GCH Loci per Cell (NOMe only)"')
c56["source"] = s.splitlines(keepends=True)
c56["outputs"] = []
c56["execution_count"] = None

json.dump(nb, open(path, "w"), indent=1)
print("updated cells 55 (load) and 56 (plot); total cells:", len(nb["cells"]))
