import json

path = "qc.ipynb"
nb = json.load(open(path))

cell = nb["cells"][49]
src = "".join(cell["source"])

old_block = '''hic_snmCseq3 = _tidy_from_cols(
    snmCseq3, "snmCseq3",
    {
        "cis_n": "cis_n",
        "cis_gt1kb_n": "cis_gt1kb_n",
        "trans_ratio": "trans_ratio",
        "cis_per_million_uniqmapq30": "cis_per_million_uniqmapq30",
    }
)'''

new_block = '''# snmCseq3 comes from MappingSummary.csv.gz and has NO pre-computed Hi-C factors;
# derive them from the raw contact columns (same as the refined cell below).
_require_cols(
    snmCseq3,
    ["CisShortContact", "CisLongContact", "TransRatio",
     "R1UniqueMappedReads", "R2UniqueMappedReads"],
    "snmCseq3",
)
hic_snmCseq3 = pd.DataFrame(index=snmCseq3.index)
hic_snmCseq3["cis_n"] = snmCseq3["CisShortContact"] + snmCseq3["CisLongContact"]
hic_snmCseq3["cis_gt1kb_n"] = snmCseq3["CisLongContact"]
hic_snmCseq3["trans_ratio"] = snmCseq3["TransRatio"]
_mapped_snmCseq3 = snmCseq3[["R1UniqueMappedReads", "R2UniqueMappedReads"]].mean(axis=1)
hic_snmCseq3["cis_per_million_uniqmapq30"] = np.where(
    _mapped_snmCseq3 > 0,
    hic_snmCseq3["cis_gt1kb_n"] / _mapped_snmCseq3 * 1e6,
    np.nan,
)
hic_snmCseq3["dataset"] = "snmCseq3"'''

assert old_block in src, "cell 49 snmCseq3 _tidy_from_cols block not found"
src = src.replace(old_block, new_block)

cell["source"] = src.splitlines(keepends=True)
# outputs / execution_count untouched

json.dump(nb, open(path, "w"), indent=1)
print("patched cell 49; preserved outputs:", len(cell.get("outputs", [])))
