import json

path = "qc.ipynb"
nb = json.load(open(path))

cell = nb["cells"][47]
src = "".join(cell["source"])

old_line = '        to_metric(snmCseq3,      "mapping_rate",           "snmCseq3"),'
new_line = '        to_metric(snmCseq3,      "mapping_rate_mapq30",    "snmCseq3"),'
assert old_line in src, "snmCseq3 to_metric line not found"

marker = "# --- Combine into tidy table (UNIFY COLUMN NAMES HERE) ---"
assert marker in src, "combine marker not found"

derive_block = (
    "# snmCseq3 is redefined upstream (cell that loads MappingSummary.csv.gz) without a\n"
    "# 'mapping_rate' column. Derive the MAPQ>=30 mapping rate (mean of R1/R2) here so\n"
    "# this cell is self-contained and consistent with the high-quality rate plotted.\n"
    'snmCseq3["mapping_rate_mapq30"] = (\n'
    '    snmCseq3["R1MappingRateMapQ30"] + snmCseq3["R2MappingRateMapQ30"]\n'
    ") / 2\n\n"
)

src = src.replace(old_line, new_line)
src = src.replace(marker, derive_block + marker)

cell["source"] = src.splitlines(keepends=True)
# outputs and execution_count left untouched

json.dump(nb, open(path, "w"), indent=1)
print("patched cell 47; outputs preserved:", len(cell.get("outputs", [])), "output items")
