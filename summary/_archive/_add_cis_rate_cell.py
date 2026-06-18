import json

path = "qc.ipynb"
nb = json.load(open(path))

new_src = '''# --- cis > 1kb rate: fraction of cis contacts that are > 1kb ---
# Uses hic_df / helpers (winsorize_to_whiskers, violin_plot) defined above.
hic_df["cis_gt1kb_rate"] = np.where(
    hic_df["cis_n"] > 0,
    hic_df["cis_gt1kb_n"] / hic_df["cis_n"] * 100.0,
    np.nan,
)

hic_df_cis_gt1kb_rate = winsorize_to_whiskers(
    hic_df.dropna(subset=["cis_gt1kb_rate"]), "cis_gt1kb_rate"
)

print("cis > 1kb rate (%) summary by dataset:")
print(
    hic_df.groupby("dataset")["cis_gt1kb_rate"]
    .describe()
    .reindex(dataset_order)
)

violin_plot(
    hic_df_cis_gt1kb_rate, "cis_gt1kb_rate",
    title="Cis > 1 kb Rate per Cell",
    ylabel="Cis > 1 kb / Cis contacts (%)",
    log=False,
    ylim=(0, 100),
    fname="figures/hic_cis_gt1kb_rate_violin.pdf",
)
'''

# cell 53 is the empty trailing code cell
cell = nb["cells"][53]
assert cell["cell_type"] == "code" and "".join(cell["source"]).strip() == "", "cell 53 not empty code cell"
cell["source"] = new_src.splitlines(keepends=True)
cell["outputs"] = []
cell["execution_count"] = None

json.dump(nb, open(path, "w"), indent=1)
print("done; cells:", len(nb["cells"]))
