import json

path = "qc.ipynb"
nb = json.load(open(path))
c = nb["cells"][55]
s = "".join(c["source"])

old = """#   HCG (all 5 methylation methods): destrand (collapse +/- per CpG) + remove GCG
#     by reference genome. Calling = Bismark methylation extraction for
#     scnome/smallwood/snmCseq2, allcools for snmCseq3 (bhmem), BisSNP-1.0.1 NOMe
#     for scnomehic (the exception). No MAPQ-filter inconsistency.
#   GCH (GpC accessibility): NOMe methods ONLY -> scnome (Bismark NOMe.GpC),
#     scnomehic (BisSNP NOMe GCH.6plus2). N/A for the bisulfite-only methods."""

new = """#   HCG (all 5 methylation methods): destrand (collapse +/- per CpG) + remove GCG
#     by reference genome. Calling = Bismark methylation extraction for
#     scnome/smallwood; YAP/allcools for snmCseq2 (mm10 mouse cells, mates merged)
#     and snmCseq3 (bhmem); BisSNP-1.0.1 NOMe for scnomehic (the exception).
#     No MAPQ-filter inconsistency (avoids BisSNP -mmq 30).
#   GCH (GpC accessibility): NOMe methods ONLY -> scnome (Bismark, destranded),
#     scnomehic (BisSNP NOMe GCH.6plus2). N/A for the bisulfite-only methods."""

assert old in s, "old comment not found"
s = s.replace(old, new)
c["source"] = s.splitlines(keepends=True)
c["outputs"] = []
c["execution_count"] = None

json.dump(nb, open(path, "w"), indent=1)
print("updated cell 55 provenance comment")
