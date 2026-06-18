import json

path = "qc.ipynb"
nb = json.load(open(path))
c = nb["cells"][56]
s = "".join(c["source"])

# palette: add snmCAT color
s = s.replace('    "scnomehic": "#FF0000",\n}',
              '    "snmCAT": "#984ea3", "scnomehic": "#FF0000",\n}')
# HCG order: include snmCAT (all 6 methylation methods)
s = s.replace('order=["scnome", "smallwood", "snmCseq2", "snmCseq3", "scnomehic"],',
              'order=["scnome", "smallwood", "snmCseq2", "snmCseq3", "snmCAT", "scnomehic"],')
# GCH order: NOMe methods incl. snmCAT
s = s.replace('order=["scnome", "scnomehic"],',
              'order=["scnome", "snmCAT", "scnomehic"],')
assert '"snmCAT": "#984ea3"' in s, "palette edit failed"
assert '"snmCseq3", "snmCAT", "scnomehic"' in s, "HCG order edit failed"
assert '["scnome", "snmCAT", "scnomehic"]' in s, "GCH order edit failed"
c["source"] = s.splitlines(keepends=True)
c["outputs"] = []
c["execution_count"] = None

# cell 55 load comment: note 6 datasets incl snmCAT (NOMe)
c55 = nb["cells"][55]
s55 = "".join(c55["source"])
s55 = s55.replace("and snmCseq3 (bhmem); BisSNP-1.0.1 NOMe",
                  "snmCseq3 and snmCAT (bhmem/YAP); BisSNP-1.0.1 NOMe")
s55 = s55.replace("NOMe methods ONLY -> scnome (Bismark, destranded),\n#     scnomehic",
                  "NOMe methods ONLY -> scnome (Bismark), snmCAT (YAP), scnomehic")
c55["source"] = s55.splitlines(keepends=True)
c55["outputs"] = []
c55["execution_count"] = None

json.dump(nb, open(path, "w"), indent=1)
print("added snmCAT to notebook cells 55/56")
