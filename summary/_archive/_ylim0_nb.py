import json

path = "qc.ipynb"
nb = json.load(open(path))
c = nb["cells"][56]
s = "".join(c["source"])

old = '''    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
'''
new = '''    ax.set_ylim(bottom=0)   # linear y-axis starting at 0 (log can't show 0)
'''
assert old in s, "log-scale block not found"
s = s.replace(old, new)
c["source"] = s.splitlines(keepends=True)
c["outputs"] = []
c["execution_count"] = None
json.dump(nb, open(path, "w"), indent=1)
print("set HCG/GCH violins to linear y-axis with ylim bottom=0")
