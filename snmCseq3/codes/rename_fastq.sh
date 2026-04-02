#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash rename_one_sra_fastq_header.sh /path/to/input.fq.gz /path/to/output.fq.gz 1
#   bash rename_one_sra_fastq_header.sh /path/to/input.fq.gz /path/to/output.fq.gz 2
#
# Example:
#   bash rename_one_sra_fastq_header.sh test_fastq/SRR123-R1.fq.gz test_fastq_renamed/SRR123-R1.fq.gz 1

INFILE="${1:?input fastq.gz required}"
OUTFILE="${2:?output fastq.gz required}"
TAG="${3:?read tag required (1 or 2)}"

if [[ "$TAG" != "1" && "$TAG" != "2" ]]; then
    echo "ERROR: TAG must be 1 or 2, got: $TAG" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTFILE")"

echo "[INFO] Input : $INFILE"
echo "[INFO] Output: $OUTFILE"
echo "[INFO] Tag   : $TAG"

if [[ ! -s "$INFILE" ]]; then
    echo "ERROR: input file missing or empty: $INFILE" >&2
    exit 1
fi

gzip -t "$INFILE"

tmp="${OUTFILE}.tmp.$$"
rm -f "$tmp"

zcat "$INFILE" | awk -v tag="$TAG" '
NR%4==1 {
    h=$0
    sub(/^@/, "", h)
    split(h, a, " ")
    print "@" a[1] "_" tag
    next
}
NR%4==3 {
    print "+"
    next
}
{ print }
' | gzip > "$tmp"

gzip -t "$tmp"

python - "$tmp" <<'PY'
import sys, gzip
p = sys.argv[1]
with gzip.open(p, "rt") as f:
    for i in range(1000):
        h = f.readline().rstrip()
        if not h:
            break
        s = f.readline().rstrip()
        pl = f.readline().rstrip()
        q = f.readline().rstrip()
        assert h.startswith("@"), f"bad header: {h}"
        assert pl.startswith("+"), f"bad plus line: {pl}"
        assert len(s) == len(q), f"seq/qual mismatch: {h}"
print("OK")
PY

mv -f "$tmp" "$OUTFILE"
echo "[INFO] Done: $OUTFILE"