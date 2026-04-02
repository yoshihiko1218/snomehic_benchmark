#!/usr/bin/env bash

# Usage:
#   bash rename_fastq_single.sh INPUT.fastq.gz OUTPUT.fq.gz 1
#   bash rename_fastq_single.sh INPUT.fastq.gz OUTPUT.fq.gz 2

INFILE="${1:?input fastq.gz required}"
OUTFILE="${2:?output fastq.gz required}"
TAG="${3:?tag required (1 or 2)}"

if [[ "$TAG" != "1" && "$TAG" != "2" ]]; then
    echo "ERROR: TAG must be 1 or 2, got: $TAG" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTFILE")"

echo "[INFO] Input : $INFILE"
echo "[INFO] Output: $OUTFILE"
echo "[INFO] Tag   : $TAG"

[[ -s "$INFILE" ]] || { echo "ERROR: missing or empty input: $INFILE" >&2; exit 1; }

gzip -t "$INFILE"

tmp="${OUTFILE}.tmp.$$"
rm -f "$tmp"

zcat "$INFILE" | awk -v tag="$TAG" '
NR%4==1 {
    h=$0
    sub(/^@/, "", h)
    split(h, a, /[[:space:]]+/)
    id=a[1]

    sub(/\/1$/, "", id)
    sub(/\/2$/, "", id)
    sub(/_1$/, "", id)
    sub(/_2$/, "", id)

    print "@" id "_" tag
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
