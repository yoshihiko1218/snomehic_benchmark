#!/usr/bin/env bash

IN_DIR="${1:?input directory required}"
OUT_DIR="${2:?output directory required}"

mkdir -p "$OUT_DIR"

echo "[INFO] Input : $IN_DIR"
echo "[INFO] Output: $OUT_DIR"

rename_one() {
    local infile="$1"
    local outfile="$2"
    local tag="$3"

    echo "[INFO] Renaming $(basename "$infile") -> $(basename "$outfile") [tag=$tag]"

    [[ -s "$infile" ]] || { echo "ERROR: missing or empty input: $infile" >&2; return 1; }

    gzip -t "$infile"

    local tmp="${outfile}.tmp.$$"
    rm -f "$tmp"

    zcat "$infile" | awk -v tag="$tag" '
    NR%4==1 {
        h=$0
        sub(/^@/, "", h)

        # keep only first token before whitespace
        split(h, a, /[[:space:]]+/)
        id=a[1]

        # strip existing pair suffix if present
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

    mv -f "$tmp" "$outfile"
}

shopt -s nullglob

n=0

for f in "$IN_DIR"/*.R1.fastq.gz; do
    [[ -e "$f" ]] || continue
    base="$(basename "$f")"
    prefix="${base%.R1.fastq.gz}"
    out="$OUT_DIR/${prefix}-R1.fq.gz"
    rename_one "$f" "$out" 1
    ((n+=1))
done

for f in "$IN_DIR"/*.R2.fastq.gz; do
    [[ -e "$f" ]] || continue
    base="$(basename "$f")"
    prefix="${base%.R2.fastq.gz}"
    out="$OUT_DIR/${prefix}-R2.fq.gz"
    rename_one "$f" "$out" 2
    ((n+=1))
done

echo "[INFO] Processed files: $n"
echo "[INFO] All done."
