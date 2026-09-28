import csv
import gzip
import os
import argparse

parser = argparse.ArgumentParser(description="Export the TaERF87 regulon AUCell score per cell.")
parser.add_argument("--input", default=os.environ.get("SCENIC_AUC_MATRIX", "auc_mtx.csv"))
parser.add_argument("--output", default="main_Fig2F_TaERF87_AUCell_per_cell.tsv.gz")
args = parser.parse_args()

SRC = args.input
OUT = args.output
TARGET = "TraesCS4A02G001300(+)"

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(SRC, "r", newline="", encoding="utf-8") as src, gzip.open(
    OUT, "wt", newline="", encoding="utf-8"
) as dst:
    reader = csv.reader(src)
    header = next(reader)
    if TARGET not in header:
        raise SystemExit(f"Missing target regulon: {TARGET}")
    idx = header.index(TARGET)
    writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
    writer.writerow(["cell_id", "TaERF87_regulon_AUCell"])
    n = 0
    for row in reader:
        writer.writerow([row[0], row[idx]])
        n += 1
print(f"wrote {n} cells to {OUT}")
