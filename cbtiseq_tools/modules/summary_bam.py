import pandas as pd
import pysam
from collections import defaultdict


def count_reads_per_barcode(bam_file, tag="XC"):
    counts = defaultdict(int)

    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for read in bam:
            if read.has_tag(tag):
                counts[read.get_tag(tag)] += 1

    return counts


def generate_summary(raw_bam, trimmed_bam, dge_summary, output_file):
    raw_counts = count_reads_per_barcode(raw_bam)
    trim_counts = count_reads_per_barcode(trimmed_bam)

    df = pd.read_csv(dge_summary, sep="\t", skiprows=6)

    data = []

    for _, row in df.iterrows():
        bc = row["CELL_BARCODE"]

        genic = row["NUM_GENIC_READS"]
        umi = row["NUM_TRANSCRIPTS"]

        sat = (1 - umi / genic) * 100 if genic > 0 else 0

        data.append({
            "CELL_BARCODE": bc,
            "RAW_READS": raw_counts.get(bc, 0),
            "TRIM_READS": trim_counts.get(bc, 0),
            "NUM_GENIC_READS": genic,
            "NUM_TRANSCRIPTS": umi,
            "NUM_GENES": row["NUM_GENES"],
            "Sequencing Saturation": sat
        })

    pd.DataFrame(data).to_csv(
        output_file,
        sep="\t",
        index=False
    )

    print(f"[DONE] Summary saved: {output_file}")
