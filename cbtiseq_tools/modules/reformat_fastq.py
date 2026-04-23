from Bio import SeqIO
from Bio.Seq import Seq
from cbtiseq_tools.utils import run_cmd


def parse_range(pos):
    if pos is None:
        return None
    start, end = map(int, pos.split("-"))
    return start - 1, end


def extract_region(seq, qual, pos):
    if pos is None:
        return "", []

    s, e = parse_range(pos)
    return seq[s:e], qual[s:e]


def remove_regions(seq, qual, regions):
    keep = [True] * len(seq)

    for region in regions:
        if region:
            s, e = parse_range(region)
            for i in range(s, min(e, len(seq))):
                keep[i] = False

    new_seq = "".join(seq[i] for i in range(len(seq)) if keep[i])
    new_qual = [qual[i] for i in range(len(seq)) if keep[i]]

    return new_seq, new_qual


def reformat_fastq(
    r1_path,
    r2_path,
    out_r1,
    out_r2,
    bcA=None,
    umiA=None,
    bcB=None,
    umiB=None
):
    count = 0

    with open(out_r1, "w") as f1_out, open(out_r2, "w") as f2_out:

        for rec1, rec2 in zip(
            SeqIO.parse(r1_path, "fastq"),
            SeqIO.parse(r2_path, "fastq")
        ):
            seq1 = str(rec1.seq)
            seq2 = str(rec2.seq)

            qual1 = rec1.letter_annotations["phred_quality"]
            qual2 = rec2.letter_annotations["phred_quality"]

            # 提取各部分 sequence + quality
            umia_seq, umia_qual = extract_region(seq1, qual1, umiA)
            umib_seq, umib_qual = extract_region(seq2, qual2, umiB)

            bca_seq, bca_qual = extract_region(seq1, qual1, bcA)
            bcb_seq, bcb_qual = extract_region(seq2, qual2, bcB)

            remain1_seq, remain1_qual = remove_regions(
                seq1, qual1, [bcA, umiA]
            )

            remain2_seq, remain2_qual = remove_regions(
                seq2, qual2, [bcB, umiB]
            )

            # 新 read1
            new_r1_seq = (
                umib_seq +
                umia_seq +
                bcb_seq +
                bca_seq +
                remain1_seq
            )

            new_r1_qual = (
                umib_qual +
                umia_qual +
                bcb_qual +
                bca_qual +
                remain1_qual
            )

            # 新 read2
            new_r2_seq = remain2_seq
            new_r2_qual = remain2_qual

            # 清空原 annotations
            rec1.letter_annotations = {}
            rec2.letter_annotations = {}

            rec1.seq = Seq(new_r1_seq)
            rec2.seq = Seq(new_r2_seq)

            rec1.letter_annotations["phred_quality"] = new_r1_qual
            rec2.letter_annotations["phred_quality"] = new_r2_qual

            SeqIO.write(rec1, f1_out, "fastq")
            SeqIO.write(rec2, f2_out, "fastq")

            count += 1

    print(f"[DONE] FASTQ reformat completed: {count} read pairs")