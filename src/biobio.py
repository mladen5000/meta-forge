from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import os

# Parameters for simulated FASTQ file
num_reads = 100000  # small sample
read_length = 150

# Nucleotides
bases = ["A", "T", "C", "G"]

# Simulate reads
records = []
for i in range(num_reads):
    seq = "".join(random.choices(bases, k=read_length))
    qual = [random.randint(30, 40) for _ in range(read_length)]  # high quality
    record = SeqRecord(
        Seq(seq),
        id=f"sim_read_{i}",
        description="",
        letter_annotations={"phred_quality": qual},
    )
    records.append(record)

# Write to FASTQ
output_path = "sample_metagenome.fastq"
with open(output_path, "w") as f:
    SeqIO.write(records, f, "fastq")

output_path
