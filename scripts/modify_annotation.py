"""
Add custom transcript entries to a GTF annotation file.

Usage:
    Edit the new_entries list below to define your custom gene/transcript/exon.
    Then run: python modify_annotation.py

After running, sort the output GTF:
    awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' \
        updated_annotation.gtf > updated_sorted_annotation.gtf
"""


def add_custom_entries(gtf_file, output_file, new_entries):
    with open(gtf_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            outfile.write(line)
        outfile.write("\n".join(new_entries) + "\n")

    print(f"Updated GTF file saved as: {output_file}")


# Define your custom entries here.
# Format: chromosome, source, feature, start, stop, score, strand, frame, attributes
new_entries = [
    'CHROM\tcustom\tgene\tSTART\tSTOP\t.\tSTRAND\t.\tgene_id "GENE_ID"; gene_name "GENE_NAME"; gene_biotype "ncRNA";',
    'CHROM\tcustom\ttranscript\tSTART\tSTOP\t.\tSTRAND\t.\tgene_id "GENE_ID"; transcript_id "TRANSCRIPT_ID"; gene_name "GENE_NAME"; gene_biotype "ncRNA";',
    'CHROM\tcustom\texon\tSTART\tSTOP\t.\tSTRAND\t.\tgene_id "GENE_ID"; transcript_id "TRANSCRIPT_ID"; exon_number "1"; gene_name "GENE_NAME"; gene_biotype "ncRNA";',
]

add_custom_entries(
    gtf_file="/path/to/annotation.gtf",
    output_file="/path/to/updated_annotation.gtf",
    new_entries=new_entries,
)
