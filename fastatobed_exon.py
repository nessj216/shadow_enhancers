from pathlib import Path
import re

fasta = Path("/Users/jillianness/Downloads/dmel-all-exon-r6.63.fasta")
bed   = Path("/Users/jillianness/Downloads/dmel-exons-r6.63.bed")



#               chr       start       end
coord_pat = re.compile(r"([A-Za-z0-9]+):[<>]?(\d+)\.\.[<>]?(\d+)")

with fasta.open() as fin, bed.open("w") as fout:
    for line in fin:
        if line.startswith(">"):
            m = coord_pat.search(line)
            if not m:                      # very rare: malformed header
                continue
            chrom, start, end = m.groups()
            start = int(start) - 1         # BED is 0-based, half-open
            exon_id = line[1:].split()[0]  # everything up to first space
            fout.write(f"{chrom}\t{start}\t{end}\t{exon_id}\n")
