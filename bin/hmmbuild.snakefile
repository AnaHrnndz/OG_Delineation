import os
from pathlib import Path


def get_sample_names(infiles):
    s = set()
    for x in infiles:
        x=str(x)
        name_f = os.path.basename(x)
        s.add(name_f)
    return sorted(list(s))


cwd = os.getcwd()


infiles = Path(cwd).glob('*_*.aln')
samples = get_sample_names(infiles)


in_fasta = cwd+"/{sample}"
out_file_hmm = cwd+"/{sample}.hmm"


rule all:
    input:
        expand(out_file_hmm, sample = samples)

rule run_mafft_auto:
    input:
        in_fasta = in_fasta
    threads: 1
    
    output:
        out_file = out_file_hmm
    shell:
        "hmmbuild {output.out_file} {input.in_fasta}"