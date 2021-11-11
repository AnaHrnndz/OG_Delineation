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


infiles = Path(cwd).glob('*_*.faa')
samples = get_sample_names(infiles)

print(samples)

in_fasta = cwd+"/{sample}"
out_file_alg = cwd+"/{sample}.aln"


rule all:
    input:
        expand(out_file_alg, sample = samples)

rule run_mafft_auto:
    input:
        in_fasta = in_fasta
    threads: 1
    
    output:
        out_file = out_file_alg
    shell:
        "mafft --auto {input.in_fasta} > {output.out_file}"