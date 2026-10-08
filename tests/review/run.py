#!/usr/bin/env python3
"""Run assembly/read-routing regressions using the pinned tool containers."""
import argparse
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--nextflow', default='nextflow')
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='genomeassembler-routing-') as temporary:
    work = Path(temporary)
    for name in ['hifi', 'short1', 'short2', 'hic1', 'hic2']:
        (work / f'{name}.fastq.gz').write_bytes(gzip.compress(b'@synthetic\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n'))
    (work / 'genome.fa').write_text('>synthetic\nACGTACGT\n')
    (work / 'nextflow.config').write_text(f'''
nextflow.enable.dsl = 2
params.outdir = '{work}/results'
params.publish_dir_mode = 'copy'
params.meryl_k = 21
params.omit_haps = false
docker.enabled = true
docker.runOptions = '--platform linux/amd64'
process {{ cpus = 1; memory = '2 GB'; time = '10 min' }}
includeConfig '{root}/conf/modules/assembly.config'
includeConfig '{root}/conf/modules/read-prep.config'
process {{ withName: HIFIASM_HIC {{ ext.when = {{ params.omit_haps.toString() != 'true' }} }} }}
''')
    (work / 'main.nf').write_text(f'''
include {{ ASSEMBLE }} from '{root}/subworkflows/local/assemble/main'
include {{ PREPARE_SHORTREADS }} from '{root}/subworkflows/local/prepare/prepare_shortreads/main'
include {{ PREPARE }} from '{root}/subworkflows/local/prepare/main'
workflow {{
    inputs = Channel.of('sampleA', 'sampleB').map {{ id ->
        [meta: [id: id, strategy: 'single', assembler_hifi: 'hifiasm', hifiasm_hic_phasing: true,
            hifireads: file('hifi.fastq.gz'), phasing_reads: [file('hic1.fastq.gz'), file('hic2.fastq.gz')],
            quast: false, busco: false, merqury: false, use_short_reads: false, lift_annotations: false]]
    }}
    ASSEMBLE(inputs, Channel.empty())
    ASSEMBLE.out.ch_main.toList().subscribe {{ rows ->
        if (params.omit_haps.toString() == 'true') return
        assert rows.collect {{ it.meta.id }}.sort() == ['sampleA-hap1', 'sampleA-hap2', 'sampleB-hap1', 'sampleB-hap2']
        assert rows.every {{ it.meta.assembly.name == "${{it.meta.id}}.fa.gz" }}
        assert rows.every {{ it.meta.hic_reads.size() == 2 }}
        println 'PASS: two phased samples produce four separately named assemblies'
    }}
    short_inputs = Channel.of(
        [id: 'raw-short', use_short_reads: true, paired: true, shortread_trim: false, merqury: true,
            shortread_F: file('short1.fastq.gz'), shortread_R: file('short2.fastq.gz'), meryl_k: 21],
        [id: 'raw-hic', use_short_reads: false, scaffold_hic: true, hic_trim: false,
            hic_F: file('hic1.fastq.gz'), hic_R: file('hic2.fastq.gz')],
        [id: 'trim-hic', use_short_reads: false, scaffold_hic: true, hic_trim: true,
            hic_F: file('hic1.fastq.gz'), hic_R: file('hic2.fastq.gz')]
    ).map {{ meta -> [meta: meta] }}
    PREPARE_SHORTREADS(short_inputs)
    PREPARE_SHORTREADS.out.main_out.toList().subscribe {{ rows ->
        assert rows.collect {{ it.meta.id }}.sort() == ['raw-hic', 'raw-short', 'trim-hic']
        assert rows.find {{ it.meta.id == 'raw-short' }}.meta.shortreads.size() == 2
        assert rows.find {{ it.meta.id == 'raw-hic' }}.meta.hic_reads*.name == ['hic1.fastq.gz', 'hic2.fastq.gz']
        assert rows.find {{ it.meta.id == 'trim-hic' }}.meta.hic_reads.size() == 2
        println 'PASS: untrimmed short reads and raw/trimmed Hi-C-only samples are retained'
    }}
    PREPARE_SHORTREADS.out.meryl_kmers.toList().subscribe {{ rows ->
        assert rows*.getAt(0) == ['raw-short']
        println 'PASS: Merqury database is emitted for untrimmed short reads'
    }}
    PREPARE(Channel.of([meta: [id: 'supplied', assembly: file('genome.fa'), qc_reads: 'ont', jellyfish: false]]))
    PREPARE.out.ch_main.toList().subscribe {{ rows ->
        assert rows*.meta*.id == ['supplied']
        assert rows[0].meta.assembly.name == 'genome.fa'
        println 'PASS: supplied assembly without long reads survives preparation'
    }}
}}
''')
    result = subprocess.run([args.nextflow, 'run', 'main.nf', '-stub-run', '-ansi-log', 'false'], cwd=work, env=os.environ,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout)
    if result.returncode:
        raise SystemExit(result.returncode)
    omitted = subprocess.run([args.nextflow, 'run', 'main.nf', '-stub-run', '-ansi-log', 'false', '--omit_haps', 'true'],
                             cwd=work, env=os.environ, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert omitted.returncode != 0 and 'mismatch' in omitted.stdout.lower(), omitted.stdout
    print('PASS: missing phased haplotypes fail the workflow')
