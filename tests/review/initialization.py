#!/usr/bin/env python3
"""Check samplesheet precedence and derived IDs through nf-schema initialization."""
import argparse
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--nextflow', default='nextflow')
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='genomeassembler-init-') as temporary:
    work = Path(temporary)
    (work / 'assets').mkdir()
    shutil.copy(root / 'assets/schema_input.json', work / 'assets/schema_input.json')
    shutil.copy(root / 'nextflow_schema.json', work / 'nextflow_schema.json')
    for name in ['hifi', 'hic1', 'hic2']:
        (work / f'{name}.fastq.gz').write_bytes(gzip.compress(b'@synthetic\nACGT\n+\nIIII\n'))
    (work / 'nextflow.config').write_text(f'''
includeConfig '{root}/nextflow.config'
params.input = '{work}/samples.csv'
params.outdir = '{work}/results'
params.hifiasm_hic_phasing = true
params.hic_trim = true
''')
    (work / 'main.nf').write_text(f'''
include {{ PIPELINE_INITIALISATION }} from '{root}/subworkflows/local/utils_nfcore_genomeassembler_pipeline/main'
workflow {{
    PIPELINE_INITIALISATION(false, false, true, [], params.outdir, params.input, false, false, false)
    PIPELINE_INITIALISATION.out.samplesheet.toList().subscribe {{ rows ->
        assert rows.size() == 3
        assert rows.find {{ it.meta.id == 'explicit' }}.meta.hifiasm_hic_phasing
        assert !rows.find {{ it.meta.id == 'disabled' }}.meta.hifiasm_hic_phasing
        assert rows.find {{ it.meta.id == 'inherited' }}.meta.hifiasm_hic_phasing
        assert rows.every {{ !it.meta.hic_trim && it.meta.qc_reads == 'hifi' }}
        println 'PASS: explicit false wins, missing phasing inherits CLI, and HiFi-only QC selects HiFi'
    }}
}}
''')
    header = 'sample,hifireads,hic_F,hic_R,hifiasm_hic_phasing,hic_trim\n'
    def row(name, phasing):
        return f'{name},{work}/hifi.fastq.gz,{work}/hic1.fastq.gz,{work}/hic2.fastq.gz,{phasing},false\n'
    (work / 'samples.csv').write_text(header + row('explicit', 'true') + row('disabled', 'false') + row('inherited', ''))
    def run():
        return subprocess.run([args.nextflow, 'run', 'main.nf', '-ansi-log', 'false'], cwd=work,
                              env=os.environ, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = run()
    if result.returncode:
        print(result.stdout)
        raise SystemExit(result.returncode)
    assert 'PASS:' in result.stdout
    print('PASS: sample boolean precedence, CLI inheritance and HiFi QC selection')
    (work / 'samples.csv').write_text(header + row('explicit', 'true') + row('explicit-hap1', 'false'))
    result = run()
    assert result.returncode != 0 and 'collide' in result.stdout, result.stdout
    print('PASS: derived haplotype/input sample ID collision rejected')
    (work / 'samples.csv').write_text(header + row('duplicate', 'false') + row('duplicate', 'false'))
    result = run()
    assert result.returncode != 0 and 'unique' in result.stdout, result.stdout
    print('PASS: duplicate samplesheet IDs rejected')
