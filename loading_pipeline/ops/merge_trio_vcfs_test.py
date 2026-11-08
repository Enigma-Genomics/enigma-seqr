# Copyright 2026 Enigma contributors
# ruff: noqa: PLR2004
import json
import subprocess
import sys

import pysam
import pytest

from loading_pipeline.ops.merge_trio_vcfs import merge

HEADER = """##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihood">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSOURCE
"""


def fixture(tmp_path, rows=None):
    if rows is None:
        rows = [
            [
                'chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:AD:PL\t0|1:8,4:10,0,20',
                'chr1\t20\t.\tG\tA\t.\tPASS\t.\tGT\t1',
            ],
            ['chr1\t10\t.\tA\tG\t.\tPASS\t.\tDP:GT:AD\t12:1/.:3,9'],
            ['chr1\t10\t.\tA\tG,C\t.\tPASS\t.\tGT:AD\t1/2:1,3,4'],
        ]
    members = []
    for i, role in enumerate(('proband', 'mother', 'father')):
        p = tmp_path / f'{role}.vcf'
        p.write_text(HEADER + '\n'.join(rows[i]) + '\n')
        members.append(
            {
                'role': role,
                'path': str(p),
                'sample_id': role,
                'vcf_sample_id': 'SOURCE',
                'genome_version': 'GRCh38',
            },
        )
    return {'genome_version': 'GRCh38', 'members': members}


def test_merge_remapping_missing_and_stale_index(tmp_path):
    manifest = fixture(tmp_path)
    for m in manifest['members']:
        with open(m['path'] + '.tbi', 'w') as f:
            f.write('stale index must not be read')
    receipt = merge(manifest, tmp_path / 'out')
    assert receipt['records'] == 2
    assert receipt['joint_genotyped'] is False
    with pysam.VariantFile(str(tmp_path / 'out/callset.vcf.gz')) as vcf:
        first, second = list(vcf)
        assert first.alleles == ('A', 'C', 'G')
        assert first.samples['proband']['GT'] == (0, 1)
        assert first.samples['proband'].phased
        assert first.samples['mother']['GT'] == (2, None)
        assert first.samples['father']['GT'] == (2, 1)
        assert second.samples['proband']['GT'] == (1,)
        assert all(a is None for a in second.samples['mother']['GT'])
    assert json.loads((tmp_path / 'out/validation.json').read_text()) == receipt
    assert len(list((tmp_path / 'out').iterdir())) == 3


@pytest.mark.parametrize(
    'change',
    ['contigs', 'reference', 'sample', 'gvcf', 'gt', 'sort', 'info', 'format', 'qual'],
)
def test_rejections_leave_no_output(tmp_path, change):
    manifest = fixture(tmp_path)
    path = tmp_path / 'mother.vcf'
    text = path.read_text()
    if change == 'contigs':
        text = text.replace('length=1000', 'length=999')
    elif change == 'reference':
        text = text.replace('GRCh38', 'GRCh37')
    elif change == 'sample':
        text = text.replace('SOURCE', 'WRONG')
    elif change == 'gvcf':
        text = text.replace('\tG\t', '\t<NON_REF>\t')
    elif change == 'gt':
        text = text.replace('12:1/.', '12:3/.')
    elif change == 'sort':
        text += text.splitlines()[-1] + '\n'
    elif change == 'info':
        text = text.replace(
            '##reference',
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Count">\n##reference',
        ).replace('PASS\t.', 'PASS\tAC=1')
    elif change == 'format':
        text = text.replace('12:1/.:3,9', '12:1/.:3')
    elif change == 'qual':
        text = text.replace('\tG\t.\tPASS', '\tG\t50\tPASS')
    path.write_text(text)
    with pytest.raises(ValueError):
        merge(manifest, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()
    assert not list(tmp_path.glob('.trio-*'))


@pytest.mark.parametrize('change', ['target', 'role', 'source', 'build'])
def test_manifest_validation(tmp_path, change):
    m = fixture(tmp_path)
    if change == 'target':
        m['members'][1]['sample_id'] = 'proband'
    elif change == 'role':
        m['members'][1]['role'] = 'proband'
    elif change == 'source':
        m['members'][1]['path'] = m['members'][0]['path']
    else:
        m.pop('genome_version')
    with pytest.raises(ValueError):
        merge(m, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()


def test_immutable_destination_and_limits(tmp_path):
    m = fixture(tmp_path)
    out = tmp_path / 'out'
    out.mkdir()
    (out / 'keep').write_text('keep')
    with pytest.raises(ValueError):
        merge(m, out)
    assert (out / 'keep').read_text() == 'keep'
    with pytest.raises(ValueError):
        merge(m, tmp_path / 'limited', max_records=1)
    assert not (tmp_path / 'limited').exists()
    with pytest.raises(subprocess.TimeoutExpired):
        merge(m, tmp_path / 'timeout', timeout=0.00001)
    assert not (tmp_path / 'timeout').exists()


def test_cli_and_haploid_pl(tmp_path):
    m = fixture(
        tmp_path,
        [
            ['chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:PL\t1:20,0'],
            ['chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT:PL\t.:.,.'],
            ['chr1\t20\t.\tG\tA\t.\tPASS\t.\tGT\t./.'],
        ],
    )
    request = tmp_path / 'manifest.json'
    request.write_text(json.dumps(m))
    result = subprocess.run(  # noqa: S603 - fixed local module
        [
            sys.executable,
            '-m',
            'loading_pipeline.ops.merge_trio_vcfs',
            str(request),
            str(tmp_path / 'out'),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'out/validation.json').exists()
