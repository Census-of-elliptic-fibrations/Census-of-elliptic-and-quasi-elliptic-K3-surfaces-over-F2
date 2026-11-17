"""Export the saved F_2 run in the published census key namespace and CSV format."""

import csv
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
MW = HERE.parent
RUN = MW / 'epyc-bundle/work/fresh/epyc-census/results/f2'
HEADER = (
    'key,a1,mw_torsion,singular_fibres,L_function,L_sign,analytic_rank,P2,'
    'zeta_sign,rho_F2,rho_geom,artin_height,br_disc,pts_F2,pts_F4,pts_F8,'
    'pts_F16,pts_F32'
)


def key_map():
    mapping = {}
    with (MW / 'verify-f2run/correspondence.csv').open() as source:
        for row in csv.DictReader(source):
            package = int(row['our_key'])
            assert package not in mapping
            mapping[package] = (int(row['upstream_key']), 'existing', '')
    assert len(mapping) == 1249569
    with (MW / 'elise/key-map.tsv').open() as source:
        for row in csv.DictReader(source, delimiter='\t'):
            package = int(row['package_key'])
            assert package not in mapping
            mapping[package] = (
                int(row['collaborator_key']), 'missing', row['unhandled']
            )
    assert len(mapping) == 1397598
    assert len({entry[0] for entry in mapping.values()}) == len(mapping)
    return mapping


def census_row(record, census_key):
    derived = record['derived']
    a1 = ' + '.join(
        term for bit, term in ((4, 's^2'), (2, 's*t'), (1, 't^2'))
        if (census_key >> 29) & bit
    ) or '0'
    torsion = ' x '.join(f'Z/{n}' for n in record['generators']['torsion'])
    fibres = sorted(
        (local['degree'], local['kodaira'], local['conductor'], local['tamagawa'])
        for local in record['local']
    )
    values = [
        census_key, a1, torsion or 'trivial',
        ';'.join(f'{degree}:{conductor}:{tamagawa}:{kodaira}'
                 for degree, kodaira, conductor, tamagawa in fibres),
        ' '.join(map(str, record['L'])), record['sign'], derived['rank'],
        ' '.join(map(str, derived['P2'])), derived['zeta_sign'],
        derived['rho'], derived['rho_geom'],
        'inf' if derived['height'] == 0 else derived['height'],
        derived['br_disc'], *derived['points'],
    ]
    assert len(values) == 18
    return ','.join(map(str, values))


def main():
    mapping = key_map()
    rows = {}
    shards = sorted((RUN / 'shards').glob('*.jsonl.gz'))
    assert len(shards) == 371
    for index, shard in enumerate(shards, 1):
        with gzip.open(shard, 'rt', encoding='ascii') as source:
            for line in source:
                record = json.loads(line)
                package = record['key']
                census_key, source_set, unhandled = mapping[package]
                assert census_key not in rows
                assert record['codes'] == [
                    (package >> 29) & 7, 0, (package >> 22) & 127,
                    (package >> 13) & 511, package & 8191,
                ]
                rows[census_key] = (package, census_row(record, census_key))
        if index % 40 == 0 or index == len(shards):
            print(f'shards={index}/{len(shards)} rows={len(rows)}', flush=True)
    assert len(rows) == len(mapping)
    with (HERE / 'census.csv').open('w', encoding='ascii', newline='') as full, \
         (HERE / 'key-map.tsv').open('w', encoding='ascii', newline='') as keyfile:
        full.write(HEADER + '\n')
        writer = csv.writer(keyfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['collaborator_key', 'package_key', 'source', 'unhandled'])
        for census_key in sorted(rows):
            package, row = rows[census_key]
            _, source_set, unhandled = mapping[package]
            full.write(row + '\n')
            writer.writerow([census_key, package, source_set, unhandled])
    print(f'EXPORT_OK rows={len(rows)}', flush=True)


if __name__ == '__main__':
    main()
