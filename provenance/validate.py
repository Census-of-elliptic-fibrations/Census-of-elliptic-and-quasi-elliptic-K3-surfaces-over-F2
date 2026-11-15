"""Check the exported CSV against both key maps, the published table and saved run data."""

import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile

HERE = Path(__file__).resolve().parent
MW = HERE.parent
RUN = MW / 'epyc-bundle/work/fresh/epyc-census/results/f2'
PUBLISHED_CENSUS = MW / 'upstream-org/census/census.csv.zip'
COLS = [
    'key', 'a1', 'mw_torsion', 'singular_fibres', 'L_function', 'L_sign',
    'analytic_rank', 'P2', 'zeta_sign', 'rho_F2', 'rho_geom', 'artin_height',
    'br_disc', 'pts_F2', 'pts_F4', 'pts_F8', 'pts_F16', 'pts_F32',
]
INT = r'(?:0|-?[1-9][0-9]*)'
INT_LIST = re.compile(INT + r'(?: ' + INT + r')*')
NAT = re.compile(r'(?:0|[1-9][0-9]*)')
FIBRE = re.compile(r'([1-9][0-9]*):([1-9][0-9]*):([1-9][0-9]*):(I[0-9]+\*?|II\*?|III\*?|IV\*?)')
A1 = {0: '0', 1: 't^2', 2: 's*t', 3: 's*t + t^2', 4: 's^2',
      5: 's^2 + t^2', 6: 's^2 + s*t', 7: 's^2 + s*t + t^2'}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def csv_lines(path):
    with path.open(encoding='ascii', newline='') as source:
        assert next(source) == ','.join(COLS) + '\n'
        for line in source:
            assert line.endswith('\n') and '\r' not in line
            fields = line[:-1].split(',')
            assert len(fields) == 18
            yield line, fields


def main():
    expected = {}
    with (MW / 'verify-f2run/correspondence.csv').open() as source:
        for row in csv.DictReader(source):
            key, package = int(row['upstream_key']), int(row['our_key'])
            assert key not in expected
            expected[key] = (package, 'existing', '')
    assert len(expected) == 1249569
    missing_keys = set()
    with (MW / 'elise/key-map.tsv').open() as source:
        for row in csv.DictReader(source, delimiter='\t'):
            key = int(row['collaborator_key'])
            assert key not in expected and row['unhandled'] in ('0', '1')
            expected[key] = (int(row['package_key']), 'missing', row['unhandled'])
            missing_keys.add(key)
    website_keys = set(map(int, (MW / 'grill/org-repos/k3-f2/data/missing_keys.txt').read_text().split()))
    assert missing_keys == website_keys and len(missing_keys) == 148029
    package_map = {item[0]: key for key, item in expected.items()}
    assert len(package_map) == len(expected) == 1397598

    previous = -1
    keymap_count = 0
    with (HERE / 'key-map.tsv').open() as source:
        reader = csv.DictReader(source, delimiter='\t')
        assert reader.fieldnames == ['collaborator_key', 'package_key', 'source', 'unhandled']
        for row in reader:
            key = int(row['collaborator_key'])
            assert key > previous
            previous = key
            assert (int(row['package_key']), row['source'], row['unhandled']) == expected[key]
            keymap_count += 1
    assert keymap_count == len(expected)

    rows = {}
    counts_a1 = collections.Counter()
    heights = collections.Counter()
    torsion_labels = collections.Counter()
    previous = -1
    for line, fields in csv_lines(HERE / 'census.csv'):
        key = int(fields[0])
        assert 0 <= key < 2**32 and key > previous and key in expected
        previous = key
        assert fields[1] == A1[key >> 29]
        assert fields[1] in ('0', 's^2', 's*t', 's^2 + s*t + t^2')
        assert fields[2] == 'trivial' or re.fullmatch(r'Z/[1-9][0-9]*(?: x Z/[1-9][0-9]*)*', fields[2])
        fibres = []
        if fields[3]:
            for token in fields[3].split(';'):
                match = FIBRE.fullmatch(token)
                assert match, (key, token)
                degree, conductor, tamagawa = map(int, match.groups()[:3])
                fibres.append((degree, match[4], conductor, tamagawa))
        assert fibres == sorted(fibres)
        for column in (4, 7):
            assert INT_LIST.fullmatch(fields[column])
            assert fields[column].split()[0] == '1'
        assert len(fields[7].split()) == 23
        assert fields[5] in ('-1', '1') and fields[8] in ('-1', '1')
        for column in (0, 6, 9, 10, 12, 13, 14, 15, 16, 17):
            assert NAT.fullmatch(fields[column])
        assert fields[11] == 'inf' or re.fullmatch(r'[1-9][0-9]*', fields[11])
        counts_a1[fields[1]] += 1
        heights[fields[11]] += 1
        torsion_labels[fields[2]] += 1
        rows[key] = line
    assert len(rows) == 1397598
    manifest = json.loads((MW / 'grill/org-repos/k3-f2/data/manifest.json').read_text())
    assert counts_a1 == manifest['expectedByA1']
    print(f'FORMAT_OK rows={len(rows)}', flush=True)

    compared = 0
    with zipfile.ZipFile(PUBLISHED_CENSUS) as archive, archive.open('census.csv') as stream:
        reader = csv.reader(io.TextIOWrapper(stream, encoding='ascii', newline=''))
        assert next(reader) == COLS
        seen = set()
        for fields in reader:
            key = int(fields[0])
            assert key not in seen and expected[key][1] == 'existing'
            seen.add(key)
            actual = rows[key].rstrip('\n').split(',')
            assert fields == actual, (key, [(COLS[i], x, y) for i, (x, y) in enumerate(zip(fields, actual)) if x != y])
            compared += 1
        assert compared == 1249569
    print(f'PUBLISHED_CENSUS_EXACT_OK rows={compared} fields={compared * 18}', flush=True)

    assert missing_keys <= rows.keys()
    nmissing = len(missing_keys)
    print(f'MISSING_KEYS_PRESENT_OK rows={nmissing}', flush=True)

    nrun = 0
    shards = sorted((RUN / 'shards').glob('*.jsonl.gz'))
    for shard in shards:
        with gzip.open(shard, 'rt', encoding='ascii') as source:
            for line in source:
                record = json.loads(line)
                package = record['key']
                fields = rows.pop(package_map[package]).rstrip('\n').split(',')
                codes = record['codes']
                assert codes == [package >> 29, 0, (package >> 22) % 128,
                                 (package >> 13) % 512, package % 8192]
                torsion = [] if fields[2] == 'trivial' else [int(x[2:]) for x in fields[2].split(' x ')]
                assert torsion == record['generators']['torsion']
                actual_fibres = sorted(tuple([int(x) for x in token.split(':')[:3]] + [token.split(':')[3]])
                                       for token in fields[3].split(';') if token)
                source_fibres = sorted((x['degree'], x['conductor'], x['tamagawa'], x['kodaira'])
                                       for x in record['local'])
                assert actual_fibres == source_fibres
                assert list(map(int, fields[4].split())) == record['L']
                assert int(fields[5]) == record['sign']
                d = record['derived']
                assert list(map(int, fields[7].split())) == d['P2']
                for index, name in ((6, 'rank'), (8, 'zeta_sign'), (9, 'rho'),
                                    (10, 'rho_geom'), (12, 'br_disc')):
                    assert int(fields[index]) == d[name]
                assert (0 if fields[11] == 'inf' else int(fields[11])) == d['height']
                assert list(map(int, fields[13:])) == d['points']
                nrun += 1
    assert nrun == 1397598 and not rows
    print(f'RUN_ROUNDTRIP_OK rows={nrun}', flush=True)

    source_paths = [
        PUBLISHED_CENSUS, MW / 'verify-f2run/correspondence.csv', MW / 'elise/key-map.tsv',
        MW / 'elise/src/functions_census.m',
        MW / 'grill/org-repos/k3-f2/data/manifest.json',
        MW / 'grill/org-repos/k3-f2/data/missing_keys.txt',
        RUN / 'summary.json', RUN / 'config.json', *shards,
    ]
    with (HERE / 'source-sha256.tsv').open('w', encoding='ascii') as out:
        out.write('sha256\tbytes\tpath_relative_to_artifacts_mw\n')
        for path in source_paths:
            out.write(f'{digest(path)}\t{path.stat().st_size}\t{path.relative_to(MW)}\n')
    result = {
        'validation': 'passed', 'rows': nrun, 'columns': COLS,
        'existing_rows_reproduced_field_for_field': compared,
        'existing_fields_identical_including_key': compared * 18,
        'missing_keys_present': nmissing, 'key_map_rows': keymap_count,
        'all_run_records_roundtrip': nrun, 'source_shards': len(shards),
        'unhandled_keys_present': sum(item[2] == '1' for item in expected.values()),
        'rows_by_a1': counts_a1, 'rows_by_artin_height': heights,
        'torsion_labels': torsion_labels,
        'format': {
            'header_exactly_matches_saved_census': True,
            'encoding': 'ASCII', 'line_ending': 'LF', 'row_order': 'numeric census key',
            'infinite_height': 'inf', 'polynomials': 'space-separated integer coefficients, constant term first',
            'fibres': 'degree:conductor:tamagawa:Kodaira; sorted by degree, Kodaira, conductor, tamagawa',
            'torsion': 'trivial or Z/n x Z/m',
        },
        'files': {name: {'sha256': digest(HERE / name), 'bytes': (HERE / name).stat().st_size}
                  for name in ('census.csv', 'key-map.tsv', 'source-sha256.tsv')},
    }
    (HERE / 'validation.json').write_text(json.dumps(result, indent=2) + '\n', encoding='ascii')
    print('VALIDATION_OK', flush=True)


if __name__ == '__main__':
    main()
