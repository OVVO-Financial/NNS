"""Regenerate the claim index and counterexamples in README.md."""
import csv
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PIN = '98bc42ea'
STATUS = {
    'PROVED-FAITHFUL': 'Proved',
    'PROVED-WITH-SIDE-CONDITION': 'Proved with conditions',
    'FALSE-AS-DOCUMENTED': 'False as documented',
    'RESIDUAL': 'Partly proved',
}


def source_link(reference):
    path, lines = reference.split(':')
    start, _, end = lines.partition('-')
    chapter = re.match(r'book/chapter-(\d+)-', path)
    label = f'Ch. {int(chapter[1])}' if chapter else Path(path).name
    anchor = f'L{start}-L{end}' if end else f'L{start}'
    return f'[{label}, {lines}](https://github.com/OVVO-Financial/NNS/blob/{PIN}/{path}#{anchor})'


def main():
    locations = {}
    for path in (ROOT / 'NNSProofs').rglob('*.lean'):
        for line, text in enumerate(path.read_text().splitlines(), 1):
            if match := re.match(r'^theorem claim_(C\d+)\b', text):
                claim = match[1]
                if claim in locations:
                    raise ValueError(f'Duplicate theorem for {claim}')
                locations[claim] = f'{path.relative_to(ROOT).as_posix()}#L{line}'

    with (ROOT / 'claims.tsv').open() as f:
        rows = sorted(csv.DictReader(f, delimiter='\t'), key=lambda row: row['claim'])
    claims = {row['claim'] for row in rows}
    if len(claims) != len(rows) or claims != locations.keys():
        raise ValueError('Claim IDs in claims.tsv and the Lean sources differ')

    readme = ROOT / 'README.md'
    heading = '\n## Claim index\n'
    intro, sep, _ = readme.read_text().partition(heading)
    if not sep:
        raise ValueError('README.md needs a claim index heading')

    counts = Counter(row['status'] for row in rows)
    partial = ', '.join(row['claim'] for row in rows if row['status'] == 'RESIDUAL')
    body = [
        f"The index covers {len(rows)} claims. Lean proves {counts['PROVED-FAITHFUL']} as stated",
        f"and {counts['PROVED-WITH-SIDE-CONDITION']} with conditions. It gives counterexamples for",
        f"{counts['FALSE-AS-DOCUMENTED']}. {partial} has a partial proof.",
        '',
        '### Counterexamples', '',
        'The links show each passage as it read before the corrections in the book.', '',
    ]
    for row in rows:
        if row['status'] == 'FALSE-AS-DOCUMENTED':
            body.append(f"- {row['claim']}. {row['note']} {source_link(row['sources'].split()[0])}")

    body.extend(['', '### Theorems', '',
                 f'| Claim | Theorem | Status | Documentation at `{PIN}` |',
                 '|---|---|---|---|'])
    for row in rows:
        claim, theorem = row['claim'], row['theorem']
        if theorem.split('.')[-1] != f'claim_{claim}':
            raise ValueError(f'{theorem} does not match {claim}')
        status = STATUS[row['status']]
        if row['status'] in {'PROVED-WITH-SIDE-CONDITION', 'RESIDUAL'}:
            status += '. ' + row['note']
        sources = ', '.join(source_link(ref) for ref in row['sources'].split())
        body.append(f"| {claim} · {row['topic']} | [`{theorem}`]({locations[claim]}) | {status} | {sources} |")
    readme.write_text(intro + heading + '\n' + '\n'.join(body) + '\n')


if __name__ == '__main__':
    main()
