"""Regenerate the claim index in README.md."""
import csv
import re
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PIN = '665e9828'
STATUS = {
    'PROVED-FAITHFUL': 'Proved',
    'PROVED-WITH-SIDE-CONDITION': 'Proved with conditions',
    'RESIDUAL': 'Partly proved',
}


def source_link(reference):
    path, lines = reference.split(':')
    start, _, end = lines.partition('-')
    chapter = re.match(r'book/chapter-(\d+)-', path)
    label = f'Ch. {int(chapter[1])}' if chapter else Path(path).name
    anchor = f'L{start}-L{end}' if end else f'L{start}'
    return f'[{label}, {lines}](https://github.com/OVVO-Financial/NNS/blob/{PIN}/{path}#{anchor})'


def check_theorems(rows, locations):
    pairs = []
    for row in rows:
        claim, name = row['claim'], row['theorem']
        if not re.fullmatch(r'(?:[A-Za-z_]\w*\.)+claim_' + re.escape(claim), name):
            raise ValueError(f'Invalid theorem name for {claim}: {name}')
        path = locations[claim].split('#')[0]
        module = '.'.join(Path(path).with_suffix('').parts)
        pairs.append(f'(`{name}, `{module})')
    checks = """import Lean
import NNSProofs
open Lean in
run_cmd
  let env ← getEnv
  let entries : List (Name × Name) := [""" + ',\n    '.join(pairs) + """]
  for (name, mod) in entries do
    unless env.find? name matches some (.thmInfo _) do
      throwError "missing theorem {name}"
    unless (env.getModuleIdxFor? name).map (env.header.moduleNames[·.toNat]!) == some mod do
      throwError "{name} is outside {mod}"
"""
    subprocess.run(['lake', 'build', 'NNSProofs'], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(['lake', 'env', 'lean', '--stdin'], cwd=ROOT, input=checks, text=True, check=True)


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

    check_theorems(rows, locations)

    readme = ROOT / 'README.md'
    heading = '\n## Claim index\n'
    intro, sep, _ = readme.read_text().partition(heading)
    if not sep:
        raise ValueError('README.md needs a claim index heading')

    unknown = {row['status'] for row in rows} - STATUS.keys()
    if unknown:
        raise ValueError(f'Unknown status: {", ".join(sorted(unknown))}')
    counts = Counter(row['status'] for row in rows)
    partial = ', '.join(row['claim'] for row in rows if row['status'] == 'RESIDUAL')
    body = [
        f"The index covers {len(rows)} claims. Lean proves {counts['PROVED-FAITHFUL']} as stated",
        f"and {counts['PROVED-WITH-SIDE-CONDITION']} with conditions."
        + (f' {partial} has a partial proof.' if partial else ''),
        '',
    ]
    body.extend(['### Theorems', '',
                 f'| Claim | Theorem | Status | Documentation at `{PIN}` |',
                 '|---|---|---|---|'])
    for row in rows:
        claim, theorem = row['claim'], row['theorem']
        status = STATUS[row['status']]
        if row['status'] in {'PROVED-WITH-SIDE-CONDITION', 'RESIDUAL'}:
            status += '. ' + row['note']
        sources = ', '.join(source_link(ref) for ref in row['sources'].split())
        body.append(f"| {claim} · {row['topic']} | [`{theorem}`]({locations[claim]}) | {status} | {sources} |")
    readme.write_text(intro + heading + '\n' + '\n'.join(body) + '\n')


if __name__ == '__main__':
    main()
