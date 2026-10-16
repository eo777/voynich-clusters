"""Read-only verification of the portable release; standard library only."""
from pathlib import Path
from urllib.parse import urlsplit, unquote
import hashlib
from html.parser import HTMLParser
import json

HERE = Path(__file__).resolve().parent


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.base = ''; self.links = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'base': self.base = a.get('href', '')
        for key in ('href', 'src'):
            if key in a and tag != 'base': self.links.append(a[key])


def main():
    manifest = json.loads((HERE/'manifest.json').read_text(encoding='utf-8'))
    checked = 0
    for name, entry in manifest['files'].items():
        path = HERE/name
        assert path.resolve().is_relative_to(HERE), name
        data = path.read_bytes()
        assert len(data) == entry['bytes'], name
        assert hashlib.sha256(data).hexdigest() == entry['sha256'], name
        if path.suffix == '.html':
            links = Links(); links.feed(data.decode('utf-8'))
            base = path.parent/links.base
            for ref in links.links:
                parts = urlsplit(ref)
                if parts.scheme or parts.netloc or not parts.path: continue
                target = (base/unquote(parts.path)).resolve()
                assert target.is_file(), (name, ref)
                checked += 1
    # Frozen scientific path keys use the original host's separator.
    root = HERE/'replication'
    study = root/'solve/botanical_enrichment/team_golf/golf_published_stability_20261016'
    frozen = 0; external = []
    for filename in ('ensemble_index.json', 'scegret_inputs.json'):
        obj = json.loads((study/filename).read_text(encoding='utf-8'))
        for name, digest in obj['source_sha256'].items():
            normalized = name.replace('\\', '/')
            if '/_deps/' in normalized:
                external.append(normalized); continue
            assert hashlib.sha256((root/normalized).read_bytes()).hexdigest() == digest, name
            frozen += 1
    print(json.dumps(dict(status='passed', files=len(manifest['files']),
        html_links=checked, frozen_source_checks=frozen,
        external_pinned_package_files=external)))


if __name__ == '__main__': main()
