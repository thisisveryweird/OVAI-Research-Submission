"""Offline PDF text evidence: python extract_pdf_evidence.py INPUT.pdf NEW_OUTPUT_DIR.
Requires existing pdfinfo and pdftotext; refuses existing output files. No OCR.
"""
import argparse, hashlib, json, pathlib, shutil, subprocess

def extract_pdf_evidence(source, output):
    source, output = pathlib.Path(source).resolve(), pathlib.Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    info_cmd = ['pdfinfo', str(source)]
    info = subprocess.run(info_cmd, check=True, capture_output=True).stdout.decode('utf-8', errors='replace')
    count = int(next(line.split(':', 1)[1] for line in info.splitlines() if line.startswith('Pages:')))
    cmd = ['pdftotext', '-layout', '-enc', 'UTF-8', str(source), '-']
    raw = subprocess.run(cmd, check=True, capture_output=True).stdout
    text = raw.decode('utf-8', errors='strict')
    pages = text.split('\f')
    if pages[-1].strip() == '': pages.pop()
    if len(pages) != count: raise ValueError(f'Page mismatch: {len(pages)} vs {count}')
    evidence = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'source_bytes': source.stat().st_size, 'pdfinfo_pages': count, 'extracted_pages': len(pages), 'raw_text_sha256': hashlib.sha256(raw).hexdigest(), 'commands': [info_cmd, cmd], 'tools': {t: shutil.which(t) for t in ['pdfinfo', 'pdftotext']}, 'pages': [{'page': i, 'characters': len(p), 'replacement_characters': p.count('\ufffd'), 'sha256': hashlib.sha256(p.encode()).hexdigest()} for i,p in enumerate(pages, 1)], 'limitations': 'Text layer only; images/visual semantics and OCR are not verified.'}
    files = {'pdfinfo.txt': info, 'extracted-pages.txt': ''.join(f'===== PDF PAGE {i:02d} =====\n{p}\n' for i,p in enumerate(pages,1)), 'extraction-evidence.json': json.dumps(evidence, indent=2)}
    for name in files:
        if (output / name).exists(): raise FileExistsError(output / name)
    for name, content in files.items():
        with (output / name).open('x', encoding='utf-8', newline='\n') as f: f.write(content)
    return evidence

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source'); parser.add_argument('output')
    args = parser.parse_args()
    print(json.dumps(extract_pdf_evidence(args.source, args.output), indent=2))
