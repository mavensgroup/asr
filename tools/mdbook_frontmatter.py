#!/usr/bin/env python3
"""
mdbook-frontmatter: Preprocessor for mdBook that strips YAML front matter
(e.g., --- title: ... ---) from markdown chapters before rendering.
"""
import sys
import json
import re

def supports_renderer(renderer):
    return True

def clean_chapter_content(content):
    # Match YAML frontmatter at the very beginning of the document
    lines = content.splitlines(keepends=True)
    if lines and lines[0].strip() == '---':
        idx = 1
        while idx < len(lines) and lines[idx].strip() != '---':
            idx += 1
        if idx < len(lines):
            # Found closing ---, strip everything up to and including it
            return ''.join(lines[idx + 1:]).lstrip('\r\n')
    return content

MATH = re.compile(r'\\\(.*?\\\)|\\\[.*?\\\]', re.S)          # \( ... \) and \[ ... \]
PUNCT_ESC = re.compile(r'\\([!-/:-@\[-`{-~])')                 # backslash before ASCII punctuation
CODE_SPAN = re.compile(r'(`+)(?:.|\n)*?(?<!`)\1(?!`)')


def protect_math(content):
    """Keep TeX math intact through Markdown.  Markdown reads a backslash before punctuation as an escape, so \\( became (
    and MathJax saw no math.  Inside every math span, outside code, such backslashes are doubled; \\sigma etc. stay."""
    out, fence = [], None
    block = []

    def flush():
        text = ''.join(block)
        parts, last = [], 0
        for m in CODE_SPAN.finditer(text):                          # leave inline code alone
            parts.append(MATH.sub(lambda s: PUNCT_ESC.sub(r'\\\\\1', s.group(0)), text[last:m.start()]))
            parts.append(m.group(0))
            last = m.end()
        parts.append(MATH.sub(lambda s: PUNCT_ESC.sub(r'\\\\\1', s.group(0)), text[last:]))
        out.append(''.join(parts))
        block.clear()

    for line in content.splitlines(keepends=True):
        mark = line.lstrip()[:3]
        if fence is None and mark in ('```', '~~~'):
            flush(); fence = mark; out.append(line)
        elif fence is not None:
            out.append(line)
            if line.lstrip().startswith(fence):
                fence = None
        else:
            block.append(line)
    flush()
    return ''.join(out)


ADMONITION_RE = re.compile(r'^!!!\s+(note|warning|tip|danger|important)(?:\s+"([^"]*)")?\s*$')
WIKILINK_RE = re.compile(r'\[\[([A-Za-z0-9_]+)(?:\([^)]*\):[A-Za-z0-9_]+)?\]\]')
PROGRAMS = {'asr_main', 'term_bench'}


def convert_admonitions(content):
    """Convert Python-Markdown/FORD '!!! kind "title"' blocks into HTML admonition divs."""
    lines = content.splitlines(keepends=True)
    out = []
    i = 0
    fence = None
    while i < len(lines):
        line = lines[i]
        mark = line.lstrip()[:3]
        if fence is None and mark in ('```', '~~~'):
            fence = mark
            out.append(line)
            i += 1
            continue
        if fence is not None:
            out.append(line)
            if line.lstrip().startswith(fence):
                fence = None
            i += 1
            continue

        m = ADMONITION_RE.match(line.rstrip('\r\n'))
        if not m:
            out.append(line)
            i += 1
            continue

        kind = m.group(1).lower()
        title = m.group(2) if m.group(2) is not None else kind.capitalize()
        i += 1
        body = []
        while i < len(lines):
            nxt = lines[i]
            if nxt.strip() == '':
                body.append('\n')
                i += 1
            elif nxt.startswith('    ') or nxt.startswith('\t'):
                body.append(nxt[4:] if nxt.startswith('    ') else nxt[1:])
                i += 1
            else:
                break
        while body and body[-1].strip() == '':
            body.pop()
        out.append(f'<div class="admonition admonition-{kind}">\n<p class="admonition-title">{title}</p>\n\n')
        out.extend(body)
        if body and not body[-1].endswith('\n'):
            out.append('\n')
        out.append('\n</div>\n')
    return ''.join(out)


def convert_wikilinks(content, chapter_path):
    """Convert FORD [[symbol]] links outside code blocks/spans into relative mdBook API links."""
    depth = chapter_path.replace('\\', '/').count('/')
    prefix = '../' * depth

    def repl(m):
        name = m.group(1)
        subdir = 'program' if name in PROGRAMS else 'module'
        return f'[{name}]({prefix}api/{subdir}/{name}.html)'

    out, fence = [], None
    block = []

    def flush():
        text = ''.join(block)
        parts, last = [], 0
        for m in CODE_SPAN.finditer(text):
            parts.append(WIKILINK_RE.sub(repl, text[last:m.start()]))
            parts.append(m.group(0))
            last = m.end()
        parts.append(WIKILINK_RE.sub(repl, text[last:]))
        out.append(''.join(parts))
        block.clear()

    for line in content.splitlines(keepends=True):
        mark = line.lstrip()[:3]
        if fence is None and mark in ('```', '~~~'):
            flush(); fence = mark; out.append(line)
        elif fence is not None:
            out.append(line)
            if line.lstrip().startswith(fence):
                fence = None
        else:
            block.append(line)
    flush()
    return ''.join(out)


def process_item(item):
    if 'Chapter' in item:
        ch = item['Chapter']
        path = ch.get('path') or ch.get('source_path') or ''
        text = clean_chapter_content(ch.get('content', ''))
        text = convert_admonitions(text)
        text = convert_wikilinks(text, path)
        ch['content'] = protect_math(text)
        for sub in ch.get('sub_items', []):
            process_item(sub)

def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'supports':
        renderer = sys.argv[2] if len(sys.argv) > 2 else ''
        sys.exit(0 if supports_renderer(renderer) else 1)

    try:
        context, book = json.load(sys.stdin)
    except Exception as e:
        sys.stderr.write(f"mdbook-frontmatter JSON decode error: {e}\n")
        sys.exit(1)

    for item in book.get('items', []):
        process_item(item)

    json.dump(book, sys.stdout)

if __name__ == '__main__':
    main()
