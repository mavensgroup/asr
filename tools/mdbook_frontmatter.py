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


def process_item(item):
    if 'Chapter' in item:
        ch = item['Chapter']
        ch['content'] = protect_math(clean_chapter_content(ch.get('content', '')))
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
