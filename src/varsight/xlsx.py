"""A small streaming .xlsx writer.

Sheets are written row by row into the zip, so a large VCF does not build up in memory. Every cell
is an inline string or number: nothing from the input can become a formula, whatever it starts with.
"""

from __future__ import annotations

import zipfile
from contextlib import closing
from pathlib import Path
from urllib.parse import urlsplit
from xml.sax.saxutils import escape, quoteattr

MAX_ROWS = 1_048_576
_MAX_LINK_LENGTH = 2000
_MAX_CELL = 32_767
_MAX_LINKS_PER_SHEET = 65_000
_INVALID_XML = {c: None for c in range(32) if c not in (9, 10, 13)}

# Style numbers for write_row.
NORMAL, HEADER, LINK, HIGH, MEDIUM, REVIEW, LOW, TITLE, BOLD, WRAP = range(10)

_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>{sheets}</Types>"""
_ROOT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>"""
_STYLES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
    '<fonts count="5">'
    '<font><sz val="11"/><name val="Calibri"/></font>'
    '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font>'
    '<font><u/><sz val="11"/><color rgb="FF0563C1"/><name val="Calibri"/></font>'
    '<font><b/><sz val="11"/><name val="Calibri"/></font>'
    '<font><b/><sz val="16"/><color rgb="FF1F3A5F"/><name val="Calibri"/></font>'
    "</fonts>"
    '<fills count="7">'
    '<fill><patternFill patternType="none"/></fill>'
    '<fill><patternFill patternType="gray125"/></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FF1F3A5F"/></patternFill></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FFF8D7DA"/></patternFill></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FFFFE8B3"/></patternFill></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FFD6E8F7"/></patternFill></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FFEDEFF2"/></patternFill></fill>'
    "</fills>"
    '<borders count="1"><border/></borders>'
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
    '<cellXfs count="10">'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
    '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
    '<xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="3" fillId="3" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
    '<xf numFmtId="0" fontId="3" fillId="4" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
    '<xf numFmtId="0" fontId="3" fillId="5" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
    '<xf numFmtId="0" fontId="0" fillId="6" borderId="0" xfId="0" applyFill="1"/>'
    '<xf numFmtId="0" fontId="4" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="3" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment wrapText="1" vertical="top"/></xf>'
    "</cellXfs>"
    '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
    "</styleSheet>"
)


def _column(index: int) -> str:
    name = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(65 + remainder) + name
    return name


def _text(value) -> str:
    text = str(value).translate(_INVALID_XML)
    return escape(text[:_MAX_CELL])


def _safe_link(value) -> bool:
    if not isinstance(value, str) or len(value) > _MAX_LINK_LENGTH:
        return False
    parts = urlsplit(value)
    return parts.scheme == "https" and bool(parts.netloc)


class XlsxWriter:
    """Write one workbook. Call add_sheet, then write_row repeatedly, and close at the end."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._zip = zipfile.ZipFile(self.path, "w", zipfile.ZIP_DEFLATED, allowZip64=True)
        self._sheets: list[str] = []
        self._stream = None
        self._row = 0
        self._columns = 0
        self._filter = False
        self._links: list[tuple[str, str]] = []

    def add_sheet(self, name: str, widths: list[int] | None = None, table: bool = False) -> None:
        """table=True freezes the first row and adds filter buttons to it."""
        self._finish_sheet()
        self._sheets.append(name)
        self._stream = self._zip.open(f"xl/worksheets/sheet{len(self._sheets)}.xml", "w", force_zip64=True)
        self._row = 0
        self._columns = 0
        self._filter = table
        self._links = []
        self._put('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
        self._put(
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        )
        if table:
            self._put(
                '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" '
                'state="frozen"/></sheetView></sheetViews>'
            )
        if widths:
            self._put("<cols>")
            for number, width in enumerate(widths, 1):
                self._put(f'<col min="{number}" max="{number}" width="{width}" customWidth="1"/>')
            self._put("</cols>")
        self._put("<sheetData>")

    def _put(self, text: str) -> None:
        self._stream.write(text.encode("utf-8"))

    def write_row(self, values: list, style: int = NORMAL, cell_styles: dict[int, int] | None = None) -> None:
        if self._row >= MAX_ROWS:
            raise ValueError(f"Excel holds at most {MAX_ROWS:,} rows; use TSV output for a file this large")
        self._row += 1
        self._columns = max(self._columns, len(values))
        cells = []
        for index, value in enumerate(values):
            if value is None or value == "":
                continue
            reference = f"{_column(index)}{self._row}"
            chosen = (cell_styles or {}).get(index, style)
            mark = f' s="{chosen}"' if chosen else ""
            if isinstance(value, bool):
                value = "yes" if value else "no"
            if isinstance(value, (int, float)):
                cells.append(f'<c r="{reference}"{mark}><v>{value}</v></c>')
            elif _safe_link(value) and style != HEADER:
                self._links.append((reference, value))
                cells.append(f'<c r="{reference}" s="{LINK}" t="inlineStr"><is><t>{_text(value)}</t></is></c>')
            else:
                cells.append(f'<c r="{reference}"{mark} t="inlineStr"><is><t xml:space="preserve">{_text(value)}</t></is></c>')
        self._put(f'<row r="{self._row}">{"".join(cells)}</row>')

    def _finish_sheet(self) -> None:
        if self._stream is None:
            return
        self._put("</sheetData>")
        if self._filter and self._row > 0 and self._columns > 0:
            self._put(f'<autoFilter ref="A1:{_column(self._columns - 1)}{self._row}"/>')
        sheet = len(self._sheets)
        kept = self._links[:_MAX_LINKS_PER_SHEET]
        if kept:
            self._put("<hyperlinks>")
            for number, (reference, _) in enumerate(kept, 1):
                self._put(f'<hyperlink ref="{reference}" r:id="rId{number}"/>')
            self._put("</hyperlinks>")
        self._put("</worksheet>")
        self._stream.close()
        self._stream = None
        if kept:
            relationships = "".join(
                f'<Relationship Id="rId{number}" '
                'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink" '
                f'Target={quoteattr(url)} TargetMode="External"/>'
                for number, (_, url) in enumerate(kept, 1)
            )
            self._zip.writestr(
                f"xl/worksheets/_rels/sheet{sheet}.xml.rels",
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                f'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">{relationships}</Relationships>',
            )
        self._links = []

    def close(self, order: list[str] | None = None) -> None:
        """order lists sheet names as they should appear as tabs; any others follow."""
        self._finish_sheet()
        numbered = {name: number for number, name in enumerate(self._sheets, 1)}
        tabs = [name for name in order or [] if name in numbered] + [name for name in self._sheets if name not in (order or [])]
        sheets = "".join(
            f'<Override PartName="/xl/worksheets/sheet{n}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            for n in numbered.values()
        )
        with closing(self._zip) as archive:
            archive.writestr("[Content_Types].xml", _CONTENT_TYPES.format(sheets=sheets))
            archive.writestr("_rels/.rels", _ROOT_RELS)
            archive.writestr("xl/styles.xml", _STYLES)
            entries = "".join(
                f'<sheet name={quoteattr(name[:31])} sheetId="{numbered[name]}" r:id="rId{numbered[name]}"/>' for name in tabs
            )
            archive.writestr(
                "xl/workbook.xml",
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
                f"<sheets>{entries}</sheets></workbook>",
            )
            relationships = "".join(
                f'<Relationship Id="rId{n}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                f'Target="worksheets/sheet{n}.xml"/>'
                for n in numbered.values()
            )
            styles_id = len(self._sheets) + 1
            archive.writestr(
                "xl/_rels/workbook.xml.rels",
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                f'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">{relationships}'
                f'<Relationship Id="rId{styles_id}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
                "</Relationships>",
            )
