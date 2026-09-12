#!/usr/bin/env python3
from __future__ import annotations

import argparse
from html import escape
from pathlib import Path
import re
import subprocess

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]


def inline(value: str) -> str:
    value = escape(value)
    value = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", value)
    value = re.sub(r"`([^`]+)`", r'<font color="#1f5b78">\1</font>', value)
    return value


def register_font(font_path: Path) -> None:
    if not font_path.exists():
        raise FileNotFoundError(
            f"Chinese TrueType font not found: {font_path}. "
            "Pass --font /path/to/a/Chinese.ttf"
        )
    pdfmetrics.registerFont(TTFont("CJK", str(font_path)))
    pdfmetrics.registerFontFamily("CJK", normal="CJK", bold="CJK")


def locate_font(explicit: Path | None) -> Path:
    if explicit is not None:
        return explicit
    result = subprocess.run(
        ["fc-match", "-f", "%{file}", "Noto Sans SC"],
        check=False,
        capture_output=True,
        text=True,
    )
    matched = Path(result.stdout.strip()) if result.stdout.strip() else None
    if matched is not None and matched.exists():
        return matched
    for candidate in (
        Path("C:/Windows/Fonts/msyh.ttc"),
        Path("C:/Windows/Fonts/simsun.ttc"),
    ):
        if candidate.exists():
            return candidate
    raise FileNotFoundError("No Chinese TrueType font found; pass --font /path/to/font.ttf")


def build(markdown: Path, output: Path, font_path: Path | None) -> None:
    font_path = locate_font(font_path)
    register_font(font_path)
    page_width, _ = A4
    content_width = page_width - 34 * mm
    base = dict(fontName="CJK", textColor=colors.HexColor("#252a30"), wordWrap="CJK")
    styles = {
        "title": ParagraphStyle(
            "title", fontName="CJK", fontSize=21, leading=29, textColor=colors.HexColor("#15364d"),
            alignment=TA_CENTER, spaceBefore=24, spaceAfter=18, wordWrap="CJK",
        ),
        "h2": ParagraphStyle(
            "h2", fontName="CJK", fontSize=14.5, leading=21, textColor=colors.HexColor("#174f6b"),
            spaceBefore=10, spaceAfter=7, keepWithNext=True, wordWrap="CJK",
        ),
        "h3": ParagraphStyle(
            "h3", fontName="CJK", fontSize=11.5, leading=17, textColor=colors.HexColor("#285d75"),
            spaceBefore=7, spaceAfter=4, keepWithNext=True, wordWrap="CJK",
        ),
        "body": ParagraphStyle("body", fontSize=9.4, leading=15, spaceAfter=6.5, **base),
        "bullet": ParagraphStyle("bullet", fontSize=9.3, leading=14.5, leftIndent=12, firstLineIndent=-8, spaceAfter=4, **base),
        "cell": ParagraphStyle("cell", fontSize=7.2, leading=10.3, spaceAfter=0, **base),
        "cell_head": ParagraphStyle(
            "cell_head", fontName="CJK", fontSize=7.3, leading=10.3, textColor=colors.white, wordWrap="CJK",
        ),
        "caption": ParagraphStyle("caption", fontSize=8.1, leading=12, alignment=TA_CENTER, textColor=colors.HexColor("#4d5962"), spaceAfter=7, fontName="CJK", wordWrap="CJK"),
        "code": ParagraphStyle(
            "code", fontName="CJK", fontSize=7.8, leading=11.5, backColor=colors.HexColor("#eef3f6"),
            borderPadding=7, spaceBefore=3, spaceAfter=7,
        ),
    }

    story = []
    lines = markdown.read_text(encoding="utf-8").splitlines()
    index = 0
    first_paragraph = True
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if line == "<!-- pagebreak -->":
            story.append(PageBreak())
        elif line.startswith("```"):
            block: list[str] = []
            index += 1
            while index < len(lines) and not lines[index].strip().startswith("```"):
                block.append(lines[index])
                index += 1
            story.append(Preformatted("\n".join(block), styles["code"]))
        elif line.startswith("!["):
            match = re.match(r"!\[(.*?)\]\((.*?)\)", line)
            if match:
                image_path = (markdown.parent / match.group(2)).resolve()
                picture = Image(str(image_path))
                scale = min(content_width / picture.imageWidth, 112 * mm / picture.imageHeight)
                picture.drawWidth = picture.imageWidth * scale
                picture.drawHeight = picture.imageHeight * scale
                picture.hAlign = "CENTER"
                story.extend([picture, Spacer(1, 3)])
        elif line.startswith("|"):
            rows: list[list[str]] = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                cells = [cell.strip() for cell in lines[index].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", cell) for cell in cells):
                    rows.append(cells)
                index += 1
            index -= 1
            column_count = len(rows[0])
            if column_count == 2:
                fractions = [0.27, 0.73]
            elif column_count == 3:
                fractions = [0.28, 0.25, 0.47]
            else:
                fractions = [0.20] + [0.80 / (column_count - 1)] * (column_count - 1)
            table_rows = [
                [Paragraph(inline(cell), styles["cell_head"] if row_index == 0 else styles["cell"]) for cell in row]
                for row_index, row in enumerate(rows)
            ]
            table = Table(
                table_rows,
                colWidths=[content_width * fraction for fraction in fractions],
                repeatRows=1,
                hAlign="LEFT",
            )
            table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#24536d")),
                        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f0f4f6")]),
                        ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#c3cdd3")),
                        ("LEFTPADDING", (0, 0), (-1, -1), 4),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                        ("TOPPADDING", (0, 0), (-1, -1), 4),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                    ]
                )
            )
            story.extend([table, Spacer(1, 8)])
        elif line.startswith("### "):
            story.append(Paragraph(inline(line[4:]), styles["h3"]))
        elif line.startswith("## "):
            story.append(Paragraph(inline(line[3:]), styles["h2"]))
        elif line.startswith("# "):
            story.append(Paragraph(inline(line[2:]), styles["title"]))
        elif re.match(r"^\d+\. ", line):
            story.append(Paragraph(inline(line), styles["bullet"]))
        elif line.startswith("- "):
            story.append(Paragraph("• " + inline(line[2:]), styles["bullet"]))
        else:
            style = styles["caption"] if line.startswith("图 ") else styles["body"]
            if first_paragraph:
                style = ParagraphStyle("subtitle", parent=styles["body"], alignment=TA_CENTER, fontSize=10.5, spaceAfter=14)
                first_paragraph = False
            story.append(Paragraph(inline(line), style))
        index += 1

    output.parent.mkdir(parents=True, exist_ok=True)
    document = SimpleDocTemplate(
        str(output), pagesize=A4, leftMargin=17 * mm, rightMargin=17 * mm,
        topMargin=16 * mm, bottomMargin=16 * mm,
        title="GvpA候选序列的多目标评价与筛选", author="GV02-03项目组",
    )

    def decorate(canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#d1d9de"))
        canvas.line(17 * mm, 12 * mm, page_width - 17 * mm, 12 * mm)
        canvas.setFont("CJK", 7.5)
        canvas.setFillColor(colors.HexColor("#64717a"))
        canvas.drawString(17 * mm, 8 * mm, "GV02-03  多目标候选评价")
        canvas.drawRightString(page_width - 17 * mm, 8 * mm, str(doc.page))
        canvas.restoreState()

    document.build(story, onFirstPage=decorate, onLaterPages=decorate)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the final Markdown report as a CJK PDF")
    parser.add_argument("--input", type=Path, default=ROOT / "reports/M5_最终报告.md")
    parser.add_argument("--output", type=Path, default=ROOT / "reports/M5_最终报告.pdf")
    parser.add_argument("--font", type=Path)
    args = parser.parse_args()
    build(
        args.input.resolve(),
        args.output.resolve(),
        args.font.resolve() if args.font is not None else None,
    )
    print(args.output.resolve())
