"""Build the M1 PDF and data figure from the Markdown report and audited inputs.

Requires reportlab and Poppler. Source FASTA files are never modified.
"""
import argparse
from collections import Counter
from html import escape
import json
import os
from pathlib import Path
import re
import subprocess

from reportlab.graphics import renderPDF
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.shapes import Drawing, String, Rect
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image, Preformatted, KeepTogether


ROOT = Path(__file__).resolve().parents[1]
REPORT_NAME = 'M1_功能需求分析与任务建模报告'


def register_fonts(font_dir):
    pdfmetrics.registerFont(TTFont('CJKBody', str(font_dir / 'simsun.ttc')))
    pdfmetrics.registerFont(TTFont('CJKHead', str(font_dir / 'simhei.ttf')))
    pdfmetrics.registerFontFamily('CJKBody', normal='CJKBody', bold='CJKHead', italic='CJKBody', boldItalic='CJKHead')


def make_length_chart(summary):
    drawing = Drawing(540, 235)
    drawing.add(String(8, 220, '输入序列长度分布', fontName='CJKHead', fontSize=13, fillColor=colors.HexColor('#17354b')))
    drawing.add(Rect(298, 218, 10, 7, fillColor=colors.HexColor('#226b8c'), strokeColor=None))
    drawing.add(String(313, 218, 'T05 候选  n=200', fontName='CJKBody', fontSize=8))
    drawing.add(Rect(411, 218, 10, 7, fillColor=colors.HexColor('#b4bac4'), strokeColor=None))
    drawing.add(String(426, 218, '名义天然并集  n=1252', fontName='CJKBody', fontSize=8))
    series = []
    for key in ['candidate_length_counts', 'nominal_gvpa_union_length_counts']:
        counts = {int(length): count for length, count in summary[key].items()}
        total = sum(counts.values())
        series.append([100 * sum(count for length, count in counts.items() if low <= length < low + 50) / total for low in range(50, 550, 50)])
    chart = VerticalBarChart()
    chart.x, chart.y, chart.width, chart.height = 42, 48, 486, 150
    chart.data = series
    chart.categoryAxis.categoryNames = ['{}-{}'.format(low, low + 49) for low in range(50, 550, 50)]
    chart.categoryAxis.labels.fontName = 'CJKBody'
    chart.categoryAxis.labels.fontSize = 7
    chart.valueAxis.labels.fontName = 'CJKBody'
    chart.valueAxis.labels.fontSize = 8
    chart.valueAxis.valueMin = 0
    chart.valueAxis.valueMax = 100
    chart.valueAxis.valueStep = 20
    chart.valueAxis.visibleGrid = True
    chart.valueAxis.gridStrokeColor = colors.HexColor('#e3e7eb')
    chart.valueAxis.gridStrokeWidth = 0.3
    chart.bars[0].fillColor = colors.HexColor('#226b8c')
    chart.bars[1].fillColor = colors.HexColor('#b4bac4')
    chart.bars.strokeColor = None
    chart.barSpacing = 1
    chart.groupSpacing = 6
    drawing.add(chart)
    drawing.add(String(5, 198, '%', fontName='CJKBody', fontSize=9))
    drawing.add(String(435, 18, '序列长度区间 / aa', fontName='CJKBody', fontSize=9))
    drawing.add(String(43, 4, '两组分别归一化；天然并集尚未完成家族和质量过滤。', fontName='CJKBody', fontSize=8, fillColor=colors.HexColor('#555f6b')))
    return drawing


def inline(text):
    text = escape(text)
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'`([^`]+)`', r'<font color="#284b63">\1</font>', text)
    # Keep citations and full URLs visible; this report uses no fetched prose.
    return text


def build(args):
    register_fonts(args.font_dir)
    report_dir = ROOT / 'reports'
    temporary = ROOT / 'tmp/report_render'
    temporary.mkdir(parents=True, exist_ok=True)
    figures = report_dir / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    summary = json.loads((ROOT / 'results/input_audit/summary.json').read_text(encoding='utf-8'))
    drawing = make_length_chart(summary)
    figure_pdf = temporary / 'input_lengths.pdf'
    renderPDF.drawToFile(drawing, str(figure_pdf))
    subprocess.run([str(args.pdftoppm), '-singlefile', '-r', '160', '-png', str(figure_pdf), str(figures / 'input_lengths')], check=True, capture_output=True)

    page_width, page_height = A4
    content_width = page_width - 36 * mm
    doc = SimpleDocTemplate(str(report_dir / (REPORT_NAME + '.pdf')), pagesize=A4,
        rightMargin=18*mm, leftMargin=18*mm, topMargin=17*mm, bottomMargin=17*mm,
        title='GV02-03 GvpA候选序列功能需求分析与任务建模报告', author='GV02-03 项目组')
    base = dict(fontName='CJKBody', textColor=colors.HexColor('#20252a'), wordWrap='CJK', splitLongWords=True)
    body = ParagraphStyle('Body', fontSize=9.6, leading=15.0, spaceAfter=7, **base)
    title = ParagraphStyle('Title', fontName='CJKHead', fontSize=21, leading=28, spaceAfter=5, textColor=colors.HexColor('#152d40'), wordWrap='CJK')
    heading = ParagraphStyle('Section', fontName='CJKHead', fontSize=14.5, leading=21, spaceBefore=9, spaceAfter=8, keepWithNext=True, textColor=colors.HexColor('#17354b'), wordWrap='CJK')
    subheading = ParagraphStyle('Subsection', fontName='CJKHead', fontSize=10.7, leading=17, spaceBefore=6, spaceAfter=5, keepWithNext=True, textColor=colors.HexColor('#17354b'), wordWrap='CJK')
    cell = ParagraphStyle('Cell', fontSize=8.0, leading=11.8, spaceAfter=0, **base)
    cell_head = ParagraphStyle('CellHead', fontName='CJKHead', fontSize=8.0, leading=11.8, textColor=colors.white, wordWrap='CJK')
    caption = ParagraphStyle('Caption', fontSize=8.3, leading=12.5, spaceAfter=8, **base)
    small = ParagraphStyle('Small', fontSize=8.1, leading=11.5, spaceAfter=4, **base)
    code = ParagraphStyle('Code', fontName='CJKBody', fontSize=8.5, leading=12.4, spaceBefore=4, spaceAfter=8, backColor=colors.HexColor('#f1f4f6'), borderPadding=7)

    lines = (report_dir / (REPORT_NAME + '.md')).read_text(encoding='utf-8').splitlines()
    story = []
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if line == '<!-- pagebreak -->':
            story.append(PageBreak())
        elif line.startswith('```'):
            block = []
            index += 1
            while index < len(lines) and not lines[index].startswith('```'):
                block.append(lines[index])
                index += 1
            story.append(Preformatted('\n'.join(block), code))
        elif line.startswith('!['):
            drawing2 = make_length_chart(summary)
            factor = content_width / drawing2.width
            drawing2.scale(factor, factor)
            drawing2.width *= factor
            drawing2.height *= factor
            story.append(drawing2)
            story.append(Spacer(1, 5))
        elif line.startswith('|'):
            rows = []
            while index < len(lines) and lines[index].strip().startswith('|'):
                pieces = [s.strip() for s in lines[index].strip().strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?', piece) for piece in pieces):
                    rows.append(pieces)
                index += 1
            index -= 1
            ncols = len(rows[0])
            if ncols == 5:
                widths = [0.42, 0.08, 0.11, 0.13, 0.26]
            elif rows[0][0] == '编号':
                widths = [0.10, 0.27, 0.63]
            elif rows[0][0] == '项目':
                widths = [0.16, 0.38, 0.46]
            else:
                widths = [0.26, 0.36, 0.38]
            cells = [[Paragraph(inline(value), cell_head if r == 0 else cell) for value in row] for r, row in enumerate(rows)]
            table = Table(cells, colWidths=[content_width*w for w in widths], repeatRows=1, hAlign='LEFT')
            table.setStyle(TableStyle([
                ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#25465e')),
                ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.HexColor('#f1f4f6'), colors.white]),
                ('VALIGN', (0,0), (-1,-1), 'TOP'),
                ('LEFTPADDING', (0,0), (-1,-1), 6), ('RIGHTPADDING', (0,0), (-1,-1), 6),
                ('TOPPADDING', (0,0), (-1,-1), 5), ('BOTTOMPADDING', (0,0), (-1,-1), 5),
                ('LINEBELOW', (0,-1), (-1,-1), 0.5, colors.HexColor('#c9d1d8')),
            ]))
            story.extend([table, Spacer(1, 9)])
        elif line.startswith('### '):
            story.append(Paragraph(inline(line[4:]), subheading))
        elif line.startswith('## '):
            story.append(Paragraph(inline(line[3:]), heading))
        elif line.startswith('# '):
            story.append(Paragraph(inline(line[2:]), title))
        else:
            chosen = caption if line.startswith('图1') else small if line.startswith(('[', '课程：', '资料核验：')) else body
            if line.startswith('- '):
                line = '• ' + line[2:]
            story.append(Paragraph(inline(line), chosen))
        index += 1

    def decorate(canvas, document):
        canvas.saveState()
        canvas.setFont('CJKBody', 7.5)
        canvas.setFillColor(colors.HexColor('#64717c'))
        canvas.drawString(18*mm, 9*mm, 'GV02-03  |  M1 问题分析')
        canvas.drawRightString(page_width-18*mm, 9*mm, str(document.page))
        canvas.restoreState()

    doc.build(story, onFirstPage=decorate, onLaterPages=decorate)
    print(report_dir / (REPORT_NAME + '.pdf'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--font-dir', type=Path, default=Path('C:/Windows/Fonts'))
    parser.add_argument('--pdftoppm', type=Path, required=True)
    build(parser.parse_args())
