from pathlib import Path
import hashlib,json

root=Path('E:/machine').resolve()
groups={
    'references/design':'另一套生成与五维mock评分工程及论文；本项目只使用T05候选，相关独有天然序列已保留在data/raw/design',
    'references/t05/src/raw_data_processing':'混合家族上游处理代码；本项目已有输入审计，正式参考清洗将独立实现',
    'references/t05/src/analysis_transformer_lm':'旧混合家族EDA，与本项目统一评价实现分开，已有基础审计替代',
}
singles={
    'data/candidates/design/candidates.fasta':'用户明确候选入口为T05，删除design已筛选的10条候选子集',
    'references/t05/T05-report.pptx':'与保留的完整来源报告重复，不作为本项目答辩材料',
    'references/t05/doc/ana_transformer_lm.md':'旧混合家族EDA说明，相关数据特征已纳入本项目报告',
    'references/t05/README_original.md':'旧完整工程说明含已删除模块及失效路径，改用当前来源说明',
}
rows=[]
for folder,reason in groups.items():
    for p in sorted((root/folder).rglob('*')):
        if p.is_file():
            rows.append({'path':p.relative_to(root).as_posix(),'reason':reason})
for path,reason in singles.items():
    rows.append({'path':path,'reason':reason})
for row in rows:
    p=(root/row['path']).resolve()
    assert root in p.parents and p.is_file()
    row.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),action='delete')
snapshot=[]
for folder in ['data','references']:
    for p in sorted((root/folder).rglob('*')):
        if p.is_file() and p.relative_to(root).as_posix() not in {r['path'] for r in rows}:
            snapshot.append({'path':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
out=root/'docs/cleanup/round2'
out.mkdir(parents=True,exist_ok=True)
(out/'deletions.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(out/'retained_before.json').write_text(json.dumps(snapshot,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'delete_files':len(rows),'delete_bytes':sum(r['bytes'] for r in rows)},indent=2))
