// Build the project deck with the Codex presentation runtime (@oai/artifact-tool).
// Inputs are published tables; charts and tables remain editable in PowerPoint.
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const runtime = process.env.ARTIFACT_TOOL_MODULE;
const { Presentation, PresentationFile, FileBlob } = await import(runtime ? pathToFileURL(runtime).href : '@oai/artifact-tool');
const skill = process.env.PRESENTATIONS_SKILL_DIR;
if (!skill) throw new Error('Set PRESENTATIONS_SKILL_DIR to the installed presentation skill');
const { applyPresentationChartFont, finalizePresentation } = await import(pathToFileURL(path.join(skill, 'container_tools/artifact_tool_utils.mjs')).href);
const buildParent = path.join(root, 'results/reproductions');
await fs.mkdir(buildParent, {recursive: true});
const build = process.env.PRESENTATION_BUILD_DIR || await fs.mkdtemp(path.join(buildParent, 'presentation-'));
const finalPath = process.env.PRESENTATION_OUTPUT || path.join(root, 'reports/PROJECT_PRESENTATION.pptx');
const reportRoot = path.join(root, 'reports');
const dataRoot = path.join(root, 'results/gv02_03_v2/d_evaluation_member_a_v1');
await fs.mkdir(build, {recursive: true});
await fs.mkdir(path.dirname(finalPath), {recursive: true});
const sha = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
const manifests = [];
for (const directory of [reportRoot, dataRoot]) {
  const bytes = await fs.readFile(path.join(directory, 'manifest.json'));
  const manifest = JSON.parse(bytes);
  for (const [name, hash] of Object.entries(manifest.output_sha256)) {
    if (sha(await fs.readFile(path.join(directory, name))) !== hash) throw new Error(`Hash mismatch: ${name}`);
  }
  manifests.push(sha(bytes));
}
const stats = JSON.parse(await fs.readFile(path.join(reportRoot, 'statistics.json'), 'utf8'));
async function tableFile(directory, name) {
  const lines = (await fs.readFile(path.join(directory, name), 'utf8')).trim().split(/\r?\n/);
  const headers = lines.shift().split('\t');
  return lines.map(line => Object.fromEntries(line.split('\t').map((v,i)=>[headers[i],v])));
}
const ranked = await tableFile(reportRoot, 'pareto_ranking.tsv');
const weights = await tableFile(dataRoot, 'weight_robustness.tsv');
const random = await tableFile(dataRoot, 'random_comparison.tsv');
const sensitivity = await tableFile(dataRoot, 'threshold_coverage_sensitivity.tsv');
const labels = ['约束满足度', '保守性', '新颖性', '候选独特性'];
const metrics = ['constraint_score','conservation_score','novelty_score','uniqueness_score'];
const strategyNames = ['等权加权','Pareto','分维度轮转'];
const palette = ['#2563A6','#14918B','#D98532'];
const font = 'Microsoft YaHei';
const ink = '#17324B', muted = '#526678';
const presentation = Presentation.create({slideSize:{width:1280,height:720}});
const slides = [], tableOwners = [], chartOwners = [];
function text(slide, value, left, top, width, height, size=27, color=ink, bold=false) {
  const shape = slide.shapes.add({geometry:'textbox',position:{left,top,width,height},fill:'none',line:{fill:'none',width:0}});
  shape.text = value;
  shape.text.style = {typeface:font,fontSize:size,color,bold,autoFit:'none'};
  return shape;
}
function slide(title, note) {
  const s = presentation.slides.add();
  s.background.fill='#FFFFFF'; slides.push(s);
  text(s,title,64,38,1152,80,39,ink,true);
  text(s,String(slides.length).padStart(2,'0'),1180,672,55,30,17,muted);
  s.speakerNotes.textFrame.setText(note + '\n来源：reports/PROJECT_REPORT.md；results/gv02_03_v2/d_evaluation_member_a_v1/。固定批次计算结论，参考暂定，未验证生物学功能。');
  return s;
}
function table(slide, values, x=64,y=160,w=1152,h=330, widths) {
  const t = slide.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,values,
    ...(widths ? {columnWidths:widths} : {})});
  for (let i=0;i<values.length;i++) for(let j=0;j<values[0].length;j++) {
    const cell = t.getCell(i,j);
    cell.fill=i===0?'#17324B':i%2?'#F1F5F8':'#FFFFFF';
    cell.text.style={typeface:font,fontSize:22,color:i===0?'#FFFFFF':ink,bold:i===0,autoFit:'none'};
  }
  tableOwners.push(slides.length);
  return t;
}
function chart(slide,type,config) {
  const c=slide.charts.add(type,{hasLegend:true,legend:{position:'bottom'},...config});
  applyPresentationChartFont(c,{fontFamily:font});
  for (const axis of [c.xAxis,c.yAxis]) {
    axis.textStyle.fontSize=18;
    if (axis.title?.textStyle) axis.title.textStyle.fontSize=19;
  }
  c.legend.textStyle.fontSize=17;
  c.dataLabels.textStyle.fontSize=17;
  applyPresentationChartFont(c,{fontFamily:font});
  chartOwners.push(slides.length);
  return c;
}
const f4 = n => Number(n).toFixed(4);
const pc = n => `${(Number(n)*100).toFixed(1)}%`;

let s=slide('GvpA 候选序列的多目标评价与筛选','研究目标是构建评价器。课程 GV02-03 不要求训练新模型，序列 VAE 是生成候选的扩展。');
text(s,'四维证据、固定预算筛选与参数敏感性分析',64,150,1130,60,31,muted);
text(s,'1,000',64,282,350,105,78,palette[0],true);
text(s,'生成候选',70,398,350,45,27);
text(s,'114',465,282,330,105,78,palette[1],true);
text(s,'证据充分的排名候选',465,398,390,45,27);
text(s,'3',942,282,210,105,78,palette[2],true);
text(s,'筛选策略',947,398,220,45,27);
text(s,'GV02-03    发布批次 seed 42\n完整报告、候选名单和复现入口随仓库发布',64,544,1120,96,25,muted);

s=slide('研究问题与评价设计','RQ1–RQ3 对应原始任务第 20–22 页。三种策略必须在相同候选池、相同预算下比较。');
table(s,[['研究问题','实验与判断依据'],['RQ1 目标之间如何冲突','Pearson / Spearman 与四维散点矩阵'],['RQ2 不同策略如何改变名单','加权、Pareto、轮转，同池比较质量与重叠'],['RQ3 名单对参数是否敏感','逐维消融、权重扰动与域/配对覆盖门槛'],['选做 同池随机筛选对照','从 114 条中抽取 20 条，重复 1,000 次']],64,148,1152,350,[365,787]);
text(s,'主分析预算为 Top-20，同时发布 Top-10 和 Top-50。',64,552,1130,60,28);

s=slide('候选池与完整审计','输入 1,000；基础 QC 1,000；家族支持 147，歧义 47，排除 806；联合门槛合格 136；最终排名 114。参考按簇隔离。');
table(s,[['筛选阶段','数量','保留依据'],['生成 / 基础 QC','1,000','格式、长度、组成与生成元数据'],['GvpA 家族支持','147','竞争性家族证据，独立于域命中'],['QC + 家族 + 域门槛','136','域分数 ≥25 bits，模型覆盖 ≥0.95'],['四维证据充分','114','天然参考已解析，可靠配对比例 ≥50%']],64,144,1152,350,[380,145,627]);
text(s,'346 条保守参考，按同源簇拆为 277 / 33 / 36。\n47 条歧义候选和所有排除原因均保留在审计表。',64,540,1130,105,27);

s=slide('四个评价维度','使用域经验百分位、参考 MSA 的 23 个位点、天然参考最近全局距离及同池已解析距离均值。');
table(s,[['目标','计算含义','注意事项'],['约束满足度','PF00741 域支持和模型覆盖','域命中不能独立证明 GvpA 身份'],['保守性','23 个参考统计保守位点一致率','候选不参与定义参考位点'],['新颖性','到可靠天然参考的最近全局距离','无可靠参考匹配时保留缺失'],['候选独特性','到联合合格池的已解析距离均值','固定 136 条参照池，记录覆盖率']],64,144,1152,370,[220,440,492]);
text(s,'四维均为 [0, 1]，越大越好。集合多样性在入选集合内另行计算。',64,562,1130,66,26);

s=slide('RQ1：当前候选池的目标冲突','数值来自两个相关表。仅解释筛选后的 114 条候选，不推广到全部生成候选或天然蛋白。');
const pairs=[];
for(let i=0;i<4;i++) for(let j=i+1;j<4;j++) pairs.push([`${labels[i]} / ${labels[j]}`,f4(stats.pearson[metrics[i]][metrics[j]]),f4(stats.spearman[metrics[i]][metrics[j]])]);
table(s,[['维度对','Pearson','Spearman'],...pairs],64,136,1152,390,[682,235,235]);
text(s,'保守性与候选独特性负相关最强。\n新颖性与候选独特性高度相关，两个排序目标在本批次有明显重叠。',64,552,1130,98,26);

s=slide('目标空间中的候选与 Pareto 前沿','散点来自 reports/pareto_ranking.tsv。投影显示两组目标关系，第一前沿仍按完整四维定义，不按二维投影重新判定。完整散点矩阵在报告。');
for (const [x,y,xlabel,ylabel,left] of [['constraint_score','novelty_score','约束满足度','新颖性',64],['novelty_score','uniqueness_score','新颖性','候选独特性',660]]) {
  const groups=[ranked.filter(r=>+r.pareto_rank!==0),ranked.filter(r=>+r.pareto_rank===0)];
  chart(s,'scatter',{position:{left,top:150,width:554,height:410},scatterOptions:{style:'marker'},
    xAxis:{title:xlabel},yAxis:{title:ylabel},series:groups.map((rows,i)=>({name:i?'第一前沿':'其他候选',xValues:rows.map(r=>+r[x]),values:rows.map(r=>+r[y]),fill:i?palette[1]:'#B3BEC8',line:{fill:'none',width:0},marker:{symbol:'circle',size:i?7:5}}))});
}
text(s,'四维第一非支配前沿有 10 条。二维图是投影，前沿归属由全部四个目标共同决定。',64,595,1130,56,24);

s=slide('RQ2：三种策略的 Top-20 质量','数据来自 strategy_summary.tsv，固定同池、固定预算；均值不是功能成功率。');
chart(s,'bar',{position:{left:64,top:145,width:810,height:410},categories:labels,
  series:stats.top20.map((r,i)=>({name:strategyNames[i],values:metrics.map(m=>r[`mean_${m}`]),fill:palette[i]})),
  barOptions:{direction:'column',grouping:'clustered'},yAxis:{title:'平均分'}});
text(s,'等权加权\n约束、保守性均值最高\n\nPareto\n新颖性、独特性最高',916,193,300,340,27);
text(s,'轮转覆盖各维度高分候选；完整九份名单同时保留预警和原始证据。',64,598,1130,58,26);

s=slide('集合多样性与证据覆盖','原 147 条主池的 10,731 对中有 3,377 对未解析。条件均值必须与有效配对比例和全部配对均值界限一起解释，界限不是置信区间。');
table(s,[['Top-20 策略','已解析对均值','配对覆盖率','全部配对均值界限'],...stats.top20.map((r,i)=>[strategyNames[i],f4(r.subset_diversity_observed_mean),pc(r.subset_diversity_resolved_fraction),`${f4(r.subset_diversity_lower_bound)}–${f4(r.subset_diversity_upper_bound)}`])],64,160,1152,280,[310,255,240,347]);
text(s,'Pareto 的条件均值较高，同时配对覆盖率较低。\n三种策略的未知距离界限重叠，不能据此断言完整集合多样性严格更高。',64,500,1130,110,29);

s=slide('同一候选池的第一 Pareto 前沿','第一前沿共 10 条，统计每个策略名单从全池第一前沿纳入多少条；不将各子集自己的前沿混为全池前沿。');
table(s,[['策略','Top-10 前沿成员','Top-20 前沿成员','Top-50 前沿成员'],...['weighted_sum','pareto','dimension_round_robin'].map((name,i)=>[strategyNames[i],...[10,20,50].map(k=>String(stats.front_by_strategy.find(r=>r.strategy===name&&r.budget===k).selected_from_pool_front))])],64,166,1152,275,[280,290,290,292]);
text(s,`114 条候选形成 ${stats.pareto_layer_count} 个非支配层。\n前沿代表目标之间的计算取舍，不代表已验证功能最优。`,64,513,1130,110,29);

s=slide('消融：约束分对入选名单影响最大','消融仅移除排序目标，硬门槛仍保持。保守性影响较小不能推广为生物学冗余。');
const ablation=stats.ablation.filter(r=>r.removed_metric);
chart(s,'bar',{position:{left:64,top:145,width:735,height:410},categories:ablation.map(r=>`移除${labels[metrics.indexOf(r.removed_metric)]}`),
  series:[{name:'Top-20 Jaccard',values:ablation.map(r=>r.jaccard_with_full),valuesFormatCode:'0.000',fill:palette[0],
    dataLabelOverrides:ablation.map((r,idx)=>({idx,text:r.jaccard_with_full.toFixed(3),position:'outEnd',showValue:false,textStyle:{fontSize:17,typeface:font}}))}],hasLegend:false,
  barOptions:{direction:'bar',grouping:'clustered'},dataLabels:{showValue:true,position:'outEnd',numberFormatCode:'0.000'},xAxis:{title:'与完整四维名单的 Jaccard'}});
text(s,'约束分移除后\n名单 Jaccard 仅 0.3333\n\n保守性在本批次接近饱和\n移除后的名单变化最小',856,190,360,350,27);
text(s,'新颖性和候选独特性高度相关，但移除任一维仍会改变名单。',64,599,1130,56,26);

s=slide('RQ3：覆盖门槛比局部域阈值更敏感','200 次独立相对权重扰动 ±20%，重新归一化；门槛共 27 组，每组会重建合格参照池。');
const bins=[[.7,.8],[.8,.9],[.9,1],[1,1.00001]];
const weightChart=chart(s,'bar',{position:{left:64,top:155,width:555,height:355},categories:['0.7–0.8','0.8–0.9','0.9–1','1.0'],
  series:[{name:'扰动次数',values:bins.map(([a,b])=>weights.filter(r=>+r.top_k_jaccard>=a&&+r.top_k_jaccard<b).length),fill:palette[0]}],
  hasLegend:false,barOptions:{direction:'column'},xAxis:{title:'Top-20 Jaccard'},yAxis:{title:'次数'}});
weightChart.xAxis.textStyle.fontSize=16;
chart(s,'line',{position:{left:670,top:155,width:546,height:355},categories:['25%','50%','75%'],lineOptions:{smooth:false},
  series:[-.1,0,.1].map((off,i)=>({name:`域覆盖 ${Math.min(.95+off,1).toFixed(2)}`,values:[.25,.5,.75].map(f=>+sensitivity.find(r=>+r.score_factor===1&&+r.coverage_offset===off&&+r.minimum_resolved_fraction===f).jaccard),line:{fill:palette[i],width:3},marker:{symbol:'circle',size:6}})),
  xAxis:{title:'最低可靠配对比例'},yAxis:{title:'Top-20 Jaccard'}});
text(s,`权重扰动 Jaccard 均值 ${f4(stats.weight_jaccard.mean)}\n门槛扰动范围 0.6000–1.0000`,64,551,540,92,26);
text(s,'50% 配对门槛下域阈值扰动不改变名单。\n25% / 75% 门槛则明显改变选择。',670,551,546,92,26);

s=slide('随机对照：优势取决于目标','同一 114 条候选池无放回抽取 20 条，重复 1,000 次。仅比较固定批次计算代理；经验比例不是生物学功能检验。');
const randomMetrics=['mean_constraint_score','mean_novelty_score','subset_diversity_lower_bound','subset_diversity_resolved_fraction'];
const randomLabels=['约束满足度','新颖性','多样性下界','配对覆盖率'];
table(s,[['指标','随机均值',...strategyNames],...randomMetrics.map((m,i)=>[randomLabels[i],f4(random.find(r=>r.metric===m).random_mean),...['weighted_sum','pareto','dimension_round_robin'].map(n=>f4(random.find(r=>r.metric===m&&r.strategy===n).observed))])],64,148,1152,335,[270,220,220,220,222]);
text(s,'三种策略提高新颖性，但 Pareto 与轮转的约束均值低于随机。\n多目标筛选没有在所有指标上同时获胜。',64,538,1130,110,29);

s=slide('适用场景与复现边界','评价及筛选链路已复现。原训练 checkpoint 未交付，不能声称原训练与原采样已重放。家族/Pfam/TM 本次校验证据而未完整重跑外部工具。');
table(s,[['使用场景','策略与检查重点'],['重视约束，验证预算有限','加权名单，预先声明权重并审阅预警'],['希望探索目标之间的取舍','Pareto 前沿，检查约束下降与未知距离'],['希望覆盖各维度高分候选','分维度轮转，固定总预算并检查均衡性']],64,145,1152,275,[415,737]);
text(s,'评价复现入口：experiments/run_full_experiment.py\n完整 QC / 比对复现：experiments/reproduce_project.py',64,463,1152,87,24);
text(s,'参考仍为暂定版；原 VAE 权重缺失。候选尚无表达、结构、组装或功能确认。',64,600,1152,58,24,muted);

const draft=path.join(build,'candidate.pptx');
await (await PresentationFile.exportPptx(presentation)).save(draft);
const result=await finalizePresentation({
  workspaceDir:root,candidatePath:draft,finalPath,
  pythonExecutable:process.env.PRESENTATION_PYTHON,
  integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',
    ...[...new Set(tableOwners)].flatMap(n=>['--require-native-table-slide',String(n)])],
  requiredNativeTableOwnerSlides:[...new Set(tableOwners)],requiredNativeChartOwnerSlides:[...new Set(chartOwners)],
  materializeLiteralChartWorkbooks:true,
  fontPolicy:{basis:'design',families:[font,'Calibri']},verifyArtifactToolImport:true,
  receiptPath:path.join(build,'validation.json'),
});
const finalPresentation=await PresentationFile.importPptx(await FileBlob.load(finalPath));
for(let i=0;i<slides.length;i++) {
  const preview=await finalPresentation.export({slide:finalPresentation.slides.items[i],format:'png',scale:1});
  await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await preview.arrayBuffer()));
}
await fs.writeFile(path.join(reportRoot,'presentation_manifest.json'),JSON.stringify({
  source_report_manifest_sha256:manifests[0],source_evaluation_manifest_sha256:manifests[1],
  builder_sha256:sha(await fs.readFile(fileURLToPath(import.meta.url))),slide_count:slides.length,
  output_sha256:{[path.basename(finalPath)]:sha(await fs.readFile(finalPath))},
},null,2)+'\n');
console.log(JSON.stringify({finalPath,slides:slides.length,build,receipt:result.receiptPath},null,2));
