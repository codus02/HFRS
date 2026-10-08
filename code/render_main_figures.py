from pathlib import Path
import ast,json,sys,re,textwrap,hashlib
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Rectangle
from matplotlib.ticker import NullLocator

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
KOREAN=ROOT/'reference_results'
OUT=ROOT/'results/reproduced_figures';OUT.mkdir(parents=True,exist_ok=True)
PREVIEW=ROOT/'results/figure_previews';PREVIEW.mkdir(parents=True,exist_ok=True)
from map_panel import draw_map_panel
layout=json.loads((ROOT/'data/figure_layout.json').read_text(encoding='utf-8'))
constants=layout['constants'];terms=layout['terms']
GROUPS=[(terms[g],ids) for g,ids in constants['GROUPS']]
ORDER=[v for _,ids in GROUPS for v in ids]
NAMES={v:terms[name] for v,name in constants['NAMES'].items()}
UNITS={v:terms.get(unit,unit).lower() for v,unit in constants['UNITS'].items()}
UNITS['V15']='% of all dwellings'
PERIODS=[(p,a,b,f'All ({a}–{b})' if p=='ALL' else f'{a}–{b}') for p,a,b,_ in constants['PERIODS']]
plt.rcParams.update({'font.family':'Times New Roman','font.size':9,'axes.unicode_minus':False,'pdf.fonttype':42})
coef=pd.read_csv(KOREAN/'irr_values.csv').set_index(['period','variable_id'])
checked=pd.read_csv(KOREAN/'irr_values.csv').set_index(['period','variable_id'])
for index,row in checked.iterrows():
    for c in ['IRR','lo','hi']:assert f'{row[c]:.3f}'==f'{coef.loc[index,c]:.3f}',(index,c)
values=pd.read_csv(KOREAN/'sir_values.csv')
geo=gpd.read_file(ROOT/'data/district_boundaries.gpkg',layer='districts')[['sigungu_cd','geometry']]
geo.sigungu_cd=pd.to_numeric(geo.sigungu_cd).astype(int)
if geo.crs is None:geo=geo.set_crs(5179)
if geo.sigungu_cd.duplicated().any():geo=geo.dissolve(by='sigungu_cd',as_index=False)
geo=geo.to_crs(5179)
omitted=geo[~geo.sigungu_cd.isin(values.sigungu_cd)].copy()
geo=geo[geo.sigungu_cd.isin(values.sigungu_cd)].copy()
assert len(geo)==228 and set(omitted.sigungu_cd)=={37430}
reports=[]

def save(fig,name):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();clipped=[];labels=[]
    for t in fig.findobj(matplotlib.text.Text):
        label=t.get_text()
        if not label:continue
        assert not re.search('[가-힣]',label),label
        labels.append(label)
    for t in fig.texts:
        b=t.get_window_extent(renderer)
        if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:clipped.append(t.get_text())
    assert not clipped,clipped
    fig.savefig(OUT/(name+'.png'),dpi=600,facecolor='white')
    fig.savefig(OUT/(name+'.pdf'),facecolor='white')
    fig.savefig(PREVIEW/(name+'.png'),dpi=140,facecolor='white')
    reports.append({'figure':name,'labels':labels,'clipped_text':clipped,'dpi':600})
    plt.close(fig)

maps={p:values[values.period.eq(p)] for p,*_ in PERIODS}
fig,_=draw_map_panel(PERIODS,geo,omitted,maps,constants['WARM'])
for t in fig.texts:
    if t.get_text().startswith('전체'):t.set_text('All (2011–2024)')
save(fig,'Figure_1_SIR')

def label(fig,x,y,s,size=9,bold=False,ha='left',color='#202020'):
    return fig.text(x,y,s,fontsize=size,fontweight='bold' if bold else 'normal',ha=ha,va='center',color=color,linespacing=1.08)
def rule(fig,y,lw=.5,color='#A0A0A0'):
    fig.add_artist(plt.Line2D([.012,.990],[y,y],transform=fig.transFigure,color=color,lw=lw))
def row_layout():
    rows=[];y=.891
    for group,ids in GROUPS:
        for j,v in enumerate(ids):
            rows.append((v,y,group if j==0 else ''))
            y-=.0375
        y-=.006
    return rows,y
def covariate_label(fig,v,y):
    name=NAMES[v]
                                                                  
    if v=='V11':name='Fiscal autonomy'
    if v=='V08':name='National Pension individually insured\npersons reporting income'
    wrapped=name if '\n' in name else textwrap.fill(name,width=35,break_long_words=False,break_on_hyphens=False)
    n=wrapped.count('\n')+1
    label(fig,.183,y+.006,wrapped,size=8.9)
    label(fig,.183,y-(.012 if n==1 else .020),'('+UNITS[v]+')',size=7.7,color='#505050')

rows,bottom=row_layout()
fig=plt.figure(figsize=(7.5,9.1),facecolor='white')
xs=[.653,.786,.919]
rule(fig,.983,.8,'black')
label(fig,.012,.958,'Domain',9.5,True)
label(fig,.183,.958,'Covariate (original unit)',9.5,True)
for x,(_,_,_,period) in zip(xs,PERIODS[1:]):label(fig,x,.958,period,9.4,True,'center')
rule(fig,.936,.8,'black')
matrix_cells=[]
for v,y,group in rows:
    if group:
        if v!=ORDER[0]:rule(fig,y+.023,.35,'#CCCCCC')
        label(fig,.012,y,textwrap.fill(group,14,break_long_words=False),8.9,True)
    covariate_label(fig,v,y)
    for x,(p,*_) in zip(xs,PERIODS[1:]):
        r=coef.loc[(p,v)];flag='higher' if r.lo>1 else 'lower' if r.hi<1 else 'uncertain'
        fill={'higher':'#F0DAD8','lower':'#D8E5EF','uncertain':'#FFFFFF'}[flag]
        fig.add_artist(Rectangle((x-.061,y-.018),.122,.036,transform=fig.transFigure,facecolor=fill,edgecolor='#D8D8D8',linewidth=.35))
        cell=f'{r.IRR:.3f}'+('*' if flag!='uncertain' else '')
        label(fig,x,y,cell,10,ha='center')
        matrix_cells.append({'period':p,'variable_id':v,'IRR':r.IRR,'lo':r.lo,'hi':r.hi,'classification':flag,'cell_text':cell})
rule(fig,bottom+.014,.8,'black')
fig.legend(handles=[Patch(facecolor='#F0DAD8',edgecolor='#D8D8D8',label='Higher incidence*'),
                    Patch(facecolor='#D8E5EF',edgecolor='#D8D8D8',label='Lower incidence*'),
                    Patch(facecolor='white',edgecolor='#D8D8D8',label='95% CrI includes 1')],
           loc='lower left',bbox_to_anchor=(.012,.055),ncol=1,frameon=False,fontsize=8.8,labelspacing=.35)
label(fig,.012,.028,'* 95% CrI excludes 1. IRRs are per 1 SD of the selected transformed covariate. Full intervals: Table S3.',8.0)
assert len(matrix_cells)==57 and sum(r['classification']!='uncertain' for r in matrix_cells)==17
pd.DataFrame(matrix_cells).to_csv(ROOT/'results/period_cells_en.csv',index=False)
save(fig,'Figure_2_Period_IRR')

sp=pd.read_csv(KOREAN/'supplement_spatial_values.csv')
assert len(sp)==684 and sp.region_id.between(1,228).all()
bounds=np.r_[geo.total_bounds.reshape(2,2),omitted.total_bounds.reshape(2,2)]
xlim=(bounds[:,0].min()-20000,bounds[:,0].max()+20000)
ylim=(bounds[:,1].min()-22000,bounds[:,1].max()+22000)
colors=['#567B9A','#A8C1D4','#F1EFE9','#D9A7A1','#AD655E']
fig=plt.figure(figsize=(7.5,4.2),facecolor='white')
for i,(p,*_,period) in enumerate(PERIODS[1:]):
    ax=fig.add_axes([.005+i*.282,.03,.278,.90])
    m=geo.merge(sp[sp.period.eq(p)],on='sigungu_cd',validate='one_to_one')
    assert len(m)==228 and m['mean'].notna().all()
    multiplier=np.exp(m['mean'])
    assert np.allclose(multiplier,m.spatial_multiplier)
    cl=np.digitize(multiplier,[.5,.8,1.25,2],right=True)
    m.plot(ax=ax,color=[colors[v] for v in cl],edgecolor='#929292',linewidth=.2)
    omitted.plot(ax=ax,facecolor='#E5E5E5',edgecolor='#777777',hatch='////',linewidth=.2)
    ax.set_xlim(xlim);ax.set_ylim(ylim);ax.set_axis_off()
    label(fig,.144+i*.282,.955,period,10,ha='center')
fig.legend(handles=[Patch(facecolor=c,edgecolor='#999999',linewidth=.3,label=l) for c,l in zip(colors,['≤0.5','>0.5–0.8','>0.8–1.25','>1.25–2','>2'])],
           loc='center left',bbox_to_anchor=(.845,.50),ncol=1,frameon=False,fontsize=8.5,
           title='Spatial'+chr(10)+'multiplier',title_fontsize=9,alignment='left')
save(fig,'Figure_3_Spatial_effects')

report={'figures':reports,'coefficient_cells_preserved':76,'matrix_cells':57,'colored_matrix_cells':17,
        'sir_district_period_rows':len(values),'spatial_district_period_rows':len(sp),'model_refitted':False,
        'sources':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [KOREAN/'irr_values.csv',KOREAN/'sir_values.csv',KOREAN/'supplement_spatial_values.csv']}}
(ROOT/'results/figure_checks_en.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ['figures','sources']},indent=2))
