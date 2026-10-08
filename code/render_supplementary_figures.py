from pathlib import Path
import json, sys
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
OUT = ROOT / 'results/reproduced_figures'; PREVIEW = ROOT / 'results/figure_previews'
KO_OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else None
EXT = ROOT / 'results/extended'
OUT.mkdir(parents=True, exist_ok=True); PREVIEW.mkdir(parents=True, exist_ok=True)
dp = pd.read_csv(EXT / 'district_observed_predicted.csv'); pit = pd.read_csv(EXT / 'pit_cpo.csv')
assert len(dp) == 228 and len(pit) == 3192
geo = gpd.read_file(ROOT / 'data/district_boundaries.gpkg')[['sigungu_cd', 'geometry']]
geo['sigungu_cd'] = pd.to_numeric(geo.sigungu_cd).astype(int)
if geo.crs is None: geo = geo.set_crs(5179)
geo = geo.to_crs(5179)
omitted = geo[~geo.sigungu_cd.isin(dp.sigungu_cd)].copy(); geo = geo[geo.sigungu_cd.isin(dp.sigungu_cd)].copy()
assert len(geo) == 228 and set(omitted.sigungu_cd) == {37430}
m = geo.merge(dp, on='sigungu_cd', validate='one_to_one')
bounds = geo.total_bounds.reshape(2, 2)                                                                                           
xlim = (bounds[:, 0].min() - 12000, bounds[:, 0].max() + 12000); ylim = (bounds[:, 1].min() - 15000, bounds[:, 1].max() + 15000)
warm = ['#FFF7BC', '#FEC44F', '#FC8D59', '#E34A33', '#990000']
div = ['#567B9A', '#A8C1D4', '#F1EFE9', '#D9A7A1', '#AD655E']
sir_labels = ['0–0.5', '>0.5–1', '>1–2', '>2–4', '>4']; res_labels = ['≤0.5', '>0.5–0.8', '>0.8–1.25', '>1.25–2', '>2']
ratio = np.exp(m.log_ratio_obs_pred)
summary = {'spearman_obs_pred_SIR': float(pd.Series(m.observed_SIR).corr(pd.Series(m.predicted_SIR), method='spearman')),
           'districts_obs_pred_ratio_within_0.8_1.25': int(((ratio > .8) & (ratio <= 1.25)).sum()),
           'districts_obs_pred_ratio_gt2_or_lt0.5': int(((ratio > 2) | (ratio <= .5)).sum()),
           'pit_mean': float(pit.pit.mean()), 'pit_sd': float(pit.pit.std()), 'pit_failures': int((pit.cpo_failure > 0).sum()),
           'pit_share_below_0.1': float((pit.pit < .1).mean()), 'pit_share_above_0.9': float((pit.pit > .9).mean()),
           'pit_ks_distance_from_uniform': float(np.max(np.abs(np.sort(pit.pit) - np.arange(1, len(pit) + 1) / len(pit))))}
(EXT / 'posterior_predictive_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')

def maps(lang, dest):
    T = {'en': dict(a='A. Observed SIR', b='B. Posterior-predicted SIR', c='C. Observed / predicted', sir='SIR', res='Ratio', font='Times New Roman'),
         'ko': dict(a='A. 관측 SIR', b='B. 사후예측 SIR', c='C. 관측/예측 비', sir='SIR', res='관측/예측', font=['Malgun Gothic', 'Times New Roman'])}[lang]
    plt.rcParams.update({'font.family': T['font'], 'font.size': 9, 'axes.unicode_minus': False, 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.5, 3.4), facecolor='white')
    panels = [(m.observed_SIR, [.5, 1, 2, 4], warm, sir_labels, T['a'], T['sir']), (m.predicted_SIR, [.5, 1, 2, 4], warm, sir_labels, T['b'], T['sir']),
              (ratio, [.5, .8, 1.25, 2], div, res_labels, T['c'], T['res'])]
    for i, (val, brk, cols, labs, title, ltitle) in enumerate(panels):
        ax = fig.add_axes([.005 + i * .27, .02, .265, .90]); cl = np.digitize(val, brk, right=True)
        m.plot(ax=ax, color=[cols[v] for v in cl], edgecolor='#929292', linewidth=.2)
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_axis_off()
        fig.text(.137 + i * .27, .965, title, ha='center', fontsize=9.5)
        if i in (0, 2):
            fig.legend(handles=[Patch(facecolor=c, edgecolor='#999999', linewidth=.3, label=l) for c, l in zip(cols, labs)], loc='center left',
                       bbox_to_anchor=(.812, .70 if i == 0 else .27), ncol=1, frameon=False, fontsize=8, title=ltitle, title_fontsize=8.5, alignment='left', labelspacing=.3)
    fig.savefig(dest / 'Figure_S1_Observed_predicted.png', dpi=600, facecolor='white'); fig.savefig(dest / 'Figure_S1_Observed_predicted.pdf', facecolor='white')
    if lang == 'en': fig.savefig(PREVIEW / 'Figure_S1_Observed_predicted.png', dpi=140, facecolor='white')
    plt.close(fig)

def pit_hist(lang, dest):
    T = {'en': dict(x='Probability integral transform (PIT) value', y='Density', font='Times New Roman'),
         'ko': dict(x='확률적분변환(PIT) 값', y='밀도', font=['Malgun Gothic', 'Times New Roman'])}[lang]
    plt.rcParams.update({'font.family': T['font'], 'font.size': 9, 'axes.unicode_minus': False, 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(4.2, 3.2), facecolor='white')
    ax.hist(pit.pit, bins=20, range=(0, 1), density=True, color='#A8C1D4', edgecolor='#567B9A', linewidth=.5)
    ax.axhline(1, color='#AD655E', lw=1, ls=(0, (3, 3)))
    ax.set_xlabel(T['x']); ax.set_ylabel(T['y']); ax.set_xlim(0, 1); ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(); fig.savefig(dest / 'Figure_S2_PIT.png', dpi=600, facecolor='white'); fig.savefig(dest / 'Figure_S2_PIT.pdf', facecolor='white')
    if lang == 'en': fig.savefig(PREVIEW / 'Figure_S2_PIT.png', dpi=140, facecolor='white')
    plt.close(fig)

maps('en', OUT); pit_hist('en', OUT)
if KO_OUT: KO_OUT.mkdir(parents=True, exist_ok=True); maps('ko', KO_OUT); pit_hist('ko', KO_OUT)
print(json.dumps(summary, indent=2))
