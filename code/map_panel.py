import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def draw_map_panel(periods, geo, excluded, maps, warm):
    fig=plt.figure(figsize=(6.8,8.2),facecolor='white')
    bounds=np.array([min(geo.total_bounds[0],excluded.total_bounds[0]),
                     min(geo.total_bounds[1],excluded.total_bounds[1]),
                     max(geo.total_bounds[2],excluded.total_bounds[2]),
                     max(geo.total_bounds[3],excluded.total_bounds[3])])
    xlim=(bounds[0]-10000,bounds[2]+12000)
    ylim=(bounds[1]-34000,bounds[3]+18000)
    panel_positions=[(.025,.585,.460,.355),(.515,.585,.460,.355),
                     (.025,.200,.460,.355),(.515,.200,.460,.355)]
    title_y=[.962,.962,.577,.577]
    axes=[]
    for i,((pid,start,stop,label),position) in enumerate(zip(periods,panel_positions)):
        ax=fig.add_axes(position)
        axes.append(ax)
        d=geo.merge(maps[pid],on='sigungu_cd',validate='one_to_one')
        assert len(d)==228 and d.SIR.notna().all()
        cl=np.digitize(d.SIR,[.5,1,2,4],right=True)
        d.plot(ax=ax,color=[warm[c] for c in cl],edgecolor='#929292',linewidth=.22)
                                                                               
        excluded.plot(ax=ax,facecolor='#E5E5E5',edgecolor='#777777',hatch='////',linewidth=.3)
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_axis_off()
        title=f'전체 ({start}–{stop})' if pid=='ALL' else f'{start}–{stop}'
        fig.text(position[0]+position[2]/2,title_y[i],title,ha='center',va='center',
                 fontsize=11,fontfamily=['Arial','Malgun Gothic'],color='#111111')

    ax=axes[0]
    ax.annotate('',xy=(.91,.90),xytext=(.91,.77),xycoords='axes fraction',
                arrowprops={'arrowstyle':'-|>','lw':.9,'color':'#333333'})
    ax.text(.91,.95,'N',transform=ax.transAxes,ha='center',va='center',
            fontsize=8,fontfamily='Arial',color='#333333')
    bx=xlim[0]+22000;by=bounds[1]-12000
    ax.plot([bx,bx+100000],[by,by],color='#333333',lw=.85)
    for x in [bx,bx+50000,bx+100000]:
        ax.plot([x,x],[by-3000,by+3000],color='#333333',lw=.6)
    ax.text(bx+50000,by-9000,'100 km',ha='center',va='top',fontsize=7,fontfamily='Arial')
    labels=['0–0.5','>0.5–1','>1–2','>2–4','>4']
    handles=[Patch(facecolor=c,edgecolor='#8B8B8B',linewidth=.35,label=l)
             for c,l in zip(warm,labels)]
    legend=fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.055,.012),
                     ncol=1,frameon=False,title='SIR',title_fontsize=9,fontsize=8.5,alignment='left',
                     handlelength=1.65,handleheight=.90,handletextpad=.55,
                     columnspacing=1.5,labelspacing=.65,borderaxespad=0)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for artist in [*fig.texts,legend]:
        b=artist.get_window_extent(renderer)
        assert b.x0>=0 and b.y0>=0 and b.x1<=fig.bbox.width and b.y1<=fig.bbox.height
    assert all(np.allclose(a.get_xlim(),axes[0].get_xlim()) and
               np.allclose(a.get_ylim(),axes[0].get_ylim()) for a in axes)
    caption=('그림 1. 전체기간 및 세 기간의 시군구별 연령보정 표준화발생비(SIR). '
             '동일한 행정경계와 2011–2024년 전국 연령군별 발생률을 기준으로 산출하였다. '
             '위쪽 왼쪽부터 행 방향으로 전체기간, 2011–2015년, 2016–2019년, 2020–2024년이다. '
             '네 지도는 동일한 축척과 색상 구간을 사용하며, 빗금은 분석에서 제외된 울릉군이다. '
             '공변량 보정 및 베이지안 평활 전의 관측 SIR이다.')
    return fig,caption
