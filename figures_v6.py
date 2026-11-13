"""V6 hydrogeochemistry diagrams. All compositions derive from meq/L ionic percentages.
The Durov square is intentionally identified as descriptive 3x3 composition fields.
"""
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Rectangle
from figures import base,style_axis,title,annotate,legend,palette,SYMBOLS,H,txy,triangle,_weight


def durov_coordinates(row):
    """Return central square coordinates: alkali cation fraction; strong-acid anion fraction."""
    x = float(row['NaK_pct'])/100.
    y = float(row['SO4NO3_pct']+row['Cl_pct'])/100.
    return x,y


def durov(df,s, color_zones=True, labels=False, grid=True):
    fig=base(s,(max(s.width,9),max(s.height,7.5)))
    ax=fig.add_subplot(111)
    ax.set_aspect('equal')
    # Central composition square (x alkalis 0-100%, y strong-acid anions 0-100%).
    if color_zones:
        colors=['#e4f2f0','#cce6df','#b5d8d0','#f9e9ca','#f4dca8','#edcd88','#e6e2f4','#d5cff0','#beb4e4']
        for row in range(3):
            for col in range(3):
                ax.add_patch(Rectangle((col/3,row/3),1/3,1/3,facecolor=colors[row*3+col],alpha=s.field_alpha,edgecolor='#535353',lw=.8,zorder=0))
                if labels:
                    ax.text((col+.5)/3,(row+.5)/3,str(row*3+col+1),fontsize=s.annotation_size,ha='center',va='center',weight='bold',alpha=.55)
    ax.add_patch(Rectangle((0,0),1,1,fill=False,edgecolor='black',lw=s.axis_width,zorder=5))
    if grid:
        for q in [.25,.5,.75]:
            ax.plot([q,q],[0,1],c='#777',lw=s.grid_width,zorder=1)
            ax.plot([0,1],[q,q],c='#777',lw=s.grid_width,zorder=1)
    # Top cation ternary and left rotated anion ternary; explicitly linked by x/y projections.
    top_offset=1.15
    # Dedicated top triangle perimeter and compositional grid
    top=np.array([[0,top_offset],[1,top_offset],[.5,top_offset+H]])
    ax.add_patch(Polygon(top,closed=True,fill=False,edgecolor='black',lw=s.axis_width,zorder=5))
    if grid:
        for q in [.2,.4,.6,.8]:
            for u,v in [((q,0,1-q),(q,1-q,0)),((0,q,1-q),(1-q,q,0)),((0,1-q,q),(1-q,0,q))]:
                x1,y1=txy(*u); x2,y2=txy(*v)
                ax.plot([x1,x2],[y1+top_offset,y2+top_offset],color='#b8c6cf',lw=s.grid_width,zorder=1)
    left_offset=-.15
    # map x,y from standard ternary to x'=left_offset-y, y'=x
    coords=np.array([(left_offset,0),(left_offset,1),(left_offset-H,.5)])
    ax.add_patch(Polygon(coords,closed=True,fill=False,edgecolor='black',lw=s.axis_width,zorder=5))
    if grid:
        for q in [.2,.4,.6,.8]:
            for u,v in [((q,0,1-q),(q,1-q,0)),((0,q,1-q),(1-q,q,0)),((0,1-q,q),(1-q,0,q))]:
                x1,y1=txy(*u); x2,y2=txy(*v)
                ax.plot([left_offset-y1,left_offset-y2],[x1,x2],color='#b8c6cf',lw=s.grid_width,zorder=1)
    cs=palette(s,len(df))
    for i,(_,r) in enumerate(df.iterrows()):
        if not np.isfinite([r['Ca_pct'],r['Mg_pct'],r['NaK_pct'],r['HCO3CO3_pct'],r['SO4NO3_pct'],r['Cl_pct']]).all():continue
        x,y=durov_coordinates(r)
        cx,cy=txy(r['Ca_pct'],r['Mg_pct'],r['NaK_pct'])
        ax.scatter(cx,cy+top_offset,c=[cs[i]],s=s.marker_size,marker=SYMBOLS[i%len(SYMBOLS)],edgecolor='black',zorder=8)
        ax.scatter(left_offset-txy(r['HCO3CO3_pct'],r['SO4NO3_pct'],r['Cl_pct'])[1],
                   txy(r['HCO3CO3_pct'],r['SO4NO3_pct'],r['Cl_pct'])[0],c=[cs[i]],s=s.marker_size,marker=SYMBOLS[i%len(SYMBOLS)],edgecolor='black',zorder=8)
        ax.scatter(x,y,c=[cs[i]],s=s.marker_size*1.15,marker=SYMBOLS[i%len(SYMBOLS)],edgecolor='black',zorder=9)
        if s.show_labels:annotate(ax,x,y,r.Sample_ID,s,box=True)
    props=dict(fontsize=s.fontsize,weight=_weight(s),family=s.family,color='black')
    ax.text(.5,top_offset+H+.10,'Mg',ha='center',**props)
    ax.text(.08,top_offset-.07,'Ca',ha='center',**props)
    ax.text(1,top_offset-.07,'Na+K',ha='center',**props)
    ax.text(left_offset-H-.08,.5,'SO₄+NO₃',rotation=90,va='center',ha='right',**props)
    ax.text(left_offset-.10,1.08,'Cl',ha='center',**props)
    ax.text(left_offset+.02,-.10,'HCO₃+CO₃',ha='center',**props)
    ax.text(.5,-.12,'Na+K fraction (meq%)',ha='center',**props)
    ax.text(1.11,.5,'Strong-acid anions (meq%)',rotation=90,va='center',**props)
    ax.set_xlim(-1.2,1.32);ax.set_ylim(-.28,2.2);ax.axis('off')
    title(ax,'Durov composition diagram',s)
    legend(ax,df.Sample_ID.tolist(),s)
    return fig


def chadha_xy(row):
    x=float(row['Ca_pct']+row['Mg_pct']-row['NaK_pct'])
    y=float(row['HCO3CO3_pct']-row['Cl_pct']-row['SO4NO3_pct'])
    return x,y


def chadha(df,s,colored=True):
    fig=base(s);ax=fig.add_subplot(111)
    if colored:
        zones=[((-100,-100),'#f7d6d2'),((0,-100),'#e3dcf5'),((-100,0),'#d8e8f4'),((0,0),'#d8eddc')]
        for (x,y),c in zones:ax.add_patch(Rectangle((x,y),100,100,facecolor=c,edgecolor='none',alpha=s.field_alpha,zorder=0))
    ax.axhline(0,c='black',lw=s.axis_width);ax.axvline(0,c='black',lw=s.axis_width)
    texts=[(-50,85,'Na–HCO₃'),(50,85,'Ca/Mg–HCO₃'),(-50,-85,'Na–Cl/SO₄'),(50,-85,'Ca/Mg–Cl/SO₄')]
    for x,y,t in texts:ax.text(x,y,t,ha='center',va='center',fontsize=max(s.annotation_size,8),weight=_weight(s),color='#333333')
    cs=palette(s,len(df))
    for i,(_,r) in enumerate(df.iterrows()):
        x,y=chadha_xy(r)
        ax.scatter(x,y,c=[cs[i]],s=s.marker_size*1.2,marker=SYMBOLS[i%len(SYMBOLS)],edgecolors='black',zorder=5)
        if s.show_labels:annotate(ax,x,y,r.Sample_ID,s,box=True)
    ax.set_xlim(-100,100);ax.set_ylim(-100,100);ax.set_aspect('equal')
    style_axis(ax,s,'(Ca + Mg) − (Na + K) (meq%)','(HCO₃ + CO₃) − (Cl + SO₄ + NO₃) (meq%)',grid=False)
    title(ax,'Chadha hydrochemical classification',s)
    legend(ax,df.Sample_ID.tolist(),s)
    return fig


def stiff_pairs(row):
    return np.array([row['Ca_meqL'],row['Mg_meqL'],row['Na_meqL']+row['K_meqL']],float),np.array([
        row['HCO3_meqL']+row['CO3_meqL'],row['SO4_meqL']+row['NO3_meqL'],row['Cl_meqL']],float)


def stiff(df,s,selected=None,overlay=False):
    target=df if overlay else df[df.Sample_ID.astype(str)==str(selected or df.Sample_ID.iloc[0])]
    fig=base(s);ax=fig.add_subplot(111)
    maximum=0
    for _,r in df.iterrows():
        a,b=stiff_pairs(r);maximum=max(maximum,np.max(a),np.max(b))
    maximum=max(1,maximum*1.15)
    cs=palette(s,len(target))
    for i,(_,r) in enumerate(target.iterrows()):
        left,right=stiff_pairs(r)
        xs=np.r_[-left, right[::-1],-left[0]]
        ys=np.r_[[2,1,0],[0,1,2],[2]]
        ax.fill(xs,ys,color=cs[i],alpha=.23,edgecolor=cs[i],lw=s.plot_width,label=r.Sample_ID,zorder=4)
        ax.plot(xs,ys,color=cs[i],lw=s.plot_width,zorder=5)
        ax.scatter(np.r_[-left,right],np.r_[[2,1,0],[2,1,0]],c=[cs[i]],s=s.marker_size*.25,zorder=6)
    ax.axvline(0,color='black',lw=s.axis_width)
    ax.set_xlim(-maximum,maximum);ax.set_ylim(-.55,2.55)
    ax.set_yticks([2,1,0],['Ca          HCO₃+CO₃','Mg          SO₄+NO₃','Na+K             Cl'])
    # y text is replaced by side labels, so their meaning is unambiguous
    ax.set_yticks([2,1,0],['','',''])
    for y,c,a in zip([2,1,0],['Ca','Mg','Na+K'],['HCO₃+CO₃','SO₄+NO₃','Cl']):
        ax.text(-maximum*1.04,y,c,ha='right',va='center',fontsize=s.fontsize,weight=_weight(s),clip_on=False)
        ax.text(maximum*1.04,y,a,ha='left',va='center',fontsize=s.fontsize,weight=_weight(s),clip_on=False)
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{abs(v):g}'))
    style_axis(ax,s,'Concentration (meq/L) — cations left, anions right',grid=False)
    title(ax,'Stiff composition polygon',s)
    if overlay and s.show_legend:ax.legend(loc='upper right',fontsize=s.legend_size,frameon=True)
    return fig
