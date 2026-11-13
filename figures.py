"""Matplotlib publication-oriented scientific figures for HydroGeoAI v5.1.
Focused on journal-style control and customer-adjustable aesthetics.
"""
from dataclasses import dataclass, replace, field
from math import sqrt
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, set_link_color_palette
from scipy.stats import zscore
from chemistry import IONS

H = sqrt(3) / 2
COLORS = ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#332288','#AA4499','#117733','#999933']
TABLEAU = ['#4E79A7','#F28E2B','#59A14F','#E15759','#76B7B2','#EDC948','#B07AA1','#FF9DA7','#9C755F','#BAB0AC']
GRAY10 = ['#111111','#303030','#4c4c4c','#666666','#808080','#999999','#b2b2b2','#cccccc','#dddddd','#eeeeee']
SYMBOLS = ['o','s','D','^','v','P','X','h','<','>']
CLUSTER_COLORS = ['#1f77b4','#2ca02c','#d62728','#9467bd','#ff7f0e','#17becf','#8c564b']

@dataclass
class Style:
    family: str='Arial'
    fontsize: int=13
    ticksize: int=11
    legend_size: int=10
    annotation_size: int=9
    bold: bool=True
    axis_width: float=2.0
    plot_width: float=1.7
    grid_width: float=.6
    marker_size: int=55
    show_grid: bool=True
    show_labels: bool=False
    show_title: bool=False
    palette: str='Colorblind'
    dpi: int=300
    width: float=10
    height: float=7
    show_legend: bool=True
    water_type_labels: bool=False
    fill_piper_fields: bool=True
    show_gibbs_domains: bool=True
    pca_mode: str='Biplot'
    color_pca_vectors: bool=True
    hca_colored_clusters: bool=True
    field_alpha: float=.70
    zoom: float=1.0
    sample_label_dx: float=6.0
    sample_label_dy: float=6.0
    sample_label_rotation: float=0.0
    loading_label_dx: float=4.0
    loading_label_dy: float=4.0
    loading_label_rotation: float=0.0
    show_domain_text: bool=True
    legend_location: str='Outside right'
    legend_columns: int=1
    legend_box: bool=False
    label_overrides: dict=field(default_factory=dict)


def with_zoom(s: Style, zoom=1.0):
    return replace(s, width=s.width*zoom, height=s.height*zoom, zoom=zoom)


def _weight(s):
    return 'bold' if s.bold else 'normal'


def palette(s, n):
    if s.palette == 'Grayscale':
        src = GRAY10
    elif s.palette == 'Tableau':
        src = TABLEAU
    else:
        src = COLORS
    return [src[i % len(src)] for i in range(max(n, 1))]


def base(s, size=None):
    fig = plt.figure(figsize=size or (s.width, s.height), dpi=120, facecolor='white', layout='constrained')
    return fig


def style_axis(ax, s, xlabel=None, ylabel=None, grid=None):
    for spine in ax.spines.values():
        spine.set_linewidth(s.axis_width)
        spine.set_color('black')
    ax.tick_params(axis='both', which='major', width=max(1, s.axis_width*.85), length=6,
                   color='black', labelcolor='black', labelsize=s.ticksize)
    ax.tick_params(axis='both', which='minor', width=max(.8, s.axis_width*.55), length=3.5, color='black')
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_fontfamily(s.family)
        t.set_fontweight(_weight(s))
        t.set_color('black')
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=s.fontsize, fontweight=_weight(s), fontfamily=s.family, color='black', labelpad=10)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=s.fontsize, fontweight=_weight(s), fontfamily=s.family, color='black', labelpad=10)
    if s.show_grid if grid is None else grid:
        ax.grid(color='#d5d9dd', linewidth=s.grid_width, alpha=.95, zorder=0)
    else:
        ax.grid(False)


def title(ax, txt, s):
    if s.show_title:
        ax.set_title(txt, fontsize=s.fontsize+2, weight=_weight(s), family=s.family, color='black', pad=14)


def legend_handles(ids, s):
    cs = palette(s, len(ids))
    return [Line2D([0],[0], marker=SYMBOLS[i%len(SYMBOLS)], linestyle='None', markersize=7,
                   markeredgecolor='black', markeredgewidth=.7, markerfacecolor=cs[i], label=str(sample))
            for i, sample in enumerate(ids)]


def legend(ax, ids, s, **kw):
    if not s.show_legend:
        return
    handles = legend_handles(ids, s)
    locations={'Inside upper right':'upper right','Inside upper left':'upper left',
               'Inside lower right':'lower right','Inside lower left':'lower left'}
    inside = s.legend_location in locations
    lg = ax.legend(handles=handles, fontsize=s.legend_size,
                   loc=locations.get(s.legend_location,'upper left'),
                   bbox_to_anchor=None if inside else (1.01,1),
                   frameon=s.legend_box, ncol=s.legend_columns, **kw)
    for t in lg.get_texts():
        t.set_fontfamily(s.family)
        t.set_fontweight(_weight(s))
        t.set_color('black')


def annotate(ax, x, y, txt, s, offset=None, color='black', box=False, rotation=None):
    if offset is None:
        offset = (s.sample_label_dx, s.sample_label_dy)
    bbox = dict(boxstyle='round,pad=0.18', fc='white', ec='none', alpha=0.72) if box else None
    ax.annotate(str(txt), xy=(x, y), xytext=offset, textcoords='offset points', fontsize=s.annotation_size,
                color=color, weight=_weight(s), family=s.family, clip_on=False, bbox=bbox,
                rotation=s.sample_label_rotation if rotation is None else rotation)


def txy(a, b, c, x0=0):
    total = a+b+c
    if total <= 0:
        return (np.nan, np.nan)
    a, b, c = a/total, b/total, c/total
    return x0 + c + .5*b, H*b


def triangle(ax, x0, s, grid=True):
    pts = np.array([(x0,0),(x0+1,0),(x0+.5,H)])
    ax.add_patch(Polygon(pts, closed=True, facecolor='none', edgecolor='black', linewidth=s.axis_width, zorder=5))
    if grid:
        for p in [0.2,0.4,0.6,0.8]:
            for points in [((p,0,1-p),(p,1-p,0)),((0,p,1-p),(1-p,p,0)),((0,1-p,p),(1-p,0,p))]:
                q1 = txy(*points[0], x0=x0); q2 = txy(*points[1], x0=x0)
                ax.plot([q1[0],q2[0]], [q1[1],q2[1]], color='#c8d2d8', lw=s.grid_width, zorder=1)


def diamond_xy(alkalis, strong, center=1.3, offset=.18):
    return center + .5*(strong-alkalis), H+offset + H*(strong+alkalis)


def diamond(ax, s, center=1.3, offset=.18, grid=True):
    coords = np.array([diamond_xy(0,0,center,offset), diamond_xy(0,1,center,offset),
                       diamond_xy(1,1,center,offset), diamond_xy(1,0,center,offset)])
    ax.add_patch(Polygon(coords, closed=True, facecolor='none', edgecolor='black', linewidth=s.axis_width, zorder=5))
    if grid:
        for p in [.2,.4,.6,.8]:
            for start,end in [((p,0),(p,1)),((0,p),(1,p))]:
                u = diamond_xy(*start, center=center, offset=offset); v = diamond_xy(*end, center=center, offset=offset)
                ax.plot([u[0],v[0]], [u[1],v[1]], color='#c8d2d8', lw=s.grid_width, zorder=1)


def _fill_triangle_facies(ax, x0, mode='cation', alpha=.75):
    centroid = (x0+0.5, H/3)
    left = (x0,0); right=(x0+1,0); top=(x0+0.5,H)
    if mode == 'cation':
        colors = ['#f2b5b0','#b8d7a5','#b9b8ef']
        verts = [[left, (x0+0.25,H/2), centroid, (x0+0.5,0)],
                 [(x0+0.25,H/2), top, (x0+0.75,H/2), centroid],
                 [(x0+0.5,0), centroid, (x0+0.75,H/2), right]]
    else:
        colors = ['#a8d9cc','#f4d35e','#b8d4ef']
        verts = [[left, (x0+0.25,H/2), centroid, (x0+0.5,0)],
                 [(x0+0.25,H/2), top, (x0+0.75,H/2), centroid],
                 [(x0+0.5,0), centroid, (x0+0.75,H/2), right]]
    for poly, c in zip(verts, colors):
        ax.add_patch(Polygon(poly, closed=True, facecolor=c, edgecolor='none', alpha=alpha, zorder=.2))


def _fill_diamond_fields(ax, center=1.3, offset=.18, alpha=.66):
    bottom = diamond_xy(0,0,center,offset)
    right = diamond_xy(0,1,center,offset)
    top = diamond_xy(1,1,center,offset)
    left = diamond_xy(1,0,center,offset)
    midl = ((bottom[0]+left[0])/2, (bottom[1]+left[1])/2)
    midr = ((bottom[0]+right[0])/2, (bottom[1]+right[1])/2)
    midtl = ((top[0]+left[0])/2, (top[1]+left[1])/2)
    midtr = ((top[0]+right[0])/2, (top[1]+right[1])/2)
    c = ((left[0]+right[0])/2, (left[1]+right[1])/2)
    polys = [([left, midtl, c, midl], '#d9bde9'), ([midtl, top, midtr, c], '#f4ef86'),
             ([c, midtr, right, midr], '#a8d89b'), ([midl, c, midr, bottom], '#f4b6ad')]
    for poly, col in polys:
        ax.add_patch(Polygon(poly, closed=True, facecolor=col, edgecolor='none', alpha=alpha, zorder=.2))


def piper(df, s):
    fig = base(s, (max(s.width,10), max(s.height,8)))
    ax = fig.add_subplot(111)
    xanion, center = 1.6, 1.3
    if s.fill_piper_fields:
        _fill_triangle_facies(ax, 0, 'cation', s.field_alpha)
        _fill_triangle_facies(ax, xanion, 'anion', s.field_alpha)
        _fill_diamond_fields(ax, center, alpha=s.field_alpha)
    triangle(ax, 0, s, s.show_grid)
    triangle(ax, xanion, s, s.show_grid)
    diamond(ax, s, center=center, grid=s.show_grid)
    cs = palette(s, len(df)); ids = df['Sample_ID'].tolist()
    for i, (_, r) in enumerate(df.iterrows()):
        vals = [r['Ca_pct'],r['Mg_pct'],r['NaK_pct'],r['HCO3CO3_pct'],r['Cl_pct'],r['SO4NO3_pct']]
        if not np.isfinite(vals).all():
            continue
        ca, mg, nak, hco, cl, so = np.array(vals)/100
        cat = txy(ca, mg, nak, 0); ani = txy(hco, so, cl, xanion); dia = diamond_xy(nak, cl+so, center=center)
        ax.scatter(*cat, c=[cs[i]], marker=SYMBOLS[i%len(SYMBOLS)], s=s.marker_size, edgecolors='black', linewidths=.8, zorder=8)
        ax.scatter(*ani, c=[cs[i]], marker=SYMBOLS[i%len(SYMBOLS)], s=s.marker_size, edgecolors='black', linewidths=.8, zorder=8)
        ax.scatter(*dia, c=[cs[i]], marker='D', s=s.marker_size*1.1, edgecolors='black', linewidths=.9, zorder=9)
        if s.show_labels:
            annotate(ax, *dia, r.Sample_ID, s, box=True)
        if s.water_type_labels:
            annotate(ax, *dia, r.Facies, s, offset=(s.sample_label_dx, s.sample_label_dy-14), color='#2f2f2f', box=True)
    labels = [(0,-.17,'Ca'),(1,-.17,'Na+K'),(.5,H+.12,'Mg'),
              (xanion,-.17,'HCO₃+CO₃'),(xanion+1,-.17,'Cl'),(xanion+.5,H+.12,'SO₄+NO₃')]
    for x, y, txt in labels:
        ax.text(x, y, txt, fontsize=s.fontsize+1, weight=_weight(s), family=s.family, ha='center', va='center', color='black')
    # Draw facies field boundaries explicitly over the color fills: three cation and
    # anion dominance sectors; diamond uses 50%-composition divisions.
    if s.fill_piper_fields:
        for x0 in [0,xanion]:
            cx=x0+.5;cy=H/3
            for vx,vy in [(x0+.5,0),(x0+.25,H/2),(x0+.75,H/2)]:
                ax.plot([cx,vx],[cy,vy],color='#585858',lw=max(1,s.plot_width*.72),zorder=3)
        # Two orthogonal mid-field boundaries in diamond (the 50% ionic sectors).
        for first,last in [((.5,0),(.5,1)),((0,.5),(1,.5))]:
            p1=diamond_xy(*first,center=center);p2=diamond_xy(*last,center=center)
            ax.plot([p1[0],p2[0]],[p1[1],p2[1]],c='#585858',lw=max(1,s.plot_width*.72),zorder=3)
    ax.set_xlim(-.28, 3.0); ax.set_ylim(-.28, 3.05); ax.set_aspect('equal'); ax.axis('off')
    title(ax, 'Piper diagram', s)
    if s.show_legend:
        legend(ax, ids, s)
    return fig


def trilinear(df, s, kind):
    fig = base(s, (s.width, s.height)); ax = fig.add_subplot(111)
    is_cat = kind == 'Cations'
    if s.fill_piper_fields:
        _fill_triangle_facies(ax, 0, 'cation' if is_cat else 'anion', s.field_alpha)
    triangle(ax, 0, s, s.show_grid)
    groups = ['Ca_pct','Mg_pct','NaK_pct'] if is_cat else ['HCO3CO3_pct','SO4NO3_pct','Cl_pct']
    labels = ['Ca','Mg','Na+K'] if is_cat else ['HCO₃+CO₃','SO₄+NO₃','Cl']
    cs = palette(s, len(df))
    for i, (_, r) in enumerate(df.iterrows()):
        q = [r[c] for c in groups]
        if not np.isfinite(q).all():
            continue
        x, y = txy(*q)
        ax.scatter(x, y, s=s.marker_size, c=[cs[i]], edgecolors='black', linewidths=.8, zorder=8, marker=SYMBOLS[i%len(SYMBOLS)])
        if s.show_labels:
            annotate(ax, x, y, r.Sample_ID, s, box=True)
    for x, y, txt in [(0,-.14,labels[0]),(.5,H+.12,labels[1]),(1,-.14,labels[2])]:
        ax.text(x, y, txt, fontsize=s.fontsize+1, weight=_weight(s), family=s.family, ha='center', va='center', color='black')
    ax.set_xlim(-.18, 1.34); ax.set_ylim(-.22, H+.2); ax.set_aspect('equal'); ax.axis('off')
    title(ax, f'{kind} trilinear plot', s)
    legend(ax, df['Sample_ID'].tolist(), s)
    return fig


def _draw_gibbs_domains(ax, s):
    # publication-like dominance domains with visible fills
    evap = np.array([[0.0, 2.8e2],[0.18,5.1e2],[0.35,9.0e2],[0.55,2.1e3],[0.78,6.2e3],[0.95,2.0e4],[1.0,4.0e4],[1.0,1e5],[0.0,1e5]])
    precip = np.array([[0.0,1.1e2],[0.15,8e1],[0.35,4.2e1],[0.55,1.9e1],[0.78,8],[0.96,3.2],[1.0,2.0],[1.0,1],[0.0,1]])
    rock = np.array([[0.0,1.05e2],[0.18,1.5e2],[0.38,2.25e2],[0.55,3.3e2],[0.68,2.15e2],[0.86,7.5e1],[1.02,1.4e1],[0.68,1.8e2],[0.55,2.8e2],[0.38,2.0e2],[0.18,1.35e2]])
    ax.add_patch(Polygon(evap, closed=True, facecolor='#ecdcc2', edgecolor='none', alpha=.52*s.field_alpha, zorder=.1))
    ax.add_patch(Polygon(precip, closed=True, facecolor='#d8dde2', edgecolor='none', alpha=.55*s.field_alpha, zorder=.1))
    ax.add_patch(Polygon(rock, closed=True, facecolor='none', edgecolor='#b84fff', linewidth=1.4, linestyle=(0,(6,3)), alpha=.95, zorder=1))
    for line in [evap[:7], precip[:7], rock[:7]]:
        ax.plot(line[:,0], line[:,1], color='#555555' if line is not rock[:7] else '#b84fff', lw=1.2, zorder=1)
    if s.show_domain_text:
        ax.text(0.17, 1.7e2, 'Rock\ndominance', fontsize=max(8,s.annotation_size), weight='bold', ha='center', color='#333333')
        ax.text(0.85, 5.5e3, 'Evaporation\ndominance', fontsize=max(8,s.annotation_size), weight='bold', ha='center', color='#333333')
        ax.text(0.72, 8.5, 'Precipitation\ndominance', fontsize=max(8,s.annotation_size), weight='bold', ha='center', color='#333333')


def gibbs(df, s, kind):
    fig = base(s); ax = fig.add_subplot(111)
    xcol = 'NaCa_ratio' if kind == 'Cation' else 'ClHCO3_ratio'
    xlabel = 'Na / (Na + Ca)' if kind == 'Cation' else 'Cl / (Cl + HCO₃)'
    measured_tds = ('TDS' in df.columns and pd.to_numeric(df['TDS'], errors='coerce').notna().all() and
                    (pd.to_numeric(df['TDS'], errors='coerce') > 0).all())
    ycol = 'TDS' if measured_tds else 'TDS_est_mgL'
    if s.show_gibbs_domains:
        _draw_gibbs_domains(ax, s)
    cs = palette(s, len(df))
    for i, (_, r) in enumerate(df.iterrows()):
        ax.scatter(r[xcol], r[ycol], s=s.marker_size*1.25, c=[cs[i]], edgecolors='black', linewidths=.9, zorder=3,
                   marker=SYMBOLS[i%len(SYMBOLS)])
        if s.show_labels:
            annotate(ax, r[xcol], r[ycol], r.Sample_ID, s, box=True)
    ydata = pd.to_numeric(df[ycol])
    ax.set_xlim(0, 1.2); ax.set_yscale('log'); ax.set_ylim(bottom=max(ydata.min()*.65, 1), top=max(ydata.max()*20.0, 1e4))
    style_axis(ax, s, xlabel, 'TDS (mg/L)' if measured_tds else 'Sum of reported ions (mg/L)')
    title(ax, f'Gibbs plot — {kind.lower()} ratio', s)
    legend(ax, df['Sample_ID'].tolist(), s)
    return fig


def correlation(df, s, annotate_values=True, method='pearson'):
    vars_ = [c for c in IONS if df[c].nunique() > 1]
    if len(vars_) < 2:
        return None
    c = df[vars_].corr(method=method)
    fig = base(s, (max(s.width,8), max(s.height,7))); ax = fig.add_subplot(111)
    im = ax.imshow(c, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal', interpolation='nearest')
    ax.set_xticks(range(len(vars_)), vars_)
    ax.set_yticks(range(len(vars_)), vars_)
    ax.tick_params(which='both', bottom=False, left=False)
    style_axis(ax, s, grid=False)
    if annotate_values:
        for i in range(len(vars_)):
            for j in range(len(vars_)):
                v = c.iloc[i,j]
                ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=s.annotation_size,
                        fontweight=_weight(s), color='white' if abs(v) > .6 else 'black', family=s.family)
    cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.03)
    cb.set_label(f'{method.capitalize()} r', fontweight=_weight(s), fontsize=s.fontsize, color='black')
    cb.ax.tick_params(labelsize=s.ticksize, colors='black', width=s.axis_width)
    cb.outline.set_linewidth(s.axis_width)
    title(ax, f'{method.capitalize()} correlation', s)
    return fig


def pca_fit(df, variables):
    X = df[variables].to_numpy(dtype=float)
    z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
    u, singular, vt = np.linalg.svd(z, full_matrices=False)
    scores = u * singular
    ev = (singular**2) / (len(X)-1)
    explained = ev / ev.sum()
    loadings = vt.T * np.sqrt(ev)
    return scores, explained, loadings, z


def pca(df, s, show_vectors=True):
    variables = [c for c in IONS if df[c].nunique() > 1]
    if len(df) < 3 or len(variables) < 2:
        return None, None
    scores, explained, loadings, _ = pca_fit(df, variables)
    mode = s.pca_mode
    fig = base(s); ax = fig.add_subplot(111)
    cs = palette(s, len(df))
    vec_colors = TABLEAU if s.color_pca_vectors else ['black']*len(variables)

    if mode in ('Biplot','Scores'):
        for i, (_, r) in enumerate(df.iterrows()):
            ax.scatter(scores[i,0], scores[i,1], s=s.marker_size, c=[cs[i]], edgecolors='black', linewidths=.8, zorder=5,
                       marker=SYMBOLS[i%len(SYMBOLS)], label=r.Sample_ID)
            if s.show_labels:
                opts=s.label_overrides.get('sample:'+str(r.Sample_ID),{})
                if opts.get('show',True):
                    annotate(ax,scores[i,0],scores[i,1],r.Sample_ID,s,
                             offset=(opts.get('dx',s.sample_label_dx),opts.get('dy',s.sample_label_dy)),
                             box=True,rotation=opts.get('rotation',s.sample_label_rotation))

    if mode in ('Biplot','Loadings') and show_vectors:
        if mode == 'Loadings':
            vec_scale = 1.0
        else:
            vec_scale = max(np.max(np.abs(scores[:,:2])), .5) * .78
        for j, var in enumerate(variables):
            x, y = loadings[j,0]*vec_scale, loadings[j,1]*vec_scale
            color = vec_colors[j % len(vec_colors)]
            ax.annotate('', xy=(x,y), xytext=(0,0), arrowprops=dict(arrowstyle='-|>', color=color, lw=s.plot_width, mutation_scale=12))
            opts=s.label_overrides.get('variable:'+str(var),{})
            if opts.get('show',True):
                annotate(ax,x,y,var,s,
                         offset=(opts.get('dx',s.loading_label_dx),opts.get('dy',s.loading_label_dy)),
                         color=color,box=False,rotation=opts.get('rotation',s.loading_label_rotation))

    ax.axhline(0, color='black', lw=s.axis_width*.9)
    ax.axvline(0, color='black', lw=s.axis_width*.9)
    if mode == 'Loadings':
        lim = max(1.0, np.nanmax(np.abs(loadings[:,:2]))*1.35)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        style_axis(ax, s, f'PC 1 ({explained[0]*100:.2f}%)', f'PC 2 ({explained[1]*100:.2f}%)', grid=False)
        title(ax, 'PCA loading plot', s)
    else:
        ax.margins(.24)
        style_axis(ax, s, f'PC 1 ({explained[0]*100:.2f}%)', f'PC 2 ({explained[1]*100:.2f}%)', grid=False)
        title(ax, 'PCA biplot' if mode == 'Biplot' else 'PCA scores plot', s)
    if mode != 'Loadings':
        legend(ax, df['Sample_ID'].tolist(), s)
    payload = {
        'scores': pd.DataFrame({'Sample_ID': df.Sample_ID, 'PC1': scores[:,0], 'PC2': scores[:,1]}),
        'loadings': pd.DataFrame(loadings[:,:2], index=variables, columns=['PC1','PC2']),
        'explained': pd.DataFrame({'PC': range(1, len(explained)+1), 'Explained_variance_ratio': explained})
    }
    return fig, payload


def hca(df, s, cut_clusters=3):
    variables = [c for c in IONS if df[c].nunique() > 1]
    if len(df) < 3 or len(variables) < 2:
        return None, None
    X = df[variables].to_numpy(float)
    Z = linkage(zscore(X, axis=0, ddof=1), method='ward')
    fig = base(s); ax = fig.add_subplot(111)
    k = min(max(2, cut_clusters), len(df))
    threshold = None
    if s.hca_colored_clusters and len(df) > 2:
        if k >= len(df):
            threshold = Z[0,2] - 1e-6
        else:
            threshold = Z[-(k-1),2] - 1e-6
        set_link_color_palette(CLUSTER_COLORS)
        dendrogram(Z, labels=df.Sample_ID.tolist(), ax=ax, leaf_rotation=50, leaf_font_size=s.ticksize,
                   color_threshold=threshold, above_threshold_color='#9b9b9b')
    else:
        dendrogram(Z, labels=df.Sample_ID.tolist(), ax=ax, above_threshold_color='black', color_threshold=0,
                   leaf_rotation=50, leaf_font_size=s.ticksize)
        for collection in ax.collections:
            collection.set_color('black')
    for collection in ax.collections:
        collection.set_linewidth(s.plot_width)
    ax.set_ylim(bottom=0)
    style_axis(ax, s, 'Samples', 'Linkage distance', grid=False)
    title(ax, 'Hierarchical cluster analysis', s)
    clusters = fcluster(Z, t=k, criterion='maxclust')
    membership = pd.DataFrame({'Sample_ID': df.Sample_ID, 'Cluster': clusters}).sort_values('Cluster')
    return fig, membership


def schoeller(df, s):
    fig = base(s); ax = fig.add_subplot(111)
    ions = ['Ca','Mg','Na','K','HCO3','Cl','SO4']
    cs = palette(s, len(df))
    for i, (_, r) in enumerate(df.iterrows()):
        ys = np.array([r[f'{c}_meqL'] for c in ions], dtype=float)
        ys[ys <= 0] = np.nan
        ax.plot(range(len(ions)), ys, color=cs[i], marker=SYMBOLS[i%len(SYMBOLS)], lw=s.plot_width,
                markersize=np.sqrt(s.marker_size), label=r.Sample_ID)
    ax.set_yscale('log'); ax.set_xticks(range(len(ions)), ions)
    style_axis(ax, s, 'Major ions', 'Concentration (meq/L)')
    title(ax, 'Schoeller plot', s)
    if s.show_legend:
        lg = ax.legend(loc='upper left', bbox_to_anchor=(1.02,1), frameon=False, fontsize=s.legend_size)
        for text in lg.get_texts():
            text.set_fontweight(_weight(s)); text.set_color('black'); text.set_fontfamily(s.family)
    return fig


def facies_hardness(df, s, kind='Hardness'):
    fig = base(s); ax = fig.add_subplot(111)
    if kind == 'Hardness':
        vals = df['Hardness_as_CaCO3_mgL']
        cs = palette(s, len(df))
        ax.bar(df.Sample_ID, vals, color=cs, edgecolor='black', linewidth=1)
        for level in [60,120,180]:
            ax.axhline(level, color='#555555', linestyle='--', linewidth=1)
        style_axis(ax, s, 'Samples', 'Hardness (mg/L as CaCO₃)', grid=False)
    else:
        vc = df.Facies.value_counts()
        ax.barh(vc.index, vc.values, color=palette(s, len(vc)), edgecolor='black')
        style_axis(ax, s, 'Number of samples', 'Hydrochemical facies', grid=False)
    title(ax, kind, s)
    ax.tick_params(axis='x', labelrotation=35 if kind == 'Hardness' else 0)
    return fig


def quality(df, s, tolerance=10):
    fig = base(s); ax = fig.add_subplot(111)
    status_colors = {'PASS':'#009E73','CHECK':'#D55E00','INSUFFICIENT IONS':'#999999'}
    cs = [status_colors.get(x,'#888888') for x in df.QA_status]
    ax.bar(df.Sample_ID, df.CBE_pct, color=cs, edgecolor='black', linewidth=.9)
    ax.axhline(tolerance, color='black', ls='--', lw=1.3)
    ax.axhline(-tolerance, color='black', ls='--', lw=1.3)
    style_axis(ax, s, 'Samples', 'Charge balance error (%)', grid=False)
    title(ax, 'Major-ion charge balance', s)
    ax.tick_params(axis='x', labelrotation=35)
    return fig


def save_bytes(fig, fmt, dpi=300):
    output = io.BytesIO()
    fig.savefig(output, format=fmt, dpi=dpi, facecolor='white', bbox_inches='tight', pad_inches=.12)
    return output.getvalue()
