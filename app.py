"""HydroGeoAI Publication Studio v6 — transparent offline scientific figure prototyping."""
from dataclasses import asdict,replace
from zipfile import ZipFile,ZIP_DEFLATED
import io,json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import streamlit as st
from PIL import Image
from chemistry import analyze,REQUIRED,IONS
from figures import Style,piper,correlation,pca,hca,schoeller,facies_hardness,quality,save_bytes
from figures_v6 import durov,chadha,stiff

st.set_page_config(page_title='HydroGeoAI | Publication Studio v6',page_icon='💧',layout='wide')
st.markdown('''<style>
.block-container{max-width:1580px;padding-top:1.3rem}
[data-testid="stSidebar"]{background:#f5f9fc}
[data-testid="stMetric"]{background:#f1f7fa;border:1px solid #dbe7ec;padding:13px;border-radius:12px}
.hero{background:linear-gradient(105deg,#092b41,#115c76);border-radius:15px;padding:25px;color:white;margin-bottom:12px}
.hero h1{color:white;margin:0;font-size:2.25rem}.hero p{color:#daf0f7;margin:5px 0 0}
</style>''',unsafe_allow_html=True)
st.markdown('<div class="hero"><h1>💧 HydroGeoAI — Publication Studio v6</h1><p>Piper · Durov · Chadha · Stiff · Schoeller · PCA · HCA · Scientific figure export</p></div>',unsafe_allow_html=True)
st.caption('Local prototype • scientifically transparent calculations • interactive image zoom/pan • vector-first publication export')

@st.cache_data
def demo_data():
    return pd.DataFrame([
    ['GW-01',72,24,14,2.4,34,15,220,1],['GW-02',54,18,10,1.6,170,18,65,2],
    ['GW-03',84,30,19,3,55,22,240,1.5],['SW-01',32,9,8,1.1,115,12,37,1],
    ['SW-02',40,12,11,1.7,135,16,58,1.4],['GW-04',65,23,13,2.1,40,25,185,2.2],
    ['SW-03',26,6,7,.9,95,10,23,1.8],['GW-05',44,17,12,1.9,100,20,120,1]],
    columns=['Sample_ID','Ca','Mg','Na','K','HCO3','Cl','SO4','NO3'])

presets={
 'High-contrast journal':dict(fontsize=13,ticksize=11,annotation_size=10,axis_width=2.2,plot_width=1.7,palette='Colorblind'),
 'Thesis / dissertation':dict(fontsize=15,ticksize=12,annotation_size=11,axis_width=2.3,plot_width=1.9,palette='Tableau'),
 'Conference / presentation':dict(fontsize=19,ticksize=15,annotation_size=14,axis_width=2.8,plot_width=2.4,palette='Tableau'),
 'Consultancy report':dict(fontsize=15,ticksize=12,annotation_size=11,axis_width=2.1,plot_width=2.0,palette='Colorblind'),
 'Minimal grayscale':dict(fontsize=13,ticksize=11,annotation_size=10,axis_width=2.2,plot_width=1.6,palette='Grayscale'),
}

with st.sidebar:
    st.header('01 · Data & QA')
    upload=st.file_uploader('Upload CSV / XLSX',type=['csv','xlsx'])
    cbe_tolerance=st.slider('Charge balance threshold (%)',2,15,10)
    st.download_button('Download synthetic example',demo_data().to_csv(index=False),file_name='HydroGeoAI_demo.csv')
    st.caption('Input concentrations in mg/L. Required: '+', '.join(REQUIRED)+'. Optional CO3, NO3, TDS.')
    st.divider()
    st.header('02 · Figure theme')
    preset=st.selectbox('Design preset',list(presets))
    pdict=presets[preset]
    font=st.selectbox('Font family',['DejaVu Sans','Arial','Times New Roman','Liberation Sans','Georgia'])
    fontsize=st.slider('Axis label font',8,28,pdict['fontsize'])
    ticksize=st.slider('Axis number font',7,26,pdict['ticksize'])
    annotationsize=st.slider('Annotation font',7,24,pdict['annotation_size'])
    legend_size=st.slider('Legend font',7,24,11)
    bold=st.checkbox('Bold text',True)
    axes=st.slider('Black axis / boundary thickness',.8,4.5,float(pdict['axis_width']),.1)
    lines=st.slider('Data / branch line thickness',.8,4.5,float(pdict['plot_width']),.1)
    markers=st.slider('Marker size',20,220,80)
    palette=st.selectbox('Color palette',['Colorblind','Tableau','Grayscale'],index=['Colorblind','Tableau','Grayscale'].index(pdict['palette']))
    legend_loc=st.selectbox('Legend location',['Outside right','Inside upper right','Inside upper left','Inside lower right','Inside lower left'])
    legend_on=st.checkbox('Show legend',True)
    legend_columns=st.slider('Legend columns',1,4,1)
    legend_box=st.checkbox('Show legend box',False)
    title_on=st.checkbox('Show figure title in export',False)
    base_w=st.slider('Figure width (inches)',6.,17.,10.,.5)
    base_h=st.slider('Figure height (inches)',5.,14.,8.,.5)
    st.divider()
    st.header('03 · Preview & export')
    interactive=st.checkbox('Interactive mouse zoom / pan preview',True)
    st.caption('Drag to zoom; use Plotly toolbar to pan or reset. Preview uses a raster image; SVG/PDF exports remain vector.')
    fmt=st.selectbox('Export type',['SVG','PDF','PNG'])
    dpi=st.selectbox('PNG DPI',[300,450,600,900],index=1)

try:
    raw=demo_data() if upload is None else (pd.read_csv(upload) if upload.name.lower().endswith('.csv') else pd.read_excel(upload))
    df=analyze(raw,cbe_tolerance)
except Exception as e:
    st.error(f'Input validation: {e}');st.stop()

style=Style(family=font,fontsize=fontsize,ticksize=ticksize,annotation_size=annotationsize,
            legend_size=legend_size,bold=bold,axis_width=axes,plot_width=lines,marker_size=markers,
            palette=palette,show_title=title_on,show_legend=legend_on,
            legend_location=legend_loc,legend_columns=legend_columns,legend_box=legend_box,
            width=base_w,height=base_h,dpi=dpi)

m=st.columns(4)
m[0].metric('Samples',len(df))
m[1].metric('CBE pass',int((df.QA_status=='PASS').sum()))
m[2].metric('QA flagged',int((df.QA_status!='PASS').sum()))
m[3].metric('Median |CBE|',f'{df.CBE_pct.abs().median():.1f}%')
if upload is None:st.info('Synthetic demonstration data loaded. Upload your own measurements to analyze real samples.')
if (df.QA_status!='PASS').any():st.warning('Some samples fail the selected charge-balance screening threshold.')

figures={}

def preview_zoom(fig,key):
    # A graphical preview with genuine mouse zoom/pan of the entire rendered figure.
    # It does not change the scientific x/y data and does not enable element-level dragging.
    png=save_bytes(fig,'png',115)
    img=Image.open(io.BytesIO(png)).convert('RGB')
    w,h=img.size
    plot=go.Figure()
    plot.add_layout_image(dict(source=img,xref='x',yref='y',x=0,y=h,sizex=w,sizey=h,
                               sizing='stretch',layer='below'))
    plot.update_layout(width=min(w,1300),height=min(h,800),margin=dict(l=0,r=0,t=0,b=0),
                       dragmode='zoom',paper_bgcolor='white',plot_bgcolor='white')
    plot.update_xaxes(visible=False,range=[0,w],constrain='domain',fixedrange=False)
    plot.update_yaxes(visible=False,range=[0,h],scaleanchor='x',scaleratio=1,fixedrange=False)
    st.plotly_chart(plot,use_container_width=True,config={'scrollZoom':True,'displaylogo':False,'modeBarButtonsToAdd':['pan2d','zoom2d','resetScale2d']},key='plot_'+key)

def render(key,fig,note=None):
    if fig is None:st.warning('Insufficient data for this figure.');return
    if interactive:preview_zoom(fig,key)
    else:st.pyplot(fig,use_container_width=True)
    if note:st.caption(note)
    ext=fmt.lower()
    st.download_button('Download '+key+' ('+fmt+')',save_bytes(fig,ext,dpi),file_name='HydroGeoAI_'+key+'.'+ext,
                       mime={'svg':'image/svg+xml','pdf':'application/pdf','png':'image/png'}[ext],key='save_'+key)
    figures[key]=fig


tabs=st.tabs(['Piper','Durov','Chadha','Stiff & Schoeller','PCA & HCA','Correlation','Facies & QA','Exports'])
with tabs[0]:
    st.subheader('Piper facies visualization')
    a,b,c=st.columns(3)
    with a:fill=a.checkbox('Colored facies fields',True,key='pf')
    with b:grid=b.checkbox('20% internal grid',True,key='pg')
    with c:alpha=c.slider('Color intensity',.1,1.,.82,.05,key='pa')
    labels=st.checkbox('Show Piper sample IDs',False)
    water_labels=st.checkbox('Show water-type text',False)
    sp=replace(style,fill_piper_fields=fill,show_grid=grid,field_alpha=alpha,show_labels=labels,water_type_labels=water_labels)
    render('Piper',piper(df,sp),'Colored boundaries show compositional regions. Facies classifications remain based on measured major-ion fractions, not on decorative region colors.')
with tabs[1]:
    st.subheader('Durov composition diagram')
    a,b,c=st.columns(3)
    with a:zones=st.checkbox('Color composition zones',True,key='dz')
    with b:zones_num=st.checkbox('Number nine square zones',False)
    with c:label=st.checkbox('Sample IDs',False,key='dl')
    render('Durov',durov(df,replace(style,show_labels=label),color_zones=zones,labels=zones_num),
           'Central nine cells are descriptive composition fields, not validated reaction pathways. Top = cations; left = anions; square axes show alkalis and strong-acid anion fractions.')
with tabs[2]:
    st.subheader('Chadha hydrochemical classification')
    a,b=st.columns(2)
    with a:zones=st.checkbox('Colored hydrochemical quadrants',True)
    with b:label=st.checkbox('Sample ID labels',False,key='chlabel')
    render('Chadha',chadha(df,replace(style,show_labels=label),colored=zones),'Quadrants use normalized ionic group differences (meq%).')
with tabs[3]:
    st.subheader('Stiff polygons')
    a,b=st.columns(2)
    with a:mode=st.radio('View',['One sample','Overlay comparison'],horizontal=True)
    with b:sample=st.selectbox('Sample',df.Sample_ID.astype(str).tolist())
    render('Stiff',stiff(df,style,selected=sample,overlay=(mode=='Overlay comparison')),
           'Left side cations; right side anions, both in meq/L. A consistent concentration range is used across samples.')
    st.subheader('Schoeller diagram')
    render('Schoeller',schoeller(df,style))
with tabs[4]:
    st.subheader('Principal component analysis')
    a,b,c=st.columns(3)
    with a:mode=st.selectbox('PCA plot',['Biplot','Scores','Loadings'])
    with b:color=st.checkbox('Colored loading vectors',True)
    with c:show_labels=st.checkbox('Show sample labels',True)
    st.caption('Edit the row for any individual sample or variable. Offset values use screen points; rotation is in degrees.')
    variables=[i for i in IONS if df[i].nunique(dropna=True)>1]
    label_names=['sample:'+str(v) for v in df.Sample_ID]+['variable:'+v for v in variables]
    defaults=pd.DataFrame([dict(Label=n,Visible=True,X_offset=7,Y_offset=6,Rotation=0) for n in label_names])
    if st.session_state.get('labels_key')!=tuple(label_names):
        st.session_state['labels_key']=tuple(label_names)
        st.session_state['labels_table']=defaults
    edited=st.data_editor(st.session_state['labels_table'],hide_index=True,use_container_width=True,
                          disabled=['Label'],key='label_editor',num_rows='fixed',
                          column_config={'Visible':st.column_config.CheckboxColumn(),
                                         'X_offset':st.column_config.NumberColumn(min_value=-150,max_value=150,step=1),
                                         'Y_offset':st.column_config.NumberColumn(min_value=-150,max_value=150,step=1),
                                         'Rotation':st.column_config.NumberColumn(min_value=-180,max_value=180,step=5)})
    st.session_state['labels_table']=edited.copy()
    overrides={str(r.Label):dict(show=bool(r.Visible),dx=float(r.X_offset),dy=float(r.Y_offset),rotation=float(r.Rotation)) for r in edited.itertuples()}
    sp=replace(style,pca_mode=mode,color_pca_vectors=color,show_labels=show_labels,label_overrides=overrides)
    pcfig,payload=pca(df,sp)
    render('PCA',pcfig,'PCA uses standardized concentrations; biplot arrows are scaled for visualization, not directly comparable numerical axes.')
    if payload:
        c1,c2=st.columns(2)
        with c1:st.dataframe(payload['explained'],hide_index=True,use_container_width=True)
        with c2:st.dataframe(payload['loadings'].round(3),use_container_width=True)
        st.download_button('Export PCA scores',payload['scores'].to_csv(index=False),file_name='PCA_scores.csv')
        st.download_button('Export PCA loadings',payload['loadings'].to_csv(),file_name='PCA_loadings.csv')
    st.divider();st.subheader('Hierarchical cluster analysis')
    k=st.slider('Cluster count',2,max(2,min(len(df)-1,12)),min(3,max(2,len(df)-1)))
    color_hca=st.checkbox('Cluster-colored branches',True)
    hc,membership=hca(df,replace(style,hca_colored_clusters=color_hca),k)
    render('HCA',hc,'Ward linkage on z-standardized ion concentrations; branch colors reflect the selected cut.')
    if membership is not None:
        st.dataframe(membership,hide_index=True,use_container_width=True)
        st.download_button('Export HCA membership',membership.to_csv(index=False),file_name='HCA_membership.csv')
with tabs[5]:
    st.subheader('Major-ion correlation')
    a,b=st.columns(2)
    with a:method=st.selectbox('Statistic',['pearson','spearman'])
    with b:values=st.checkbox('Show coefficients',True)
    render('Correlation',correlation(df,style,values,method))
with tabs[6]:
    a,b=st.columns(2)
    with a:render('Facies',facies_hardness(df,style,'Facies'))
    with b:render('Hardness',facies_hardness(df,style,'Hardness'))
    render('Charge_balance',quality(df,style,cbe_tolerance))
    st.dataframe(df[['Sample_ID','CBE_pct','QA_status','Facies']].round(2),hide_index=True,use_container_width=True)
with tabs[7]:
    st.subheader('Batch publication export')
    st.download_button('Calculated chemistry CSV',df.to_csv(index=False),file_name='HydroGeoAI_chemistry.csv')
    settings=asdict(style)
    settings.pop('label_overrides',None)
    st.download_button('Download theme settings JSON',json.dumps(settings,indent=2),file_name='HydroGeoAI_theme.json')
    workbook=io.BytesIO()
    with pd.ExcelWriter(workbook,engine='openpyxl') as w:
        df.to_excel(w,'Results',index=False)
        pd.DataFrame(list(settings.items()),columns=['Setting','Value']).to_excel(w,'Figure_settings',index=False)
    st.download_button('Results + settings XLSX',workbook.getvalue(),file_name='HydroGeoAI_results.xlsx')
    archive=io.BytesIO()
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for name,fig in figures.items():z.writestr(name+'.'+fmt.lower(),save_bytes(fig,fmt.lower(),dpi))
        z.writestr('results.csv',df.to_csv(index=False))
    st.download_button('All diagrams and results ZIP',archive.getvalue(),file_name='HydroGeoAI_v6_exports.zip')
    st.caption('Scientific verification, user accounts, paid billing, and server-side project storage are NOT included in this local prototype.')
for fig in figures.values():plt.close(fig)
