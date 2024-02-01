# figure imports
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd

def prf_violins_plot(data, subject, ecc_th=[None,None], size_th=[None,None], rsq_th=[None,None]) :
    """
    Make violins plots for pRF r2/loo_r2, ecc and size

    Parameters
    ----------
    data : A data frame with prf_rsq, prf_ecc, prf_size, prf_loo_r2, rois and subject columns
    
    Returns
    -------
    fig : the figure 
    """
    # Replace all data outer threshold with NaN data
    data.loc[(data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) | 
             (data.prf_size < size_th[0]) | (data.prf_size > size_th[1]) | 
             (data.prf_loo_r2 <=rsq_th[0])] = np.nan
    
    data = data.dropna()
    rois = pd.unique(data.rois)

    
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    
    
    rows, cols = 2,2
    fig_height, fig_width = 1080,1920
    

    fig = make_subplots(rows=rows, cols=cols, 
                        print_grid=False, 
                        vertical_spacing=0.08, 
                        horizontal_spacing=0.1)
    
    for j, roi in enumerate(rois):
        
        df = data.loc[(data.subject == subject) & (data.rois == roi)]
        
        # df = df.sort_values('prf_rsq_loo', ascending=False)
        # df = df.head(250)

        # pRF loo r2
        fig.add_trace(go.Violin(x=df.rois[df.rois==roi], 
                                y=df.prf_loo_r2, 
                                name=roi, 
                                showlegend=True, 
                                legendgroup='loo', 
                                points=False, 
                                scalemode='width', 
                                width=0.75, 
                                side='negative', 
                                line_color = roi_colors[j], 
                                meanline_visible=True), 
                      row=1, col=1)
        
        
        # pRF r2
        fig.add_trace(go.Violin(x=df.rois[df.rois==roi], 
                                y=df.prf_rsq, 
                                name=roi, 
                                showlegend=False, 
                                legendgroup='avg', 
                                points=False, 
                                scalemode='width', 
                                width=0.75, 
                                side='positive', 
                                line_color = roi_colors[j], 
                                meanline_visible=True, 
                                fillcolor='rgb(255,255,255)'), 
                      row=1, col=1)
        
        
        # pRF size
        fig.add_trace(go.Violin(x=df.rois[df.rois==roi], 
                                y=df.prf_size, 
                                name=roi, 
                                showlegend=False, 
                                legendgroup='avg', 
                                points=False, 
                                scalemode='width', 
                                width=0.75, 
                                line_color = roi_colors[j], 
                                meanline_visible=True), 
                      row=1, col=2)
        
        # pRF ecc
        fig.add_trace(go.Violin(x=df.rois[df.rois==roi], 
                                y=df.prf_ecc, 
                                name=roi, 
                                showlegend=False, 
                                legendgroup='avg', 
                                points=False, 
                                scalemode='width', 
                                width=0.75,  
                                line_color = roi_colors[j], 
                                meanline_visible=True), 
                      row=2, col=1)
 
        
        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(range=[0,1], 
                         title_text='r2', 
                         row=1, col=1)
        
        fig.update_yaxes(range=[0,20], 
                         title_text='pRF size (dva)', 
                         row=1, col=2)
        
        fig.update_yaxes(range=[0,15], 
                         nticks=4, 
                         title_text='pRF eccentricity (dva)', 
                         row=2, col=1)
        
        fig.update_xaxes(showline=True, 
                         ticklen=0, 
                         linecolor=('rgba(255,255,255,0)'), 
                         tickfont=dict(size=18))
        
        fig.update_traces(spanmode='manual', 
                          span=[0,1], 
                          row=1, col=1)  
        
        fig.update_traces(spanmode='manual', 
                          span=[0.1,20], 
                          row=1, col=2)
        
        fig.update_traces(spanmode='manual', 
                          span=[0,15], 
                          row=2, col=1)
        
        
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=True, 
                      legend=dict(orientation="h", 
                                  yanchor='top', 
                                  y=1.15, 
                                  xanchor='left', 
                                  x=0.22, 
                                  traceorder='normal', 
                                  itemwidth=50), 
                      template='simple_white', 
                      font=dict(size=16))
    
    return fig 

def prf_ecc_size_plot(data, subject, ecc_th=[None,None], size_th=[None,None], rsq_th=[None,None]) :
    """
    Make violins plots for pRF r2/loo_r2, ecc and size

    Parameters
    ----------
    data : A data frame with prf_rsq, prf_ecc, prf_size, prf_loo_r2, rois and subject columns
    
    Returns
    -------
    fig : the figure 
    """
    
    from maths_utils import weighted_regression
    import scipy
    from scipy import stats
    
    
    fig_height, fig_width = 400, 800
    rows, cols = 1,4
    
    # Replace all data outer threshold with NaN data
    data.loc[(data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) | 
             (data.prf_size < size_th[0]) | (data.prf_size > size_th[1]) | 
             (data.prf_loo_r2 <=rsq_th[0])] = np.nan
    
    data = data.dropna()

    # Define colors
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    
    lines = [['V1', 'V2', 'V3'],['V3AB', 'LO', 'VO'],['hMT+', 'iIPS', 'sIPS'],['iPCS', 'sPCS', 'mPCS']]

    fig = make_subplots(rows=rows, cols=cols, print_grid=False)
    for l, line_label in enumerate(lines):
        for j, roi in enumerate(line_label):
            
            # Sorting best datas
            df = data.loc[(data.subject == subject) & (data.rois == roi)]
            
            # Parametring colors
            roi_color = roi_colors[j + l * 3]
            roi_color_opac = f"rgba{roi_color[3:-1]}, 0.15)"
            
            # Grouping by eccentricities
            df_grouped = df.groupby(pd.cut(df['prf_ecc'], bins=np.arange(0, 17.5, 2.5)))
            df_sorted = df.sort_values('prf_ecc')
            
            ecc_mean = np.array(df_grouped['prf_ecc'].mean())
            sd_mean = np.array(df_grouped['prf_size'].mean())
            r2_mean = np.array(df_grouped['prf_loo_r2'].mean())
            
            # CI95 for each group of ecc
            ci = df_grouped['prf_size'].apply(lambda x: stats.t.interval(0.95, len(x)-1, loc=np.nanmean(x), scale=scipy.stats.sem(x, nan_policy='omit')))
            upper_bound = np.array(ci.apply(lambda x: x[1]))
            lower_bound = np.array(ci.apply(lambda x: x[0]))
            
            # Linear regression
            slope, intercept = weighted_regression(ecc_mean, sd_mean, r2_mean)
            slope_upper, intercept_upper = weighted_regression(ecc_mean[np.where(~np.isnan(upper_bound))], upper_bound[~np.isnan(upper_bound)], r2_mean[np.where(~np.isnan(upper_bound))])
            slope_lower, intercept_lower = weighted_regression(ecc_mean[np.where(~np.isnan(lower_bound))], lower_bound[~np.isnan(lower_bound)], r2_mean[np.where(~np.isnan(lower_bound))])
            line = slope[0][0] * np.array(df_sorted.prf_ecc) + intercept[0]
            line_upper = slope_upper[0][0] * np.array(df_sorted.prf_ecc) + intercept_upper[0]
            line_lower = slope_lower[0][0] * np.array(df_sorted.prf_ecc) + intercept_lower[0]

            fig.add_trace(go.Scatter(x=np.array(df_sorted.prf_ecc), y=line, mode='lines', name=roi, legendgroup=roi, 
                                     line=dict(color=roi_color, width=3), showlegend=False), 
                          row=1, col=l+1)

            # Error area
            fig.add_trace(go.Scatter(x=np.concatenate([df_sorted.prf_ecc, df_sorted.prf_ecc[::-1]]), 
                                     y=np.concatenate([list(line_upper), list(line_lower[::-1])]), 
                                     mode='lines', fill='toself', fillcolor=roi_color_opac, 
                                     line=dict(color=roi_color_opac, width=0), showlegend=False), 
                          row=1, col=l+1)

            # Markers
            fig.add_trace(go.Scatter(x=ecc_mean, y=sd_mean, mode='markers', 
                                     error_y=dict(type='data', array=ci.apply(lambda x: (x[1] - x[0]) / 2).tolist(), visible=True, thickness=3, width=0, color =roi_color),
                                     marker=dict(color='white', size=8, 
                                                 line=dict(color=roi_color,width=3)
                                                ), showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=1, y=15-j*1.5, text=roi, xanchor='left',
                                              showarrow=False, font=dict(color=roi_color, size=12))
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(title_text='pRF size (dva)', row=1, col=1)
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[0,15], row=1, col=l+1)
        fig.update_yaxes(range=[0,15])
        fig.update_layout(height=fig_height, width=fig_width, showlegend=False, template='simple_white')
        
    return fig
    
        
        