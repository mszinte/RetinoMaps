# figure imports
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd


def prf_violins_plot(data, subject, fig_height, fig_width, ecc_th=[None,None], size_th=[None,None], rsq_th=[None,None], pcm_th=[None,None]) :
    """
    Make violins plots for pRF r2/loo_r2, ecc and size

    Parameters
    ----------
    data : A data frame with prf_rsq, prf_size, prf_loo_r2, pcm, rois and subject columns
    
    Returns
    -------
    fig : the figure 
    """
    data = data.copy()
    # Replace all data outer threshold with NaN data
    data.loc[(data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) | 
             (data.prf_size < size_th[0]) | (data.prf_size > size_th[1]) | 
             (data.pcm < pcm_th[0]) | (data.pcm > pcm_th[1]) |
             (data.prf_loo_r2 <=rsq_th[0])] = np.nan
    
    data = data.dropna()
    rois = pd.unique(data.rois)

    
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    
    
    rows, cols = 2,2
    # fig_height, fig_width = 1080,1920
    

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
        
        # pRF n
        fig.add_trace(go.Violin(x=df.rois[df.rois==roi], 
                                y=df.prf_n, 
                                name=roi, 
                                showlegend=False, 
                                legendgroup='avg', 
                                points=False, 
                                scalemode='width', 
                                width=0.75,  
                                line_color = roi_colors[j], 
                                meanline_visible=True), 
                      row=2, col=1)
        
        # pcm
        fig.add_trace(go.Violin(x=df.rois[df.rois==roi], 
                                y=df.pcm, 
                                name=roi, 
                                showlegend=False, 
                                legendgroup='avg', 
                                points=False, 
                                scalemode='width', 
                                width=0.75,  
                                line_color = roi_colors[j], 
                                meanline_visible=True), 
                      row=2, col=2)
 
        
        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(range=[0,1],
                         nticks=5, 
                         title_text='R<sup>2</sup>', 
                         row=1, col=1)
        
        fig.update_yaxes(range=[0,20], 
                         nticks=5, 
                         title_text='pRF size (dva)', 
                         row=1, col=2)
        
        fig.update_yaxes(range=[0,2], 
                         nticks=5, 
                         title_text='pRF n', 
                         row=2, col=1)
        
        fig.update_yaxes(range=[0,20], 
                         nticks=5, 
                         title_text='pCM', 
                         row=2, col=2)
        
        fig.update_xaxes(showline=True, 
                         ticklen=0, 
                         linecolor=('rgba(255,255,255,0)'), 
                         tickfont=dict(size=18))
        
        # fig.update_traces(spanmode='manual', 
        #                   span=[0,1], 
        #                   row=1, col=1)  
        
        # fig.update_traces(spanmode='manual', 
        #                   span=[0.1,20], 
        #                   row=1, col=2)
        
        # fig.update_traces(spanmode='manual', 
        #                   span=[0,2], 
        #                   row=2, col=1)
        
        
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




def prf_ecc_size_plot(data, subject, fig_height, fig_width, ecc_th=[None,None], size_th=[None,None], rsq_th=[None,None]) :
    """
    Make violins plots for pRF r2/loo_r2, ecc and size

    Parameters
    ----------
    data : A data frame with prf_rsq, prf_ecc, prf_size, prf_loo_r2, rois and subject columns
    
    Returns
    -------
    fig : the figure 
    """
    
    from maths_utils import weighted_regression, bootstrap_ci_mean
    data = data.copy()
    
    
    # fig_height, fig_width = 1080, 190
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
            

            ci = df_grouped['prf_size'].apply(lambda x: bootstrap_ci_mean(x))
            upper_bound = np.array(ci.apply(lambda x: x[1] if not np.isnan(x[1]) else np.nan))
            lower_bound = np.array(ci.apply(lambda x: x[0] if not np.isnan(x[0]) else np.nan))
            
            # Linear regression
            slope, intercept = weighted_regression(ecc_mean, 
                                                    sd_mean, 
                                                    r2_mean,
                                                    model='linear')
            
            slope_upper, intercept_upper = weighted_regression(ecc_mean[np.where(~np.isnan(upper_bound))], 
                                                                upper_bound[~np.isnan(upper_bound)], 
                                                                r2_mean[np.where(~np.isnan(upper_bound))], 
                                                                model='linear')
            
            slope_lower, intercept_lower = weighted_regression(ecc_mean[np.where(~np.isnan(lower_bound))], 
                                                                lower_bound[~np.isnan(lower_bound)], 
                                                                r2_mean[np.where(~np.isnan(lower_bound))], 
                                                                model='linear')
            
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
            fig.add_trace(go.Scatter(x=ecc_mean, 
                                     y=sd_mean, mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=ci.apply(lambda x: (x[1] - x[0]) / 2).tolist(), 
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                      marker=dict(color='white', 
                                                  size=8, 
                                                  line=dict(color=roi_color, 
                                                            width=3)), 
                                      showlegend=False), 
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

def prf_ecc_pcm_plot(data, subject, fig_height, fig_width, ecc_th=[None,None], pcm_th=[None,None], rsq_th=[None,None]) :
    """
    Make figure of ecc pcm relation 

    Parameters
    ----------
    data : A data frame with prf_rsq, prf_ecc, prf_size, prf_loo_r2, pcm, rois and subject columns
    
    Returns
    -------
    fig : the figure 
    """

    from maths_utils import weighted_regression, bootstrap_ci_mean
    data = data.copy()

    # fig_height, fig_width = 1080, 1920
    rows, cols = 1,4
    
    # Replace all data outer threshold with NaN data
    data.loc[(data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) | 
              (data.pcm < pcm_th[0]) | (data.pcm > pcm_th[1]) | 
              (data.prf_loo_r2 <= rsq_th[0])] = np.nan
    
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
            sd_mean  = np.array(df_grouped['pcm'].mean())
            r2_mean  = np.array(df_grouped['prf_loo_r2'].mean())
            
            # CI95 for each group of ecc
            ci = df_grouped['pcm'].apply(lambda x: bootstrap_ci_mean(x))
            upper_bound = np.array(ci.apply(lambda x: x[1] if not np.isnan(x[1]) else np.nan))
            lower_bound = np.array(ci.apply(lambda x: x[0] if not np.isnan(x[0]) else np.nan))
            
            # Linear regression
            slope, intercept = weighted_regression(ecc_mean, 
                                                    sd_mean, 
                                                    r2_mean, 
                                                    model='pcm')
            
            slope_upper, intercept_upper = weighted_regression(ecc_mean[~np.isnan(upper_bound)], 
                                                                upper_bound[~np.isnan(upper_bound)], 
                                                                r2_mean[~np.isnan(upper_bound)], 
                                                                model='pcm')
            
            slope_lower, intercept_lower = weighted_regression(ecc_mean[~np.isnan(lower_bound)], 
                                                                lower_bound[~np.isnan(lower_bound)], 
                                                                r2_mean[~np.isnan(lower_bound)], 
                                                                model='pcm')
            
            line = 1 / (slope * np.array(df_sorted.prf_ecc)) + intercept
            line_upper = 1 / (slope_upper * np.array(df_sorted.prf_ecc)) + intercept_upper
            line_lower = 1 / (slope_lower * np.array(df_sorted.prf_ecc)) + intercept_lower

            fig.add_trace(go.Scatter(x=np.array(df_sorted.prf_ecc), 
                                      y=line, 
                                      mode='lines', 
                                      name=roi, 
                                      legendgroup=roi, 
                                      line=dict(color=roi_color, width=3), 
                                      showlegend=False), 
                          row=1, col=l+1)

            # Error area
            fig.add_trace(go.Scatter(x=np.concatenate([df_sorted.prf_ecc, df_sorted.prf_ecc[::-1]]), 
                                      y=np.concatenate([list(line_upper), list(line_lower[::-1])]), 
                                      mode='lines', fill='toself', fillcolor=roi_color_opac, 
                                      line=dict(color=roi_color_opac, width=0), showlegend=False), 
                          row=1, col=l+1)

            # Markers
            fig.add_trace(go.Scatter(x=ecc_mean, 
                                      y=sd_mean, 
                                      mode='markers', 
                                      error_y=dict(type='data', array=ci.apply(lambda x: (x[1] - x[0]) / 2).tolist(), visible=True, thickness=3, width=0, color=roi_color),
                                      marker=dict(color='white', size=8, line=dict(color=roi_color,width=3)), 
                                      showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=10, 
                                              y=10-j*1.5, 
                                              text=roi, 
                                              xanchor='left',
                                              showarrow=False, 
                                              font=dict(color=roi_color, size=12))
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(title_text='pCM (mm/dva)', row=1, col=1)
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[0,15], row=1, col=l+1)
        fig.update_yaxes(range=[0,10])
        fig.update_layout(height=fig_height, width=fig_width, showlegend=False, template='simple_white')
        
    return fig





def prf_polar_plot(data, subject, fig_height, fig_width, ecc_th=[None,None], size_th=[None,None], rsq_th=[None,None]) :    
    """
     Make polar plots
    
     Parameters
     ----------
     data : A data frame with prf_rsq, prf_ecc, prf_size, prf_loo_r2, rois and subject columns
     
     Returns
     -------
     figs : a list of three figures
     hemispheres : a list of corresponding hemispheres
    """
    data = data.copy()
    data['prf_angle'] = np.angle(data.polar_real + 1j * data.polar_imag)
    
    # Replace all data outer threshold with NaN data
    data.loc[(data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) | 
             (data.prf_size < size_th[0]) | (data.prf_size > size_th[1]) | 
             (data.prf_loo_r2 <=rsq_th[0])] = np.nan
    data = data.dropna()
    
    rois = pd.unique(data.rois)
    hemis = np.append(pd.unique(data.hemi), 'brain')
    
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    
    rows, cols = 1, 12
    # fig_height, fig_width = 300, 1920
    specs = [[{'type': 'polar'}] * cols]
    num_slices = 12
    
    figs = []
    hemispheres = []
    for i, hemi in enumerate(hemis):
        fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs)
        
        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'brain' else [hemi]
    
        for j, roi in enumerate(rois):
    
            df = data.loc[(data.subject == subject) & (data.rois==roi) & (data.hemi.isin(hemi_values))]
            df = df.sort_values('prf_loo_r2', ascending=False)
            df = df.head(250)
    
            #Conversion
            df.prf_angle = np.degrees(df.prf_angle)
            df.prf_angle = np.mod(df.prf_angle, 360)
    
            # Parts of polar angles and number of voxels in each part
            theta_slices = np.linspace(0, 360, num_slices+1, endpoint=True)
            voxels_counts, _ = np.histogram(df.prf_angle, bins=theta_slices)
    
            # barpolar
            fig.add_trace(go.Barpolar(r=voxels_counts, 
                                      theta=theta_slices, 
                                      width=30, 
                                      marker_color=roi_colors[j], 
                                      marker_line_color='black', 
                                      marker_line_width=1, 
                                      opacity=0.8), 
                          row=1, col=j+1)
    
        # Define parameters
        fig.update_polars(angularaxis=dict(visible=False), 
                          radialaxis=dict(visible=False))
        
        fig.update_layout(title='{}'.format(hemi), 
                          height=fig_height, 
                          width=fig_width, 
                          showlegend=False, 
                          template='simple_white')
        figs.append(fig)
        hemispheres.append(hemi)
        
    return figs, hemispheres



def prf_contralaterality_plot(data, subject, fig_height, fig_width, ecc_th=[None,None], size_th=[None,None], rsq_th=[None,None]) :    
    """
     Make polar plots
    
     Parameters
     ----------
     data : A data frame with prf_rsq, prf_ecc, prf_size, prf_loo_r2, rois and subject columns
     
     Returns
     -------
     figs : a list of three figures
     hemispheres : a list of corresponding hemispheres
    """
    data = data.copy()
    # Replace all data outer threshold with NaN data
    data.loc[(data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) | 
             (data.prf_size < size_th[0]) | (data.prf_size > size_th[1]) | 
             (data.prf_loo_r2 <=rsq_th[0])] = np.nan
    data = data.dropna()
    
    rois = pd.unique(data.rois)
    
    rows, cols = 1, 12
    # fig_height, fig_width = 300, 1920
    specs = [[{'type': 'pie'}] * cols]
    
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    

    fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs)
    for j, roi in enumerate(rois):
        
        df_rh = data.loc[(data.subject == subject) & (data.rois == roi) & (data.hemi == 'hemi-R')]
        df_lh = data.loc[(data.subject == subject) & (data.rois == roi) & (data.hemi == 'hemi-L')]
        
        # # Calculer le pourcentage de contralatéralité pour chaque hémisphère pondéré par RSQ
        # percentage_right = sum(df_rh.loc[df_rh.prf_x < 0].prf_loo_r2) / sum(df_rh.prf_loo_r2) * 100
        # percentage_left = sum(df_lh.loc[df_lh.prf_x > 0].prf_loo_r2) / sum(df_lh.prf_loo_r2) * 100

        # Calculer le pourcentage d'excentricité total
        percentage_total = (sum(df_rh.loc[df_rh.prf_x < 0].prf_loo_r2) + sum(df_lh.loc[df_lh.prf_x > 0].prf_loo_r2))/ (sum(df_rh.prf_loo_r2)+sum(df_lh.prf_loo_r2)) *100
    
        # print("{} - Contralateralité in {}: {:.1f}%".format(subject, roi, percentage_total))
        percentage_rest = 100 - percentage_total        
        
        values = [percentage_total, percentage_rest]
        
        fig.add_trace(go.Pie(values=values,
                             hole=0.5,
                             marker=dict(colors=[roi_colors[j], 'white'], line=dict(color=['black', 'white'], width=[1,0])),
                             showlegend=False,
                             pull=[0, 0.01],
                             rotation=percentage_total*3.6 if percentage_total < percentage_rest else 0,
                             ), row=1, col=j+1)

    # Define parameters
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=False, 
                      template='simple_white')
    
    return fig 
    
def categories_proportions_roi_plot(data, subject, fig_height, fig_width):
    data = data.copy()
    filtered_data = data[data['stats_final'] != 'non_responding']
    
    # Sort categories
    categories_order = ['vision', 'vision_and_pursuit_and_saccade', 'pursuit_and_saccade', 'vision_and_saccade', 'vision_and_pursuit', 'saccade', 'pursuit']
    filtered_data['stats_final'] = pd.Categorical(filtered_data['stats_final'], categories=categories_order, ordered=True)
    filtered_data = filtered_data.sort_values(['rois', 'stats_final'])

    
    #  Defines colors settings 
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    # To write make the percent visible only for choose categories 
    percent_color =  {'pursuit': 'rgba(255,255,255,0)', 
                      'saccade': 'rgba(255,255,255,0)', 
                      'pursuit_and_saccade': 'rgba(0, 0, 0, 1)', 
                      'vision': 'rgba(0, 0, 0, 1)', 
                      'vision_and_pursuit': 'rgba(255,255,255,0)', 
                      'vision_and_saccade': 'rgba(255,255,255,0)', 
                      'vision_and_pursuit_and_saccade': 'rgba(0, 0, 0, 1)'}
    
    categorie_color_map = {'pursuit': 'rgba(255,255,255,0)', 
                           'saccade': 'rgba(255,255,255,0)', 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': 'rgba(255,255,255,0)', 
                           'vision_and_saccade': 'rgba(255,255,255,0)', 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    rois = pd.unique(data.rois)
    #  Make the subplot
    # fig_height, fig_width = 300, 1920
    rows = 2 
    cols =len(rois)
    specs = [[{'type': 'domain'}] * cols,  [{'type': 'xy'}] * cols]
    
    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs, row_heights=[1,0.05])
    
    
    for i, roi in enumerate(rois):
        df_rois = filtered_data.loc[filtered_data.rois == roi]
        #  Colors for categories 
        categorie_colors = [categorie_color_map[label] for label in df_rois.stats_final]
        #  Colors for the percentages 
        percentage_colors = [percent_color[label] for label in df_rois.stats_final]
        
        
        fig.add_trace(go.Pie(labels=df_rois.stats_final, 
                             values=df_rois.vertex_surf, 
                             showlegend=False, 
                             sort=False,
                             textinfo='percent',
                             textposition='inside',
                             direction='clockwise',
                             name= roi,    
                             marker=dict(colors=categorie_colors),
                             insidetextfont=dict(color=percentage_colors),
                             hole=0.3), 
                      row=1, col=i+1)
        
        fig.add_annotation(text=roi, 
                           yshift =10,
                           showarrow=False, 
                           font=dict(size=13,color=roi_colors[i]), 
                           row=2, col=i+1)
    
    
        
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      template='simple_white')  
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(color='rgba(255,255,255,0)'))
    
    fig.update_yaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(color='rgba(255,255,255,0)'))
    
    return fig

def surface_rois_categories_plot(data, subject, fig_height, fig_width):   
    data = data.copy()
    #  Defines colors settings 
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    
    
    
    categorie_color_map = {'non_responding': stats_categories_colors[0], 
                           'pursuit': stats_categories_colors[1], 
                           'saccade': stats_categories_colors[2], 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': stats_categories_colors[5], 
                           'vision_and_saccade': stats_categories_colors[6], 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    
    #  grpup df 
    group_df_rois = data.groupby(['rois'], sort=False)['vertex_surf'].sum().reset_index()
    group_rois_categories = data.groupby(['rois', 'stats_final'], sort=False)['vertex_surf'].sum().reset_index()
    
    #  Make subplot 
    # fig_height, fig_width = 1080, 1920
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    
    
    fig.add_trace(go.Bar(x=group_df_rois.rois, 
                         y=group_df_rois.vertex_surf, 
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.1)), 
                  secondary_y=False)
    
    #Choose categories to plot
    stats_categories = ['vision',  'vision_and_pursuit_and_saccade','pursuit_and_saccade']
    for categorie in stats_categories:
        df = group_rois_categories.loc[group_rois_categories.stats_final == categorie]
    
        fig.add_trace(go.Bar(x=df.rois, 
                             y=df.vertex_surf, 
                             name=categorie,  
                             legendgroup=categorie, 
                             marker_color=categorie_color_map[categorie]), 
                      secondary_y=True) 
    
    fig.update_layout(yaxis2=dict(overlaying='y',
                                  side='right',
                                  range=[0, 5000],  
                                  showticklabels=False, 
                                  ticklen=0, 
                                  linecolor=('rgba(255,255,255,0)')))
    
    
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(size=12))      
    
    fig.update_yaxes(range=[0,5000], 
                     nticks=5, 
                     title_text='Surface in mm<sup>2</sup>',secondary_y=False)
    
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      barmode='stack',
                      showlegend=True, 
                      template='simple_white')  
    
    
    return fig 


def surface_rois_all_categories_plot(data, subject, fig_height, fig_width):  
    data = data.copy()
    #  Defines colors settings 
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    
    
    
    categorie_color_map = {'non_responding': stats_categories_colors[0], 
                           'pursuit': stats_categories_colors[1], 
                           'saccade': stats_categories_colors[2], 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': stats_categories_colors[5], 
                           'vision_and_saccade': stats_categories_colors[6], 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    #  grpup df 
    group_df = data.groupby(['rois', 'stats_final'], sort=False)['vertex_surf'].sum().reset_index()
    
    #  Figure settings
    # fig_height, fig_width = 1080, 1920
    fig = go.Figure()
    
    
    stats_categories= ['non_responding', 'vision','vision_and_pursuit_and_saccade', 'pursuit_and_saccade', 'pursuit', 'saccade', 'vision_and_pursuit', 'vision_and_saccade']
    for categorie in stats_categories:
        df = group_df.loc[group_df.stats_final == categorie]
    
        fig.add_trace(go.Bar(x=df.rois, 
                             y=df.vertex_surf, 
                             name=categorie,  
                             legendgroup=categorie, 
                             marker_color=categorie_color_map[categorie])) 
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(size=12))      
    
    fig.update_yaxes(range=[0,5000], 
                     nticks=5, 
                     title_text='Surface in mm<sup>2</sup>')
    
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      barmode='stack',
                      showlegend=True, 
                      template='simple_white')  
    
    return fig 
    
    
    

    
        
    
    
    
    
        
           
            