def eventsMatrix(design_file, task, tr):
    """
    Returns the events matrix for the GLM. Works for the Sac/Pur Localisers and Sac/Pur Visual/Endogenous Localisers
    Parameters
    ----------
    design_file         : path to the tsv/csv file containing the events
    task                : task name (SacLoc / PurLoc / SacVELoc / PurVELoc)
    tr                  : time repetition in seconds (e.g. 1.2)
    Returns
    -------
    new_events_glm - pandas DataFrame containing the events for the GLM 
    """
    import pandas as pd
    import numpy as np
    import ipdb
    deb = ipdb.set_trace
    tr_dur = tr
    events = pd.read_table(design_file)


    if 'pMF' in task:
        events_glm = events[['onset','duration','trial_type']].copy(deep=True)
        events_glm.replace({'trial_type': {3: 'Fix', 1: 'PurSac', 2: 'PurSac'}},inplace=True)
        events_glm.duration = tr_dur
        events_glm.onset = 0
        
        events_glm_groups = events_glm.groupby((events_glm.trial_type!=events_glm.trial_type.shift()).cumsum())
        
        new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
        
        for idx, group in enumerate(events_glm_groups):
            onset = np.round(group[1]['onset'][group[1].index[0]],2)
            dur = np.round(sum(group[1]['duration']),2)
            ttype = group[1]['trial_type'][group[1].index[0]]
            

            new_events_glm = new_events_glm.append(pd.Series([onset,dur,ttype],index =['onset', 'duration', 'trial_type']), ignore_index=True)

    else:
        if 'VE' not in task: # Sac/Pur Loc
            events_glm = events[['onset','duration','trial_type']].copy(deep=True)
            events_glm.replace({'trial_type': {3: 'Fix', 1: 'Sac', 2: 'Pur'}},inplace=True)
            events_glm.duration = tr_dur
            events_glm.onset = 0
            events_glm_groups = events_glm.groupby((events_glm.trial_type!=events_glm.trial_type.shift()).cumsum())

            new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
            for idx, group in enumerate(events_glm_groups):
                onset = np.round(group[1]['onset'][group[1].index[0]],2)
                dur = np.round(sum(group[1]['duration']),2)
                ttype = group[1]['trial_type'][group[1].index[0]]

                new_events_glm = new_events_glm.append(pd.Series([onset,dur,ttype],index =['onset', 'duration', 'trial_type']), ignore_index=True)


        elif 'SacVE' in task: # Sac/Pur Loc: # Visual-Endogenous Sac/Pur Loc 
            
            events_glm = events[['onset','duration','trial_type', 'eyemov_vis_end']].copy(deep=True)
            events_glm.replace({'trial_type': {3: 'Fix', 1: 'Sac', 2: 'Pur'}},inplace=True)
            events_glm.replace({'eyemov_vis_end': {3: 'Fix', 1: 'SacExo', 2: 'SacEndo'}},inplace=True)
            events_glm.duration = tr_dur
            events_glm.onset = 0
            events_glm_groups = events_glm.groupby((events_glm.eyemov_vis_end!=events_glm.eyemov_vis_end.shift()).cumsum())

            new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type'])
            
            for idx, group in enumerate(events_glm_groups):
                onset = np.round(group[1]['onset'][group[1].index[0]],2)
                dur = np.round(sum(group[1]['duration']),2)
                vis_end = group[1]['eyemov_vis_end'][group[1].index[0]]

                new_events_glm = new_events_glm.append(pd.Series([onset,dur,vis_end],index =['onset', 'duration', 'trial_type']), ignore_index=True)

        elif 'PurVE' in task: # Sac/Pur Loc: # Visual-Endogenous Sac/Pur Loc 
            events_glm = events[['onset','duration','trial_type', 'eyemov_vis_end']].copy(deep=True)
            events_glm.replace({'trial_type': {3: 'Fix', 1: 'Sac', 2: 'Pur'}},inplace=True)
            events_glm.replace({'eyemov_vis_end': {3: 'Fix', 1: 'PurExo', 2: 'PurEndo'}},inplace=True)
            events_glm.duration = tr_dur
            events_glm.onset = 0
            events_glm_groups = events_glm.groupby((events_glm.eyemov_vis_end!=events_glm.eyemov_vis_end.shift()).cumsum())

            new_events_glm = pd.DataFrame([], columns=['onset', 'duration', 'trial_type', 'vis_end'])
            for idx, group in enumerate(events_glm_groups):
                onset = np.round(group[1]['onset'][group[1].index[0]],2)
                dur = np.round(sum(group[1]['duration']),2)
                ttype = group[1]['trial_type'][group[1].index[0]]
                vis_end = group[1]['eyemov_vis_end'][group[1].index[0]]

                new_events_glm = new_events_glm.append(pd.Series([onset,dur,vis_end],index =['onset', 'duration', 'trial_type']), ignore_index=True)

    for idx in new_events_glm.index:
        if idx==0:
            new_events_glm.at[idx, 'onset'] = 0
        else:
            new_events_glm.at[idx, 'onset'] = np.round(new_events_glm.at[idx-1, 'onset'] + new_events_glm.at[idx-1, 'duration'],2)

    return new_events_glm