def fit2deriv(fit_array, data_array, pred_array):
    """
    Compute pRF derivatives out of fitting output and predictions

    Parameters
    ----------
    fit_array: fit parameters 2D array
    data_array: data timeseries 2D array
    pred_array: prediction timeseries 2D array
    
    Returns
    -------
    deriv_array: 2D array with pRF derivatives

    stucture output:
    columns: 1->size of input
    rows: derivatives parameters
    
    """

    # Imports
    # -------
    # General imports
    import os
    import numpy as np
    import ipdb
    from sklearn.metrics import r2_score
    deb = ipdb.set_trace

    # Compute derived measures from prfs/pmfs
    # ---------------------------------------
    # get data index
    x_idx, y_idx, sigma_idx, beta_idx, baseline_idx, rsq_idx = 0, 1, 2, 3, 4, 5

    # change to nan empty voxels
    fit_array[fit_array[...,rsq_idx]==0] = np.nan
    
    # r-square
    rsq = fit_array[...,rsq_idx]

    # eccentricity
    ecc = np.nan_to_num(np.sqrt(fit_array[...,x_idx]**2 + fit_array[...,y_idx]**2))

    # polar angle
    complex_polar = fit_array[...,x_idx] + 1j * fit_array[...,y_idx]
    normed_polar = complex_polar / np.abs(complex_polar)
    polar_real = np.real(normed_polar)
    polar_imag = np.imag(normed_polar)
    
    # size
    size_ = fit_array[...,sigma_idx].astype(np.float64)
    size_[size_<1e-4] = 1e-4

    # amplitude
    amp = fit_array[...,beta_idx]
    
    # baseline
    baseline = fit_array[...,baseline_idx]

    # x
    x = fit_array[...,x_idx]

    # y
    y = fit_array[...,y_idx]
    
    # r-square between data and model prediction
    num_elmt = data_array.shape[0]*data_array.shape[1]*data_array.shape[2]
    data_array_flat = data_array.reshape((num_elmt,data_array.shape[-1]))
    pred_array_flat = pred_array.reshape((num_elmt,pred_array.shape[-1]))
    pred_array_flat[np.isnan(pred_array_flat)]=0
    rsq_pred = np.power(r2_score(data_array_flat.T, pred_array_flat.T, multioutput='raw_values'),2)
    rsq_pred = rsq_pred.reshape(data_array.shape[:-1])
    
    # Save results
    if np.ndim(fit_array) == 4:
        deriv_array = np.zeros((fit_array.shape[0],fit_array.shape[1],fit_array.shape[2],10))*np.nan
    elif np.ndim(fit_array) == 2:
        deriv_array = np.zeros((fit_array.shape[0],10))*np.nan

    deriv_array[...,0] = rsq
    deriv_array[...,1] = rsq_pred
    deriv_array[...,2] = ecc
    deriv_array[...,3] = polar_real
    deriv_array[...,4] = polar_imag
    deriv_array[...,5] = size_
    deriv_array[...,6] = amp
    deriv_array[...,7] = baseline
    deriv_array[...,8] = x
    deriv_array[...,9] = y

    deriv_array = deriv_array.astype(np.float32)

    return deriv_array

