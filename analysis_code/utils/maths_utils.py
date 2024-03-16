def weighted_regression(x_reg, y_reg, weight_reg, model):
    """
    Function to compute regression parameter weighted by a matrix (e.g. r2 value),
    where the regression model is y = 1/(cx) + d.

    Parameters
    ----------
    x_reg : array (1D)
        x values to regress
    y_reg : array
        y values to regress
    weight_reg : array (1D) 
        weight values (0 to 1) for weighted regression
    model : str
        Type of regression model, either 'pcm' for the original model or 'linear' for a linear model.

    Returns
    -------
    coef_reg : float or array
        regression coefficient(s)
    intercept_reg : float or str
        regression intercept or a string indicating no intercept (for linear model)
    """

    import numpy as np
    from scipy.optimize import curve_fit
    from sklearn import linear_model
    
    x_reg = np.array(x_reg)
    y_reg = np.array(y_reg)
    weight_reg = np.array(weight_reg)

    # Filter out NaN values
    x_reg_nan = x_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
    y_reg_nan = y_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
    weight_reg_nan = weight_reg[~np.isnan(weight_reg)]
    
    if model == 'pcm':
        # Define the model function
        def model_function(x, c, d):
            return 1 / (c * x) + d
       
        # Perform curve fitting
        params, _ = curve_fit(model_function, x_reg_nan, y_reg_nan, sigma=weight_reg_nan)
        
        # Extract parameters
        c, d = params
        
        return c, d
    
    elif model == 'linear':
        regr = linear_model.LinearRegression()
        
        # Filter out NaN values
        x_reg_nan = x_reg_nan.reshape(-1, 1)
        y_reg_nan = y_reg_nan.reshape(-1, 1)
        
        regr.fit(x_reg_nan, y_reg_nan, sample_weight=weight_reg_nan)
        coef_reg, intercept_reg = regr.coef_, regr.intercept_

        return coef_reg, intercept_reg
    
    else:
        raise ValueError("Invalid model type. Supported models are 'pcm' and 'linear'.")


# def r2_score_surf(bold_signal, model_prediction):
#     """
#     Compute r2 between bold signal and model. The gestion of nan values 
#     is down with created a non nan mask on the model prediction 

#     Parameters
#     ----------
#     bold_signal: bold signal in 2-dimensional np.array (time, vertex)
#     model_prediction: model prediction in 2-dimensional np.array (time, vertex)
    
#     Returns
#     -------
#     r2_scores: the R2 score for each vertex
#     """
#     import numpy as np
#     from sklearn.metrics import r2_score
    
#     # Check for NaN values in bold_signal
#     nan_mask = np.isnan(model_prediction).any(axis=0)
#     valid_vertices = ~nan_mask
    
#     # Set R2 scores for vertices with NaN values to NaN
#     r2_scores = np.full_like(nan_mask, np.nan, dtype=float)
    
#     # Compute R2 scores for vertices without NaN values
#     r2_scores[valid_vertices] = r2_score(bold_signal[:, valid_vertices], model_prediction[:, valid_vertices], multioutput='raw_values')
    
#     return r2_scores

def r2_score_surf(bold_signal, model_prediction):
    """
    Compute r2 between bold signal and model. The gestion of nan values 
    is down with created a non nan mask on the model prediction 

    Parameters
    ----------
    bold_signal: bold signal in 2-dimensional np.array (time, vertex)
    model_prediction: model prediction in 2-dimensional np.array (time, vertex)
    
    Returns
    -------
    r2_scores: the R2 score for each vertex
    """
    import numpy as np
    from sklearn.metrics import r2_score
    
    # Check for NaN values in both bold_signal and model_prediction
    nan_mask = np.isnan(model_prediction).any(axis=0) | np.isnan(bold_signal).any(axis=0)
    valid_vertices = ~nan_mask
    
    # Set R2 scores for vertices with NaN values to NaN
    r2_scores = np.full_like(nan_mask, np.nan, dtype=float)
    
    # Compute R2 scores for vertices without NaN values
    r2_scores[valid_vertices] = r2_score(bold_signal[:, valid_vertices], model_prediction[:, valid_vertices], multioutput='raw_values')
    
    return r2_scores


# def weighted_regression(x_reg, y_reg, weight_reg):
#     """
#     Function to compute regression parameter weighted by a matrix (e.g. r2 value).

#     Parameters
#     ----------
#     x_reg : array (1D)
#         x values to regress
#     y_reg : array
#         y values to regress
#     weight_reg : array (1D) 
#         weight values (0 to 1) for weighted regression

#     Returns
#     -------
#     coef_reg : array
#         regression coefficient
#     intercept_reg : str
#         regression intercept
#     """

#     from sklearn import linear_model
#     import numpy as np
    
#     regr = linear_model.LinearRegression()
    
#     x_reg = np.array(x_reg)
#     y_reg = np.array(y_reg)
#     weight_reg = np.array(weight_reg)
    
#     def m(x, w):
#         return np.sum(x * w) / np.sum(w)

#     def cov(x, y, w):
#         # see https://www2.microstrategy.com/producthelp/archive/10.8/FunctionsRef/Content/FuncRef/WeightedCov__weighted_covariance_.htm
#         return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

#     def weighted_corr(x, y, w):
#         # see https://www2.microstrategy.com/producthelp/10.4/FunctionsRef/Content/FuncRef/WeightedCorr__weighted_correlation_.htm
#         return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))

    
#     x_reg_nan = x_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
#     y_reg_nan = y_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
#     weight_reg_nan = weight_reg[~np.isnan(weight_reg)]

#     regr.fit(x_reg_nan.reshape(-1, 1), y_reg_nan.reshape(-1, 1),weight_reg_nan)
#     coef_reg, intercept_reg = regr.coef_, regr.intercept_

#     return coef_reg, intercept_reg



# def weighted_regression(x_reg, y_reg, weight_reg, model):
#     """
#     Function to compute regression parameter weighted by a matrix (e.g. r2 value),
#     where the regression model is y = 1/(cx) + d.

#     Parameters
#     ----------
#     x_reg : array (1D)
#         x values to regress
#     y_reg : array
#         y values to regress
#     weight_reg : array (1D) 
#         weight values (0 to 1) for weighted regression

#     Returns
#     -------
#     coef_reg : float
#         regression coefficient
#     intercept_reg : float
#         regression intercept
#     """

#     import numpy as np
#     from scipy.optimize import curve_fit

#     x_reg = np.array(x_reg)
#     y_reg = np.array(y_reg)
#     weight_reg = np.array(weight_reg)

#     # Filter out NaN values
#     x_reg_nan = x_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
#     y_reg_nan = y_reg[(~np.isnan(x_reg) & ~np.isnan(y_reg))]
#     weight_reg_nan = weight_reg[~np.isnan(weight_reg)]
    
#     if model == 'pcm':
#         # Define the model function
#         def model_function(x, c, d):
#             return 1 / (c * x) + d
#     if model == 'linear':
#         def model_function(x, a, b):
#             return a*x + b

   
#     # Perform curve fitting
#     # params, _ = curve_fit(model_function, x_reg_nan, y_reg_nan, p0=initial_guess, sigma=1/np.sqrt(weight_reg_nan))
#     params, _ = curve_fit(model_function, x_reg_nan, y_reg_nan, sigma=weight_reg_nan)

#     # Extract parameters
#     c, d = params

#     return c, d
