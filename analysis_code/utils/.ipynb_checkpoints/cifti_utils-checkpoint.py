def from_170k_to_59k(img, data, return_concat_hemis=False, return_59k_mask=False, return_170k_mask=False):
    """
    Transform 170k data into 59k data by retaining only the left and right cortex and medial wall vertices.

    Parameters
    ----------
    img : the cifti image of your 170k 
    data : your 170k data
    return_concat_hemis : if False, will return two arrays for the two hemispheres 
                          (2 arrays with dimensions (time x 59k vertices))
                          
                          if True, will return the concatenation of hemi-L and hemi-R 
                          (1 array with dimensions (time x 118 vertices))
                          
    return_59k_mask : if True, will return a mask where True corresponds to cortex vertices and
                      False to medial wall vertices.
    
    return_170k_mask : if True, will return a mask where True corresponds to cortex vertices 
                       (independent of left or right) and False corresponds to non-cortex vertices.
                       
                       WARNING: Applying 170k_mask will yield a 54k surface (without medial wall vertices) 
                       and not a 59k surface (i.e., cortex + medial wall vertices)

    Returns
    -------
    If only one result is requested, it is returned directly. 
    If multiple results are requested, they are returned as a tuple.

    """
    import numpy as np 
    
    brain_model = img.header.get_axis(1)
    surf_name_L = 'CIFTI_STRUCTURE_CORTEX_LEFT'
    surf_name_R = 'CIFTI_STRUCTURE_CORTEX_RIGHT'

    for structure_name, data_indices, model in brain_model.iter_structures(): 
        


        if structure_name == surf_name_L: 
            data_L = data.T[data_indices]
            vtx_indices_L = model.vertex 

            # include inter hemi vertex
            surf_data_L = np.zeros((vtx_indices_L.max() + 1,) + data_L.shape[1:], dtype=data_L.dtype)
            surf_data_L[vtx_indices_L] = data_L    
            surf_data_L = surf_data_L.T 

            # Know where are inter hemi vertex
            mask_L = np.any(surf_data_L != 0, axis=0)






        elif structure_name == surf_name_R: 
            data_R = data.T[data_indices]
            vtx_indices_R = model.vertex 
        
            # include inter hemi vertex
            surf_data_R = np.zeros((vtx_indices_R.max() + 1,) + data_R.shape[1:], dtype=data_R.dtype)
            surf_data_R[vtx_indices_R] = data_R
            surf_data_R = surf_data_R.T 
            
            # Know where are inter hemi vertex
            mask_R = np.any(surf_data_R != 0, axis=0)
            
            
            

    brain_mask_59k = np.concatenate((mask_L, mask_R))
    brain_mask_170k = np.concatenate((brain_mask_59k, np.zeros(brain_model.size - len(brain_mask_59k), dtype=bool)))
    
    result = []

    if return_concat_hemis:
        result.append(np.concatenate((surf_data_L, surf_data_R),axis=1))

    if return_59k_mask:
        result.append(brain_mask_59k)

    if return_170k_mask:
        result.append(brain_mask_170k)

    if len(result) == 1:  
        return result[0]  
    elif result:
        return tuple(result)
    else:
        return surf_data_L, surf_data_R


def from_59k_to_170k(data_59k, data_170k_orig, brain_mask_59k, brain_mask_170k):
    """
    Transform 59k data into 170k data by filling non-59k vertices with numpy.nan.

    Parameters
    ----------
    
    data_59k : The 59k data you want to transform into 170k.
    data_170k_orig : The original 170k data from which your 59k data originated.
    brain_mask_59k : 59k brain mask output from from_170k_to_59k.
    brain_mask_170k : 170k brain mask output from from_170k_to_59k.
    
    Returns
    -------
    The transformed data in 170k format with non-59k vertices filled with numpy.nan.
    """
    
    # mask 59k data to optain only cortex vertex (and not medial wall vertices ) 
    data_59k = data_59k[:,brain_mask_59k]
    
    # create an 170k full nan array
    data_170k_final = np.full(data_170k_orig.shape,np.nan)
    
    # fill the 170k array with the cortex data
    data_170k_final[:,brain_mask_170k] = data_59k
    
    return data_170k_final

