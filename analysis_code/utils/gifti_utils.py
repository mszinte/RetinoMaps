def make_giti_image(source_img, data_to_write) : 
    """
    make a gifti image with data

    Parameters
    ----------
    source_img : image from with new data derives
    data_to_write : data you want to write on your image 
    
    Returns
    -------
    The gifti image of your data with the same strucure the source image 
    
    """
    # import 
    import nibabel as nb
    
    # gete source img informations
    header = source_img.header
    meta = source_img.meta
    file_map = source_img.file_map
    labeltable = source_img.labeltable
    
    # initialise the final image 
    final_img = nb.gifti.GiftiImage(header = header, meta = meta,file_map = file_map,labeltable = labeltable)
    
    # fill final image  
    for i, data in enumerate(data_to_write):
        darray = nb.gifti.GiftiDataArray(data,datatype = 'NIFTI_TYPE_FLOAT32',intent=source_img.darrays[i].intent,meta=source_img.darrays[i].meta,coordsys=source_img.darrays[i].coordsys)
        final_img.add_gifti_data_array(darray)
    
    return final_img