#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import SimpleITK as sitk
from scipy.fft import fftfreq, fft
from sklearn.decomposition import PCA
import matplotlib.backends.backend_pdf
import math


plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 
plt.rc('axes', labelsize=15) 

# # Take the arguments from bash
input_epi=os.path.abspath(sys.argv[1])
input_roi=os.path.abspath(sys.argv[2])
desired_slice=os.path.abspath(sys.argv[3])
TR=os.path.abspath(sys.argv[4])


# # Functions

def extract_an_roi(slices, PE_matrix_size, FE_matrix_size):
    #this function creates an one-slice 10x10 roi in the center of the EPI image, if no manually defined ROI is specified
    
    #create empty matrix with the same dimensions as EPI image
    roi_matrix = np.zeros((slices, PE_matrix_size, FE_matrix_size))
    
    #decide where the center of the ROI will be within the matrix
    slice_of_roi = int(slices/2)
    center_of_roi_PE = int(PE_matrix_size/2)
    center_of_roi_FE = int(FE_matrix_size/2)
    
    #populate the 10x10 voxels around the center with ones to define the ROI
    for i in range(0,10):
        for j in range(-5,5):
            roi_matrix[slice_of_roi, center_of_roi_PE + i, center_of_roi_FE + j] = 1

    return roi_matrix


def extract_residuals(phantom_epi, roi, time, slices, PE_matrix_size, FE_matrix_size, num_rep):
    
    if roi is None:
        #fit the second order polynomial to the data
        phantom_epi_flat = phantom_epi.transpose(1,2,3,0).reshape(-1,phantom_epi.shape[0])
        model = np.polyfit(time, phantom_epi_flat.T, 2)

        #generate the predicted polynomial curve based on the fitted model(for each voxel)
        tot_vox = slices*PE_matrix_size*FE_matrix_size
        predicted = np.zeros([tot_vox, len(time)])
        for i in range(0,tot_vox):
            predicted[i,:] = np.polyval(model[:,i],time)

        #detrend the data by removing second order polynomial
        phantom_epi_flat_detrended = phantom_epi_flat - predicted
        
        phantom_epi_roi_mean = None #this variable doesn't apply if no ROI, but we want to return it if there is an ROI
        
    else:    
        #extract only the voxels in the ROI
        l=[]
        for i in range(phantom_epi.shape[0]): #iterate over timepoints
            l.append(phantom_epi[i,roi.astype(bool)])
        phantom_epi_roi = np.array(l)

        ###############################################################################################
        # Obtain residuals within each roi

        #first fit the second order polynomial to the data
        phantom_epi_roi_mean = np.mean(phantom_epi_roi, axis = 1) #gives mean across voxels, at each timepoint
        model_roi = np.polyfit(time, phantom_epi_roi_mean, 2)
        predicted = np.polyval(model_roi,time)
        phantom_epi_flat_detrended = phantom_epi_roi_mean - predicted

    
    return phantom_epi_flat_detrended, predicted, phantom_epi_roi_mean


def voxelwise_wholephantom_analysis(phantom_epi, roi_to_plot, time, slice_num, slices, PE_matrix_size, FE_matrix_size, num_rep):
    #calculate signal image (average across the timepoints, voxel-wise)
    signal_image = np.mean(phantom_epi, axis = 0)
    
    ###############################################################################################
    #calculate the temporal fluctuation noise image (std of residuals after detrending timeseries with 2nd order polynomial)
    phantom_epi_flat_detrended, predicted, phantom_epi_roi_mean = extract_residuals(phantom_epi, None, time, slices,
                                                                               PE_matrix_size, FE_matrix_size, num_rep)
    phantom_epi_detrended = phantom_epi_flat_detrended.reshape(slices, PE_matrix_size, FE_matrix_size,len(time))
    temp_fluc_noise_image = np.std(phantom_epi_detrended, axis = 3)
    
    ############################################################################################
    #compute signal to fluctuation noise ratio (SFNR)
    sfnr_image = signal_image/temp_fluc_noise_image
    
    ###########################################################################################
    #compute static spatial noise image
    sumeven = 0
    sumodd = 0
    for i in range(0,len(time)-1,2):
        sumeven = sumeven + phantom_epi[i,:,:,:]
        sumodd = sumodd + phantom_epi[i+1,:,:,:]
    static_spatial_noise_im = sumodd - sumeven

    ################################### PLOT ##########################################################
    #find which slice the ROI is drawn in and plot that slice in the last subplot
    for i in range(0,slices):
        if np.sum(roi_to_plot[i, :,:]) > 0:
            slice_to_plot = i
            
    #if a desired slice number is given for the other images, plot that slice. Otherwise, plot same slice as ROI.
    if slice_num is None:
        slice_num = slice_to_plot
        
    fig, axs = plt.subplots(1, 5, figsize = (15,6), sharey = True)
    fig.suptitle('Whole Phantom Signal Analysis', y = 0.7, fontsize = 15) 
    
    axs[0].set_title('Mean Signal Image', fontsize = 15)
    s = axs[0].imshow(signal_image[slice_num,:,:], origin = 'lower')
    cbar = plt.colorbar(s, ax = axs[0], orientation = 'horizontal')
    
    axs[1].set_title('Temporal Fluctuation Noise Image', fontsize = 15)
    t = axs[1].imshow(temp_fluc_noise_image[slice_num,:,:], origin = 'lower')
    cbar = plt.colorbar(t, ax = axs[1], orientation = 'horizontal')
    
    axs[2].set_title('SFNR Image', fontsize = 15)
    sf = axs[2].imshow(sfnr_image[slice_num,:,:], origin = 'lower')
    cbar = plt.colorbar(sf, ax = axs[2], orientation = 'horizontal')
    
    axs[3].set_title('Static Spatial Noise Image', fontsize = 15)
    g = axs[3].imshow(static_spatial_noise_im[slice_num,:,:], origin = 'lower')
    cbar = plt.colorbar(g, ax = axs[3], orientation = 'horizontal')      
    
    axs[4].set_title('Location of ROI', fontsize = 15)
    roi_ontop_of_signal = 500*roi_to_plot + signal_image
    r = axs[4].imshow(roi_ontop_of_signal[slice_to_plot,:,:], origin = 'lower')
    cbar = plt.colorbar(r, ax = axs[4], orientation = 'horizontal')
    
    fig.tight_layout()
    
    
    return fig, signal_image, sfnr_image, static_spatial_noise_im, phantom_epi_flat_detrended


def roi_residuals_analysis(phantom_epi, roi, time, signal_image, sfnr_image, static_spatial_noise_im, TR, num_rep):
    residuals_in_roi, predicted_roi, phantom_epi_roi_mean = extract_residuals(phantom_epi, roi, time, 0,0,0, num_rep)
    
    ###################################### CALC METRICS WITHIN ROI #######################################3
    signal_summary_value = np.mean(signal_image[roi.astype(bool)])
    sfnr_summary_value = np.mean(sfnr_image[roi.astype(bool)])
    intrinsic_noise = np.var(static_spatial_noise_im[roi.astype(bool)])
    snr = signal_summary_value/np.sqrt(intrinsic_noise/355)
    
    percent_fluc = 100*np.std(residuals_in_roi)/np.mean(phantom_epi_roi_mean)
    diff = max(predicted_roi)- min(predicted_roi)
    drift = 100*diff/np.mean(phantom_epi_roi_mean) #not sure if this is right
    drift_alt = 100*diff/signal_summary_value
    
    ####################################### Fourier Analysis ####################################################
    N = phantom_epi.shape[0] #length

    yf = scipy.fft.fft(residuals_in_roi)
    yf_half = np.abs(yf[1:(N+1)//2])
    xf = scipy.fft.fftfreq(N, TR)[1:(N+1)//2]
    location_of_peak = np.argwhere(yf_half > max(yf_half) - 0.1)
    value_of_peak = xf[location_of_peak]
    
    ########################################difference with Gaussian (qq correlation)#############################
    (osm, osr),(slope, intercept, r) = stats.probplot(residuals_in_roi)
    
    ##################################### PLOT ##################################################################
    fig0, axs = plt.subplots(1, 5, figsize = (20,4))
    fig0.suptitle('Analysis of Residuals within an ROI', y = 1, fontsize = 15) 
    axs[0].set_title('Polynomial Fit (ROI average)', fontsize=15)
    axs[0].plot(time, phantom_epi_roi_mean)
    axs[0].plot(time, predicted_roi)
    axs[0].set_xlabel('Time (s)')
    axs[0].set_ylabel('Signal Intensity')
    axs[1].set_title('Residuals', fontsize=15)
    axs[1].plot(time, residuals_in_roi)
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('Signal Intensity')
    axs[2].set_title('FFT Spectrum', fontsize=15)
    axs[2].plot(xf, np.abs(yf[1:(N+1)//2]))
    axs[2].set_xlabel('Frequency (Hz)')
    axs[2].set_ylabel('FFT Magnitude')
    axs[3].set_title('Histogram of Residuals', fontsize=15)
    axs[3].hist(residuals_in_roi, bins = 30)
    axs[3].set_xlabel('Residual intensity')
    axs[3].set_ylabel('Frequency')
    stats.probplot(residuals_in_roi, plot=axs[4])
    fig0.text(0,-0.1,'QQ Correlation of residuals: '+str(r), fontsize = 15)
    fig0.text(0,-0.2, "The SFNR summary value (inside roi) is: " + str(sfnr_summary_value), fontsize = 15)
    fig0.text(0,-0.3, "The SNR summary value (inside roi) is: " + str(snr), fontsize = 15)
    fig0.text(0,-0.4, "The percent fluctuation (inside roi) is: " + str(percent_fluc), fontsize = 15)
    fig0.text(0,-0.5, "The drift (inside roi) is: " + str(drift_alt), fontsize = 15)
    fig0.text(0,-0.6, "The strongest frequency in the FFT spectrum is: " + str(value_of_peak[0]) + ' Hz', fontsize = 15)
    fig0.tight_layout()
    
    return fig0, sfnr_summary_value, snr, percent_fluc, drift_alt, value_of_peak


def pca_analysis(agar_epi_flat_detrended, time, slices, PE_matrix_size, FE_matrix_size, num_rep, TR):
    
    pca = PCA()
    pc_space = pca.fit_transform(agar_epi_flat_detrended)
    pc_time = pca.components_
    pc_exp_var = pca.explained_variance_
    
    num_components = len(pc_exp_var)
    xf = scipy.fft.fftfreq(num_rep, TR)[1:(num_rep+1)//2]
    
    fig1, axs = plt.subplots(2, 6, figsize = (26,8))
    fig1.suptitle('Temporal PCA Across All (Detrended) Voxels', fontsize = 15)
    
    for i in range(0,6):
        #plot timecourses of each pc_time_sliceselectdir in first row
        axs[0,i].plot(time, pc_time[i,:])
        axs[0,i].set_title('Component ' + str(i) + (' (') + str(round(pc_exp_var[i],2)) + ('%)'), fontsize = 15)
        axs[0,i].set_xlabel('Time (s)', fontsize = 15)
        axs[0,0].set_ylabel('Amplitude (a.u.)', fontsize = 15)
        
        #plot fourier transform of each pc_time in second row
        axs[1,i].plot(xf, np.abs(scipy.fft.fft(pc_time[i,:])[1:(num_rep+1)//2]))
        axs[1,i].set_xlabel('Frequency (Hz)', fontsize = 15)
        axs[1,0].set_ylabel('Amplitude (a.u.)', fontsize = 15)
    #plt.show()
    
    #################################### Plot the spatial pattern of the first 3 components ############################
    
    #reshape the 1d array into the original image dimensions
    pc_space_im = pc_space.reshape(slices,PE_matrix_size, FE_matrix_size,num_components)
    
    #decide how many subplots are necessary (based on the number of slices)
    root = np.sqrt(slices)+1
    subplot_dim1 = math.ceil(root)
    subplot_dim2 = math.ceil(slices/subplot_dim1)
    fig_dim1 = 2*subplot_dim1
    fig_dim2 = 2*subplot_dim2
    
    fig2, axs2 = plt.subplots(subplot_dim2, subplot_dim1, figsize = (fig_dim1,fig_dim2), sharex = True, sharey = True)
    fig2.suptitle('Spatial Pattern of PC 0', fontsize = 15)
    fig3, axs3 = plt.subplots(subplot_dim2, subplot_dim1, figsize = (fig_dim1,fig_dim2), sharex = True, sharey = True)
    fig3.suptitle('Spatial Pattern of PC 1', fontsize = 15)
    fig4, axs4 = plt.subplots(subplot_dim2, subplot_dim1, figsize = (fig_dim1,fig_dim2), sharex = True, sharey = True)
    fig4.suptitle('Spatial Pattern of PC 2', fontsize = 15)
    fig5, axs5 = plt.subplots(subplot_dim2, subplot_dim1, figsize = (fig_dim1,fig_dim2), sharex = True, sharey = True)
    fig5.suptitle('Spatial Pattern of PC 3', fontsize = 15)
    fig6, axs6 = plt.subplots(subplot_dim2, subplot_dim1, figsize = (fig_dim1,fig_dim2), sharex = True, sharey = True)
    fig6.suptitle('Spatial Pattern of PC 4', fontsize = 15)
    fig7, axs7 = plt.subplots(subplot_dim2, subplot_dim1, figsize = (fig_dim1,fig_dim2), sharex = True, sharey = True)
    fig7.suptitle('Spatial Pattern of PC 5', fontsize = 15)
    
    
    slice_num = 0
    max_val = 5
    min_val = -5
    for j in range(0,subplot_dim2):
        for k in range(0,subplot_dim1):
            
            if slice_num >= slices:
                break
                
            im2 = axs2[j,k].imshow(pc_space_im[slice_num,:,:,0], vmax = max_val, vmin = min_val)
            axs2[j,k].set_title('Slice #' + str(slice_num), fontsize = 15)
            
            im3 = axs3[j,k].imshow(pc_space_im[slice_num,:,:,1], vmax = max_val, vmin = min_val)
            axs3[j,k].set_title('Slice #' + str(slice_num), fontsize = 15)
            
            im4 = axs4[j,k].imshow(pc_space_im[slice_num,:,:,2], vmax = max_val, vmin = min_val)
            axs4[j,k].set_title('Slice #' + str(slice_num), fontsize = 15)
            
            im5 = axs5[j,k].imshow(pc_space_im[slice_num,:,:,3], vmax = max_val, vmin = min_val)
            axs5[j,k].set_title('Slice #' + str(slice_num), fontsize = 15)
            
            im6 = axs6[j,k].imshow(pc_space_im[slice_num,:,:,4], vmax = max_val, vmin = min_val)
            axs6[j,k].set_title('Slice #' + str(slice_num), fontsize = 15)
            
            im7 = axs7[j,k].imshow(pc_space_im[slice_num,:,:,5], vmax = max_val, vmin = min_val)
            axs7[j,k].set_title('Slice #' + str(slice_num), fontsize = 15)
            
            slice_num = slice_num + 1
    fig2.colorbar(im2, ax = axs2, orientation = 'horizontal')
    fig3.colorbar(im3, ax = axs3, orientation = 'horizontal')
    fig4.colorbar(im4, ax = axs4, orientation = 'horizontal')
    fig5.colorbar(im5, ax = axs5, orientation = 'horizontal')
    fig6.colorbar(im6, ax = axs6, orientation = 'horizontal')
    fig7.colorbar(im7, ax = axs7, orientation = 'horizontal')
    plt.show
            
    return fig1, fig2, fig3, fig4, fig5, fig6, fig7


def full_analysis(phantom_epi_filepath, roi_filepath, output_filepath, slice_to_plot, TR):
    
    #load the images, then convert them to arrays
    agar_epi_image = sitk.ReadImage(phantom_epi_filepath)
    agar_epi_full = sitk.GetArrayFromImage(agar_epi_image)

    #extract dimensions
    num_rep = agar_epi_full.shape[0]
    slices = agar_epi_full.shape[1]
    PE_matrix_size = agar_epi_full.shape[2]
    FE_matrix_size = agar_epi_full.shape[3]
    
    #remove dummy scans from EPI
    num_dummy_scans = int(round(num_rep*0.013))
    num_rep_no_dummy = num_rep - num_dummy_scans
    agar_epi = agar_epi_full[num_dummy_scans:num_rep,:,:,:]
    
    #define a time array that corresponds with EPI (without dummy scans)
    time_arr = np.linspace(0, num_rep_no_dummy-1, num_rep_no_dummy)
    
    #if there is no manually drawn roi provided, extract a 10x10, one-slice roi from the middle slice
    if roi_filepath is None:
        roi = extract_an_roi(slices, PE_matrix_size, FE_matrix_size)
    else:
        #Deals with MINC ROIS
        roi_image = sitk.ReadImage(roi_filepath)
        roi = sitk.GetArrayFromImage(roi_image)
        roi = roi.swapaxes(0,1)
        
    #perform whole phantom analysis
    [figure_voxelwise_wholephantom, signal_image, sfnr_image, 
     static_spatial_noise_im, agar_epi_flat_detrended] = voxelwise_wholephantom_analysis(agar_epi, roi, time_arr, 
                                                                                         slice_to_plot, slices, PE_matrix_size, 
                                                                                         FE_matrix_size, num_rep)
        
    #perform within roi analysis
    [figure_roi_analysis, sfnr_summary_value, snr, percent_fluc,
     drift_alt, value_of_peak] = roi_residuals_analysis(agar_epi, roi, time_arr, signal_image, sfnr_image, 
                                                      static_spatial_noise_im, TR, num_rep)
    #also perform PCA
    [figure_pca_time, figure_pca_space0, 
     figure_pca_space1, figure_pca_space2,
    figure_pca_space3, figure_pca_space4,figure_pca_space5] = pca_analysis(agar_epi_flat_detrended, time_arr, slices, 
                                                                           PE_matrix_size, FE_matrix_size, num_rep, TR)
    
    #export all figures to pdf
    pdf_multiplot = matplotlib.backends.backend_pdf.PdfPages(output_filepath)
    pdf_multiplot.savefig(figure_voxelwise_wholephantom)
    pdf_multiplot.savefig(figure_roi_analysis, bbox_inches="tight")
    pdf_multiplot.savefig(figure_pca_time)
    pdf_multiplot.savefig(figure_pca_space0)
    pdf_multiplot.savefig(figure_pca_space1)
    pdf_multiplot.savefig(figure_pca_space2)
    pdf_multiplot.savefig(figure_pca_space3)
    pdf_multiplot.savefig(figure_pca_space4)
    pdf_multiplot.savefig(figure_pca_space5)
    pdf_multiplot.close()


# # Call the function

full_analysis(input_epi, input_roi, output_path, desired_slice, TR)

