#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import SimpleITK as sitk
from scipy.fft import fftfreq, fft
from sklearn.decomposition import PCA
import matplotlib.backends.backend_pdf
import math
import ast
from skimage.filters import threshold_li
import csv
import pandas as pd


plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)

# # Take the arguments from bash
input_epi=os.path.abspath(sys.argv[1])
output_path=os.path.abspath(sys.argv[2])
TR=ast.literal_eval(sys.argv[3])
input_roi=sys.argv[4]
desired_slice=ast.literal_eval(sys.argv[5])
weisskoff_max_roi_width=ast.literal_eval(sys.argv[6])
longitudinal_csv=os.path.abspath(sys.argv[7])

if input_roi == 'None':
    input_roi = None
if 'None' in longitudinal_csv:
    longitudinal_csv = None
    
# # Functions

def extract_an_roi(slices, PE_matrix_size, FE_matrix_size, width):
    #this function creates an one-slice 10x10 roi in the center of the EPI image, if no manually defined ROI is specified

    #create empty matrix with the same dimensions as EPI image
    roi_matrix = np.zeros((slices, PE_matrix_size, FE_matrix_size))

    #decide where the center of the ROI will be within the matrix
    slice_of_roi = int(slices/2)
    center_of_roi_PE = int(PE_matrix_size/2)
    center_of_roi_FE = int(FE_matrix_size/2)

    #populate the 10x10 voxels around the center with ones to define the ROI
    for i in range(math.ceil(-width/2),math.ceil(width/2)):
        for j in range(math.ceil(-width/2),math.ceil(width/2)):
            roi_matrix[slice_of_roi, center_of_roi_PE + i, center_of_roi_FE + j] = 1

    return roi_matrix

def extract_an_roi_background(slices, PE_matrix_size, FE_matrix_size, width):
    #this function creates an one-slice widthxwidth roi in the background of the EPI (for weisskoff analysis)

    #create empty matrix with the same dimensions as EPI image
    roi_matrix = np.zeros((slices, PE_matrix_size, FE_matrix_size))

    #decide which corner to use as roi location
    slice_of_roi = int(slices/2)
    corner_of_roi_PE = 0
    corner_of_roi_FE = 0

    #populate the width x width voxels around the center with ones to define the ROI
    for i in range(0,width):
        for j in range(0,width):
            roi_matrix[slice_of_roi, corner_of_roi_PE + i, corner_of_roi_FE + j] = 1

    return roi_matrix


def extract_residuals(phantom_epi, roi, time, slices, PE_matrix_size, FE_matrix_size):

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

        phantom_epi_mean_timeseries_in_roi = None #this variable doesn't apply if no ROI, but we want to return it if there is
        phantom_epi_spatial_std_acrosstime_in_roi = None

    else:
        #extract only the voxels in the ROI
        l=[]
        for i in range(phantom_epi.shape[0]): #iterate over timepoints
            l.append(phantom_epi[i,roi.astype(bool)])
        phantom_epi_roi = np.array(l)

        #obtain the mean timeseries within roi, prior to detrending
        phantom_epi_mean_timeseries_in_roi = np.mean(phantom_epi_roi, axis = 1) #gives mean across voxels, at each timepoint

        #obtain the spatial std within roi, at each timepoint
        phantom_epi_spatial_std_acrosstime_in_roi = np.std(phantom_epi_roi, axis = 1)

        ###############################################################################################
        # Obtain residuals within each roi

        #first fit the second order polynomial to the data
        model_roi = np.polyfit(time, phantom_epi_mean_timeseries_in_roi, 2)
        predicted = np.polyval(model_roi,time)
        phantom_epi_flat_detrended = phantom_epi_mean_timeseries_in_roi - predicted


    return phantom_epi_flat_detrended, predicted, phantom_epi_mean_timeseries_in_roi, phantom_epi_spatial_std_acrosstime_in_roi


def voxelwise_wholephantom_analysis(phantom_epi, roi_to_plot, time, slice_num, slices, PE_matrix_size, FE_matrix_size, num_rep):
    #calculate signal image (average across the timepoints, voxel-wise)
    signal_image = np.mean(phantom_epi, axis = 0)

    ###############################################################################################
    #calculate the temporal fluctuation noise image (std of residuals after detrending timeseries with 2nd order polynomial)
    phantom_epi_flat_detrended, predicted, phantom_epi_roi_mean, tmp = extract_residuals(phantom_epi, None, time, slices,
                                                                               PE_matrix_size, FE_matrix_size)
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

    fig, axs = plt.subplots(1, 5, figsize = (15,7), sharey = True)
    fig.suptitle('Standard Stability Metrics - voxelwise', y = 0.7, fontsize = 15)

    axs[0].set_title('Mean Signal Map', fontsize = 15)
    s = axs[0].imshow(signal_image[slice_num,:,:], origin = 'lower', vmax = 80, vmin=0)
    cbar = plt.colorbar(s, ax = axs[0], orientation = 'horizontal')
    cbar.set_label('Intensity (a.u.)')
    axs[0].axis('off')

    axs[1].set_title('Standard Deviation Map', fontsize = 15)
    t = axs[1].imshow(temp_fluc_noise_image[slice_num,:,:], origin = 'lower', vmax = 1.5, vmin = 0)
    cbar = plt.colorbar(t, ax = axs[1], orientation = 'horizontal')
    cbar.set_label('Intensity (a.u.)')
    axs[1].axis('off')

    axs[2].set_title('tSNR Map', fontsize = 15)
    sf = axs[2].imshow(sfnr_image[slice_num,:,:], origin = 'lower', vmax = 115, vmin = 0)
    cbar = plt.colorbar(sf, ax = axs[2], orientation = 'horizontal')
    cbar.set_label('Intensity (a.u.)')
    axs[2].axis('off')

    axs[3].set_title('Static Spatial Noise Map', fontsize = 15)
    g = axs[3].imshow(static_spatial_noise_im[slice_num,:,:], origin = 'lower', vmax = 70, vmin = -70)
    cbar = plt.colorbar(g, ax = axs[3], orientation = 'horizontal')
    cbar.set_label('Intensity (a.u.)')
    axs[3].axis('off')

    axs[4].set_title('Location of ROI', fontsize = 15)
    roi_ontop_of_signal = 500*roi_to_plot + signal_image
    r = axs[4].imshow(roi_ontop_of_signal[slice_to_plot,:,:], origin = 'lower')
    cbar = plt.colorbar(r, ax = axs[4], orientation = 'horizontal')
    axs[4].axis('off')

    fig.tight_layout()


    return fig, signal_image, temp_fluc_noise_image, sfnr_image, static_spatial_noise_im, phantom_epi_flat_detrended, slice_num


def roi_residuals_analysis(phantom_epi, roi, time, signal_image, temp_fluc_noise_image, sfnr_image, static_spatial_noise_im, TR, num_rep_no_dummy):
    residuals_in_roi, predicted_roi, phantom_epi_roi_mean,tmp = extract_residuals(phantom_epi, roi, time, 0,0,0)

    ###################################### CALC METRICS WITHIN ROI #######################################3
    signal_summary_value = np.mean(signal_image[roi.astype(bool)])
    std_summary_value = np.mean(temp_fluc_noise_image[roi.astype(bool)])
    sfnr_summary_value = np.mean(sfnr_image[roi.astype(bool)])
    intrinsic_noise = np.var(static_spatial_noise_im[roi.astype(bool)])
    snr = signal_summary_value/np.sqrt(intrinsic_noise/num_rep_no_dummy)

    percent_fluc = 100*np.std(residuals_in_roi)/np.mean(phantom_epi_roi_mean)
    diff = max(predicted_roi)- min(predicted_roi)
    driftfit = 100*diff/np.mean(phantom_epi_roi_mean)
    drift = 100*(max(phantom_epi_roi_mean) - min(phantom_epi_roi_mean))/np.mean(phantom_epi_roi_mean)

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
    fig0, axs = plt.subplots(1, 5, figsize = (15,4))
    fig0.suptitle('Standard Stability Metrics - Average within ROI', y = 1, fontsize = 15)
    axs[0].set_title('Polynomial Fit (ROI average)', fontsize = 15)
    axs[0].plot(time, phantom_epi_roi_mean)
    axs[0].plot(time, predicted_roi)
    axs[0].set_xlabel('Time (s)', fontsize = 15)
    axs[0].set_ylabel('Signal Intensity', fontsize = 15)
    axs[1].set_title('Residuals', fontsize = 15)
    axs[1].plot(time, residuals_in_roi)
    axs[1].set_xlabel('Time (s)', fontsize = 15)
    axs[1].set_ylabel('Signal Intensity', fontsize = 15)
    axs[2].set_title('FFT Spectrum', fontsize = 15)
    axs[2].plot(xf, np.abs(yf[1:(N+1)//2]))
    axs[2].set_xlabel('Frequency (Hz)', fontsize = 15)
    axs[2].set_ylabel('FFT Magnitude', fontsize = 15)
    axs[3].set_title('Histogram of Residuals', fontsize = 15)
    axs[3].hist(residuals_in_roi, bins = 30)
    axs[3].set_xlabel('Residual intensity', fontsize = 15)
    axs[3].set_ylabel('Frequency', fontsize = 15)
    stats.probplot(residuals_in_roi, plot=axs[4])
    fig0.text(0,-0.1, "The mean signal (inside roi) is: " + str(signal_summary_value), fontsize = 15)
    fig0.text(0,-0.2, "The temporal std (inside roi) is: " + str(std_summary_value), fontsize = 15)
    fig0.text(0,-0.3, "The tSNR summary value (inside roi) is: " + str(sfnr_summary_value), fontsize = 15)
    fig0.text(0,-0.4, "The SNR summary value (inside roi) is: " + str(snr), fontsize = 15)
    fig0.text(0,-0.5, "The drift (inside roi) is: " + str(drift), fontsize = 15)
    fig0.text(0,-0.6, "The drift of the fit (inside roi) is: " + str(driftfit), fontsize = 15)
    fig0.text(0,-0.7, "The percent fluctuation (inside roi) is: " + str(percent_fluc), fontsize = 15)
    fig0.text(0,-0.8, "The strongest frequency in the FFT spectrum is: " + str(value_of_peak[0][0]) + ' Hz', fontsize = 15)
    fig0.text(0,-0.9,'QQ Correlation of residuals: '+str(r), fontsize = 15)
    fig0.tight_layout()

    return fig0, signal_summary_value, std_summary_value, sfnr_summary_value, snr, drift, driftfit, percent_fluc, value_of_peak, r

def ghosting_analysis(phantom_epi, time_arr, PE_matrix_size, num_rep_no_dummy, slice_to_plot):
    
    #create a mask around the phantom
    threshold_for_binarizing = threshold_li(phantom_epi.flatten())
    mask_per_time = phantom_epi > threshold_for_binarizing
    phantom_mask = np.mean(mask_per_time, axis = 0) > 0   
    nophantom_mask = np.mean(mask_per_time, axis = 0) == 0
    
    #shift the phantom mask to create a mask for the Nyquist ghost
    pixel_shift = round(PE_matrix_size/2) #N/2 ghosts are shifted by N/2 pixels
    phantom_mask_rolled = np.roll(phantom_mask, pixel_shift, axis=1)
    ghost_mask_tmp = np.subtract(phantom_mask_rolled.astype(int), phantom_mask.astype(int))>0
    ghost_mask = ghost_mask_tmp.astype(int)
    
    #finally, create a mask for the background (no ghosts or phantom)
    nobackground_mask = phantom_mask.astype(int) + ghost_mask
    background_mask = 1 - nobackground_mask
    background_mask_eroded = scipy.ndimage.morphology.binary_erosion(background_mask, iterations = 3, border_value = 1)
    
    #apply masks at each timepoint
    mean_ghosting_per_rep = np.zeros(num_rep_no_dummy)
    ratio_ghosting_per_rep = np.zeros(num_rep_no_dummy)
    for i in range(0, num_rep_no_dummy):
        epi_ghosting_image = np.multiply(phantom_epi[i,:,:,:], ghost_mask)
        epi_background_image = np.multiply(phantom_epi[i,:,:,:], background_mask_eroded)
        ratio_ghosting_per_rep[i] = np.divide(epi_ghosting_image[epi_ghosting_image != 0.00].mean(), 
                                              epi_background_image[epi_background_image != 0.00].mean())
    mean_ghost = np.mean(ratio_ghosting_per_rep)
    fig1, axs = plt.subplots(1, 4, figsize = (15,4))
    fig1.suptitle('Ghosting Analysis', y = 1, fontsize = 15)
    plot = axs[0].plot(time_arr, ratio_ghosting_per_rep)
    axs[0].set_xlabel('Repetition (#)', fontsize = 12)
    axs[0].set_ylabel('Ghost to Background Intensities', fontsize = 12)
    axs[0].set_title('Ghosting Level across time, (mean = ' + str(mean_ghost) + ')', fontsize = 15)
    
    toplot = phantom_epi[int(len(time_arr)/2), slice_to_plot,:,:] + 500*ghost_mask[slice_to_plot,:,:]
    im1 = axs[1].imshow(toplot, origin = 'lower')
    cbar = plt.colorbar(im1, ax = axs[1], orientation = 'horizontal')
    axs[1].set_title('Ghost Mask Location', fontsize = 15)
    axs[1].axis('off')

    im2 = axs[2].imshow(epi_ghosting_image[slice_to_plot,:,:], origin = 'lower')
    cbar = plt.colorbar(im2, ax = axs[2], orientation = 'horizontal')
    cbar.set_label('Intensity (a.u.)')
    axs[2].set_title('Ghosting image', fontsize = 15)
    axs[2].axis('off')
    
    im3 = axs[3].imshow(epi_background_image[slice_to_plot,:,:], origin = 'lower')
    cbar = plt.colorbar(im3, ax = axs[3], orientation = 'horizontal')
    cbar.set_label('Intensity (a.u.)')
    axs[3].set_title('Background image', fontsize = 15)
    axs[3].axis('off')
  
    
    return fig1, mean_ghost

def weisskoff_analysis(phantom_epi, time, slices, PE_matrix_size, FE_matrix_size, num_rep_no_dummy, max_roi_width):
    
    ############################################ calculate SNR0 ########################################33
    #extract equal ROIs in background and center
    background_roi = extract_an_roi_background(slices, PE_matrix_size, FE_matrix_size, 15)
    center_roi = extract_an_roi(slices, PE_matrix_size, FE_matrix_size, 15)
    
    #extract the average timeseries in each roi
    tmp, tmp, tmp, spatial_std_background = extract_residuals(phantom_epi, background_roi, time, slices,
                                                                             PE_matrix_size, FE_matrix_size)
    tmp, tmp, spatial_mean_center, tmp = extract_residuals(phantom_epi, center_roi, time, slices,
                                                                             PE_matrix_size, FE_matrix_size)
    #calculate SNR0
    SNR0 = np.mean(spatial_mean_center)/(1.53*np.mean(spatial_std_background))
    expected_deviation_of_one_pixel = 100/SNR0
    
    ######################################### Find CV across ROI widths ###################################
    cv_of_mean_timeseries = np.zeros(max_roi_width)
    
    #plot the location of all the rois (to check that they make sense and don't go outside phantom)
    num_rows_in_fig = int(max_roi_width/5)
    fig1, axs = plt.subplots(num_rows_in_fig, 5, figsize = (15, 7))
    fig1.suptitle('ROI Locations for Weisskoff analysis', fontsize = 15)
    
    #iterate over various roi widths
    for roi_width in range(1,max_roi_width+1):
        #create rois of various sizes, plot each one
        new_roi = extract_an_roi(slices, PE_matrix_size, FE_matrix_size, roi_width)
        image_to_plot = 100*(new_roi[int(slices/2), :, :] + background_roi[int(slices/2), :, :])+ np.mean(phantom_epi, axis = 0)[int(slices/2),:,:]
        current_row = int((roi_width-1)/5)
        current_col = np.mod(roi_width-1, 5)
        axs[current_row, current_col].imshow(image_to_plot, origin='lower')
        axs[current_row, current_col].axis('off')
        axs[current_row, current_col].set_title(roi_width, fontsize = 15)
        
        #within each roi, extract the mean timeseries (both detrended and non-detrended)
        epi_mean_detrended_timeseries_in_roi,tmp, epi_mean_timeseries_in_roi, tmp = extract_residuals(phantom_epi, new_roi, 
                                                                                                       time, slices,
                                                                                                       PE_matrix_size, 
                                                                                                       FE_matrix_size)
        cv_of_mean_timeseries[roi_width-1] = np.std(epi_mean_detrended_timeseries_in_roi)/np.mean(epi_mean_timeseries_in_roi)
    fig1.tight_layout()
    
    ############################### Compute 'radius of decorrelation' ####################
    rdc = expected_deviation_of_one_pixel/(100*cv_of_mean_timeseries[max_roi_width-1])    
        
    fig2 = plt.figure(figsize = (15,7))
    roi_widths_arr = np.arange(1,max_roi_width+1)
    theoretical = expected_deviation_of_one_pixel/roi_widths_arr
    plt.plot(roi_widths_arr, 100*cv_of_mean_timeseries, 'o', label = 'Measured')
    plt.plot(roi_widths_arr, theoretical, label = 'Theoretical')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Weisskoff Analysis', fontsize = 15)
    plt.xlabel('ROI width (# of voxels)', fontsize = 15)
    plt.ylabel('Coefficient of variation', fontsize = 15)
    fig2.text(0,-0.2, "Radius of decorrelation: " + str(round(rdc,2)) + ' pixels', fontsize = 15)
    fig2.legend(fontsize = 15)
    
    return rdc, fig1, fig2


def pca_analysis(agar_epi_flat_detrended, time, slices, PE_matrix_size, FE_matrix_size, num_rep, TR):

    pca = PCA()
    pc_space = pca.fit_transform(agar_epi_flat_detrended)
    pc_time = pca.components_
    pc_exp_var = pca.explained_variance_

    num_components = len(pc_exp_var)
    xf = scipy.fft.fftfreq(num_rep, TR)[1:(num_rep+1)//2]

    fig1, axs = plt.subplots(2, 6, figsize = (15,7))
    fig1.suptitle('Principal Component Analysis - timecourse and FT of top 6 PCs', fontsize = 15)

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
    fig1.tight_layout()

    #################################### Plot the spatial pattern of the first 5 components ############################

    #reshape the 1d array into the original image dimensions
    pc_space_im = pc_space.reshape(slices,PE_matrix_size, FE_matrix_size,num_components)
    
    #decide how many subplots are necessary (based on the number of slices)
    root = np.sqrt(slices)+1
    subplot_dim1 = math.floor(root)
    subplot_dim2 = math.ceil(slices/subplot_dim1)
    fig_dim1 = 18
    fig_dim2 = 11
    
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
    
    if (subplot_dim1 == 1) or (subplot_dim2 == 1):
        axs2 = np.expand_dims(axs2, axis=0)
        axs3 = np.expand_dims(axs3, axis=0)
        axs4 = np.expand_dims(axs4, axis=0)
        axs5 = np.expand_dims(axs5, axis=0)
        axs6 = np.expand_dims(axs6, axis=0)
        axs7 = np.expand_dims(axs7, axis=0)
        
    slice_num = 0
    max_val = 5
    min_val = -5
    for j in range(0,subplot_dim2):
        for k in range(0,subplot_dim1):

            if slice_num >= slices:
                break
            im2 = axs2[j,k].imshow(pc_space_im[slice_num,:,:,0], vmax = max_val, vmin = min_val, origin = 'lower')
            axs2[j,k].set_title('Slice #' + str(slice_num), fontsize = 12)
            axs2[j,k].axis('off')

            im3 = axs3[j,k].imshow(pc_space_im[slice_num,:,:,1], vmax = max_val, vmin = min_val, origin = 'lower')
            axs3[j,k].set_title('Slice #' + str(slice_num), fontsize = 12)
            axs3[j,k].axis('off')

            im4 = axs4[j,k].imshow(pc_space_im[slice_num,:,:,2], vmax = max_val, vmin = min_val, origin = 'lower')
            axs4[j,k].set_title('Slice #' + str(slice_num), fontsize = 12)
            axs4[j,k].axis('off')

            im5 = axs5[j,k].imshow(pc_space_im[slice_num,:,:,3], vmax = max_val, vmin = min_val, origin = 'lower')
            axs5[j,k].set_title('Slice #' + str(slice_num), fontsize = 12)
            axs5[j,k].axis('off')

            im6 = axs6[j,k].imshow(pc_space_im[slice_num,:,:,4], vmax = max_val, vmin = min_val, origin = 'lower')
            axs6[j,k].set_title('Slice #' + str(slice_num), fontsize = 12)
            axs6[j,k].axis('off')

            im7 = axs7[j,k].imshow(pc_space_im[slice_num,:,:,5], vmax = max_val, vmin = min_val, origin = 'lower')
            axs7[j,k].set_title('Slice #' + str(slice_num), fontsize = 12)
            axs7[j,k].axis('off')

            slice_num = slice_num + 1
    cbar1 = fig2.colorbar(im2, ax = axs2, orientation = 'horizontal', fraction = 0.05)
    cbar2 = fig3.colorbar(im3, ax = axs3, orientation = 'horizontal', fraction = 0.05)
    cbar3 = fig4.colorbar(im4, ax = axs4, orientation = 'horizontal', fraction = 0.05)
    cbar4 = fig5.colorbar(im5, ax = axs5, orientation = 'horizontal', fraction = 0.05)
    cbar5 = fig6.colorbar(im6, ax = axs6, orientation = 'horizontal', fraction = 0.05)
    cbar6 = fig7.colorbar(im7, ax = axs7, orientation = 'horizontal', fraction = 0.05)
    cbar1.set_label('PC loading strength', fontsize = 12)
    cbar2.set_label('PC loading strength', fontsize = 12)
    cbar3.set_label('PC loading strength', fontsize = 12)
    cbar4.set_label('PC loading strength', fontsize = 12)
    cbar5.set_label('PC loading strength', fontsize = 12)
    cbar6.set_label('PC loading strength', fontsize = 12)
    cbar1.ax.tick_params(labelsize=12)
    cbar2.ax.tick_params(labelsize=12)
    cbar3.ax.tick_params(labelsize=12)
    cbar4.ax.tick_params(labelsize=12)
    cbar5.ax.tick_params(labelsize=12)
    cbar6.ax.tick_params(labelsize=12)

    return fig1, fig2, fig3, fig4, fig5, fig6, fig7

def export_csv(csv_metrics, output_filepath, longitudinal_csv):
    
    #define header
    csv_header = ["Mean signal", "Temporal standard deviation", "tSNR", "SNR", "Drift", "Drift of fit", "Percent fluctuation", "Peak Fourier frequency",
                   "QQ correlation", "Mean ghost", "Weisskoff radius of decorrelation"]

    #create a clean csv with just this session's data
    with open(output_filepath + ".csv", 'w', newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(csv_header)
        csv_writer.writerow(csv_metrics)
            
    #if a longitudinal csv is provided, append to that one
    if longitudinal_csv is not None:
        
        #read provided csv
        longitudinal_df = pd.read_csv(longitudinal_csv).reset_index(drop=True)
        
        #check that the provided csv has the proper headings - if not, raise error
        if longitudinal_df.axes[1].values.tolist() != csv_header:
             raise ValueError(f""" The provided longitudinal_csv does not have the correct format. Please provide a csv
             outputed by this script from a previous session with no modifications. """)
                
        #append the new metrics from the current session to the df
        longitudinal_df.loc[len(longitudinal_df)] = csv_metrics
        
        #convert_to_csv
        num_ses = len(longitudinal_df)
        longitudinal_df = longitudinal_df.reset_index(drop=True)
        longitudinal_df.to_csv(output_filepath + "-" + str(num_ses)  + "_multisession.csv", index = False)
                               
        #plot the results across sessions
        fig, axs = plt.subplots(len(csv_header),1, figsize = (15,19), sharex = True)
        fig.suptitle('Stability metrics across sessions', fontsize = 15)

        for i in range(0, len(csv_header)):                           
            axs[i].set_title(csv_header[i], fontsize = 15)
            axs[i].plot(longitudinal_df[csv_header[i]], 'o-')
            axs[i].set_xlabel('Session number', fontsize = 15)
            axs[i].set_title(csv_header[i], fontsize = 15)
        fig.tight_layout()

        return fig

def full_analysis(phantom_epi_filepath, roi_filepath, output_filepath, slice_to_plot, TR, weisskoff_max_roi_width,
                 longitudinal_csv):

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
    rep_arr = np.linspace(0, num_rep_no_dummy-1, num_rep_no_dummy)
    time_arr = np.divide(rep_arr, TR)

    #if there is no manually drawn roi provided, extract a 10x10, one-slice roi from the middle slice
    if roi_filepath is None:
        roi = extract_an_roi(slices, PE_matrix_size, FE_matrix_size,10)
    else:
        #Deals with MINC ROIS
        roi_image = sitk.ReadImage(roi_filepath)
        roi = sitk.GetArrayFromImage(roi_image)
        roi = roi.swapaxes(0,1)

    #perform whole phantom analysis
    [figure_voxelwise_wholephantom, signal_image, std_image, sfnr_image,
     static_spatial_noise_im, agar_epi_flat_detrended, slice_num] = voxelwise_wholephantom_analysis(agar_epi, roi, time_arr,
                                                                                         slice_to_plot, slices, PE_matrix_size,
                                                                                         FE_matrix_size, num_rep)

    #perform within roi analysis
    [figure_roi_analysis, signal_summary_value, std_summary_value, sfnr_summary_value, snr,
     drift, driftfit, percent_fluc, value_of_peak, r] = roi_residuals_analysis(agar_epi, roi, time_arr, signal_image, std_image, sfnr_image,
                                                      static_spatial_noise_im, TR, num_rep)
    #perform ghosting analysis
    figure_ghosting_analysis, mean_ghost = ghosting_analysis(agar_epi, time_arr, PE_matrix_size, num_rep_no_dummy, slice_num)
    
    #perform weisskoff analysis
    [rdc, figure_weisskoff_roi_positions, figure_weisskoff_rdc] = weisskoff_analysis(agar_epi, time_arr, slices,PE_matrix_size,
                                                                          FE_matrix_size, num_rep, weisskoff_max_roi_width)

    #also perform PCA
    [figure_pca_time, figure_pca_space0,
     figure_pca_space1, figure_pca_space2,
    figure_pca_space3, figure_pca_space4,figure_pca_space5] = pca_analysis(agar_epi_flat_detrended, time_arr, slices,
                                                                           PE_matrix_size, FE_matrix_size, num_rep, TR)
       
    #export csv
    csv_metrics = [signal_summary_value, std_summary_value, sfnr_summary_value, snr,
                    drift, driftfit, percent_fluc, value_of_peak[0][0], r, mean_ghost, rdc]
    figure_longitudinal_metrics = export_csv(csv_metrics, output_filepath, longitudinal_csv)

    #export all figures to pdf
    pdf_multiplot = matplotlib.backends.backend_pdf.PdfPages(output_filepath + ".pdf")
    pdf_multiplot.savefig(figure_pca_time, bbox_inches="tight")
    pdf_multiplot.savefig(figure_pca_space0, bbox_inches="tight")
    pdf_multiplot.savefig(figure_pca_space1, bbox_inches="tight")
    pdf_multiplot.savefig(figure_pca_space2, bbox_inches="tight")
    pdf_multiplot.savefig(figure_pca_space3, bbox_inches="tight")
    pdf_multiplot.savefig(figure_pca_space4, bbox_inches="tight")
    pdf_multiplot.savefig(figure_pca_space5, bbox_inches="tight")
    pdf_multiplot.savefig(figure_ghosting_analysis, bbox_inches="tight")
    pdf_multiplot.savefig(figure_voxelwise_wholephantom, bbox_inches="tight")
    pdf_multiplot.savefig(figure_roi_analysis, bbox_inches="tight")
    pdf_multiplot.savefig(figure_weisskoff_roi_positions, bbox_inches="tight")
    pdf_multiplot.savefig(figure_weisskoff_rdc, bbox_inches="tight")
    if longitudinal_csv is not None:
        pdf_multiplot.savefig(figure_longitudinal_metrics, bbox_inches="tight")
    pdf_multiplot.close()


# # Call the function

full_analysis(input_epi, input_roi, output_path, desired_slice, TR, weisskoff_max_roi_width, longitudinal_csv)
