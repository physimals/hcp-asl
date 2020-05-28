"""
Set of functions for estimating the MT effect
"""
from hcpasl.m0_mt_correction import load_json, update_json
from hcpasl.initial_bookkeeping import create_dirs
from fsl.data.image import Image
import numpy as np
from sklearn.linear_model import LinearRegression
from pathlib import Path
import multiprocessing
import matplotlib.pyplot as plt

T1_VALS = {
    'wm': 1.0,
    'gm': 1.3,
    'csf': 4.3
}
BAND_RANGE = {
    'wm': range(1, 5),
    'gm': range(1, 5),
    'csf': [3, ],
    'combined': range(1, 5)
}
PLOT_LIMS = {
    'wm': 1000,
    'gm': 1000,
    'csf': 1500,
    'combined': 1000
}
# scan parameters
slice_in_band = np.tile(np.arange(0, 10), 6).reshape(1, 1, -1)
slicedt = 0.059
def slicetime_correction(image, tissue, tr):
    """
    Rescale data to the given TR to account for T1 relaxation.
    """
    slice_times = tr + (slice_in_band * slicedt)
    denominator = 1 - np.exp(-slice_times/T1_VALS[tissue])
    numerator = 1 - np.exp(-tr/T1_VALS[tissue])
    print(numerator / denominator)
    rescaled_image = image * (numerator / denominator)
    return rescaled_image

def undo_st_correction(rescaled_image, tissue, ti):
    slice_times = 8 + (slice_in_band * slicedt)
    numerator = 1 - np.exp(-slice_times/T1_VALS[tissue])
    denominator = 1 - np.exp(-ti/T1_VALS[tissue])
    descaled_image = rescaled_image * (numerator / denominator)
    return descaled_image

def estimate_mt(subject_dirs, rois=['wm', ], tr=8):
    """
    Estimates the slice-dependent MT effect on the given subject's 
    calibration images. Performs the estimation using a linear 
    model and calculates scaling factors which can be used to 
    correct the effect.
    """
    for tissue in rois:
        # initialise array to store image-level means
        mean_array = np.zeros((60, 2*len(subject_dirs)))
        # iterate over subjects
        for n1, subject_dir in enumerate(subject_dirs):
            print(subject_dir)
            # load subject's json
            json_dict = load_json(subject_dir)
            # calculate mean per slice of masked in both calib images
            masked_names = (
                json_dict[f'calib0_{tissue}_masked'],
                json_dict[f'calib1_{tissue}_masked']
            )
            for n2, masked_name in enumerate(masked_names):
                if tissue == 'combined':
                    gm_masked, wm_masked = masked_name
                    gm_masked_data = slicetime_correction(
                        image=Image(gm_masked).data, 
                        tissue='gm',
                        tr=tr
                    )
                    wm_masked_data = slicetime_correction(
                        image=Image(wm_masked).data, 
                        tissue='wm',
                        tr=tr
                    )
                    masked_data = gm_masked_data + wm_masked_data
                else:
                    # load masked calibration data
                    masked_data = slicetime_correction(
                        image=Image(masked_name).data,
                        tissue=tissue,
                        tr=tr
                    )
                # find zero indices
                masked_data[masked_data==0] = np.nan
                # calculate slicewise summary stats
                slicewise_mean = np.nanmean(masked_data, axis=(0, 1))
                mean_array[:, 2*n1 + n2] = slicewise_mean
        
        # calculate non-zero slicewise mean of mean_array
        slice_means = np.nanmean(mean_array, axis=1)

        # fit linear models to central 4 bands
        # estimate scaling factors using these models
        X = np.arange(0, 10, 1).reshape(-1, 1)
        scaling_factors = np.ones((86, 86, 60))
        y_pred = np.zeros((40000, 1))
        for band in BAND_RANGE[tissue]:
            y = slice_means[10*band: 10*band+10]
            model = LinearRegression()
            model.fit(X, y)
            scaling_factors[:, :, 10*band: 10*band+10] = model.intercept_ / model.predict(X)
            y_pred[10000*(band-1): 10000*(band-1)+10000, 0] = model.predict(np.arange(0, 10, 0.001).reshape(-1, 1))
        
        # plot
        slice_numbers = np.arange(0, 60, 1)
        x_coords = np.arange(0, 60, 10)
        plt.figure(figsize=(8, 4.5))
        plt.scatter(slice_numbers, slice_means)
        plt.scatter(np.arange(10, 50, 0.001), y_pred.flatten(), color='k', s=0.1)
        plt.ylim([0, PLOT_LIMS[tissue]])
        plt.xlim([0, 60])
        plt.title(f'Mean signal per slice in {tissue} across 47 subjects.')
        plt.xlabel('Slice number')
        plt.ylabel('Mean signal')
        for x_coord in x_coords:
            plt.axvline(x_coord, linestyle='-', linewidth=0.1, color='k')
        # save plot
        plt_name = Path().cwd() / f'{tissue}_mean_per_slice.png'
        plt.savefig(plt_name)

        # # the scaling factors have been estimated on images which have been 
        # # slice-timing corrected - the scaling factors should hence be 
        # # adjusted to account for this, as they will be applied to images 
        # # which haven't had this correction
        # scaling_factors = undo_st_correction(scaling_factors, tissue, tr)

        for subject_dir in subject_dirs:
            json_dict = load_json(subject_dir)
            # load calibration image
            calib_img = Image(json_dict['calib0_img'])
            # create and save scaling factors image
            scaling_img = Image(scaling_factors, header=calib_img.header)
            scaling_dir = Path(json_dict['calib_dir']) / 'MTEstimation'
            create_dirs([scaling_dir, ])
            scaling_name = scaling_dir / f'MTcorr_SFs_{tissue}.nii.gz'
            scaling_img.save(scaling_name)
            
