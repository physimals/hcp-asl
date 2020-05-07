"""
Set of functions for estimating the MT effect
"""
from m0_mt_correction import load_json, update_json
from initial_bookkeeping import create_dirs
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

def estimate_mt(subject_dirs, rois=['wm', ]):
    # scan parameters
    slice_in_band = np.tile(np.arange(0, 10), 6).reshape(1, 1, -1)
    slicedt = 0.059
    delayed_times = 8 + (slice_in_band * slicedt)
    for tissue in rois:
        t1 = T1_VALS[tissue]
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
                # load masked calibration data
                masked_data = Image(masked_name).data
                # correct for timing
                denominator = 1 - np.exp(-delayed_times / t1)
                numerator = 1 - np.exp(-8/t1)
                masked_data = masked_data * (numerator / denominator)
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
        for band in range(1, 5):
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
        plt.ylim([0, 1000])
        for x_coord in x_coords:
            plt.axvline(x_coord, linestyle='-', linewidth=0.1, color='k')

        # save scaling factors
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
            figure_name = scaling_dir / f'MTcorr_SFs_{tissue}.png'
            plt.savefig(figure_name)