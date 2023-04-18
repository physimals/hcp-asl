"""
Set of functions for estimating the MT effect
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from fsl.data.image import Image
from sklearn.linear_model import LinearRegression

T1_VALS = {"wm": 1.0, "gm": 1.3, "csf": 4.3}
BAND_RANGE = {
    "wm": range(1, 5),
    "gm": range(1, 5),
    "csf": [
        3,
    ],
    "combined": range(1, 5),
}
PLOT_LIMS = {"wm": 1000, "gm": 1000, "csf": 1500, "combined": 1000}
# scan parameters
slice_in_band = np.tile(np.arange(0, 10), 6).reshape(1, 1, -1)
slicedt = 0.059


def slicetime_correction(image, tissue, tr):
    """
    Rescale data to the given TR to account for T1 relaxation.
    """
    slice_times = tr + (slice_in_band * slicedt)
    denominator = 1 - np.exp(-slice_times / T1_VALS[tissue])
    numerator = 1 - np.exp(-tr / T1_VALS[tissue])
    rescaled_image = image * (numerator / denominator)
    return rescaled_image


def undo_st_correction(rescaled_image, tissue, ti):
    slice_times = 8 + (slice_in_band * slicedt)
    numerator = 1 - np.exp(-slice_times / T1_VALS[tissue])
    denominator = 1 - np.exp(-ti / T1_VALS[tissue])
    descaled_image = rescaled_image * (numerator / denominator)
    return descaled_image


def fit_linear_model(slice_means, method="separate", resolution=10000):
    X = np.arange(0, 10, 1).reshape(-1, 1)
    scaling_factors = np.ones((6, 10))
    y_pred = np.zeros(shape=(resolution * 6, 1))
    if method == "separate":
        X_pred = np.arange(0, 10, 10 / resolution).reshape(-1, 1)
        for band in range(1, 5):
            y = slice_means[10 * band : 10 * band + 10]
            if np.isnan(y).any():
                scaling_factors[band, :] = 1
                y_pred[resolution * band : resolution * (band + 1), 0] = 0
                continue
            model = LinearRegression()
            model.fit(X, y)
            scaling_factors[band, :] = model.intercept_ / model.predict(X)
            y_pred[resolution * band : resolution * (band + 1), 0] = model.predict(
                X_pred
            )
    elif method == "together":
        X_pred = np.tile(np.arange(0, 10, 10 / resolution), 4)[..., np.newaxis]
        y_train = np.vstack(np.split(slice_means, 6)[1:5])
        y_train = np.nanmean(y_train, axis=0).reshape(-1, 1)
        model = LinearRegression()
        model.fit(X, y_train)
        sfs = model.intercept_ / model.predict(np.tile(X, (4, 1))).flatten()
        scaling_factors[1:5, :] = sfs.reshape(4, 10)
        y_pred[resolution : resolution * 5] = model.predict(X_pred)
    scaling_factors[[0, 5], :] = scaling_factors[1:5, :].mean(axis=0)
    scaling_factors = scaling_factors.flatten()
    return scaling_factors, X_pred, y_pred


def estimate_mt(
    subject_dirs,
    rois=[
        "wm",
    ],
    tr=8,
    method="separate",
    outdir=None,
    ignore_dropouts=False,
):
    """
    Estimates the slice-dependent MT effect on the given subject's
    calibration images. Performs the estimation using a linear
    model and calculates scaling factors which can be used to
    correct the effect.
    """
    outdir = Path(outdir).resolve(strict=True) if outdir else Path.cwd()
    errors = []
    suf = "_ignoredropouts" if ignore_dropouts else ""
    error_free_subs = []
    for tissue in rois:
        # initialise array to store image-level means
        mean_array = np.zeros((60, 2 * len(subject_dirs)))
        count_array = np.zeros((60, 2 * len(subject_dirs), 2))  # wm and gm
        # iterate over subjects
        for n1, subject_dir in enumerate(subject_dirs):
            try:
                print(subject_dir)
                mask_dirs = [
                    subject_dir
                    / "ASL/Calib"
                    / c
                    / f"SEbased_MT_t1mask{suf}/DistCorr/masks"
                    for c in ("Calib0", "Calib1")
                ]
                tissues = ("gm", "wm") if tissue == "combined" else (tissue,)
                masked_names = [
                    [mask_dir / tissue / f"calib{n}_{t}_masked.nii.gz" for t in tissues]
                    for n, mask_dir in enumerate(mask_dirs)
                ]
                for n2, masked_name in enumerate(masked_names):
                    if tissue == "combined":
                        gm_masked, wm_masked = masked_name
                        gm_masked_data = slicetime_correction(
                            image=Image(str(gm_masked)).data, tissue="gm", tr=tr
                        )
                        wm_masked_data = slicetime_correction(
                            image=Image(str(wm_masked)).data, tissue="wm", tr=tr
                        )
                        masked_data = gm_masked_data + wm_masked_data
                        gm_bin = np.where(gm_masked_data > 0, 1, 0)
                        gm_count = np.sum(gm_bin, axis=(0, 1))[..., np.newaxis]
                        wm_bin = np.where(wm_masked_data > 0, 1, 0)
                        wm_count = np.sum(wm_bin, axis=(0, 1))[..., np.newaxis]
                        count_array[:, 2 * n1 + n2, :] = np.hstack((wm_count, gm_count))
                    else:
                        # load masked calibration data
                        masked_data = slicetime_correction(
                            image=Image(str(*masked_name)).data, tissue=tissue, tr=tr
                        )
                    # find zero indices
                    masked_data[masked_data == 0] = np.nan
                    # calculate slicewise summary stats
                    slicewise_mean = np.nanmean(masked_data, axis=(0, 1))
                    mean_array[:, 2 * n1 + n2] = slicewise_mean
                error_free_subs.append(subject_dir)
            except:
                errors.append(tissue + " " + str(subject_dir))
        # calculate non-zero slicewise mean of mean_array
        slice_means = np.nanmean(mean_array, axis=1)
        slice_std = np.nanstd(mean_array, axis=1)

        # calculate slicewise mean of tissue type counts
        count_means = np.nanmean(count_array, axis=1)

        # fit linear models to central 4 bands
        # estimate scaling factors using these models
        scaling_factors, X_pred, y_pred = fit_linear_model(slice_means, method=method)
        # plot slicewise mean signal
        slice_numbers = np.arange(0, 60, 1)
        x_coords = np.arange(0, 60, 10)
        plt.figure(figsize=(8, 4.5))
        plt.scatter(slice_numbers, slice_means)
        plt.errorbar(slice_numbers, slice_means, slice_std, linestyle="None", capsize=3)
        plt.ylim([0, PLOT_LIMS[tissue]])
        plt.xlim([0, 60])
        if tissue == "combined":
            plt.title(
                f"Mean signal per slice in GM and WM across {int(len(error_free_subs)/len(rois))} subjects."
            )
        else:
            plt.title(
                f"Mean signal per slice in {tissue} ({method}) across {int(len(error_free_subs)/len(rois))} subjects."
            )
        plt.xlabel("Slice number")
        plt.ylabel("Mean signal")
        for x_coord in x_coords:
            plt.axvline(x_coord, linestyle="-", linewidth=0.1, color="k")
        # save plot
        plt_name = outdir / f"{method}_{tissue}_mean_per_slice_t1.png"
        plt.savefig(plt_name)
        # add linear models on top
        plt.scatter(
            np.arange(10, 50, 0.001), y_pred.flatten()[10000:50000], color="k", s=0.1
        )
        plt_name = outdir / f"{method}_{tissue}_mean_per_slice_with_lin_sebased.png"
        plt.savefig(plt_name)

        # plot rescaled slice-means
        fig, ax = plt.subplots(figsize=(8, 4.5))
        rescaled_means = slice_means * scaling_factors
        plt.scatter(slice_numbers, rescaled_means)
        plt.ylim([0, PLOT_LIMS[tissue]])
        plt.xlim([0, 60])
        if tissue == "combined":
            plt.title(
                f"Rescaled mean signal per slice in GM and WM across {int(len(error_free_subs)/len(rois))} subjects."
            )
        else:
            plt.title(
                f"Rescaled mean signal per slice in {tissue} across {int(len(error_free_subs)/len(rois))} subjects."
            )
        plt.xlabel("Slice number")
        plt.ylabel("Rescaled mean signal")
        for x_coord in x_coords:
            plt.axvline(x_coord, linestyle="-", linewidth=0.1, color="k")
        # save plot
        plt_name = outdir / f"{method}_{tissue}_mean_per_slice_rescaled_sebased.png"
        plt.savefig(plt_name)

        # plot slicewise mean tissue count for WM and GM
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.scatter(
            slice_numbers, count_means[:, 0], c="c", label="WM pve > 70%"
        )  # wm mean
        ax.scatter(
            slice_numbers, count_means[:, 1], c="g", label="GM pve > 70%"
        )  # gm mean
        ax.legend()
        for x_coord in x_coords:
            ax.axvline(x_coord, linestyle="-", linewidth=0.1, color="k")
        plt.title(
            "Mean number of voxels per slice with"
            + f" PVE $\geqslant$ 70% across {int(len(error_free_subs)/len(rois))} subjects."
        )
        plt.xlabel("Slice number")
        plt.ylabel("Mean number of voxels with PVE $\geqslant$ 70% in a given tissue")
        plt_name = outdir / f"mean_voxel_count_sebased.png"
        plt.savefig(plt_name)

        # # the scaling factors have been estimated on images which have been
        # # slice-timing corrected - the scaling factors should hence be
        # # adjusted to account for this, as they will be applied to images
        # # which haven't had this correction
        # scaling_factors = undo_st_correction(scaling_factors, tissue, tr)

        # save scaling factors as a .txt file
        sfs_savename = outdir / f"{method}_{tissue}_scaling_factors_sebased.txt"
        np.savetxt(sfs_savename, scaling_factors, fmt="%.5f")
        # create array from scaling_factors
        scaling_factors = np.tile(scaling_factors, (86, 86, 1))
        for subject_dir in subject_dirs:
            # load bias (and possibly distortion) corrected calibration image
            method_dir = (
                subject_dir / f"ASL/Calib/Calib0/SEbased_MT_t1mask{suf}/DistCorr"
            )
            calib_name = method_dir / "calib0_restore.nii.gz"
            calib_img = Image(str(calib_name))
            # create and save scaling factors image
            scaling_img = Image(scaling_factors, header=calib_img.header)
            mtcorr_dir = method_dir / "MTCorr"
            mtcorr_dir.mkdir(exist_ok=True)
            scaling_name = mtcorr_dir / f"MTcorr_SFs_{method}_{tissue}_sebased.nii.gz"
            scaling_img.save(scaling_name)

            # apply scaling factors to image to perform MT correction
            mtcorr_name = mtcorr_dir / f"calib0_mtcorr_{method}_{tissue}_sebased.nii.gz"
            mtcorr_img = Image(
                calib_img.data * scaling_factors, header=calib_img.header
            )
            mtcorr_img.save(str(mtcorr_name))
    return errors
