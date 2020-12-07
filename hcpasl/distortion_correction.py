import os
import os.path as op
from pathlib import Path
import subprocess as sp
import regtricks as rt
import numpy as np
from scipy.ndimage import binary_fill_holes
import nibabel as nb
from fsl.wrappers import bet

def generate_gdc_warp(vol, coeffs_path, distcorr_dir, interpolation=1):
    """
    Generate distortion correction warp via gradient_unwarp. 

    Args: 
        vol: path to first volume we wish to correct
        coeffs_path: path to coefficients file for the scanner (.grad)
        distcorr_dir: directory in which to put output
        interpolation: order of interpolation to be used, default is 1 
                        (this is the gradient_unwarp.py default also)
    
    Returns: 
        n/a, file 'fullWarp_abs.nii.gz' will be created in output dir
    """

    # Need to run in the output directory to make sure files end up in the
    # right place
    pwd = os.getcwd()
    os.chdir(distcorr_dir)
    cmd = ("gradient_unwarp.py {} gdc_corr_vol1.nii.gz siemens -g {} --interp_order {}"
            .format(vol, coeffs_path, interpolation))
    sp.run(cmd, shell=True)
    os.chdir(pwd)

def generate_topup_params(pars_filepath):
    """
    Generate a file containing the parameters used by topup
    """
    if os.path.isfile(pars_filepath):
        os.remove(pars_filepath)
    with open(pars_filepath, "a") as t_pars:
        t_pars.write("0 1 0 0.04845" + "\n")
        t_pars.write("0 -1 0 0.04845")

def stack_fmaps(pa_sefm, ap_sefm, savename):
    rt.ImageSpace.save_like(
        ref=str(pa_sefm), 
        data=np.stack((nb.load(pa_sefm).get_fdata(), nb.load(ap_sefm).get_fdata()), axis=-1), 
        path=str(savename)
    )

def generate_fmaps(pa_ap_sefms, params, config, distcorr_dir): 
    """
    Generate fieldmaps via topup for use with asl_reg. 

    Args: 
        asl_vol0: path to image of stacked blipped images (ie, PEdir as vol0,
            (oPEdir as vol1), in this case stacked as pa then ap)
        params: path to text file for topup --datain, PE directions/times
        config: path to text file for topup --config, other args 
        distcorr_dir: directory in which to put output
    
    Returns: 
        n/a, files 'fmap, fmapmag, fmapmagbrain.nii.gz' will be created in output dir
    """

    pwd = os.getcwd()
    os.chdir(distcorr_dir)

    # Run topup to get fmap in Hz 
    topup_fmap = op.join(distcorr_dir, 'topup_fmap_hz.nii.gz')        
    cmd = (("topup --imain={} --datain={}".format(pa_ap_sefms, params)
            + " --config={} --out=topup".format(config))
            + " --fout={} --iout={}".format(topup_fmap,
                op.join(distcorr_dir, 'corrected_sefms.nii.gz')))
    sp.run(cmd, shell=True)

    fmap, fmapmag, fmapmagbrain = [ 
        op.join(distcorr_dir, '{}.nii.gz'.format(s)) 
        for s in [ 'fmap', 'fmapmag', 'fmapmagbrain' ]
    ]    

    # Convert fmap from Hz to rad/s
    fmap_spc = rt.ImageSpace(topup_fmap)
    fmap_arr_hz = nb.load(topup_fmap).get_data()
    fmap_arr = fmap_arr_hz * 2 * np.pi
    fmap_spc.save_image(fmap_arr, fmap)

    # Mean across volumes of corrected sefms to get fmapmag
    fmapmag_arr = nb.load(op.join(
                    distcorr_dir, "corrected_sefms.nii.gz")).get_data()
    fmapmag_arr = fmapmag_arr.mean(-1)
    fmap_spc.save_image(fmapmag_arr, fmapmag)

    # Run BET on fmapmag to get brain only version 
    bet(fmap_spc.make_nifti(fmapmag_arr), output=fmapmagbrain)

    os.chdir(pwd)
    
def register_fmap(fmapmag, fmapmagbrain, s, sbet, out_dir, wm_tissseg):
    # create output directory
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    # get schedule
    fsldir = os.environ.get('FSLDIR')
    schedule = Path(fsldir)/'etc/flirtsch/bbr.sch'

    # set up commands
    init_xform = out_dir/'fmapmag2struct_init.mat'
    sec_xform = out_dir/'fmapmag2struct_sec.mat'
    bbr_xform = out_dir/'fmapmag2struct_bbr.mat'
    init_cmd = [
        'flirt',
        '-in', fmapmagbrain,
        '-ref', sbet,
        '-dof', '6',
        '-omat', init_xform
    ]
    sec_cmd = [
        'flirt',
        '-in', fmapmag,
        '-ref', s,
        '-dof', '6',
        '-init', init_xform,
        '-omat', sec_xform,
        '-nosearch'
    ]
    bbr_cmd = [
        'flirt',
        '-ref', s,
        '-in', fmapmag,
        '-dof', '6',
        '-cost', 'bbr',
        '-wmseg', wm_tissseg,
        '-init', sec_xform,
        '-omat', bbr_xform,
        '-schedule', schedule
    ]
    for cmd in (init_cmd, sec_cmd, bbr_cmd):
        sp.run(cmd, check=True)
    return str(bbr_xform)

def gradunwarp_and_topup(vol, coeffs_path, distcorr_dir, pa_sefm, ap_sefm, 
                         interpolation=1, force_refresh=True):
    """
    Run gradient_unwarp and topup.

    Parameters
    ----------
    vol: path to volume to correct
    coeffs_path: path to gradient coefficients
    distcorr_dir: pathlib.Path to base results directory
    pa_sefm: path to PA spin-echo fieldmap image
    ap_sefm: path to AP spin-echo fieldmap image
    interpolation: integer order for image interpolation, default 1
    force_refresh: Boolean whether to refresh already existing files, default True

    Returns
    -------
    n/a: Saves outputs to file in subdirectories within distcorr_dir.
    """
    # run gradient_unwarp
    gdc_dir = distcorr_dir/"gradient_unwarp"
    gdc_dir.mkdir(exist_ok=True)
    gdc_warp = gdc_dir/"fullWarp_abs.nii.gz"
    if not gdc_warp.exists() or force_refresh:
        generate_gdc_warp(vol, coeffs_path, gdc_dir, interpolation)

    # create topup results directory
    topup_dir = distcorr_dir/"topup"
    topup_dir.mkdir(exist_ok=True)
    # stack epi images together for use with topup
    pa_ap_sefms = topup_dir/"merged_sefms.nii.gz"
    if not pa_ap_sefms.exists() or force_refresh:
        stack_fmaps(pa_sefm, ap_sefm, pa_ap_sefms)
    # generate topup params
    topup_params = topup_dir/"topup_params.txt"
    if not topup_params.exists() or force_refresh:
        generate_topup_params(topup_params)
    # run topup
    topup_config = "b02b0.cnf"
    fmap, fmapmag, fmapmagbrain = [topup_dir/f"fmap{ext}.nii.gz" for ext in ('', 'mag', 'magbrain')]
    if not all([f.exists() for f in (fmap, fmapmag, fmapmagbrain)]) or force_refresh:
        generate_fmaps(pa_ap_sefms, topup_params, topup_config, topup_dir)

def generate_epidc_warp(asl_vol0_brain, struct, struct_brain, asl_mask,
                       wmmask, asl2struct, fmap, fmapmag, fmapmagbrain, 
                       distcorr_dir, interpolation=3):
    """
    Generate EPI distortion correction warp via asl_reg. 

    Args: 
        asl_vol0_brain: path to first volume of ASL series, brain-extracted
        struct: path to T1 image, ac_dc_restore
        struct_brain: path to brain-extracted T1 image, ac_dc_restore_brain
        asl_mask: path to brain mask in ASL space 
        wmmask: path to WM mask in T1 space 
        asl2struct: regtricks.Registration for asl to structural
        fmap: path to topup's field map in rad/s
        fmapmag: path to topup's field map magnitude 
        fmapmagbrain: path to topup's field map magnitude, brain only
        distcorr_dir: path to directory in which to place output 

    Returns: 
        n/a, file 'asl2struct_warp.nii.gz' is created in output directory 
    """

    a2s_fsl = op.join(distcorr_dir, 'asl2struct.mat')
    asl2struct.save_fsl(a2s_fsl, asl_vol0_brain, struct)

    # get linear registration from fieldmaps to structural
    fmap_struct_dir = op.join(distcorr_dir, "fmap_struct_reg")
    bbr_fmap2struct = register_fmap(fmapmag, fmapmagbrain, struct, struct_brain, fmap_struct_dir, wmmask)

    # apply linear registration to fmap, fmapmag and fmapmagbrain
    bbr_fmap2struct = rt.Registration.from_flirt(bbr_fmap2struct, src=str(fmapmag), ref=str(struct))
    fmap_struct, fmapmag_struct, fmapmagbrain_struct = [
        op.join(fmap_struct_dir, f"fmap{ext}_struct.nii.gz") for ext in ("", "mag", "magbrain")
    ]
    for fmap_name, fmapstruct_name in zip((fmap, fmapmag, fmapmagbrain),
                                          (fmap_struct, fmapmag_struct, fmapmagbrain_struct)):
        fmapstruct_img = bbr_fmap2struct.apply_to_image(
            str(fmap_name), str(struct), order=interpolation
        )
        nb.save(fmapstruct_img, fmapstruct_name)

    # run asl_reg using pre-registered fieldmap images
    cmd = ("asl_reg -i {} -o {} ".format(asl_vol0_brain, distcorr_dir)
           + "-s {} --sbet={} -m {} ".format(struct, struct_brain, asl_mask)
           + "--tissseg={} --imat={} --finalonly ".format(wmmask, a2s_fsl)
           + "--fmap={} --fmapmag={} ".format(fmap_struct, fmapmag_struct)
           + "--fmapmagbrain={} --nofmapreg ".format(fmapmagbrain_struct)
           + "--echospacing=0.00057 --pedir=y")
    sp.run(cmd, shell=True)

def generate_asl_mask(struct_brain, asl, asl2struct):
    """
    Generate brain mask in ASL space 

    Args: 
        struct_brain: path to T1 brain-extracted, ac_dc_restore_brain
        asl: path to ASL image 
        asl2struct: regtricks.Registration for asl to structural 

    Returns: 
        np.array, logical mask. 
    """

    brain_mask = (nb.load(struct_brain).get_data() > 0).astype(np.float32)
    asl_mask = asl2struct.inverse().apply_to_array(brain_mask, struct_brain, asl)
    asl_mask = binary_fill_holes(asl_mask > 0.25)
    return asl_mask