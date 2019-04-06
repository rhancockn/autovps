#!/usr/bin/env python
import argparse
import numpy as np
import autovps.dataset.siemens as siemens
import os
parser = argparse.ArgumentParser(
    description='Process Siemens SVS data from DICOMs')
parser.add_argument('dicom', type=str, nargs='+',
                    help="Path to DICOM directory")
parser.add_argument('--prefix', type=str,
                    help="Output prefix")
parser.add_argument('--t1', type=str, required=False,
                    help="Path to co-registered T1 NIFTI (saves VOI mask)")
parser.add_argument('--load', dest='load', action='store_true',
                    help="Load FID data (required for analysis)")
parser.add_argument('--no-load', dest='load', action='store_false',
                    help="Do not load FID data")
parser.add_argument('--fida', dest='fida', action='store_true',
                    help="Run FID-A preprocessing (requires matlab)")
parser.add_argument('--no-fida', dest='fida', action='store_false',
                    help="Do not run FID-A processing")
parser.add_argument('--tarquin', dest='tarquin', action='store_true',
                    help="Run Tarquin fit")
parser.add_argument('--no-tarquin', dest='tarquin', action='store_false',
                    help="Do not run Tarquin fitting")

parser.set_defaults(load=True)
parser.set_defaults(fida=True)
parser.set_defaults(tarquin=True)

args = parser.parse_args(['--prefix', '/Users/roh17004/Downloads/TD919/autotest3',
                         '--t1','/Users/roh17004/Downloads/TD919/5/TD919-T1w.nii.gz',
                         '/Users/roh17004/Downloads/TD919/39', '/Users/roh17004/Downloads/TD919/42',
                         '/Users/roh17004/Downloads/TD919/45'])

args = parser.parse_args()

# Parse the DICOMs
runs = []
for run,fname in enumerate(args.dicom):
    print('Reading metadata from %s...\n' % fname)
    dcm = siemens.Siemens(fname)
    dcm.calculate_transform()
    # save transform
    np.savetxt('%s_run-%02d_from-device_to-orig_mode-image_xfm.mat' % (args.prefix, run), dcm.qform.get_matrix())

    # make voxel
    if args.t1 is not None:
        print('Creating VOI mask...\n')
        from autovps import make_voi
        import nibabel as nib
        tform = dcm.calculate_transform()
        t1 = nib.load(args.t1)
        img = make_voi.make_voi(t1, tform)
        nib.save(img, '%s_run-%02d_space-orig_roi.nii.gz' % (args.prefix, run))

    # Read data
    if args.load:
        print('Loading data...\n')
        svs = dcm.get_svsdata()

        # Save to FID-A .mat
        fida_fname = '%s_run-%02d_fida.mat' % (args.prefix, run)
        svs.save_fida(fida_fname)

        # TODO: save NIFTI, LCM

    runs.append(svs)

# merge runs

fids = runs[0].fid
specs = runs[0].spec

if len(runs) > 1:
    print('Merging multiple runs...')
    for idx in range(1, len(runs)):
        # check the sequence is the same
        if runs[0].sequence_name != runs[idx].sequence_name:
            raise Exception('Not all runs are the same sequence!')
        if runs[0].te != runs[idx].te:
            raise Exception('Not all runs have the same TE!')
        if runs[0].tr != runs[idx].tr:
            raise Exception('Not all runs have the same TR!')
        if runs[0].sw != runs[idx].sw:
            raise Exception('Not all runs have the same spectral width!')
        # TODO: add equality operator to Transform class
        if not np.all(runs[0].transform.get_matrix() == runs[idx].transform.get_matrix()):
            raise Exception('Not all runs are from the same location!')
        # TODO: check channels match

        fids = np.concatenate((fids, runs[idx].fid), 1)
        specs = np.concatenate((specs, runs[idx].spec), 1)
    svs = runs[0]
    svs.fid = fids
    svs.spec = specs
    fida_fname = '%s_acq-merged_fida.mat' % (args.prefix)
    svs.save_fida(fida_fname)

# Preprocess using FID-A

# TODO: Mark STEAM and sLASER
if args.fida:
    print('Preprocessing using FID-A...')

    import subprocess
    from subprocess import Popen
    import platform

    # check for macOS and try to wake up the display
    # so matlab figures will work
    p = None
    if platform.system() == 'Darwin':
        print('macOS might not work!')
        p = Popen(['caffeinate', '-u'])

    if svs.sequence_type == 'PRESS':
        cmd = "preproc_press('%s_report', [], '%s'); quit" % (args.prefix, fida_fname)
        subprocess.run(['matlab', '-nodesktop', '-r', cmd])
        if args.tarquin:
            print('Processing PRESS spectra using Tarquin')
            cmd = """tarquin --input %s_report/spectra_jmrui_ave.txt --format jmrui_txt \
                     --lipid_filter true --auto_phase true --auto_ref true \
                     --crlb_optim false --ref_signals 1h_naa_cr_cho --fs %d \
                     --ft %d --echo %f --pul_seq press \
                     --output_txt %s_report/tarquin.txt --output_pdf %s_report/tarquin.pdf \
                     --output_image %s_report/tarquin_img.pdf \
                     --ext_pdf true --output_xml %s_report/tarquin.xml --int_basis 1h_brain_glth \
                     --svs_only true --stack_pdf true --te1 .014""" % (args.prefix, svs.sw, svs.larmor, svs.te, args.prefix, args.prefix, args.prefix, args.prefix)
            print(cmd)
            subprocess.run(cmd, shell=True)

    if svs.sequence_type == 'MEGAPRESS':
        cmd = "preproc_megapress('%s_report', [], '%s'); quit" % (args.prefix, fida_fname)
        subprocess.run(['matlab', '-nodesktop', '-r', cmd])

        if args.tarquin:
            print('Processing MEGA-PRESS difference spectra using Tarquin')
            cmd = """tarquin --input %s_report/spectra_jmrui_diff_ave.txt --format jmrui_txt \
                     --lipid_filter true --auto_phase true --auto_ref true \
                     --crlb_optim false --ref_signals 1h_naa --fs %d \
                     --ft %d --echo %f --pul_seq mega_press \
                     --output_txt %s_report/tarquin.txt --output_pdf %s_report/tarquin.pdf \
                     --output_image %s_report/tarquin_img.pdf \
                     --ext_pdf true --output_xml %s_report/tarquin.xml --int_basis megapress_gaba \
                     --svs_only true --stack_pdf true --te1 .014""" % (args.prefix, svs.sw, svs.larmor, svs.te, args.prefix, args.prefix, args.prefix, args.prefix)
            print(cmd)
            subprocess.run(cmd, shell=True)

    if p is not None:
        p.terminate()


# TODO: use water for coil combinations
# TODO: merge after preprocessing
