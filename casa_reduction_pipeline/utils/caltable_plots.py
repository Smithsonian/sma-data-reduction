
import os
from glob import glob
import numpy as np

from casatools import logsink

casalog = logsink()

CALTABLE_MAPPING = {'bandpass_amp': {'output_folder': 'final_caltable_txt',
                                'search_string': '*.bandpass.solnorm_true.bcal',
                                'x': 'freq',
                                'y': 'amp',
                                'iter': 'spw',
                                'colorby': 'ant'},
                    'bandpass_phase': {'output_folder': 'final_caltable_txt',
                                'search_string': '*.bandpass.solnorm_true.bcal',
                                'x': 'freq',
                                'y': 'phase',
                                'iter': 'spw',
                                'colorby': 'ant'},
                    'BPinitialgain': {'output_folder': 'final_caltable_txt',
                                    'search_string': '*.bpself.ap.gcal',
                                    'x': 'time',
                                    'y': 'phase',
                                    'iter': 'ant',
                                    'colorby': 'spw'},
                    'BLcal': {'output_folder': 'final_caltable_txt',
                                    'search_string': '*.blcal',
                                    'x': 'freq',
                                    'y': 'amp',
                                    'iter': 'ant',
                                    'colorby': 'spw'},
                    'phaseshortgaincal': {'output_folder': 'final_caltable_txt',
                                    'search_string': '*.intphase_combinespw.gcal',
                                    'x': 'time',
                                    'y': 'phase',
                                    'iter': 'ant',
                                    'colorby': 'spw'},
                    'ampgaincal_time': {'output_folder': 'final_caltable_txt',
                                    'search_string': '*.amp.gcal',
                                    'x': 'time',
                                    'y': 'amp',
                                    'iter': 'ant',
                                    'colorby': 'spw'},
                    'ampgaincal_freq': {'output_folder': 'final_caltable_txt',
                                    'search_string': '*.amp.gcal',
                                    'x': 'freq',
                                    'y': 'amp',
                                    'iter': 'ant',
                                    'colorby': 'spw'},
                    'phasegaincal': {'output_folder': 'final_caltable_txt',
                                    'search_string': '*.scanphase_combinespw.gcal',
                                    'x': 'time',
                                    'y': 'phase',
                                    'iter': 'ant',
                                    'colorby': 'spw'}}


def make_caltable_txt(ms_active, caltable_type):
    '''
    Output txt files using plotms to make plots of various calibration tables.
    See definitions in `CALTABLE_MAPPING`.
    The naming convention follows the VLA pipeline table names from `hifv_finalcals`

    '''

    caltable_values = CALTABLE_MAPPING[caltable_type]

    # from taskinit import tb, casalog
    from casatools import table

    tb = table()

    from casaplotms import plotms

    casalog.post(f"Running make_caltable_txt on {caltable_type} to export txt files for QA.")
    print(f"Running make_caltable_txt on {caltable_type} to export txt files for QA.")

    mySDM = ms_active.rstrip(".ms")

    if not os.path.exists(caltable_values['output_folder']):
        os.mkdir(caltable_values['output_folder'])

    # Final BP cal table now includes the stage number and step
    caltable_name = glob(mySDM + caltable_values['search_string'])
    if len(caltable_name) == 0:
        raise ValueError("Cannot find {} table name.".format(caltable_values['search_string']))
    # Blindly assume we want the first name
    caltable_name = caltable_name[0]

    tb.open(caltable_name)
    spw_vals = np.unique(tb.getcol("SPECTRAL_WINDOW_ID"))
    ant_vals = np.unique(tb.getcol("ANTENNA1"))
    tb.close()

    # Make txt files per SPW.
    iteraxis = spw_vals if caltable_values['iter'] == 'spw' else ant_vals

    for ii in iteraxis:

        print("On spectral window: {}".format(ii))
        casalog.post("On spectral window: {}".format(ii))

        # Output text names
        # name_xaxis_yaxis_iter num
        out_filename = '{0}_{1}_{2}_{3}{4}.txt'.format(caltable_name,
                                                       caltable_values['x'],
                                                       caltable_values['y'],
                                                       caltable_values['iter'],
                                                       ii)

        thisplotfile = os.path.join(caltable_values['output_folder'], out_filename)

        if not os.path.exists(thisplotfile):

            print(caltable_name)

            plotms(vis=caltable_name,
                   xaxis=caltable_values['x'],
                   yaxis=caltable_values['y'],
                   field='',
                   antenna=str(ii) if caltable_values['iter'] == 'ant' else "",
                   spw=str(ii) if caltable_values['iter'] == 'spw' else "",
                   timerange='',
                   showgui=False,
                   # avgtime='1e8',
                   averagedata=True,
                   plotfile=thisplotfile)
        else:
            casalog.post("File {} already exists. Skipping".format(thisplotfile))


def make_all_caltable_txt(msname):

    for key in CALTABLE_MAPPING:
        make_caltable_txt(msname, key)

