
from scipy.stats import percentileofscore
import numpy as np


def find_short_solint(myvis, bpcal, bpscan, refant,
                      int_time=7, max_solint=60, nsteps=5,
                      target_flag_perc=20, target_snr=5,
                      minsnr=0.0, bpchans='0~11', minblperant=3,
                      priorcals=[],
                      removezerosnr=True,
                      cleanup_tables=True):
    '''
    Optimize the solution interval time.
    '''

    raise NotImplementedError("This function does not currently work.")

    solints = np.linspace(int_time, max_solint, nsteps)

    # log.info(f"Trying solution intervals of {solints} sec")
    print(f"Trying solution intervals of {solints} sec")

    table_names = []

    flag_percs = []
    median_snrs = []

    for this_solint in solints:

        phaseshortgaincal_table = f"{myvis[:-3]}_solint_test_{this_solint}"

        gaincal(vis=myvis,caltable=phaseshortgaincal_table,
                field=bpcal,spw=bpchans,refant=refant, scan=bpscan,
                calmode='p',solint=f'{this_solint}s',
                minsnr=minsnr,
                minblperant=minblperant,
                gaintable=priorcals)

        table_names.append(phaseshortgaincal_table)

        tb.open(phaseshortgaincal_table)
        snr = tb.getcol('SNR').squeeze()
        tb.close()

        # The table will write out for at least 2 correlations, even if one isn't
        # present in the data (i.e. a split tuning or a single rx in the MS)
        # Check for all zeros in one of these correlations:

        snr_valid = []
        for corr_idx in range(snr.shape[0]):
            snr_corr0 = snr[corr_idx]
            if np.unique(snr_corr0).size > 1:
                snr_valid.append(snr_corr0)

        # Reshape into an array:
        snr_valid = np.dstack(snr_valid).squeeze()

        # Optionally ignore all 0.0 missing data solutions.
        if removezerosnr:
            snr_valid = snr_valid[np.nonzero(snr_valid)]

        flag_percs.append(percentileofscore(snr_valid, target_snr))
        median_snrs.append(np.median(snr_valid))

    flag_percs = np.array(flag_percs)
    median_snrs = np.array(median_snrs)

    # Find solint that best matches target:

    satisfy_flag = flag_percs <= target_flag_perc
    satisfy_snr = median_snrs >= target_snr

    satisfy_both = np.logical_and(satisfy_flag, satisfy_snr)

    print(f"solint flagging percentages: {flag_percs}")
    print(f"solint median SNRs: {median_snrs}")

    print(f"solints satisfying flag fraction: {solints[satisfy_flag]}")
    print(f"solints satisfying flag fraction: {solints[satisfy_snr]}")
    print(f"solints satisfying both: {solints[satisfy_both]}")

    if cleanup_tables:
        for table_name in table_names:
            rmtables(table_name)

    target_solints = solints[satisfy_both]
    if target_solints.size == 0:
        print("Found no acceptable solint.")
        target_solint = None
    else:
        target_solint = target_solints[0]

    return target_solint
