
import numpy as np
import os

# Inspect the fit solution and visualize.
def fluxscale_to_tables(fluxresults):
    '''
    Make a fit solution table and flux per SPW table with measured +
    predicted fluxes.
    '''

    from astropy.table import Table, Column

    freqs = Column(fluxresults['freq'], name='freq')
    spws = Column(fluxresults['spwID'], name='spw')

    num_spws = len(spws)

    # numfield = num_keys - (freq, spwID, spwName)
    num_fields = len(fluxresults.keys()) - 3

    field_ids = [key for key in fluxresults.keys() if key.isdigit()]

    # Make a column of table names for the data table and fit table
    field_ids_fit = Column(field_ids, name='field_idx')
    field_names_fit = Column([fluxresults[idx]['fieldName'] for idx in field_ids_fit],
                              name='field')

    all_fields_ids = []
    all_field_spws = []
    all_field_freqs = []
    for field_id in field_ids:
        all_fields_ids.extend([field_id] * num_spws)
        all_field_spws.extend(list(fluxresults['spwID']))
        all_field_freqs.extend(list(fluxresults['freq']))
    field_ids_data = Column(all_fields_ids, name='field_idx')
    field_spws_data = Column(all_field_spws, name='spw')
    field_freqs_data = Column(all_field_freqs, name='freq')
    field_names_data = Column([fluxresults[idx]['fieldName'] for idx in all_fields_ids],
                              name='field')

    for nn, field_id in enumerate(field_ids):

        # field_name = fluxresults[field_id]['fieldName']

        fit_flux = fluxresults[field_id]['fitFluxd']
        fit_fluxerr = fluxresults[field_id]['fitFluxdErr']
        fit_flux_reffreq = fluxresults[field_id]['fitRefFreq']

        # 𝑙𝑜𝑔(𝑆𝜈)=𝑎𝑜+𝑎1∗(𝑙𝑜𝑔(𝜈/𝜈0))+𝑎2∗(𝑙𝑜𝑔(𝜈/𝜈0))∗∗2
        fit_pars = fluxresults[field_id]['spidx']
        fit_errs = fluxresults[field_id]['spidxerr']

        num_params = len(fit_pars)

        fit_vals = np.array([fit_flux, fit_fluxerr, fit_flux_reffreq])
        fit_vals = np.append(fit_vals, fit_pars)
        fit_vals = np.append(fit_vals, fit_errs)

        # 4 * [IQUV], 4 * [IQUV err]
        data_vals = np.zeros((len(spws), 8))

        for ii, this_spw in enumerate(spws):
            data_vals[ii, :] = np.append(fluxresults[field_id][str(this_spw)]['fluxd'],
                                         fluxresults[field_id][str(this_spw)]['fluxdErr'])

        if nn == 0:
            all_data_vals = data_vals
            all_fit_vals = fit_vals
        else:
            all_data_vals = np.vstack([all_data_vals, data_vals])
            all_fit_vals = np.vstack([all_fit_vals, fit_vals])

    # Now construct the final tables:
    # Fit table
    fit_cols = ['fitFluxd', 'fitFluxdErr', 'fitRefFreq']
    fit_cols.extend([f'a{num}' for num in range(num_params)])
    fit_cols.extend([f'a{num}_err' for num in range(num_params)])

    fit_table = Table(data=all_fit_vals, names=fit_cols)
    fit_table.add_column(field_ids_fit)
    fit_table.add_column(field_names_fit)

    # Data table
    data_cols = ['I', 'Q', 'U', 'V', 'I_err', 'Q_err', 'U_err', 'V_err']
    data_table = Table(data=all_data_vals, names=data_cols)
    data_table.add_column(field_ids_data)
    data_table.add_column(field_names_data)
    data_table.add_column(field_spws_data)
    data_table.add_column(field_freqs_data)

    return data_table, fit_table


def log_fit_func(freq, a0, a1, a2, ref_freq):
    return a0 + a1 * np.log10(freq / ref_freq) + a2 * np.log10(freq / ref_freq)**2


def return_flux_models(fit_table):
    '''
    𝑙𝑜𝑔(𝑆𝜈)=𝑎𝑜+𝑎1∗(𝑙𝑜𝑔(𝜈/𝜈0))+𝑎2∗(𝑙𝑜𝑔(𝜈/𝜈0))∗∗2
    '''

    from functools import partial

    flux_models = {}

    for ii, field in enumerate(fit_table['field']):

        a0 = fit_table['a0'][ii]
        a1 = fit_table['a1'][ii]
        ref_freq = fit_table['fitRefFreq'][ii]
        flux0 = fit_table['fitFluxd'][ii]

        if 'a2' in fit_table.colnames:
            a2 = fit_table['a2'][ii]
        else:
            a2 = 0.0

        flux_models[field] = partial(log_fit_func, a0=a0, a1=a1,
                                     a2=a2, ref_freq=ref_freq)

    return flux_models


def return_flux_fit_residuals(fit_table, data_table, flux_models=None):

    from astropy.table import Table, Column

    residuals = {}

    if flux_models is None:
        flux_models = return_flux_models(fit_table)

    for ii, field in enumerate(fit_table['field']):

        this_model_func = flux_models[field]

        mask = data_table['field'] == field

        freqs = data_table['freq'][mask]
        spws = data_table['spw'][mask]

        # Could generalize this at some point.
        stokes_I = data_table['I'][mask]
        stokes_I_err = data_table['I_err'][mask]

        resids = stokes_I - 10**this_model_func(freqs)

        # -1.0 values indicates that SPW was not used in the fit
        resids[stokes_I == -1.] = 0.

        # Calc the N-sigma deviation from the model
        sigma_resids = resids / stokes_I_err

        residual_table = Table([spws, freqs,
                                Column(resids, name="I_resid"),
                                Column(sigma_resids, name="sigma_resid")])

        residuals[field] = residual_table

    return residuals


def plot_flux_fits(fit_table, data_table,
                   flux_models=None,
                   resids_dict=None,
                   output_folder="fluxfit_plots"):
    '''
    Makes a plot per field.
    '''

    import matplotlib.pyplot as plt

    if flux_models is None:
        flux_models = return_flux_models(fit_table)

    if resids_dict is None:
        resids_dict = return_flux_fit_residuals(fit_table, data_table,
                                                flux_models=flux_models)

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for this_field in fit_table['field']:

        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True,
                                figsize=(12, 12))

        freqs = resids_dict[this_field]['freq']
        stokes_I = data_table['I'][data_table['field'] == this_field]
        stokes_I_err = data_table['I_err'][data_table['field'] == this_field]

        valid_data = stokes_I != -1.

        stokes_I_resid = resids_dict[this_field]['I_resid']

        this_model_func = flux_models[this_field]
        stokes_I_model = 10**this_model_func(freqs)

        axs[0].errorbar(freqs[valid_data] / 1e9,
                        stokes_I[valid_data],
                        yerr=stokes_I_err[valid_data])

        a0 = fit_table['a0'][fit_table['field'] == this_field][0]
        a1 = fit_table['a1'][fit_table['field'] == this_field][0]
        ref_freq_GHz = fit_table['fitRefFreq'][fit_table['field'] == this_field][0] / 1.e9
        model_label = f"{a0:.2f} + {a1:.2f} * (freq / {ref_freq_GHz:.1f})"

        axs[0].plot(freqs / 1e9,
                    stokes_I_model,
                    label=model_label)

        axs[0].legend(frameon=True, fontsize=10)

        # Residual
        axs[1].errorbar(freqs[valid_data] / 1e9,
                        stokes_I_resid[valid_data],
                        yerr=stokes_I_err[valid_data])

        axs[1].axhline(0., color='k', linewidth=2, linestyle='--', zorder=1)

        # Set symmetric y axis for residual
        resid_max = np.abs(stokes_I_resid).max()
        axs[1].set_ylim([-1.25 * resid_max, 1.25 * resid_max])

        axs[1].set_xlabel("Freq (GHz)", fontsize=10)

        axs[0].set_ylabel("Amplitude (Jy)", fontsize=10)
        axs[1].set_ylabel("Residual (Jy)", fontsize=10)

        fig.savefig(f"{output_folder}/{this_field}_fluxscale_fit.png", dpi=300)

        plt.close()


def flag_on_large_fit_residuals(fit_table, data_table):
    '''
    Pass to the logger when large fit residuals are present.
    '''
    pass
