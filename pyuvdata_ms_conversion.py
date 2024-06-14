
'''
Convert the MIR file format to a CASA measurement set (MS).

Convert a single MIR file to an MS:
>>> python pyuvdata_ms_conversion.py -f [mir_filename] -r [rechunk_factor] -o [output_folder]

Convert multiple MIR files to a single MS
>>> python pyuvdata_ms_conversion.py -f [mir_filename1,mir_filename2,...] -r [rechunk_factor] -o [output_folder]

'''

import time
import tracemalloc

def convert_mir_to_ms(filenames,
                      rechunk=1,
                      outputdir='./',
                      antfile=None,
                      verbose=True,
                      run_check=True,
                      track_mem_usage=True):


    from pyuvdata import UVData
    from pyuvdata.uvdata.mir import generate_sma_antpos_dict

    # Check if multiple filenames were given.

    if len(filenames) > 1:
        print("Given multiple filenames. Will concatenate together.")

    uv_data = UVData()

    if track_mem_usage:
        tracemalloc.start()

    t0 = time.perf_counter()

    if verbose:
        print(f"Reading in MIR: {[str(file) for file in filenames]}")

    uv_data.read(filenames, rechunk=rechunk, run_check=run_check)

    if antfile is not None:
        if verbose:
            print("Updating antenna positions.")
        antpos = generate_sma_antpos_dict(antfile)
        uv_data.update_antenna_positions(new_positions=antpos, delta_antpos=True)

    t1 = time.perf_counter()
    if verbose:
        print(f"Total read in time: {(t1 - t0) / 60:.4f} min")

    if track_mem_usage:
        current, peak = tracemalloc.get_traced_memory()
        print(f"MIR read: Current memory during usage is {current / 10**6:.1f}MB; Peak was {peak / 10**6:.1f}MB")
        tracemalloc.stop()

    # Now write out the file! Note that setting clobber=True will
    # overwrite any identically named measurement set.
    if verbose:
        print("Writing out MS.")

    if track_mem_usage:
        tracemalloc.start()

    # Grab the basename

    if len(filenames) == 1:
        filename_base = filenames[0].name
    else:
        # I'm assuming this is always 1 track and the date is therefore the same.
        date_str = filenames[0].name.split("_")[0]
        time_str = filenames[0].name.split("_")[1]
        for filename in filenames[1:]:
            this_time_str = filename.name.split("_")[1]
            time_str += f"_{this_time_str}"
        filename_base = f"{date_str}_{time_str}"

    out_filename = outputdir / filename_base
    uv_data.write_ms(f"{out_filename}_bin{rechunk}.ms", clobber=True,
                     run_check=run_check)
    t2 = time.perf_counter()
    if verbose:
        print(f"Done write in: {(t2 - t1) / 60:.4f} min")

    if verbose:
        print(f"Read time: {(t1 - t0) / 60:.4f} min")
        print(f"Write time: {(t2 - t1) / 60:.4f} min")
        print(f"Total time: {(t2 - t0) / 60:.4f} min")

    if track_mem_usage:
        current, peak = tracemalloc.get_traced_memory()
        print(f"MS write: Current memory during usage is {current / 10**6:.1f}MB; Peak was {peak / 10**6:.1f}MB")
        tracemalloc.stop()


if __name__ == "__main__":

    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description='Convert MIR to MS.')

    parser.add_argument('-f', '--filename',
                        type=str,
                        help='MIR filename (comma separated when reading multiple files)')

    parser.add_argument('-r', '--rechunk',
                        type=int,
                        help='rechunk factor (integer)')

    parser.add_argument('-o', '--outpath',
                        type=str,
                        help='Output directory for MS')

    parser.add_argument('-n', '--nocheck', action='store_false')

    parser.add_argument('-a', '--antfile',
                        type=str,
                        help="Updated antennas file (for baseline corrections)",
                        default=None)


    args = parser.parse_args()

    print(f"Is check enabled? {args.nocheck}")

    convert_mir_to_ms([Path(filename.strip()) for filename in args.filename.split(",")],
                      rechunk=args.rechunk,
                      outputdir=Path(args.outpath),
                      run_check=args.nocheck,
                      antfile=args.antfile if args.antfile is None else Path(args.antfile),
                      )

