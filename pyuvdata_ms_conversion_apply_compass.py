
import time
import tracemalloc

def convert_mir_to_ms(filenames,
                      compass_filename,
                      rechunk=1,
                      outputdir='./',
                      verbose=True,
                      run_check=True,
                      track_mem_usage=True):


    from pyuvdata import UVData

    import pyuvdata

    print(pyuvdata.__version__)

    # Check if multiple filenames were given.

    if len(filenames) > 1:
        print("Given multiple filenames. Will concatenate together.")

    uv_data = UVData()

    if track_mem_usage:
        tracemalloc.start()

    t0 = time.perf_counter()

    for ii, filename in enumerate(filenames):

        if verbose:
            print(f"Reading in MIR: {filename}")

        t1 = time.perf_counter()
        if ii == 0:

            uv_data.read(filename, rechunk=rechunk,
                         run_check=run_check,
                         compass_soln=compass_filename)
        else:
            uv_data_part = UVData()
            uv_data_part.read(filename, rechunk=rechunk,
                              run_check=run_check,
                              compass_soln=compass_filename)

            uv_data = uv_data + uv_data_part

        t2 = time.perf_counter()
        if verbose:
            print(f"Done read in: {(t2 - t1) / 60:.4f} min")

    t3 = time.perf_counter()
    if verbose:
        print(f"Total read in time: {(t3 - t0) / 60:.4f} min")

    if track_mem_usage:
        current, peak = tracemalloc.get_traced_memory()
        print(f"MIR read: Current memory during usage is {current / 10**6:.1f}MB; Peak was {peak / 10**6:.1f}MB")
        tracemalloc.stop()

    # Note that this step is needed temporarily to work around a limitation
    # in CASA w/ coordinate frames, and will not be required in the future
    # for entry in uv_data.phase_center_catalog.keys():
    #     uv_data.phase_center_catalog[entry]["cat_epoch"] = 2000.0

    # TODO: ask to propagate this back in uvdata. Remove spacing in antenna names
    uv_data.antenna_names = [ant.replace(" ", "") for ant in uv_data.antenna_names]

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
    t4 = time.perf_counter()
    if verbose:
        print(f"Done write in: {(t4 - t3) / 60:.4f} min")

    if verbose:
        print(f"Read time: {(t3 - t0) / 60:.4f} min")
        print(f"Write time: {(t4 - t3) / 60:.4f} min")
        print(f"Total time: {(t4 - t0) / 60:.4f} min")

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
                        help='MIR filename')

    parser.add_argument('-s', '--compsolution',
                        type=str,
                        help='COMPASS solutions')

    parser.add_argument('-r', '--rechunk',
                        type=int,
                        help='rechunk factor (integer)')

    parser.add_argument('-o', '--outpath',
                        type=str,
                        help='Output directory for MS')

    parser.add_argument('-n', '--nocheck', action='store_false')

    args = parser.parse_args()

    print(f"Is check enabled? {args.nocheck}")

    print(f"Running on {args.filename}")

    print(args)

    if args.filename is None:
        raise FileNotFoundError("No filename given.")

    convert_mir_to_ms([Path(filename.strip()) for filename in args.filename.split(",")],
                      Path(args.compsolution),
                      rechunk=args.rechunk,
                      outputdir=Path(args.outpath),
                      run_check=args.nocheck,
                      )

