"""
Simulate a population of ULX pulsars in NGC253, with an aperiodically
variable thermal spectrum and a hard cutoff PL pulsed spectrum.

Possible improvements:
+ Better physical spectral parameters, allowing various classes of sources
  (e.g. magnetars with thermal spectra)
+ Better code design: isolate work on different sources, allow for manual
  modification before merging simput file


"""

import sys
import os
import shutil
import subprocess as sp

import tqdm
import random

import numpy as np
from scipy.stats import vonmises
from astropy.io import fits
from astropy.table import Table

import matplotlib.pyplot as plt

from xspec import AllModels, Model, Xset


def generate_list_of_sources(fname):
    """Read a Chandra catalogue and return a list of dictionaries."""
    srcs = []
    if fname == "test":
        srcs.append({"Name": "0000+0000", "Flux": 1e-12, "RA": 0, "Dec": 0, "nH": 1})
        return srcs

    table = Table.read(fname)
    table.sort("FH8", reverse=True)
    for row in table:
        src = {}
        src["Name"] = row["NAME"]
        src["RA"] = row["RA"]
        src["Dec"] = row["DEC"]
        src["nH"] = row["NH"]
        src["Flux"] = row["FH8"]
        srcs.append(src)

    return srcs


def generate_pulsed_spectrum(fname, nH=None):
    """Write an Xspec file with a cutoff powerlaw spectrum."""
    model_str = "tbabs * cutoffpl"
    if nH is None:
        nH = 10**np.random.uniform(-1, 2)  # times 10^22
    gamma = np.random.uniform(-2, 0)
    ecut = np.random.uniform(1, 20)

    model_par = [nH, gamma, ecut, 1]
    AllModels.clear()
    m1 = Model(model_str)
    for i, par in enumerate(model_par):
        m1(i + 1).values = par

    # outfname = fname.replace(".xcm", "")
    if os.path.exists(fname):
        os.remove(fname)
    Xset.save(fname, "m")


def generate_thermal_spectrum(fname, nH=None):
    """Write an Xspec file with a thermal-like spectrum."""
    model_str = "tbabs * diskpbb"
    if nH is None:
        nH = 10**np.random.uniform(-1, 2)  # times 10^22

    tin = np.random.uniform(0.1, 0.7)
    p = np.random.uniform(0.5, 1)

    model_par = [nH, tin, p, 1]
    AllModels.clear()
    m1 = Model(model_str)
    for i, par in enumerate(model_par):
        m1(i + 1).values = par

    # outfname = fname.replace(".xcm", "")
    if os.path.exists(fname):
        os.remove(fname)
    Xset.save(fname, "m")


def generate_pulsar_lc(fname, force_pars=None, nbins=None):
    """Add a pulse profile to a simput file."""
    dummy_file = f"f_{np.random.random()}_{np.random.random()}.txt"

    pars = {
        "extname": "TIMING",
        "LCFile": dummy_file,
        "MJDREF": 55000.0
    }
    if force_pars is not None:
        pars.update(force_pars)

    # This is kappa parameter of the Von Mises distribution
    kappa = np.random.uniform(1, 20)

    if nbins is None:
        nbins = max(20, int(kappa * 8))

    pf = np.random.uniform(0, 1)

    rv = vonmises(kappa)

    phases = np.linspace(0, 1, nbins + 1)
    profile = rv.pdf((phases - 0.5) * 2 * np.pi)
    # Normalize the profile from 0 to 1
    profile -= profile.min()
    profile /= profile.max()

    profile = profile * pf + (1 - pf)

    np.savetxt(dummy_file, np.array([phases, profile]).T)

    cmd = f"simputlc {fname} " + " ".join([f"{key}={val}" for key, val in pars.items()])
    sp.check_call(cmd.split())

    # Update file with pulse ephemeris
    if os.path.exists(fname):
        with fits.open(fname) as hdul:
            hdu = hdul[pars["extname"]]
            cols = hdu.data.columns
            cols.change_name("TIME", "PHASE")
            hdu.header["PHASE0"] = np.random.uniform(0, 1)
            hdu.header["PERIOD"] = 10**np.random.uniform(-3, 2)
            # hdu.header["DPERIOD"] = 10**np.random.uniform(-16, -7)

            hdul.writeto(fname, overwrite=True)

    os.remove(dummy_file)


def generate_psd(fname, force_pars=None):
    """Generate s SIMPUT extension with a complicated PSD"""
    pars = {
        "extname": "TIMING",
        "PSDnpt": 10000,
        "PSDfmin": 1e-8,
        "PSDfmax": 1e4,
        "LFQ": np.random.uniform(0.1, 10.),
        "LFrms": np.random.uniform(0.01, 0.4),
        "HBOf": np.random.uniform(0.1, 10.),
        "HBOQ": np.random.uniform(2, 10.),
        "HBOrms": np.random.uniform(0.01, 0.1),
        "Q1f": np.random.uniform(10, 100.),
        "Q1Q": np.random.uniform(2, 10.),
        "Q1rms": np.random.uniform(0.01, 0.1),
        "Q2f": np.random.uniform(100, 500.),
        "Q2Q": np.random.uniform(2, 10.),
        "Q2rms": np.random.uniform(0.01, 0.1),
        "Q3f": np.random.uniform(500, 1300.),
        "Q3Q": np.random.uniform(2, 10.),
        "Q3rms": np.random.uniform(0.01, 0.1),
    }
    if force_pars is not None:
        pars = pars.update(force_pars)

    cmd = f"simputpsd {fname} " + " ".join([f"{key}={val}" for key, val in pars.items()])
    sp.check_call(cmd.split())


def _generate_pulsed_simput(pulse_file, simpar, src_id, pulse_xspecfile, pulsed_flux, pul_ext, nH=None):
    generate_pulsed_spectrum(pulse_xspecfile, nH=nH)
    pulse_pars = {"Simput":pulse_file, "Src_ID": src_id, "Src_Name": name, "srcFlux":f"{pulsed_flux:.2e}", "XSPECFile": pulse_xspecfile}
    pulse_pars.update(simpar)

    cmd = f"simputfile " + " ".join(f"{key}={val}" for key, val in pulse_pars.items())
    sp.check_call(cmd.split())
    os.remove(pulse_xspecfile)
    generate_pulsar_lc(pulse_file, force_pars={"extname": pul_ext})


def _generate_unpulsed_simput(var_file, simpar, src_id, var_xspecfile, var_flux, var_ext, nH=None, add_var=False):
    generate_thermal_spectrum(var_xspecfile, nH=nH)
    var_pars = {"Simput":var_file, "Src_ID": src_id + 1, "Src_Name": name, "srcFlux":f"{var_flux:.2e}", "XSPECFile": var_xspecfile}
    var_pars.update(simpar)
    cmd = f"simputfile " + " ".join(f"{key}={val}" for key, val in var_pars.items())
    sp.check_call(cmd.split())
    os.remove(var_xspecfile)
    if add_var:
        generate_psd(var_file, force_pars={"extname": var_ext})
        # os.remove(var_file)


def _generate_new_source(flux, nH, ra, dec, name, source_file, src_id=0, add_var=False, emin=0.5, emax=7.0, elow=0.01, eup=300, lower_flux_limit=1e-16):
    pul_ext = f"PUL{i}"
    pulse_file = f"{name}_pulse.fits"
    pulse_xspecfile = f"{name}_pulse.xcm"
    balance_pulsations_to_var = np.random.uniform(0, 1)
    var_xspecfile = f"{name}_var.xcm"
    var_file = f"{name}_var.fits"
    var_ext = f"VAR{i}"

    pulsed_flux = flux * balance_pulsations_to_var
    var_flux = flux * (1 - balance_pulsations_to_var)
    if flux < lower_flux_limit:
        return None

    if pulsed_flux < lower_flux_limit:
        # If the pulsed flux is too low, push it a little to reach the lower flux limit
        pulsed_flux = lower_flux_limit
        var_flux = var_flux - lower_flux_limit

    simpar = {"RA": ra, "Dec": dec, "Emin":emin, "Emax":emax, "Elow": elow, "Eup": eup, "clobber": "yes", "chatter": 1}

    valid_sources = []
    _generate_pulsed_simput(pulse_file, simpar, src_id, pulse_xspecfile, pulsed_flux, pul_ext, nH=nH)
    valid_sources.append(pulse_file)

    if var_flux >= lower_flux_limit:
        _generate_unpulsed_simput(var_file, simpar, src_id, var_xspecfile, var_flux, var_ext, nH=nH, add_var=add_var)
        valid_sources.append(var_file)

    if len(valid_sources) == 2:
        sp.check_call(f"simputmerge {valid_sources[0]},{valid_sources[1]} {source_file} FetchExtensions=yes".split())
    else:
        shutil.copyfile(valid_sources[0], source_file)

    return source_file


def merge_all_files(all_files, output_fname):
   sp.check_call(f"simputmerge {','.join(all_files)} {output_fname} FetchExtensions=yes".split())


input_fname = sys.argv[1]
output_fname = sys.argv[2]
diffuse_file = diffuse_xcm_file = None
diffuse_source_flux = None
n_max_sources = 70

if len(sys.argv) > 3:
    diffuse_file = sys.argv[3]
    diffuse_xcm_file = sys.argv[4]
    diffuse_source_flux = float(sys.argv[5])
    with fits.open(diffuse_file) as hdul:
        for suffix in ["_PNT", "_NOM", "TARG", "_OBJ", ""]:
            if "RA" + suffix in hdul[0].header:
                break
        else:
            raise ValueError("Pointing information not found in the header")

        ra = hdul[0].header["RA" + suffix]
        dec = hdul[0].header["DEC" + suffix]
        elow = hdul[0].header["ENERGYLO"]
        ehigh = hdul[0].header["ENERGYHI"]

    sp.check_call(f"simputfile RA={ra} Dec={dec} srcFlux={diffuse_source_flux} XSPECFile={diffuse_xcm_file} Simput=diffuse.fits Image={diffuse_file} Emin={elow} Emax={ehigh} Elow=0.1 Eup=100 clobber=yes chatter=0".split())

# If set to True, add a PSD to the spectrum
add_var = False
sources = generate_list_of_sources(input_fname)

all_files = ["diffuse.fits"]
for i, src in enumerate(sources):
    # try:
    n_src = i * 2
    flux = src["Flux"]
    nH = src["nH"]
    if src["nH"] is np.ma.masked:
        nH = 10**np.random.uniform(-1, 2)
    ra = src["RA"]
    dec = src["Dec"]
    name = "J" + src["Name"]

    source_file = _generate_new_source(flux, nH, ra, dec, name, source_file=f"{name}.fits", src_id=n_src, emin=0.5, emax=7.0, elow=0.01, eup=300, add_var=add_var)

    if source_file is not None:
        all_files.append(source_file)
    if i >= n_max_sources:
        break

merge_all_files(all_files, output_fname)