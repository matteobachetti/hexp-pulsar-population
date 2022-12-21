import sys
import numpy as np
from astropy.io import fits
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astropy.table import Table
from regions import CircleSkyRegion, Regions
import matplotlib.pyplot as plt

from stingray import EventList, Powerspectrum
from stingray.pulse.pulsar import htest
from stingray.pulse.search import plot_profile
from stingray.stats import z2_n_logprobability, equivalent_gaussian_Nsigma_from_logp


def process_data(evt_file, catalog_file, region_size, src_ids):
    root = evt_file.replace(".fits", "")
    all_regions = None
    with (fits.open(catalog_file) as catalog_hdul, fits.open(evt_file) as evt_hdul):
        data = catalog_hdul["SRC_CAT"].data
        # Access the data.
        Table(data)
        event_times = evt_hdul["EVENTS"].data["TIME"]
        event_ra = evt_hdul["EVENTS"].data["RA"]
        event_dec = evt_hdul["EVENTS"].data["DEC"]
        # event_src_ids = evt_hdul["EVENTS"].data["SRC_ID"][:,0]
        # print(event_src_ids)
        event_skycoords = SkyCoord(evt_hdul["EVENTS"].data["RA"], evt_hdul["EVENTS"].data["DEC"], unit=(u.deg, u.deg))

        plt.figure("Events")
        plt.scatter(event_ra, event_dec, s=0.01)
        plt.gca().set_aspect('equal', 'box')
        # plt.xlim([(sky_center.ra - delta_ra / 2).value, (sky_center.ra + delta_ra / 2).value])
        # plt.ylim([(sky_center.dec - delta_dec / 2).value, (sky_center.dec + delta_dec / 2).value])
        plt.xlabel("RA")
        plt.ylabel("Dec")

        for src_id in src_ids:
            good = np.where(data["SRC_ID"] == int(src_id))
            if not np.any(good):
                continue
            row = data[good][0]
            source_ra = row["RA"]
            source_dec = row["DEC"]
            timing_ext = "".join(row["TIMING"]).rstrip("1] ").lstrip("[ ").rstrip(",")

            try:
                timing_header = catalog_hdul[timing_ext].header
            except KeyError:
                print(f"Source {src_id} has no timing information")
                continue

            sky_center = SkyCoord(source_ra * u.deg, source_dec * u.deg)
            sky_radius = Angle(region_size, 'arcsec')

            good = np.abs(event_skycoords.separation(sky_center)) < sky_radius
            times = event_times[good]
            src_ra = event_ra[good]
            src_dec = event_dec[good]

            # all_source_events = event_src_ids == src_id
            # these_source_events = good & all_source_events

            Nev = np.count_nonzero(good)

            period = timing_header["PERIOD"]

            print(f"Source {src_id}: {Nev} events, {period:g} s ({1/period:g} Hz)")
            print(f"    RA={source_ra}, Dec={source_dec}")

            if period < 0.005:
                continue

            if Nev <= 10:
                continue

            sample_time = max(min(period / 20, 0.1), 0.0005)

            phases = times / period
            phases -= np.floor(phases)

            M, H = htest(phases, datatype="events")

            logp = z2_n_logprobability(H - 4 + 4 * M, M, ntrial=int((times[1] - times[0]) / sample_time))
            nsigma = equivalent_gaussian_Nsigma_from_logp(logp)
            print(f"    H={H} (M={M}), {nsigma:.2f} sigma")

            # ev = EventList(times)
            # lc = ev.to_lc(dt=sample_time)   #Lightcurve.make_lightcurve(times, dt=sample_time, skip_checks=True)
            # pds = Powerspectrum(lc, norm="leahy")

            # print(f"     Max Power={np.max(pds.power)}")
            # good = pds.power > 35
            # great = pds.power > 55
            label = ""

            # if np.count_nonzero(good) > 1 or np.any(great):
            if nsigma > 5:
                plt.figure("Events")
                plt.scatter(src_ra, src_dec, s=1, label=f"Src {src_id}")
                label = "_DETECTION"

            fig = plt.figure()
            plt.title(f"Source {src_id}")
            # plt.plot(pds.freq, pds.power)
            profile, bins=np.histogram(phases, bins=np.linspace(0, 1, max(M * 4, 16) + 1))
            phase_center_bins = (bins[1:] + bins[:-1]) / 2
            plot_profile(phase_center_bins, profile, ax=plt.gca())
            plt.savefig(f"{root}_src{src_id}{label}.jpg")
            plt.close(fig)
            res = Table(
                {"Phase": np.concatenate((phase_center_bins, phase_center_bins + 1)),
                 "Profile": np.concatenate((profile, profile)) / np.mean(profile)})
            res.meta.update({"H": H, "M": M, "sigma": nsigma, "period": period})
            res.write(f"{root}_src{src_id}{label}_profile.hdf5", overwrite=True)

            # pds.to_astropy_table().write(f"{root}_src{src_id}{label}_pds.hdf5", overwrite=True)
            region = CircleSkyRegion(sky_center, sky_radius)
            region.meta["label"] = f"{src_id}"
            region.write(f"{root}_src{src_id}.reg", overwrite=True)
            if all_regions is None:
                all_regions = Regions([region])
            else:
                all_regions.append(region)

        plt.figure("Events")
        plt.legend()
        all_regions.write(f"{root}_all_sources.reg", overwrite=True)
        plt.savefig(f"{root}_all_detections.jpg")
        plt.show()


catalog_file = sys.argv[1]
evt_file = sys.argv[2]
region_size = float(sys.argv[3])
if len(sys.argv) > 4:
    src_ids = sys.argv[4:]
else:
    src_ids = np.arange(76)

process_data(evt_file, catalog_file, region_size, src_ids)
