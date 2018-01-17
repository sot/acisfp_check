#!/usr/bin/env python

"""
========================
dpa_check
========================

This code generates backstop load review outputs for checking the ACIS
focal plane temperature: FP_TEMP11. It also generates FP_TEMP11 model 
validation plots comparing predicted values to telemetry for the 
previous three weeks.
"""
from __future__ import print_function

# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import glob
import re
from Ska.Matplotlib import pointpair, \
    cxctime2plotdate
import Ska.engarchive.fetch_sci as fetch
from Chandra.Time import DateTime, date2secs
from astropy.io import ascii
from collections import defaultdict
import numpy as np
import xija
from acis_thermal_check import \
    ACISThermalCheck, \
    calc_off_nom_rolls, \
    get_options, \
    make_state_builder, \
    mylog
from acis_thermal_check.utils import \
    plot_two
import os
import sys
from kadi import events

#
# Import ACIS-specific observation extraction, filtering 
# and attribute support routines.
#
from .acis_obs import ObsidFindFilter

model_path = os.path.abspath(os.path.dirname(__file__))
default_nopref_list = os.path.join(model_path, "FPS_NoPref.txt")

#
# INIT
#
MSID = {"acisfp": "FPTEMP"}
YELLOW = {"acisfp": -50.0}
MARGIN = {"acisfp": 2.5}

# This is the cutoff temperature for any FPTEMP sensitive observation
# if the FP temp goes above this number, and the obswervation is sensitive to
# the focal plane temperature, it has to be flagged
FP_TEMP_SENSITIVE = {"acisfp": -118.7}
FP_TEMP_MINIMUM = {"acisfp": -118.7}

# This is the new maximum temperature for all ACIS-S observations (4/26/16)
ACIS_S_RED = {"acisfp": -112.0}

# ACIS-I max temperatures remain at -114 deg. C
ACIS_I_RED = {"acisfp": -114.0}

VALIDATION_LIMITS = {'PITCH': [(1, 3.0), (99, 3.0)],
                     'TSCPOS': [(1, 2.5), (99, 2.5)]
                     }

HIST_LIMIT = [(-120.0, -112.0)]

URL = "file:///home/gregg/git/xija/models/acisfp"


def calc_model(model_spec, states, start, stop, T_acisfp=None,
               T_acisfp_times=None, dh_heater=None, dh_heater_times=None):
    """
    Create and run the Thermal Model for the Focal Plane temperature.

    Given: Model name (some string)
           Commanded States collected by make_week_predict
           start time
           Stop Time
           T_acisfp
           T_acisfp_times
    """
    # Start by creating the basic modeling framework for a XIJA Thermal Model
    # Give it some name, start and stop time and the name of the JSON file
    # --------------
    # THERMAL_MODEL
    # --------------
    model = xija.ThermalModel('acisfp', start=start, stop=stop, model_spec=model_spec)

    # create a numpy array of the start and stop times in the 
    # commanded states array fetched by make_week_predict
    times = np.array([states['tstart'], states['tstop']])

    # Now set any data values for the components of your model
    # What you have to push in manually are:
    #       any states information like vid_board or ccd count
    #       any pseudo-MSID's such as 1cbat (because the node does not reflect the MSID)
    #       any single value initializations you think ought to be made.
    #         - e.g. fptemp in this case since it's what you are looking for.

    # For each item in the Commanded States data structure which matters to us,
    # insert the values in the commanded states data structure into the model:
    # 
    # Telemetry doesn't have to be pushed in - the model handles that.But items in the states
    # array have to be manually shoved in.
    #
    # pitch comes from the telemetry
    model.comp['eclipse'].set_data(False)
    model.comp['sim_z'].set_data(states['simpos'], times)
    model.comp['fptemp'].set_data(T_acisfp, T_acisfp_times)

    if "roll" in model.comp:
        model.comp['roll'].set_data(calc_off_nom_rolls(states), times)

    for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
        model.comp[name].set_data(states[name], times)

    # model.comp Not in xija documentation
    model.comp['dh_heater'].set_data(dh_heater, dh_heater_times)

    #  "orbitephem0_x","orbitephem0_y","orbitephem0_z" are not in Commanded
    #  states  but they are in telemetry
    # We have to manually insert the aoattqt<x> valued because some items 
    # are sampled on 5 minute intervals and some are not.
    for i in range(1, 5):
        name = 'aoattqt{}'.format(i)
        state_name = 'q{}'.format(i)
        model.comp[name].set_data(states[state_name], times)

    # Get ephemeris from eng archive
    for axis in "xyz":
        name = 'orbitephem0_{}'.format(axis)
        msid = fetch.Msid(name, model.tstart - 2000, model.tstop + 2000)
        model.comp[name].set_data(msid.vals, msid.times)

    # Set some initial values. Some of these are superfluous. You do this because some
    # of these values may not be set at the actual start time. The telemetry might not have
    # reached that point
    model.comp['dpa_power'].set_data(0.0)
    model.comp['1cbat'].set_data(-53.0)
    model.comp['sim_px'].set_data(-120.0)

    # Create the model
    model.make()
    model.calc()

    return model

class ACISFPCheck(ACISThermalCheck):

    def calc_model_wrapper(self, model_spec, states, tstart, tstop, state0=None):
        if state0 is None:
            start_msid = None
            dh_heater = None
            dh_heater_times = None
        else:
            start_msid = state0[self.msid]
            htrbfn = os.path.join(self.bsdir, 'dahtbon_history.rdb')
            mylog.info('Reading file of dahtrb commands from file %s' % htrbfn)
            htrb = ascii.read(htrbfn, format='rdb')
            dh_heater_times = date2secs(htrb['time'])
            dh_heater = htrb['dahtbon'].astype(bool)
        return self.calc_model(model_spec, states, tstart, tstop, T_acisfp=start_msid,
                               T_acisfp_times=None, dh_heater=dh_heater, 
                               dh_heater_times=dh_heater_times)

    def make_week_predict(self, tstart, tstop, tlm, T_init, model_spec,
                          outdir):

        mylog.info('Calculating %s thermal model' % self.name.upper())

        # The first step is to build a list of all the perigee passages.
        # We will get those from the relevant CRM pad time file 
        # (e.g. DO12143_CRM_Pad.txt) inside the bsdir directory
        # Each line is either an inbound  or outbound CTI
        #
        # The reason we are doing this is because we want to draw vertical 
        # lines denoting each perigee passage on the plots
        #
        # Open the file
        crm_file_path = glob.glob(self.bsdir + "/*CRM*")[0]
        crm_file = open(crm_file_path, 'r')

        alines = crm_file.readlines()

        idx = None
        # Keep reading until you hit the last header line which is all "*"'s
        for i, aline in enumerate(alines):
            if len(aline) > 0 and aline[0] == "*":
                idx = i+1
                break

        if idx is None:
            raise RuntimeError("Couldn't find the end of the CRM Pad Time file header!")

        # Found the last line of the header. Start processing Perigee Passages
        # initialize the resultant Perigee Passage list to empty
        # This is the product of this section of code.
        perigee_passages = []

        # While there are still lines to be read
        for aline in alines[idx:]:
            # create an empty Peri. Passage instance location
            passage = []

            # split the CRM Pad Time file line read in and extract the
            # relevant information
            splitline = aline.split()
            passage.append(splitline[6])  # Radzone entry/exit
            passage.append(splitline[9])  # Perigee Passage time

            # append this passage to the passages list
            perigee_passages.append(passage)

        # Done with the CRM Pad Time file - close it
        crm_file.close()

        # Read in the FP Sensitive Nopref file and form nopref array from it.
        nopref_array = process_nopref_list(self.fps_nopref)

        # Get commanded states and set initial temperature
        states, state0 = self.get_states(tlm, T_init)

        # calc_model_wrapper actually does the model calculation by running
        # model-specific code.
        model = self.calc_model_wrapper(model_spec, states, state0['tstart'],
                                        tstop, state0=state0)

        # Now that we have run the model, we gather the perigee passages
        # that occur from the beginning of the model run up to the start
        # of the load from kadi
        rzs = events.rad_zones.filter(state0['tstart'], tstart)
        perigee_passages = [[rz.start, rz.perigee] for rz in rzs] + \
                           [[rz.stop, rz.perigee] for rz in rzs] + \
                           perigee_passages

        # Make the limit check plots and data files
        plt.rc("axes", labelsize=10, titlesize=12)
        plt.rc("xtick", labelsize=10)
        plt.rc("ytick", labelsize=10)
        temps = {self.name: model.comp[self.msid].mvals}

        # obs_with_sensitivity contains all ACIS and all CTI observations 
        # and has had the sensitivity boolean added.
        plots, obs_with_sensitivity = self.make_prediction_plots(outdir, states, 
                                                                 model.times, temps, 
                                                                 tstart, perigee_passages, 
                                                                 nopref_array)

        viols = self.make_prediction_viols(states, model.times, temps, tstart, 
                                           obs_with_sensitivity, nopref_array)

        # write_states writes the commanded states to states.dat
        self.write_states(outdir, states)
        # write_temps writes the temperatures to temperatures.dat
        self.write_temps(outdir, model.times, temps)

        return dict(states=states, times=model.times, temps=temps,
                    plots=plots, ACIS_I_viols=viols[0], ACIS_S_viols=viols[1],
                    cti_viols=viols[2], fp_sens_viols=viols[3])

    def driver(self, args, state_builder):
        """
        The main interface to all of ACISFPCheck's functions.
        This method must be called by the particular thermal model
        implementation to actually run the code and make the webpage.

        Parameters
        ----------
        args : ArgumentParser arguments
            The command-line options object, which has the options
            attached to it as attributes
        state_builder : StateBuilder object
            The StateBuilder object used to construct commanded states
        """
        self.state_builder = state_builder
        self.fps_nopref = args.fps_nopref

        proc = self._setup_proc_and_logger(args)

        is_weekly_load = args.backstop_file is not None
        tstart, tstop, tnow = self._determine_times(args.run_start,
                                                    is_weekly_load)

        # Get the telemetry values which will be used
        # for prediction and validation
        tlm = self.get_telem_values(min(tstart, tnow), days=args.days)

        # make predictions on a backstop file if defined
        if args.backstop_file is not None:
            pred = self.make_week_predict(tstart, tstop, tlm, args.T_init,
                                          args.model_spec, args.outdir)
        else:
            pred = defaultdict(lambda: None)

        # Validation
        # Make the validation plots
        plots_validation = self.make_validation_plots(tlm, args.model_spec,
                                                      args.outdir, args.run_start)

        # Determine violations of temperature validation
        valid_viols = self.make_validation_viols(plots_validation)

        # if you found some violations....
        if len(valid_viols) > 0:
            # generate daily plot url if outdir in expected year/day format
            daymatch = re.match('.*(\d{4})/(\d{3})', args.outdir)
            if self.bsdir is None and daymatch:
                url = os.path.join(URL, daymatch.group(1), daymatch.group(2))
                mylog.info('validation warning(s) at %s' % url)
            else:
                mylog.info('validation warning(s) in output at %s' % args.outdir)

        # Write everything to the web page.
        # First, write the reStructuredText file.

        # Set up the context for the reST file
        context = {'bsdir': self.bsdir,
                   'plots': pred["plots"],
                   'valid_viols': valid_viols,
                   'proc': proc,
                   'plots_validation': plots_validation,
                   'ACIS_I_viols': pred["ACIS_I_viols"],
                   'ACIS_S_viols': pred["ACIS_S_viols"],
                   'fp_sens_viols': pred["fp_sens_viols"],
                   'cti_viols': pred["cti_viols"]}

        template_path = os.path.join(model_path, 'templates')

        self.write_index_rst(self.bsdir, args.outdir, context, 
                             template_path=template_path)

        # Second, convert reST to HTML
        self.rst_to_html(args.outdir, proc)

        return

    def make_prediction_viols(self, states, times, temps, load_start,
                              obs_with_sensitivity, nopref_array):
        """
        Find limit violations where predicted temperature is above the
        red minus margin.

        MSID is a global 

        obs_with_sensitivity contains all ACIS and CTI observations 
        and they have had FP sensitivity boolean added. In other words it's
        All ACIS and ECS runs.

        We will create a list of CTI-ONLY runs, and a list of all
        ACIS science runs without CTI runs. These two lists will
        be used to assess the categories of violations:

            1) Any ACIS-I observation that violates the -114 red limit 
               is a violation and a load killer
                 - science_viols

            2) Any ACIS-S observation that violates the -112 red limit 
               is a violation and a load killer
                 - science_viols

            3) Any ACIS FP TEMP sensitive obs that gets warmer than -118.7 
               results in a "Preferences Not Met" indicator.
                 - fp_sense_viols

            3) Any CTI run that violates the -114 RED limit needs to be
               tracked and is NOT a load killer
                 - cti_viols

        """
        mylog.info('\nMAKE VIOLS Checking for limit violations in ' + 
                   str(len(states)) + ' states and\n ' + 
                   str(len(obs_with_sensitivity)) + 
                   " total science observations")

        # create an instance of ObsidFindFilter()
        eandf = ObsidFindFilter()

        # ------------------------------------------------------
        #   Create subsets of all the observations
        # ------------------------------------------------------
        # Just the CTI runs
        cti_only_obs = eandf.cti_only_filter(obs_with_sensitivity)
        # Now divide out observations by ACIS-S and ACIS-I
        ACIS_S_obs = eandf.get_all_specific_instrument(obs_with_sensitivity, "ACIS-S")
        ACIS_I_obs = eandf.get_all_specific_instrument(obs_with_sensitivity, "ACIS-I")

        # ACIS SCIENCE observations only  - no HRC; no CTI
        non_cti_obs = eandf.cti_filter(obs_with_sensitivity)

        # ACIS SCIENCE OBS which are sensitive to FP TEMP
        fp_sens_only_obs = eandf.fp_sens_filter(non_cti_obs)

        # Now if there is an Obsid in the FP Sense list that is ALSO in the
        # No Pref list - remove that Obsid from the FP Sense list:
        fp_sense_without_noprefs = []

        nopref_list = nopref_array['obsid'].tolist()

        for eachobs in fp_sens_only_obs:
            if str(eandf.get_obsid(eachobs)) not in nopref_list:
                fp_sense_without_noprefs.append(eachobs)

        temp = temps[self.name]

        # ------------------------------------
        #  CTI-ONLY, -118.7 violation check
        # ------------------------------------
        mylog.info('\n\nFP SENSITIVE -118.7 CTI ONLY violations')
        # Collect any -118.7C violations of CTI runs. These are not
        # load killers but need to be reported

        plan_limit = FP_TEMP_SENSITIVE[self.name]
        cti_viols = search_obsids_for_viols(self.msid, self.name, plan_limit, cti_only_obs, 
                                            temp, times, load_start)

        # ------------------------------------------------------------
        #  FP TEMP sensitive observations; -118.7 violation check
        #     These are not load killers
        # ------------------------------------------------------------
        mylog.info('\n\nFP SENSITIVE -118.7 SCIENCE ONLY violations')
        # Set the limit for those observations that are sensitive to the FP Temp
        plan_limit = FP_TEMP_SENSITIVE[self.name]

        fp_sens_viols = search_obsids_for_viols(self.msid, self.name, plan_limit, 
                                                fp_sense_without_noprefs, temp, times,
                                                load_start)

        # --------------------------------------------------------------
        #  ACIS-S - Collect any -112C violations of any non-CTI ACIS-S science run. 
        #  These are load killers
        # --------------------------------------------------------------
        # 
        mylog.info('\n\n ACIS-S -112 SCIENCE ONLY violations')

        # Set the limit 
        plan_limit = ACIS_S_RED[self.name]
        ACIS_S_viols = search_obsids_for_viols(self.msid, self.name, plan_limit, 
                                               ACIS_S_obs, temp, times, load_start)

        # --------------------------------------------------------------
        #  ACIS-I - Collect any -114C violations of any non-CTI ACIS science run. 
        #  These are load killers
        # --------------------------------------------------------------
        # 
        mylog.info('\n\n ACIS-I -114 SCIENCE ONLY violations')

        # set the planning limit to the -114 C Red limit for ACIS-I observations
        plan_limit = ACIS_I_RED[self.name]

        # Create the violation data structure.
        ACIS_I_viols = search_obsids_for_viols(self.msid, self.name, plan_limit, 
                                               ACIS_I_obs, temp, times, load_start)

        return ACIS_I_viols, ACIS_S_viols, cti_viols, fp_sens_viols

    def get_histogram_mask(self, tlm, limit):
        """
        This method determines which values of telemetry
        should be used to construct the temperature 
        histogram plots, using limits provided by the 
        calling program to mask the array via a logical
        operation.

        The implementation here in ACISFPCheck is to plot
        values which fall between a lower and an upper
        limit.

        Parameters
        ----------
        tlm : NumPy record array
            NumPy record array of telemetry
        limit : array of floats
            The limit or limits to use in the masking.
        """
        return (tlm[self.msid] >= limit[0]) & (tlm[self.msid] <= limit[1])

    def make_prediction_plots(self, outdir, states, times, temps, tstart, 
                              perigee_passages, nopref_array):
        """
        Make output plots.

        :param outdir: the directory to which the products are written
        :param states: commanded states
        :param times: time stamps (sec) for temperature arrays
        :param temps: dict of temperatures
        :param tstart: load start time
        :rtype: dict of review information including plot file names
    
        This function assumes that ACIS Ops LR has been run and that the directory 
        is populated with
        """
    
        # Next we need to find all the ACIS-S observations within the start/stop
        # times so that we can paint those on the plots as well. We will get
        # those from the commanded states data structure called "states" 
        # 
        # Create an instance of the ObsidFindFilter class. This class provides
        # methods to extract obsid intervals from the commanded states based 
        # upon ACIS definitions and considerations. It also provides
        # various methods to filter the interval set based upon pitch range, 
        # number of ccd's, filter out CTI observations, and a range of exposure 
        # times.
        extract_and_filter = ObsidFindFilter()
    
        # extract the OBSID's from the commanded states. NOTE: this contains all
        # observations including CTI runs and HRC observations
        observation_intervals = extract_and_filter.find_obsid_intervals(states, None)
    
        # Filter out any HRC science observations BUT keep ACIS CTI observations
        acis_and_cti_obs = extract_and_filter.hrc_science_obs_filter(observation_intervals)
    
        # Ok so now you have all the ACIS observations collected. Also,
        # they have been identified by ObsidFindFilter as to who is in the focal plane.
        # Some apps, like this one, care about FP_TEMP sensitivity. Some do not. 
        # Since we do, then checking that and assigning a sensitivity must be done
        # 
        # Open the sensitive observation list file, which is found in the LR 
        # directory,
        # read each line, extract the OBSID and add that to a list.
        sensefile = open(self.bsdir + '/fp_sensitive.txt', 'r')
    
        # The list_of_sensitive_obs is the list of all FP TEMP sensitive 
        # observations extracted from the file in the load review directory
        list_of_sensitive_obs = []
    
        # Get the list of FP_TEMP sensitive observations
        for eachline in sensefile.readlines()[1:]:
            # Extract the OBSID from each line; the obsid is in the second
            # column of this line. Append it to the list of FP_TEMP sensitive
            # observations
            #
            # NOTE: The obsid here is a STRING
            list_of_sensitive_obs.append(eachline.split()[1])
        # Done with the file - close it
        sensefile.close()
    
        # Now that you have the latest list of temperature sensitive OBSID's, 
        # run through each observation and append either "*FP SENS*" or
        # "NOT FP SENS" to the end of each observation. 
        #
        # Create an empty observation list which will hold the results. This
        # list contains all ACIS and all CTI observations and will have the 
        # sensitivity boolean added.
        obs_with_sensitivity = []
    
        # Now run through the observation list attribute of the ObsidFindFilter class
        for eachobservation in acis_and_cti_obs:
            # Pull the obsid from the observation and turn it into a string
    
            obsid = str(extract_and_filter.get_obsid(eachobservation))
            # See if it's in the sensitive list. If so, indicate whether or
            # not this observation is FP Senstive in the new list. This will be
            # used later in make_prediction_viols to catch violations.
            if obsid in list_of_sensitive_obs:
                eachobservation.append(True)
            else:
                eachobservation.append(False)
    
            obs_with_sensitivity.append(eachobservation)
    
        #
        # create an empty dictionary called plots to contain the returned 
        # figures, axes 1  and axes 2 of the plot_two call
        plots = {}
    
        # Start time of loads being reviewed expressed in units for plotdate()
        load_start = cxctime2plotdate([tstart])[0]
        # Value for left side of plots
        plot_start = max(load_start-2.0, cxctime2plotdate([times[0]])[0])

        #
        # Make  plots of FPTEMP and pitch vs time 
        #
        #
        # For each MSID in the loop statement, make a plot from -120 to 20
        #   - kind of a superfluous loop in this case.
        #
        for fig_id, msid in enumerate((self.name,)):
            # -----------------------------------------------------
            #
            #   PLOT 1 -  fptemp_11 plt with ylim from -120 to -90
            # 
            # -----------------------------------------------------
            plots[msid+"_1"] = plot_two(fig_id=fig_id + 1,
                                        x=times,
                                        y=temps[msid],
                                        x2=pointpair(states['tstart'], states['tstop']),
                                        y2=pointpair(states['pitch']),
                                        title=MSID[msid] + " (ACIS-I obs. in red; ACIS-S in green)",
                                        xlabel='Date',
                                        ylabel='Temperature (C)',
                                        ylabel2='Pitch (deg)',
                                        ylim=(-120, -90),
                                        ylim2=(40, 180),
                                        figsize=(14, 7),
                                        )
            plots[msid+"_1"]['ax'].axhline(ACIS_I_RED[msid], linestyle='--', color='red',
                                           linewidth=2.0)
            # Add a vertical line to mark the start time of the load
            plots[msid+"_1"]['ax'].axvline(load_start, linestyle='-', color='g',
                                           linewidth=2.0)

            # The next line is to ensure that the width of the axes
            # of all the weekly prediction plots are the same--used
            # later down in the function
            w1, h1 = plots[msid+"_1"]['fig'].get_size_inches()

            #
            # Now plot any perigee passages that occur between xmin and xmax
            # for eachpassage in perigee_passages:
            paint_perigee(perigee_passages, states, plots, msid+"_1")
    
            # Now draw horizontal lines on the plot running from start to stop
            # and label them with the Obsid
            ypos = -110.0
            endcapstart = -111.0
            endcapstop = -109.0
            textypos = -108.0
            fontsize = 12
            draw_obsids(extract_and_filter, obs_with_sensitivity, nopref_array,
                        plots, msid+"_1", ypos, endcapstart, endcapstop, textypos, fontsize)
            # Set the left limit of the plot to be -2 days before the load start
            plots[msid+"_1"]['ax'].set_xlim(plot_start, None)
            # Build the file name and output the plot to a file
            filename = MSID[msid].lower() + 'M120toM90.png'
            outfile = os.path.join(outdir, filename)
            mylog.info('Writing plot file %s' % outfile)
            plots[msid+"_1"]['fig'].savefig(outfile)
            plots[msid+"_1"]['filename'] = filename
    
            # ------------------------------------------------------
            #
            #   PLOT 2 -  fptemp_11 plt with ylim from -120 to -119
            # 
            # ------------------------------------------------------
            #
            plots[msid+"_2"] = plot_two(fig_id=fig_id + 1,
                                        x=times,
                                        y=temps[msid],
                                        x2=pointpair(states['tstart'], states['tstop']),
                                        y2=pointpair(states['pitch']),
                                        title=MSID[msid] + " (ACIS-I obs. in red; ACIS-S in green)",
                                        xlabel='Date',
                                        ylabel='Temperature (C)',
                                        ylabel2='Pitch (deg)',
                                        ylim=(-120, -119),
                                        ylim2=(40, 180),
                                        figsize=(14, 7),
                                        )
            # Add a vertical line to mark the start time of the load
            plots[msid+"_2"]['ax'].axvline(load_start, linestyle='-', color='g',
                                           linewidth=2.0)
            # The next several lines ensure that the width of the axes
            # of all the weekly prediction plots are the same.
            w2, h2 = plots[msid+"_2"]['fig'].get_size_inches()
            lm = plots[msid+"_1"]['fig'].subplotpars.left*w1/w2
            rm = plots[msid+"_1"]['fig'].subplotpars.right*w1/w2
            plots[msid+"_2"]['fig'].subplots_adjust(left=lm, right=rm)

            #
            # Now plot any perigee passages that occur between xmin and xmax
            #
            paint_perigee(perigee_passages, states, plots, msid+"_2")
    
            # Now draw horizontal lines on the plot running form start to stop
            # and label them with the Obsid
            ypos = -119.35
            endcapstart = ypos + 0.05
            endcapstop = ypos - 0.05
            textypos = ypos + 0.05
            fontsize = 9
            draw_obsids(extract_and_filter, obs_with_sensitivity, nopref_array,
                        plots, msid+"_2", ypos, endcapstart, endcapstop, textypos, fontsize)
            # Set the left limit of the plot to be -2 days before the load start
            plots[msid+"_2"]['ax'].set_xlim(plot_start, None)
            # Build the file name and output the file
            filename = MSID[msid].lower() + 'M120toM119.png'
            outfile = os.path.join(outdir, filename)
            mylog.info('   Writing plot file %s' % outfile)
            plots[msid+"_2"]['fig'].savefig(outfile)
            plots[msid+"_2"]['filename'] = filename
    
            # ------------------------------------------------------
            #
            #   PLOT 3 -  fptemp_11 plt with ylim from -120 to -112
            # 
            # ------------------------------------------------------
            #
            plots[msid+"_3"] = plot_two(fig_id=fig_id + 1,
                                        x=times,
                                        y=temps[msid],
                                        x2=pointpair(states['tstart'], states['tstop']),
                                        y2=pointpair(states['pitch']),
                                        title=MSID[msid] + " (ACIS-I obs. in red; ACIS-S in green)",
                                        xlabel='Date',
                                        ylabel='Temperature (C)',
                                        ylabel2='Pitch (deg)',
                                        ylim=(-120, -111.5),
                                        ylim2=(40, 180),
                                        figsize=(14, 7),
                                        )
            # Add a vertical line to mark the start time of the load
            plots[msid+"_3"]['ax'].axvline(load_start, linestyle='-', color='g',
                                           linewidth=2.0)
            # The next several lines ensure that the width of the axes                                             
            # of all the weekly prediction plots are the same.                                                      
            w2, h2 = plots[msid+"_3"]['fig'].get_size_inches()
            lm = plots[msid+"_1"]['fig'].subplotpars.left*w1/w2
            rm = plots[msid+"_1"]['fig'].subplotpars.right*w1/w2
            plots[msid+"_3"]['fig'].subplots_adjust(left=lm, right=rm)

            #
            # Now plot any perigee passages that occur between xmin and xmax
            #
            paint_perigee(perigee_passages, states, plots, msid+"_3")
    
            # Now draw horizontal lines on the plot running from start to stop
            # and label them with the Obsid
            ypos = -116
            endcapstart = -116.2
            endcapstop = -115.8
            textypos = -115.7
            fontsize = 9
    
            draw_obsids(extract_and_filter, obs_with_sensitivity, nopref_array,
                        plots, msid+"_3", ypos, endcapstart, endcapstop, textypos, fontsize)
    
            # Draw a horizontal line indicating the FP Sensitive Observation Cut off
            plots[msid+"_3"]['ax'].axhline(FP_TEMP_SENSITIVE[msid], linestyle='--', color='red', linewidth=2.0)
            # Draw a horizontal line showing the ACIS-I -114 deg. C cutoff
            plots[msid+"_3"]['ax'].axhline(ACIS_I_RED[msid], linestyle='--', color='purple', linewidth=1.0)
            # Draw a horizontal line showing the ACIS-S -112 deg. C cutoff
            plots[msid+"_3"]['ax'].axhline(ACIS_S_RED[msid], linestyle='--', color='blue', linewidth=1.0)
            # Set the left limit of the plot to be -2 days before the load start
            plots[msid+"_3"]['ax'].set_xlim(plot_start, None)

            # The next several lines ensure that the width of the axes                                                
            # of all the weekly prediction plots are the same.                                                       
            w2, h2 = plots[msid+"_3"]['fig'].get_size_inches()
            lm = plots[msid+"_1"]['fig'].subplotpars.left*w1/w2
            rm = plots[msid+"_1"]['fig'].subplotpars.right*w1/w2
            plots[msid+"_3"]['fig'].subplots_adjust(left=lm, right=rm)

            # Build the file name and output the file
            filename = MSID[msid].lower() + 'M120toM112.png'
            outfile = os.path.join(outdir, filename)
            mylog.info('   Writing plot file %s' % outfile)
            plots[msid+"_3"]['fig'].savefig(outfile)
            plots[msid+"_3"]['filename'] = filename
    
        # end of for fig_id, msid in enumerate(('fptemp',)):
    
        # Now create the plot of ACIS CCD number and Sim-Z position
        plots['pow_sim'] = plot_two(
            fig_id=fig_id + 1,
            title='ACIS CCDs and SIM-Z position',
            xlabel='Date',
            x=pointpair(states['tstart'], states['tstop']),
            y=pointpair(states['ccd_count']),
            ylabel='CCD_COUNT',
            ylim=(-0.1, 6.1),
            x2=pointpair(states['tstart'], states['tstop']),
            y2=pointpair(states['simpos']),
            ylabel2='SIM-Z (steps)',
            ylim2=(-105000, 105000),
        )
        # Add a vertical line to mark the start time of the load
        plots['pow_sim']['ax'].axvline(load_start, linestyle='-', color='g',
                                       linewidth=2.0)
        # Set the left limit of the plot to be -2 days before the load start
        plots['pow_sim']['ax'].set_xlim(plot_start, None)
        # The next several lines ensure that the width of the axes                                              
        # of all the weekly prediction plots are the same.                                                   
        w2, h2 = plots["pow_sim"]['fig'].get_size_inches()
        lm = plots[msid+"_1"]['fig'].subplotpars.left*w1/w2
        rm = plots[msid+"_1"]['fig'].subplotpars.right*w1/w2
        plots["pow_sim"]['fig'].subplots_adjust(left=lm, right=rm)

        filename = 'pow_sim.png'
        outfile = os.path.join(outdir, filename)
        mylog.info('   Writing plot file %s' % outfile)
        plots['pow_sim']['fig'].savefig(outfile)
        plots['pow_sim']['filename'] = filename

        return plots, obs_with_sensitivity


def search_obsids_for_viols(msid, name, plan_limit, observations, temp, times,
                            load_start):
    """
    Given a planning limit and a list of observations, find those time intervals
    where the temp gets warmer than the planning limit and identify which 
    observations (if any) include part or all of those intervals.
    """

    # create an instance of ObsidFindFilter()
    eandf = ObsidFindFilter()

    viols_list = defaultdict(list)

    bad = np.concatenate([[False], temp >= plan_limit, [False]])

    # changes is a list of lists. Each sublist is a tuple which
    # contains indices into the times list. 0 = start times
    # and 1 = stop time
    changes = np.flatnonzero(bad[1:] != bad[:-1]).reshape(-1, 2)

    # Add any obs violations to the viols_list
    for change in changes:
        tstart = times[change[0]]
        tstop = times[change[1] - 1]

        # Only report violations which occur after the load being
        # reviewed starts.
        in_load = tstart > load_start or \
                  (tstart < load_start < tstop)

        if not in_load:
            continue

        # find the observations that contains all or part of this time interval
        #  add this to the violations list "viols[msid]"
        #
        # First create an empty obsid list. This list represents all the
        # observations that contain the violation represented by this change
        # (if any)
        obsid_list = ''
 
        # See if any  observation that contains this time interval
        # If it is, add this to the violations list "viols_list[msid]"
        for eachobs in observations:
            # Get the observation tstart and tstop times, and obsid
            obs_tstart = eandf.get_tstart(eachobs)
            obs_tstop = eandf.get_tstop(eachobs)

            # If either tstart is inside the obs tstart/tstop
            # OR
            #    tstop is inside the obs tstart/tstop
            if obs_tstop >= tstart >= obs_tstart or \
               obs_tstop >= tstop >= obs_tstart or \
               (tstart <= obs_tstart and tstop >= obs_tstop):
                # Fetch the obsid for this observation and append to list
                obsid_list = obsid_list + ' '+ str(eandf.get_obsid(eachobs))

        # If obsid_list is not empty, then create the violation
        if obsid_list != '':
            # Figure out the max temp for the
            # Then create the violation
            viol = {'datestart': DateTime(times[change[0]]).date,
                    'datestop': DateTime(times[change[1] - 1]).date,
                    'maxtemp': temp[change[0]:change[1]].max(),
                    'obsid': obsid_list
                   }

            # .........and then add it to the violations list
            viols_list[msid].append(viol)

            mylog.info('   VIOLATION: %s  exceeds planning limit of %.2f '
                        'degC from %s to %s'
                        % (MSID[name], plan_limit, viol['datestart'],
                        viol['datestop']))

    # Finished - return the violations list
    return viols_list

#----------------------------------------------------------------------
#
#   paint_perigee
#
#----------------------------------------------------------------------
def paint_perigee(perigee_passages, states, plots, msid):
    """
    This function draws vertical dahsed lines for EEF, Perigee and XEF
    events in the load.EEF and XEF lines are black; Perigee is red.

    You supply the list of perigee passage events which are:
        Radzone Start/Stop time
        Perigee Passage time

        The states you created in main

        The dictionary of plots you created

        The MSID (in this case FP_TEMP) used to access the dictionary
    """
    #
    # Now plot any perigee passages that occur between xmin and xmax
    for eachpassage in perigee_passages:
        # The index [1] item is always the Perigee Passage time. Draw that line in red
        # If this line is between tstart and tstop then it needs to be drawn 
        # on the plot. otherwise ignore
        if states['tstop'][-1] >= DateTime(eachpassage[0]).secs >= states['tstart'][0]:
            # Have to convert this time into the new x axis time scale necessitated by SKA
            xpos = cxctime2plotdate([DateTime(eachpassage[0]).secs])

            # now plot the line.
            plots[msid]['ax'].vlines(xpos, -120, 20, linestyle=':', color='red', linewidth=2.0)

            # Plot the perigee passage time so long as it was specified in the CTI_report file
            if eachpassage[1] != "Not-within-load":
                perigee_time = cxctime2plotdate([DateTime(eachpassage[1]).secs])
                plots[msid]['ax'].vlines(perigee_time, -120, 20, linestyle=':', 
                                         color='black', linewidth=2.0)

def draw_obsids(extract_and_filter, 
                obs_with_sensitivity, 
                nopref_array,
                plots,
                msid, 
                ypos, 
                endcapstart, 
                endcapstop, 
                textypos, 
                fontsize):
    """
    This functiion draws visual indicators across the top of the plot showing
    which observations are ACIS; whether they are ACIS-I (red) or ACIS-S (green)
    when they start and stop, and whether or not any observation is sensitive to the
    focal plane temperature.  The list of observations sensitive to the focal plane
    is found by reading the fp_sensitive.dat file that is located in each LR
    directory and is created by the LR script.

    No CTI measurements are indicated - only science runs.

    The caller supplies:
               Options from the Command line supplied by the user at runtime
               The instance of the ObsidFindFilter() class created 
               nopref rec array
               The plot dictionary
               The MSID used to index into the plot dictinary (superfluous but required)
               The position on the Y axis you'd like these indicators to appear
               The Y position of the bottom of the end caps
               The Y position of the top of the end caps
               The starting position of the OBSID number text
               The font size
    """
    # Now run through the observation list attribute of the ObsidFindFilter class
    for eachobservation in obs_with_sensitivity:
        # extract the obsid

        obsid = str(extract_and_filter.get_obsid(eachobservation))

        # Color all ACIS-S observations green; all ACIS-I 
        # observations red
        if eachobservation[extract_and_filter.in_focal_plane] == "ACIS-I":
            color = 'red'
        else:
            color = 'green'

        # Add the sensitivity text if the observation was found to be FP TEMP
        # sensitive
        #
        # If the observation is FP sensitive in the first place............
        if eachobservation[extract_and_filter.is_fp_sensitive]:
            # extract the obsid for convenience
            this_obsid = extract_and_filter.get_obsid(eachobservation)

            # But if it's also in the nopref list AND the upcased CandS_status entry is "NO PREF"
            where_words = np.where(nopref_array['obsid'] == str(this_obsid))
            if str(this_obsid) in nopref_array['obsid'] and \
               nopref_array['CandS_status'][where_words[0][0]].upper()[:7] == 'NO_PREF':
                color = 'purple'
                obsid = obsid + ' * NO PREF *'
            else:
                obsid = obsid + ' * FP SENS *'

        # Convert the start and stop times into the Ska-required format
        obs_start = cxctime2plotdate([extract_and_filter.get_tstart(eachobservation)])
        obs_stop = cxctime2plotdate([extract_and_filter.get_tstop(eachobservation)])

        if eachobservation[extract_and_filter.in_focal_plane].startswith("ACIS-"):
            # For each ACIS Obsid, draw a horizontal line to show 
            # its start and stop
            plots[msid]['ax'].hlines(ypos, 
                                     obs_start, 
                                     obs_stop, 
                                     linestyle='-', 
                                     color=color, 
                                     linewidth=2.0)

            # Plot vertical end caps for each obsid to visually show start/stop
            plots[msid]['ax'].vlines(obs_start, 
                                     endcapstart, 
                                     endcapstop, 
                                     color=color, 
                                     linewidth=2.0)
            plots[msid]['ax'].vlines(obs_stop, 
                                     endcapstart, 
                                     endcapstop, 
                                     color=color, 
                                     linewidth=2.0)

            # Now print the obsid in the middle of the time span, 
            # above the line, and rotate 90 degrees. 

            # Now plot the obsid.
            plots[msid]['ax'].text(obs_start + ((obs_stop - obs_start)/2), 
                                   textypos, 
                                   obsid,  
                                   color = color, 
                                   va='bottom', 
                                   ma='left', 
                                   rotation = 90, 
                                   fontsize = fontsize )


#------------------------------------------------------------------------------
#
#   process_nopref_list - read and store the list of observations which
#                         prefer a Cold and Stable focal plane, but
#                         which have had that desire waived.
#
#                          Input:  nopref file specification
#                         Output:  nopref_array (numpy array)
#
#------------------------------------------------------------------------------
def process_nopref_list(filespec=default_nopref_list):
    # Create the dtype for the nopref list rec array
    nopref_dtype = [('obsid', '|S10'), ('Seq_no', '|S10'), 
                    ('prop', '|S10'), ('CandS_status', '|S15')]

    # Create an empty nopref array
    nopref_array = np.array([], dtype=nopref_dtype)

    # Open the nopref list file
    nopreflist = open(filespec, "r")

    # For each line in the file......
    for line in nopreflist.readlines():

        # split the line out into it's constituent parts
        splitline = line.encode().split()

        # If the line is NOT a comment (i.e. does not start with a "#")
        if splitline[0][0] != '#':

            # Create the row
            one_entry = np.array(splitline, dtype=nopref_dtype)
            # Append this row onto the array
            nopref_array = np.append(nopref_array, one_entry, axis=0)

    # Done with the nopref list file
    nopreflist.close()

    # Now return the nopref array
    return nopref_array 

def main():
    opts = [("fps_nopref", {"default": default_nopref_list,
             "help": "Full path to the FP sensitive nopref file"})]
    args = get_options("acisfp", model_path, opts=opts)
    state_builder = make_state_builder(args.state_builder, args)
    acisfp_check = ACISFPCheck("fptemp", "acisfp", MSID, YELLOW,
                               MARGIN, VALIDATION_LIMITS, HIST_LIMIT, 
                               calc_model, other_telem=['1dahtbon'],
                               other_map={'1dahtbon': 'dh_heater', 
                                          "fptemp_11": "fptemp"})
    try:
        acisfp_check.driver(args, state_builder)
    except Exception as msg:
        if args.traceback:
            raise
        else:
            print("ERROR:", msg)
            sys.exit(1)

if __name__ == '__main__':
    main()
