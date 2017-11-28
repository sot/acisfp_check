#!/usr/bin/env python

"""
========================
acisfp_check
========================

This code generates backstop load review outputs for checking the ACIS
Focal Plane Temperature: FP_TEMP.  It also generates FP_TEMP model validation
plots comparing predicted values to telemetry for the previous three weeks.
"""

import sys
import os
import glob
import logging
from pprint import pformat
import re
import time
import shutil
import pickle

import numpy as np
import Ska.DBI
import Ska.Table
import Ska.Numpy
import Ska.engarchive.fetch_sci as fetch
from Chandra.Time import DateTime
from Chandra.Time import date2secs

import Chandra.cmd_states as cmd_states
# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import Ska.Matplotlib

import xija

#
# Import ACIS-specific observation extraction, filtering 
# and attribute support routines.
#
import ACISobs

#
# INIT
#
MSID = dict(fptemp = 'FPTEMP')

YELLOW = dict(fptemp = -50.0)

MARGIN = dict(fptemp = 2.5)

# This is the cutoff temperature for any FPTEMP sensitive observation
# if the FP temp goes above this number, and the obswervation is sensitive to
# the focal plane temperature, it has to be flagged
FP_TEMP_SENSITIVE = dict(fptemp = -118.7)
FP_TEMP_MINIMUM = dict(fptemp = -118.7)

# This is the new maximum temperature for all ACIS-S observations (4/26/16)
ACIS_S_RED = dict(fptemp = -112.0)

# ACIS-I max temperatures remain at -114 deg. C
ACIS_I_RED = dict(fptemp = -114.0)

VALIDATION_LIMITS = {'PITCH': [(1, 3.0),
                                  (99, 3.0)],
                     'TSCPOS': [(1, 2.5),
                                (99, 2.5)]
                     }

TASK_DATA = os.path.dirname(__file__)

#
# Initialize the logger. Set up a file for log output
#
logger = logging.getLogger('acisfp_check')

#
URL = "file:///home/gregg/git/xija/models/acisfp"

_versionfile = os.path.join(os.path.dirname(__file__), 'VERSION')
VERSION = open(_versionfile).read().strip()
#------------------------------------------------------------------------------
#
#   get_option - extract the command line arguments
#
#------------------------------------------------------------------------------
def get_options():
    """
    Set command line argument defaults, then extract any command line arguments
    specified by the user.

    NOTE: This parser used optparse. According to the PYTHON documents, This is
          a deprecated module:

          http://docs.python.org/library/optparse.html

          Deprecated since version 2.7: The optparse module is deprecated 
          and will not be developed further; development will continue 
          with the argparse module.

    """
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--outdir",
                      default="out",
                      help="Output directory")

    parser.add_option("--oflsdir",
                       help="Load products OFLS directory")

    parser.add_option("--model-spec",
                      default=os.path.join(TASK_DATA, 'acisfp_spec.json'),
                       help="ACIS FOCAL PLANE model specification file")

    parser.add_option("--days",
                      type='float',
                      default=21.0,
                      help="Days of validation data (days)")

    parser.add_option("--run-start",
                      help="Reference time to replace run start time "
                           "for regression testing")

    parser.add_option("--traceback",
                      default=True,
                      help='Enable tracebacks')

    parser.add_option("--verbose",
                      type='int',
                      default=1,
                      help="Verbosity (0=quiet, 1=normal, 2=debug)")

    parser.add_option("--ccd-count",
                      type='int',
                      default=6,
                      help="Initial number of CCDs (default=6)")

    parser.add_option("--fep-count",
                      type='int',
                      default=6,
                      help="Initial number of FEPs (default=6)")

    parser.add_option("--vid-board",
                      type='int',
                      default=1,
                      help="Initial state of ACIS vid_board (default=1)")

    parser.add_option("--clocking",
                      type='int',
                      default=1,
                      help="Initial state of ACIS clocking (default=1)")

    parser.add_option("--simpos",
                      default=75616,
                      type='float',
                      help="Starting SIM-Z position (steps)")

    parser.add_option("--pitch",
                      default=150.0,
                      type='float',
                      help="Starting pitch (deg)")

    parser.add_option("--T_acisfp",
                      type='float',
                      help="Starting FPTEMP temperature (degC)")

    # adding dh_heater
    parser.add_option("--dh_heater",
                      type='int',
                      default=0,
                      help="Starting Detector Housing Heater state")

    # adding FP Sensitive NoPref File
    parser.add_option("--fps_nopref",
                      default="/data/acis/LoadReviews/script/fp_temp_predictor/FPS_NoPref.txt",
                      help="full path to the FP sensitive nopref file")

    parser.add_option("--version",
                      action='store_true',
                      help="Print version")

    opt, args = parser.parse_args()
    return opt, args

###############################################################################
#
#   MAIN
#
###############################################################################
"""
Main execution routine - input is opt - which is the list of user-specified 
command line options.  Main calls are: get_bs_cmds 
                                       get_tlm_values
                                       make_week_predict
                                       make_validation_plots
                                       make_validation_viols
                                       write_index_rst
                                       rst_to_html
"""
def main(opt):

    if not os.path.exists(opt.outdir):
        os.mkdir(opt.outdir)

    config_logging(opt.outdir, opt.verbose)

    # Store info relevant to processing for use in outputs
    # FP_TEMP limits set in globals above.
    # ctime is the present time

    proc = dict(run_user=os.environ['USER'],
                run_time=time.ctime(),
                errors=[],
                )

    logger.info('##############################'
                '#######################################')
    logger.info('# acisfp_check.py run at %s by %s'
                % (proc['run_time'], proc['run_user']))
    logger.info('# acisfp_check version = {}'.format(VERSION))
    logger.info('# model_spec file = %s' % os.path.abspath(opt.model_spec))
    logger.info('# Focal Plane Sensitivity NoPref file Path: %s', opt.fps_nopref)
    logger.info('###############################'
                '######################################\n')

    logger.info('Command line options:\n%s\n' % pformat(opt.__dict__))

    # Connect to the SKA.DBI database (NEED TO USE aca_read)
    logger.info('Connecting to SKA DBI database to get cmd_states')
    db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read',
                     database='aca')

    tnow = DateTime(opt.run_start).secs


    #
    # GET_BS_CMDS 
    # If you have an ofls directory, get tstart, tstop, 
    # and the commands from the backstop file in that directory.  
    # Note that tstart and tstop are taken from the first and last 
    # element in the list. This DOES NOT necessarily coincide with the 
    # 
    if opt.oflsdir is not None:
        # Get tstart, tstop, commands from backstop file in opt.oflsdir
        #------------
        # GET_BS_CMDS - Get the backstop file commands for this week
        #------------
        logger.info('\nObtaining the commands from the prediction load backstop file')
        bs_cmds = get_bs_cmds(opt.oflsdir)
        tstart = bs_cmds[0]['time']
        tstop = bs_cmds[-1]['time']

        logger.info('   BACSTOP FILE TSTART is: %s' % DateTime(tstart).date)
        logger.info('   BACSTOP FILE TSTOP is:  %s\n' % DateTime(tstop).date)

        proc.update(dict(datestart=DateTime(tstart).date,
                         datestop=DateTime(tstop).date))
    else:
        tstart = tnow


    #------------------
    # GET_TELEM_VALUES Get temperature telemetry for 3 weeks prior to 
    #                  min(tstart, NOW)    
    #-----------------
    # Fetch the telemetry data for the specified MSID's. Note that the 
    # undocumented purpose of name_map is to allow the user to set a 
    # numpy structured array column header to a user specified value RATHER
    # than the original msid name.
    # Changed 10/2015
    tlm = get_telem_values(min(tstart, tnow),
                           ['fptemp_11',
                            'sim_z',                # dpa json
                            'aosares1',
                            'dp_dpa_power','1dahtbon'],
                           days=opt.days,
                           name_map={'sim_z': 'tscpos',
                                     'aosares1': 'pitch',
                                     'fptemp_11': 'fptemp',
                                     '1dahtbon': 'dh_heater'})

    # Convert the returned SIM position to the numbers used by
    # us ( as in -99616. for HRC-S)
    tlm['tscpos'] = tlm['tscpos'] * -397.7225924607

    #------------------
    # MAKE_WEEK_PREDICT - make predictions on oflsdir if defined
    #------------------
    # 
    if opt.oflsdir is not None:
        pred = make_week_predict(opt, tstart, tstop, bs_cmds, tlm, db)
    else:
        pred = dict(plots=None, ACIS_I_viols=None, ACIS_S_viols=None, cti_viols=None, fp_sens_viols=None, times=None, states=None, temps=None)

    #----------------------
    # MAKE_VALIDATION_PLOTS
    #----------------------

    plots_validation = make_validation_plots(opt, tlm, db)
    
    #----------------------
    # MAKE_VALIDATION_VIOLS
    #----------------------
    valid_viols = make_validation_viols(plots_validation)

    # if you found some violations....
    if len(valid_viols) > 0:
        # generate daily plot url if outdir in expected year/day format
        daymatch = re.match('.*(\d{4})/(\d{3})', opt.outdir)
        if opt.oflsdir is None and daymatch:
            url = os.path.join(URL, daymatch.group(1), daymatch.group(2))
            logger.info('validation warning(s) at %s' % url)
        else:
            logger.info('validation warning(s) in output at %s' % opt.outdir)

    #------------------
    # WRITE INDEX RST
    #------------------
    write_index_rst(opt, 
                    proc, 
                    plots_validation, 
                    valid_viols=valid_viols,
                    plots=pred['plots'], 
                    ACIS_I_viols=pred['ACIS_I_viols'], 
                    ACIS_S_viols=pred['ACIS_S_viols'], 
                    fp_sens_viols=pred['fp_sens_viols'],
                    cti_viols = pred['cti_viols'])

    rst_to_html(opt, proc)

    print "ALL DONE!"
    return dict(opt=opt, states=pred['states'], times=pred['times'],
                temps=pred['temps'], plots=pred['plots'],
                ACIS_I_viols=pred['ACIS_I_viols'], 
                ACIS_S_viols=pred['ACIS_S_viols'], 
                cti_viols=pred['cti_viols'], 
                fp_sens_viols=pred['fp_sens_viols'], 

                proc=proc,
                plots_validation=plots_validation)

    

#------------------------------------------------------------------------------
#
#   calc_model
#
#------------------------------------------------------------------------------
def calc_model(model_spec, states, start, stop, T_acisfp=None, T_acisfp_times=None, dh_heater=None, dh_heater_times=None):
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
    #--------------
    # THERMAL_MODEL
    #--------------
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

    #  
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
    for axis in ('x', 'y', 'z'):
        name = 'orbitephem0_{}'.format(axis)
        msid = fetch.Msid(name, model.tstart - 2000, model.tstop + 2000)
        model.comp[name].set_data(msid.vals, msid.times)

    # Set some initial values. Some of these are superfluous. You do this because some
    # of these values may not be set at the actual start time. The telemetry might not have
    # reached that point
    model.comp['dpa_power'].set_data(0.0)
    model.comp['1cbat'].set_data(-53.0)
    model.comp['sim_px'].set_data(-120.0)
    model.comp['fptemp'].set_data(-120.0)

    # Create the model
    model.make()
    model.calc()

    return model


#------------------------------------------------------------------------------
#
#   make_week_predict
#
#------------------------------------------------------------------------------
def make_week_predict(opt, tstart, tstop, bs_cmds, tlm, db):
    """
    Given:  opt - the options entered by the user on the command
                  line
                  tstart - from bs_cmds[0]['time']
		  tstop  - from bs_cmds[-1]['time']
		  bs_cmds - output from get_bs_cmds
		  tlm    - output from get_telem_values numpy structured array
		  db     - Ska.DBI.DBI
  
      Find the start and end times that incorporate the initial state0 
      calculated time, through the first backstop file time (if different)
      and to the end of the backstop file. Get the commanded states from 
      state0 through the end of backstop commands.

    Return: A dictionary containing:
            opt=opt, 
            states=states, 
            times=model.times, 
            temps=temps,
            plots=plots, 
            viols=viols

    Calls: calc_model
           make_check_plots
           make_viols
    """
    #
    # Ok the first step is to build a list of all the perigee passages.
    # We will get those from the relevant CRM pad time file 
    # (e.g. DO12143_CRM_Pad.txt) inside the opt.oflsdir directory
    # Each line is either an inbound  or outbound CTI
    #
    # The reason we are doing this is because we want to draw vertical 
    # lines denoting each perigee passage on the plots
    #
    # Open the file
    crm_file_path = glob.glob(opt.oflsdir+"/*CRM*")[0]
    crm_file = open(crm_file_path, 'r')

    #read the first line. This is sure to be a header line
    aline = crm_file.readline()

    # Keep reading until you hit the last header line which is all "*"'s
    while (len(aline) > 0) and (aline[0] != '*'):
        aline = crm_file.readline()

    # Found the last line of the header. Start processing Perigee Passages
    # initialize the resultant Perigee Passage list to empty
    # This is the product of this section of code.
    perigee_passages = []

    #
    # Get a Perigee Passage line
    aline = crm_file.readline()

    # While there are still lines to be read
    while (aline):
        # create an empty Peri. Passage instance location
        passage = []

        # split the CRM Pad Time file line read in and extract the 
        #relevant information
        splitline = aline.split()
        passage.append(splitline[0])   #  Event Type (EEF or XEF)
        passage.append(splitline[6])   #  CTI Start time
        passage.append(splitline[7])   #  CTI Stop time
        passage.append(splitline[9])   #  Perigee Passage time

        # append this passage to the passages list
        perigee_passages.append(passage)

        # read the next line (if there is one
        aline = crm_file.readline()

    # Done with the CRM Pad Time file - close it
    crm_file.close()

    # Read in the FP Sensitive Nopref file and form nopref array from it.
    nopref_array = process_nopref_list(opt.fps_nopref)

    # Make initial state0 from cmd line options. The CL options all 
    # have defaults.
    # NOTE - a state looks like one of the entries in TA's states.dat file
    state0 = dict((x, getattr(opt, x)) for x in ('pitch', 'simpos', 'ccd_count', 'fep_count',
                                                 'vid_board', 'clocking', 'T_acisfp', 'dh_heater'))

    # Populate the state with info from the input paramters like tstart
    # and some defaults (e.g. quarternions)
    state0.update({'tstart': tstart - 30,
                   'tstop': tstart,
                   'datestart': DateTime(tstart - 30).date,
                   'datestop': DateTime(tstart).date,
                   'q1': 0.0, 'q2': 0.0, 'q3': 0.0, 'q4': 1.0,
                   }
                  )

    # If cmd lines options were not fully specified then get state0 as last
    # cmd_state that starts within available telemetry.  Update with the
    # mean temperatures at the start of state0.
    if None in state0.values():
        state0 = cmd_states.get_state0(tlm['date'][-5], db,
                                       datepar='datestart')
        ok = ((tlm['date'] >= state0['tstart'] - 700) &
              (tlm['date'] <= state0['tstart'] + 700))

        state0.update({'T_acisfp': np.mean(tlm['fptemp'][ok])})

    # TEMPORARY HACK: core model doesn't actually support predictive
    # active heater yet.  Initial temperature determines active heater
    # state for predictions now.
    #
    # So clamp T_acisfp floor to 15.0
    if state0['T_acisfp'] < 15:
        state0['T_acisfp'] = 15.0

    logger.debug('\n    MWP - state0 at %s is\n%s' % (DateTime(state0['tstart']).date,
                                           pformat(state0)))

    # Get commanded states after end of state0 through first backstop command time
    cmds_datestart = state0['datestop']
    cmds_datestop = bs_cmds[0]['date']

    logger.info("    \n    fetchall date start and stop: %s %s " %(cmds_datestart, cmds_datestop))

    # Get timeline load segments including state0 and beyond.
    timeline_loads = db.fetchall("""SELECT * from timeline_loads
                                 WHERE datestop > '%s'
                                 and datestart < '%s'"""
                                 % (cmds_datestart, cmds_datestop))

    logger.info('    Found {} timeline_loads  after {}'.format(
                len(timeline_loads), cmds_datestart))

    # Get cmds since datestart within timeline_loads
    db_cmds = cmd_states.get_cmds(cmds_datestart, db=db, update_db=False,
                                  timeline_loads=timeline_loads)

    # Delete non-load cmds that are within the backstop time span
    # => Keep if timeline_id is not None or date < bs_cmds[0]['time']
    db_cmds = [x for x in db_cmds if (x['timeline_id'] is not None or
                                      x['time'] < bs_cmds[0]['time'])]


    logger.info('    Got %d cmds from database between %s and %s' %
                  (len(db_cmds), cmds_datestart, cmds_datestop))

    # Get the commanded states from state0 through the end of backstop commands
    states = cmd_states.get_states(state0, db_cmds + bs_cmds)

    states[-1].datestop = bs_cmds[-1]['date']
    states[-1].tstop = bs_cmds[-1]['time']


    logger.info('    Found %d commanded states from %s to %s' %
                 (len(states), states[0]['datestart'], states[-1]['datestop']))

    # October 2015 - DH Heater addition
    htrbfn='dahtbon_history.rdb'
    logger.info('    Reading file of dahtrb commands from file %s' % htrbfn)
    htrb=Ska.Table.read_ascii_table(htrbfn,headerrow=2,headertype='rdb')
    dh_heater_times=date2secs(htrb['time'])

    dh_heater=htrb['dahtbon'].astype(bool)

    #
    # CALC_MODEL - Calculate the model
    #
    model = calc_model(opt.model_spec, states, state0['tstart'], tstop,
                       state0['T_acisfp'], None, dh_heater, dh_heater_times)

    # Make the DPA limit check plots and data files
    plt.rc("axes", labelsize=10, titlesize=12)
    plt.rc("xtick", labelsize=10)
    plt.rc("ytick", labelsize=10)

    temps = {'fptemp': model.comp['fptemp'].mvals}

    #------------------------
    #
    #  call make_check_plots
    #
    #-------------------------
    #
    # obs_with_sensitivity contains all ACIS and all CTI observations 
    # and has had the sensitivity boolean added.
    plots, obs_with_sensitivity = make_check_plots(opt, states, model.times, temps, tstart, perigee_passages, nopref_array)

    #-------------------
    #
    #  call make_viols
    #
    #-------------------
    ACIS_I_viols, cti_viols, fp_sens_viols, ACIS_S_viols = make_viols(opt, states, model.times, temps, obs_with_sensitivity, nopref_array)

    # Write out the states.dat and temperatures.dat files
    write_states(opt, states)
    write_temps(opt, model.times, temps)

    return dict(opt=opt, states=states, times=model.times, temps=temps,
               plots=plots, cti_viols=cti_viols, ACIS_I_viols=ACIS_I_viols, ACIS_S_viols = ACIS_S_viols, fp_sens_viols=fp_sens_viols)
#------------------------------------------------------------------------------
#
#   make_validation_viols
#
#------------------------------------------------------------------------------
def make_validation_viols(plots_validation):
    """
    Find limit violations where MSID quantile values are outside the
    allowed range.
    """
    logger.info('\nChecking for validation violations')

    viols = []

    for plot in plots_validation:
        # 'plot' is actually a structure with plot info and stats about the
        #  plotted data for a particular MSID.  'msid' can be a real MSID
        #  (1DPAMZT) or pseudo like 'POWER'
        msid = plot['msid']

        # Make sure validation limits exist for this MSID
        if msid not in VALIDATION_LIMITS:
            continue

        # Cycle through defined quantiles (e.g. 99 for 99%) and corresponding
        # limit values for this MSID.
        for quantile, limit in VALIDATION_LIMITS[msid]:
            # Get the quantile statistic as calculated when making plots
            msid_quantile_value = float(plot['quant%02d' % quantile])

            # Check for a violation and take appropriate action
            if abs(msid_quantile_value) > limit:
                viol = {'msid': msid,
                        'value': msid_quantile_value,
                        'limit': limit,
                        'quant': quantile,
                        }
                viols.append(viol)
                logger.info('   WARNING: %s %d%% quantile value of %s exceeds '
                            'limit of %.2f' %
                            (msid, quantile, msid_quantile_value, limit))

    return viols

#------------------------------------------------------------------------------
#
#   get_bs_cmds - 
#
#------------------------------------------------------------------------------
def get_bs_cmds(oflsdir):
    """Return commands from the backstop file in opt.oflsdir.
    """
    import Ska.ParseCM
    backstop_file = globfile(os.path.join(oflsdir, 'CR*.backstop'))
    logger.info('   log - Using backstop file %s' % backstop_file)
    bs_cmds = Ska.ParseCM.read_backstop(backstop_file)
    logger.info('   log - Found %d backstop commands between %s and %s' %
                  (len(bs_cmds), bs_cmds[0]['date'], bs_cmds[-1]['date']))

    return bs_cmds
#------------------------------------------------------------------------------
#
#   get_telem_values
#
#------------------------------------------------------------------------------
def get_telem_values(tstart, msids, days=14, name_map={}):
    """
    Given: tstart: start time for telemetry (secs)
           msids: fetch msids list
           days: length of telemetry request before ``tstart``
           dt: sample time (secs)
           name_map: dict mapping msid to recarray col name

    NOTE: that the undocumented purpose of name_map is to allow the 
          user to set a numpy structured array column header to 
          a user specified value RATHER than the original msid name.
          Not every msid needs to have a mapped name.

    Fetch last ``days`` of available ``msids`` telemetry values before
    time ``tstart``.

    :returns: np recarray of requested telemetry values from fetch
     
    Using the list of msid's in msid, fetch the data from the telemetry db
    in an fetch.MSIDset call.
    """
    tstart = DateTime(tstart).secs

    #Calculate the days PRIOR to tstart to use for the start time in
    # the fetch.MSID call
    start = DateTime(tstart - days * 86400).date
    stop = DateTime(tstart).date
    logger.info('\nFetching telemetry between %s and %s' % (start, stop))

    # msidset is the results of a fetch.MSID set on the "msid's" list
    # of items to get from telemetry. It is an array, containing one list
    # which contains one Ska.engarchive.fetch.MSID object for each 
    # msid. e.g. 
    # MSIDset([('1dpamzt', <Ska.engarchive.fetch.MSID object at 0x84fc3d0>),
    #          ('fptemp_11', <Ska.engarchive.fetch.MSID object at 0x7051450>),
    #          ..........])
    #  
    #    each msid object containing "times" and "vals" just like any other
    #      e.g.  msidset['1dpamzt'].times

    msidset = fetch.MSIDset(msids, start, stop, stat='5min')
    start = max(x.times[0] for x in msidset.values())
    stop = min(x.times[-1] for x in msidset.values())
    msidset.interpolate(328.0, start, stop + 1)  # 328 for '5min' stat

    # Finished when we found at least 4 good records (20 mins)
    if len(msidset.times) < 4:
        raise ValueError('   Found no telemetry within %d days of %s'
                         % (days, str(tstart)))

    # x is each element in the msid list, 
    # name_map.get(x, x) returns name_map y contents for that key OR,
    # the x value itself if it is not a key in the dictionary    
    # So for example, msids has "fptemp_11"; name_map maps that to
    # "fptemp". So either the original msids name goes into outnames,
    # Or the name mapped to the msid name in name_map.
    #
    # outnames is a list of the values - used as headers for the array
    outnames = ['date'] + [name_map.get(x, x) for x in msids]

    # vals is a dictionary whose keys are either the original msid name OR
    # the name the user decided to map it to in name_map
    # and whose values are the list of each fetch results for that msid
    #    =      keys:             the values
    vals = {name_map.get(x, x): msidset[x].vals for x in msids}
    #   e.g. vals['pitch']
   
    # Now add a date column by taking the times out of msidset, and 
    # appending as a n array whose dictionary key is "date"
    vals['date'] = msidset.times


    # Now turn vals into a structured array, and use outnames 
    # column headers.
    #
    # NOTE: at this point, outnames and vals keys are the same list
    out = Ska.Numpy.structured_array(vals, colnames=outnames)

    return out


#------------------------------------------------------------------------------
#
#   rst_to_html
#
#------------------------------------------------------------------------------
def rst_to_html(opt, proc):
    """Run rst2html.py to render index.rst as HTML"""

    # First copy CSS files to outdir
    import Ska.Shell
    import docutils.writers.html4css1
    dirname = os.path.dirname(docutils.writers.html4css1.__file__)
    shutil.copy2(os.path.join(dirname, 'html4css1.css'), opt.outdir)

    shutil.copy2(os.path.join(TASK_DATA, 'acisfp_check.css'), opt.outdir)

    spawn = Ska.Shell.Spawn(stdout=None)
    infile = os.path.join(opt.outdir, 'index.rst')
    outfile = os.path.join(opt.outdir, 'index.html')
    status = spawn.run(['rst2html.py',
                        '--stylesheet-path={}'
                        .format(os.path.join(opt.outdir, 'acisfp_check.css')),
                        infile, outfile])
    if status != 0:
        proc['errors'].append('rst2html.py failed with status {}: see run log'
                              .format(status))
        logger.error('rst2html.py failed')
        logger.error(''.join(spawn.outlines) + '\n')

    # Remove the stupid <colgroup> field that docbook inserts.  This
    # <colgroup> prevents HTML table auto-sizing.
    del_colgroup = re.compile(r'<colgroup>.*?</colgroup>', re.DOTALL)
    outtext = del_colgroup.sub('', open(outfile).read())
    open(outfile, 'w').write(outtext)


#------------------------------------------------------------------------------
#
#   config_logging
#
#------------------------------------------------------------------------------
def config_logging(outdir, verbose):
    """Set up file and console logger.
    See http://docs.python.org/library/logging.html
              #logging-to-multiple-destinations
    """
    # Disable auto-configuration of root logger by adding a null handler.
    # This prevents other modules (e.g. Chandra.cmd_states) from generating
    # a streamhandler by just calling logging.info(..).
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass
    rootlogger = logging.getLogger()
    rootlogger.addHandler(NullHandler())

    loglevel = {0: logging.CRITICAL,
                1: logging.INFO,
                2: logging.DEBUG}.get(verbose, logging.INFO)

    logger = logging.getLogger('acisfp_check')
    logger.setLevel(loglevel)

    formatter = logging.Formatter('%(message)s')

    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(loglevel);
    logger.addHandler(console)

    filehandler = logging.FileHandler(
        filename=os.path.join(outdir, 'run.dat'), mode='w')
    filehandler.setFormatter(formatter)

    # Set the file loglevel to be at least INFO,
    # but override to DEBUG if that is requested at the
    # command line
    filehandler.setLevel(logging.INFO)
    if loglevel == logging.DEBUG:
        filehandler.setLevel(logging.DEBUG)

    logger.addHandler(filehandler)


#------------------------------------------------------------------------------
#
#   write_states
#
#------------------------------------------------------------------------------
def write_states(opt, states):
    """Write states recarray to file states.dat"""
    outfile = os.path.join(opt.outdir, 'states.dat')
    logger.info('Writing states to %s' % outfile)
    out = open(outfile, 'w')
    fmt = {'power': '%.1f',
           'pitch': '%.2f',
           'tstart': '%.2f',
           'tstop': '%.2f',
           }
    newcols = list(states.dtype.names)
    newcols.remove('T_acisfp')
    newstates = np.rec.fromarrays([states[x] for x in newcols], names=newcols)
    Ska.Numpy.pprint(newstates, fmt, out)
    out.close()


#------------------------------------------------------------------------------
#
#   write_temps
#
#------------------------------------------------------------------------------
def write_temps(opt, times, temps):
    """Write temperature predictions to file temperatures.dat"""
    outfile = os.path.join(opt.outdir, 'temperatures.dat')
    logger.info('Writing temperatures to %s' % outfile)
    T_acisfp = temps['fptemp']
    temp_recs = [(times[i], DateTime(times[i]).date, T_acisfp[i])
                 for i in xrange(len(times))]
    temp_array = np.rec.fromrecords(
        temp_recs, names=('time', 'date', 'fptemp'))

    fmt = {'fptemp': '%.2f',
           'time': '%.2f'}
    out = open(outfile, 'w')
    Ska.Numpy.pprint(temp_array, fmt, out)
    out.close()


#------------------------------------------------------------------------------
#
#   write_index_rst
#
#------------------------------------------------------------------------------
def write_index_rst(opt, 
                    proc, 
                    plots_validation, 
                    valid_viols=None,
                    plots=None, 
                    ACIS_I_viols=None, 
                    ACIS_S_viols=None, 
                    fp_sens_viols=None,
                    cti_viols=None):
    """
    Make output text (in ReST format) in opt.outdir.

    opt: a dictionary containing user supplied command line options
    proc: Process info such as username, start/stop times, OS environment
    plots_validation: returned list from make_validation_plots()
    valid_viols:      returned list from make_validation_viols()
    plots: dictionary containing prediction temp plots
    viols: dictionary containing violation temp plots
    """

    # Django setup (used for template rendering)
    import django.template
    import django.conf
    try:
        django.conf.settings.configure()
    except RuntimeError, msg:
        print msg

    outfile = os.path.join(opt.outdir, 'index.rst')
    logger.info('Writing report file %s' % outfile)

    # Right now, the method to get the fptemp plot of interest 
    # (fptempM120toM114.png) to be the main interest is pretty kludgy. 
    # I happen to make 3 different plots of FP_TEMP vs Time at 3 different Y limit
    # ranges. I placed the plot of interest last in the creation
    # order so it handily appears in the "plots" dictionary used below
    django_context = django.template.Context(
        {'opt': opt,
         'plots': plots,
         'ACIS_I_viols': ACIS_I_viols,
         'ACIS_S_viols': ACIS_S_viols,
         'fp_sens_viols': fp_sens_viols,
         'cti_viols': cti_viols,
         'valid_viols': valid_viols,
         'proc': proc,
         'plots_validation': plots_validation,
         })
    index_template_file = ('index_template.rst'
                           if opt.oflsdir else
                           'index_template_val_only.rst')
    index_template = open(os.path.join(TASK_DATA, index_template_file)).read()
    index_template = re.sub(r' %}\n', ' %}', index_template)
    template = django.template.Template(index_template)
    open(outfile, 'w').write(template.render(django_context))

#------------------------------------------------------------------------------
#
#   search_obsids_for_viols
#
#------------------------------------------------------------------------------
def search_obsids_for_viols(msid, plan_limit, observations, temp, times):
    """
    Given a planning limit and a list of observations, find those time intervals
    where the temp gets warmer than the planning limit and identify which 
    observations (if any) include part or all of those intervals.
    """

    # create an instance of ACISobs.ObsidFindFilter()
    eandf = ACISobs.ObsidFindFilter()

    viols_list = dict((x, []) for x in MSID)

    bad = np.concatenate(([False],
                          temp >= plan_limit,
                          [False]))
    # changes is a list of lists. Each sublist is a 2-ple which
    # contains indices into the times list. 0 = start times 
    # and 1 = stop time
    changes = np.flatnonzero(bad[1:] != bad[:-1]).reshape(-1, 2)

    # Add any obs violations to the viols_list
    for change in changes:
        datestart = DateTime(times[change[0]]).date
        tstart = times[change[0]]
        datestop  = DateTime(times[change[1] - 1]).date
        tstop = times[change[1] - 1]
        maxtemp = temp[change[0]:change[1]].max()

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
            obs_tstop  = eandf.get_tstop(eachobs)
            obsid = eandf.get_obsid(eachobs)

            # If either tstart is inside the obs tstart/tstop
            # OR 
            #    tstop is inside the obs tstart/tstop
            # 
            if  ( (tstart >= obs_tstart) and \
                  (tstart <= obs_tstop )) or \
                ( (tstop >= obs_tstart) and \
                  (tstop <= obs_tstop)) or \
                ( (tstart <= obs_tstart) and \
                  (tstop >= obs_tstop)):

                # Fetch the obsid for this observation and append to list
                obsid_list = obsid_list + ' '+ str(eandf.get_obsid(eachobs))

        # If obsid_list is not empty, then create the violation
        if (obsid_list != ''):
            # Figure out the max temp for the 
            # Then create the violation
            viol = {'datestart': DateTime(times[change[0]]).date,
                    'datestop': DateTime(times[change[1] - 1]).date,
                    'maxtemp': temp[change[0]:change[1]].max(),
                    'obsid': obsid_list
                   }

            # .........and then add it to the violations list
            viols_list[msid].append(viol)

            logger.info('   VIOLATION: %s  exceeds planning limit of %.2f '
                        'degC from %s to %s'
                        % (MSID[msid], plan_limit, viol['datestart'],
                        viol['datestop']))

    # Finished - return the violations list
    return viols_list




#------------------------------------------------------------------------------
#
#   make_viols
#
#------------------------------------------------------------------------------

def make_viols(opt, states, times, temps, obs_with_sensitivity, nopref_array):
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
    logger.info('\nMAKE VIOLS Checking for limit violations in '+str(len(states))+' states and\n '+ str(len(obs_with_sensitivity))+ " total science observations")

#    viols = dict((x, []) for x in MSID)

    # create an instance of ACISobs.ObsidFindFilter()
    eandf = ACISobs.ObsidFindFilter()

    #------------------------------------------------------
    #   Create subsets of all the observations
    #------------------------------------------------------
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

    # superfluous loop in this case since there is only one MSID
    for msid in MSID:
        temp = temps[msid]

        #------------------------------------
        #  CTI-ONLY, -118.7 violation check
        #------------------------------------
        logger.info('\n\nCTI ONLY -118.7 FP SENSE violations')
        # Collect any -118.7C violations of CTI runs. These are not
        # load killers but need to be reported

        plan_limit = FP_TEMP_SENSITIVE[msid]
        cti_viols = search_obsids_for_viols(msid, plan_limit, cti_only_obs, temp, times)

        #------------------------------------------------------------
        #  FP TEMP sensitive observations; -118.7 violation check
        #     These are not load killers
        #------------------------------------------------------------
        logger.info('\n\nFP SENSITIVE SCIENCE ONLY -118.7 violations')
        # Set the limit for  Thos eObservations that are sensitive to the FP Temp
        plan_limit = FP_TEMP_SENSITIVE[msid]

        #fp_sens_viols = search_obsids_for_viols(msid, plan_limit, fp_sens_only_obs, temp, times)
        fp_sens_viols = search_obsids_for_viols(msid, plan_limit, fp_sense_without_noprefs, temp, times)

        #--------------------------------------------------------------
        #  ACIS-S - Collect any -112C violations of any non-CTI ACIS-S science run. 
        #  These are load killers
        #--------------------------------------------------------------
        # 
        logger.info('\n\n-112 ACIS-S SCIENCE ONLY violations')
        # Set the limit 
        plan_limit = ACIS_S_RED[msid]
        ACIS_S_viols = search_obsids_for_viols(msid,  plan_limit, ACIS_S_obs, temp, times)
     
        #--------------------------------------------------------------
        #  ACIS-I - Collect any -114C violations of any non-CTI ACIS science run. 
        #  These are load killers
        #--------------------------------------------------------------
        # 
        logger.info('\n\n ACIS-I -114 SCIENCE ONLY violations')
        # Create the violation data structure.
        ACIS_I_viols = dict((x, []) for x in MSID)

        # set the planning limit to the -114 C Red limit for ACIS-I observations
        plan_limit = ACIS_I_RED[msid]

        ACIS_I_viols = search_obsids_for_viols(msid, plan_limit,ACIS_I_obs , temp, times)

    return(ACIS_I_viols, cti_viols, fp_sens_viols, ACIS_S_viols)

#------------------------------------------------------------------------------
#
#   plot_two
#
#------------------------------------------------------------------------------
def plot_two(fig_id, x, y, x2, y2,
             linestyle='-', linestyle2='-',
             color='blue', color2='magenta',
             ylim=None, ylim2=None,
             xlabel='', ylabel='', ylabel2='', title='',
             figsize=(7, 3.5),
             ):
    """Plot two quantities with a date x-axis"""

    xt = Ska.Matplotlib.cxctime2plotdate(x)
    fig = plt.figure(fig_id, figsize=figsize)
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot_date(xt, y, fmt='-', linestyle=linestyle, color=color)


    
    ax.set_xlim(min(xt), max(xt))
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid()

    ax2 = ax.twinx()

    xt2 = Ska.Matplotlib.cxctime2plotdate(x2)
    ax2.plot_date(xt2, y2, fmt='-', linestyle=linestyle2, color=color2)
    ax2.set_xlim(min(xt), max(xt))
    if ylim2:
        ax2.set_ylim(*ylim2)
    ax2.set_ylabel(ylabel2, color=color2)
    ax2.xaxis.set_visible(False)

    Ska.Matplotlib.set_time_ticks(ax)
    [label.set_rotation(30) for label in ax.xaxis.get_ticklabels()]
    [label.set_color(color2) for label in ax2.yaxis.get_ticklabels()]

    fig.subplots_adjust(bottom=0.22)

    return {'fig': fig, 'ax': ax, 'ax2': ax2}


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
        Event Type (EEF or XEF)
        CTI Start time
        CTI Stop time
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
        if(DateTime(eachpassage[1]).secs >= states['tstart'][0]) and \
          (DateTime(eachpassage[1]).secs <= states['tstop'][-1]):
          # Have to convert this time into the new x axis time scale necessitated by SKA
          xpos = Ska.Matplotlib.cxctime2plotdate([DateTime(eachpassage[1]).secs])

          # now plot the line.
          plots[msid]['ax'].vlines(xpos, -120, 20, linestyle=':', color='red', linewidth=2.0)

          # Plot the perigee passage time so long as it was specified in the CTI_report file
          if eachpassage[3] != "Not-within-load":
              perigee_time = Ska.Matplotlib.cxctime2plotdate([DateTime(eachpassage[3]).secs])
              plots[msid]['ax'].vlines(perigee_time, -120, 20, linestyle=':', color='black', linewidth=2.0)
          

#----------------------------------------------------------------------
#
#   draw_obsids
#
#----------------------------------------------------------------------
def draw_obsids(opt, 
                extract_and_filter, 
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
               The instance of the  ACISobs.ObsidFindFilter() class created 
               nopref rec array
               The plot dictionary
               The MSID used to index into the plot dictinary (superfluous but required)
               The position on the Y axis you'd like these indicators to appear
               The Y position of the bottom of the end caps
               The Y position of the top of the end caps
               The starting position of the OBSID number text
               The font size
    """

    # Set the color for ACIS_S observations - green
    color = 'green'


    # Now run through the observation list attribute of the ACISobs class
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
        if (eachobservation[extract_and_filter.is_fp_sensitive] == True):
            # extract the obsid for convenience
            this_obsid = extract_and_filter.get_obsid(eachobservation)

            # But if it's also in the nopref list AND the upcased CandS_status entry is " NO PREF "
            if (str(this_obsid) in nopref_array['obsid']) and \
               ( nopref_array['CandS_status'][ np.where(nopref_array['obsid'] == str(this_obsid) )[0][0] ].upper()[:7]  == 'NO_PREF' ):
                color = 'purple'
                obsid = obsid + ' * NO PREF *'
            else:
                obsid = obsid + ' * FP SENS *'


        # Convert the start and stop times into the Ska-required format
        obs_start = Ska.Matplotlib.cxctime2plotdate([extract_and_filter.get_tstart(eachobservation)])
        obs_stop = Ska.Matplotlib.cxctime2plotdate([extract_and_filter.get_tstop(eachobservation)])

        if eachobservation[extract_and_filter.in_focal_plane] == "ACIS-I" or \
           eachobservation[extract_and_filter.in_focal_plane] == "ACIS-S":
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
def process_nopref_list(filespec = "/data/acis/LoadReviews/script/fp_temp_predictor/FPS_NoPref.txt"):
    # Create the dtype for the nopref list rec array
    nopref_dtype = [('obsid', '|S10'), ('Seq_no', '|S10'), ('prop', '|S10'), ('CandS_status', '|S15')]

    # Create an empty nopref array
    nopref_array = np.array([], dtype = nopref_dtype)

    # Open the nopref list file
    nopreflist = open(filespec, "r")
    
    # For each line in the file......
    for line in nopreflist.readlines():

        # split the line out into it's constiutent parts
        splitline = line.split()

        # If the line is NOT a comment (i.e. does not start with a "#")
        if splitline[0][0] != '#':

            # Create the row
            one_entry = np.array( [ (splitline[0], splitline[1], splitline[2], splitline[3] )], dtype=nopref_dtype)
            # Append this row onto the array
            nopref_array = np.append(nopref_array, one_entry, axis = 0)

    # Done with the nopref list file
    nopreflist.close()
    
    # Now return the nopref array
    return(nopref_array)
        

#------------------------------------------------------------------------------
#
#   make_check_plots   
#
#------------------------------------------------------------------------------
def make_check_plots(opt, states, times, temps, tstart, perigee_passages, nopref_array):
    """
    Make output plots.

    :param opt: options
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
    extract_and_filter = ACISobs.ObsidFindFilter()

    # extract the OBSID's from the commanded states. NOTE: this contains all
    # observations including CTI runs and HRC observations
    observation_intervals = extract_and_filter.find_obsid_intervals(states, None)

    # Filter out any HRC science observations BUT keep ACIS CTI observations
    acis_and_cti_obs = extract_and_filter.hrc_science_obs_filter(observation_intervals)

    # Ok so now you have all the  ACIS observations collected. Also,
    # they have been identified by ACISobs as to who is in the focal plane.
    # Some apps, like this one, care about FP_TEMP sensitivity. Some do not. 
    # Since we do, then checking that and assigning a sensitivity must be done
    # 
    # Open the sensitive observation list file, which is found in the LR 
    # directory,
    # read each line, extract the OBSID and add that to a list.
    sensefile = open(opt.oflsdir+'/fp_sensitive.txt', 'r')

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
        list_of_sensitive_obs.append(eachline.split()[1] )
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

    # Now run through the observation list attribute of the ACISobs class
    for eachobservation in acis_and_cti_obs:
        # Pull the obsid from the observation and turn it into a string

        obsid = str(extract_and_filter.get_obsid(eachobservation))
        # See if it's in the sensitive list. If so, indicate whether or
        # not this observation is FP Senstive in the new list. This will be
        # used later in make_viols to catch violations.
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
    load_start = Ska.Matplotlib.cxctime2plotdate([tstart])[0]

    #
    # Make  plots of FPTEMP and pitch vs time 
    #
    #
    # For each MSID in the loop statement, make a plot from -120 to 20
    #   - kind of a superfluous loop in this case.
    #
    for fig_id, msid in enumerate(('fptemp',)):
        #-----------------------------------------------------
        #
        #   PLOT 1 -  fptemp_11 plt with ylim from -120 to -90
        # 
        #-----------------------------------------------------
        plots[msid] = plot_two(fig_id=fig_id + 1,
                               x=times,
                               y=temps[msid],
                               x2=pointpair(states['tstart'], states['tstop']),
                               y2=pointpair(states['pitch']),
                               title=MSID[msid]+" (ACIS-I obs. in red; ACIS-S in green)",
                               xlabel='Date',
                               ylabel='Temperature (C)',
                               ylabel2='Pitch (deg)',
                               ylim = (-120, -90),
                               ylim2=(40, 180),
                               figsize=(14, 7),
                               )
        plots[msid]['ax'].axhline(ACIS_I_RED[msid], linestyle='--', color='red',
                                  linewidth=2.0)

        plots[msid]['ax'].axvline(load_start, linestyle=':', color='g',
                                  linewidth=2.0)

        #
        # Now plot any perigee passages that occur between xmin and xmax
        #for eachpassage in perigee_passages:
        paint_perigee(perigee_passages, states, plots, msid)

        # Now draw horizontal lines on the plot running from start to stop
        # and label them with the Obsid
        ypos = -110.0
        endcapstart =  -111.0
        endcapstop = -109.0
        textypos = -108.0
        fontsize = 12
        draw_obsids(opt, extract_and_filter, obs_with_sensitivity, nopref_array, plots, msid, ypos, endcapstart, endcapstop, textypos, fontsize)

        # Build the file name and output the plot to a file
        filename = MSID[msid].lower() + 'M120toM90.png'
        outfile = os.path.join(opt.outdir, filename)
        logger.info('   Writing plot file %s' % outfile)
        plots[msid]['fig'].savefig(outfile)
        plots[msid]['filename'] = filename


        #------------------------------------------------------
        #
        #   PLOT 2 -  fptemp_11 plt with ylim from -120 to -119
        # 
        #------------------------------------------------------
        #
        plots[msid] = plot_two(fig_id=fig_id + 1,
                               x=times,
                               y=temps[msid],
                               x2=pointpair(states['tstart'], states['tstop']),
                               y2=pointpair(states['pitch']),
                               title=MSID[msid]+" (ACIS-I obs. in red; ACIS-S in green)",
                               xlabel='Date',
                               ylabel='Temperature (C)',
                               ylabel2='Pitch (deg)',
                               ylim = (-120, -119),
                               ylim2=(40, 180),
                               figsize=(14, 7),
                               )

        #
        # Now plot any perigee passages that occur between xmin and xmax
        #
        paint_perigee(perigee_passages, states, plots, msid)

        # Now draw horizontal lines on the plot running form start to stop
        # and label them with the Obsid
        ypos = -119.35
        endcapstart = ypos + 0.05
        endcapstop = ypos - 0.05
        textypos = ypos + 0.05
        fontsize = 9
        draw_obsids(opt, extract_and_filter, obs_with_sensitivity, nopref_array, plots, msid, ypos, endcapstart, endcapstop, textypos, fontsize)

        # Build the file name and output the file
        filename = MSID[msid].lower() + 'M120toM119.png'
        outfile = os.path.join(opt.outdir, filename)
        logger.info('   Writing plot file %s' % outfile)
        plots[msid]['fig'].savefig(outfile)
        plots[msid]['filename'] = filename

        #------------------------------------------------------
        #
        #   PLOT 3 -  fptemp_11 plt with ylim from -120 to -112
        # 
        #------------------------------------------------------
        #
        plots[msid] = plot_two(fig_id=fig_id + 1,
                               x=times,
                               y=temps[msid],
                               x2=pointpair(states['tstart'], states['tstop']),
                               y2=pointpair(states['pitch']),
                               title=MSID[msid]+" (ACIS-I obs. in red; ACIS-S in green)",
                               xlabel='Date',
                               ylabel='Temperature (C)',
                               ylabel2='Pitch (deg)',
                               ylim = (-120, -111.5),
                               ylim2=(40, 180),
                               figsize=(14, 7),
                               )
        #
        # Now plot any perigee passages that occur between xmin and xmax
        #
        paint_perigee(perigee_passages, states, plots, msid)

        # Now draw horizontal lines on the plot running from start to stop
        # and label them with the Obsid
        ypos = -116
        endcapstart =  -116.2
        endcapstop = -115.8
        textypos = -115.7
        fontsize = 9

        draw_obsids(opt, extract_and_filter, obs_with_sensitivity, nopref_array, plots, msid, ypos, endcapstart, endcapstop, textypos, fontsize)

        # Draw a horizontal line indicating the FP Sensitive Observation Cut off
        plots[msid]['ax'].axhline(FP_TEMP_SENSITIVE[msid], linestyle='--', color='red', linewidth=2.0)
        # Draw a horizontal line showing the ACIS-I -114 deg. C cutoff
        plots[msid]['ax'].axhline(ACIS_I_RED[msid], linestyle='--', color='purple', linewidth=1.0)
        # Draw a horizontal line showing the ACIS-S -112 deg. C cutoff
        plots[msid]['ax'].axhline(ACIS_S_RED[msid], linestyle='--', color='blue', linewidth=1.0)
 
        # Build the file name and output the file
        filename = MSID[msid].lower() + 'M120toM112.png'
        outfile = os.path.join(opt.outdir, filename)
        logger.info('   Writing plot file %s' % outfile)
        plots[msid]['fig'].savefig(outfile)
        plots[msid]['filename'] = filename

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
    plots['pow_sim']['ax'].axvline(load_start, linestyle=':', color='g',
                                   linewidth=1.0)
    plots['pow_sim']['fig'].subplots_adjust(right=0.85)
    filename = 'pow_sim.png'
    outfile = os.path.join(opt.outdir, filename)
    logger.info('   Writing plot file %s' % outfile)
    plots['pow_sim']['fig'].savefig(outfile)
    plots['pow_sim']['filename'] = filename

    return plots, obs_with_sensitivity

#------------------------------------------------------------------------------
#
#   get_states  
#
#------------------------------------------------------------------------------
def get_states(datestart, datestop, db):
    """Get states exactly covering date range

    :param datestart: start date
    :param datestop: stop date
    :param db: database handle
    :returns: np recarry of states
    """
    datestart = DateTime(datestart).date
    datestop = DateTime(datestop).date
    logger.info('Getting commanded states between %s - %s' %
                 (datestart, datestop))

    # Get all states that intersect specified date range
    cmd = """SELECT * FROM cmd_states
             WHERE datestop > '%s' AND datestart < '%s'
             ORDER BY datestart""" % (datestart, datestop)
    logger.debug('   Query command: %s' % cmd)
    states = db.fetchall(cmd)
    logger.info('   Found %d commanded states' % len(states))

    # Set start and end state date/times to match telemetry span.  Extend the
    # state durations by a small amount because of a precision issue converting
    # to date and back to secs.  (The reference tstop could be just over the
    # 0.001 precision of date and thus cause an out-of-bounds error when
    # interpolating state values).
    states[0].tstart = DateTime(datestart).secs - 0.01
    states[0].datestart = DateTime(states[0].tstart).date
    states[-1].tstop = DateTime(datestop).secs + 0.01
    states[-1].datestop = DateTime(states[-1].tstop).date

    return states
#------------------------------------------------------------------------------
#
#   make_validation_plots
#
#------------------------------------------------------------------------------
def make_validation_plots(opt, tlm, db):
    """
    Make validation output plots.

    :param outdir: output directory
    :param tlm: telemetry obtained from MAIN; get_telem_values
                '1dpamzt',              
                'fptemp_11',
                '1cbat',
                'sim_z',              
                'aosares1',           
                'dp_dpa_power
    :param db: database handle
    :returns: list of plot info including plot file names
    """
    outdir = opt.outdir
    start = tlm['date'][0]
    stop = tlm['date'][-1]
    states = get_states(start, stop, db)

    # Create array of times at which to calculate FP temperatures, then do it
    logger.info('\nCalculating Focal Plane thermal model for validation')

    model = calc_model(opt.model_spec, states, start, stop)

    # Interpolate states onto the tlm.date grid
    # TA REMOVED state_vals = cmd_states.interpolate_states(states, model.times)

    #? make a list of items for which predictions need be made
    pred = {'fptemp': model.comp['fptemp'].mvals,
            'pitch': model.comp['pitch'].mvals,
            'tscpos': model.comp['sim_z'].mvals
            }

    idxs = Ska.Numpy.interpolate(np.arange(len(tlm)), tlm['date'], model.times,
                                 method='nearest')

    tlm = tlm[idxs]


    labels = {'fptemp': 'Degrees (C)',
              'pitch': 'Pitch (degrees)',
              'tscpos': 'SIM-Z (steps/1000)',
              }

    scales = {'tscpos': 1000.}

    fmts = {'fptemp': '%.2f',
            'pitch': '%.3f',
            'tscpos': '%d'}

    # 10/2015 DH Heater acisfp_check.py Addition
    good_mask = np.ones(len(tlm),dtype='bool')
    for interval in model.bad_times:
        bad = ((tlm['date'] >= DateTime(interval[0]).secs)
            & (tlm['date'] < DateTime(interval[1]).secs))
        good_mask[bad] = False


    plots = []
    logger.info('   Making FPTEMP model validation plots and quantile table')
    quantiles = (1, 5, 16, 50, 84, 95, 99)
    # store lines of quantile table in a string and write out later
    quant_table = ''
    quant_head = ",".join(['MSID'] + ["quant%d" % x for x in quantiles])
    quant_table += quant_head + "\n"
    for fig_id, msid in enumerate(sorted(pred)): 
        plot = dict(msid=msid.upper())
        fig = plt.figure(10 + fig_id, figsize=(7, 3.5))
        fig.clf()
        scale = scales.get(msid, 1.0)
        ticklocs, fig, ax = plot_cxctime(model.times, tlm[msid] / scale,
                                         fig=fig, fmt='-r')
        ticklocs, fig, ax = plot_cxctime(model.times, pred[msid] / scale,
                                         fig=fig, fmt='-b')

        # Added 10/2015 DH Heater Addition
        if  np.any(~good_mask) :
            ticklocs, fig, ax = plot_cxctime(model.times[~good_mask], tlm[msid][~good_mask] / scale,
                                         fig=fig, fmt='.c')


        ax.set_title(msid.upper() + ' validation')
        ax.set_ylabel(labels[msid])
        ax.grid()
        filename = msid + '_valid.png'
        outfile = os.path.join(outdir, filename)
        logger.info('   Writing plot file %s' % outfile)
        fig.savefig(outfile)
        plot['lines'] = filename

        # Make quantiles
        #
        # For FPTEMP, the only quantiles we want are those where the temperature is 
        # -120.0 <= fp temp <= -112.0
        if msid == 'fptemp':
            ok = (tlm[msid] >= -120.0) & (tlm[msid] <= -112.0) # TA
        else:
            ok = np.ones(len(tlm[msid]), dtype=bool)


        diff = np.sort(tlm[msid][ok] - pred[msid][ok])

        quant_line = "%s" % msid
        for quant in quantiles:
            quant_val = diff[(len(diff) * quant) // 100]
            plot['quant%02d' % quant] = fmts[msid] % quant_val
            quant_line += (',' + fmts[msid] % quant_val)
        quant_table += quant_line + "\n"

        for histscale in ('log', 'lin'):
            fig = plt.figure(20 + fig_id, figsize=(4, 3))
            fig.clf()
            ax = fig.gca()
            ax.hist(diff / scale, bins=50, log=(histscale == 'log'))
            ax.set_title(msid.upper() + ' residuals: data - model')
            ax.set_xlabel(labels[msid])
            fig.subplots_adjust(bottom=0.18)
            filename = '%s_valid_hist_%s.png' % (msid, histscale)
            outfile = os.path.join(outdir, filename)
            logger.info('   Writing plot file %s' % outfile)
            fig.savefig(outfile)
            plot['hist' + histscale] = filename

        plots.append(plot)

    filename = os.path.join(outdir, 'validation_quant.csv')
    logger.info('   Writing quantile table %s' % filename)
    f = open(filename, 'w')
    f.write(quant_table)
    f.close()

    # If run_start is specified this is likely for regression testing
    # or other debugging.  In this case write out the full predicted and
    # telemetered dataset as a pickle.
    if opt.run_start:
        filename = os.path.join(outdir, 'validation_data.pkl')
        logger.info('   Writing validation data %s' % filename)
        f = open(filename, 'w')
        pickle.dump({'pred': pred, 'tlm': tlm}, f, protocol=-1)
        f.close()

    return plots
#------------------------------------------------------------------------------
#
#   plot_cxctime
#
#------------------------------------------------------------------------------
def plot_cxctime(times, y, fig=None, **kwargs):
    """Make a date plot where the X-axis values are in CXC time.  If no ``fig``
    value is supplied then the current figure will be used (and created
    automatically if needed).  Any additional keyword arguments
    (e.g. ``fmt='b-'``) are passed through to the ``plot_date()`` function.

    :param times: CXC time values for x-axis (date)
    :param y: y values
    :param fig: pyplot figure object (optional)
    :param **kwargs: keyword args passed through to ``plot_date()``

    :rtype: ticklocs, fig, ax = tick locations, figure, and axes object.
    """
    # if no Matplotlib figure has been supplied, use the current one or
    # make a new one
    if fig is None:
        fig = plt.gcf()

    # set ax to the current axes
    ax = fig.gca()

    # This importation should be moved to the top
    import Ska.Matplotlib

    # Plot yoru y value vs time
    ax.plot_date(Ska.Matplotlib.cxctime2plotdate(times), y, **kwargs)

    # Obviously this uses some sort of x tick mark labeler in the
    # Ska environment but the author hasn't commented what it is or
    # how it's used. i will backfil once I look at the function
    ticklocs = Ska.Matplotlib.set_time_ticks(ax)
    fig.autofmt_xdate()

    return ticklocs, fig, ax






#------------------------------------------------------------------------------
#
#   pointpair
#
#------------------------------------------------------------------------------
def pointpair(x, y=None):
    """
    I have no idea what this function does because it isn't commented
    and I haven't had the time to try and figure it out from the code
    """
    if y is None:
        y = x
    return np.array([x, y]).reshape(-1, order='F')


#------------------------------------------------------------------------------
#
#   globfile
#
#------------------------------------------------------------------------------
def globfile(pathglob):
    """Return the one file name matching ``pathglob``.  Zero or multiple
    matches raises an IOError exception."""

    files = glob.glob(pathglob)
    if len(files) == 0:
        raise IOError('No files matching %s' % pathglob)
    elif len(files) > 1:
        raise IOError('Multiple files matching %s' % pathglob)
    else:
        return files[0]


if __name__ == '__main__':
    opt, args = get_options()
    if opt.version:
        print VERSION
        sys.exit(0)

    try:
        main(opt)
    except Exception, msg:
        if opt.traceback:
            raise
        else:
            print "ERROR:", msg
            sys.exit(1)
