import ska_helpers

__version__ = ska_helpers.get_version(__package__)

from .acisfp_check import \
    ACISFPCheck, opts, main, \
    model_path

from .acis_obs import \
    ObsidFindFilter
