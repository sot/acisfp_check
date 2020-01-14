import ska_helpers

__version__ = ska_helpers.get_version(__package__)

from .acisfp_check import \
    ACISFPCheck, calc_model

from .acis_obs import \
    ObsidFindFilter
