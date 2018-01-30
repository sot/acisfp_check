from ..acisfp_check import VALIDATION_LIMITS, \
    HIST_LIMIT, calc_model, model_path, ACISFPCheck
from acis_thermal_check.regression_testing import \
    RegressionTester

atc_kwargs = {"other_telem": ['1dahtbon'],
              "other_map": {'1dahtbon': 'dh_heater',
                            "fptemp_11": "fptemp"}}

acisfp_rt = RegressionTester("fptemp", "acisfp", model_path, VALIDATION_LIMITS,
                             HIST_LIMIT, calc_model, atc_class=ACISFPCheck,
                             atc_kwargs=atc_kwargs)

def test_acisfp_loads(answer_store):
    acisfp_rt.run_test_arrays(answer_store)

