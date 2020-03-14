from ..acisfp_check import model_path, ACISFPCheck
from acis_thermal_check.regression_testing import \
    RegressionTester, all_loads
import pytest

acisfp_rt = RegressionTester(ACISFPCheck, model_path, "acisfp_test_spec.json")

# ACIS state builder tests

acisfp_rt.run_models(state_builder='acis')

# Prediction tests


@pytest.mark.parametrize('load', all_loads)
def test_prediction(answer_store, load):
    acisfp_rt.run_test("prediction", load, answer_store=answer_store)

# Validation tests


@pytest.mark.parametrize('load', all_loads)
def test_validation(answer_store, load):
    acisfp_rt.run_test("validation", load, answer_store=answer_store)
