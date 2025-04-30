from Modules.PAM_parametrizer.PAM_data_classes import ParametrizationResults

def test_if_ParametrizationResults_is_initialized_properly():
    obj1 = ParametrizationResults(substrate_uptake_reactions=['EX_glc__D_e'])
    obj2 = ParametrizationResults(substrate_uptake_reactions=['EX_glc__D_e'])
    for attr in obj1.__dataclass_fields__:
        assert id(getattr(obj1, attr)) == id(getattr(obj2, attr),f"Warning: Field '{attr}' is shared between instances!"