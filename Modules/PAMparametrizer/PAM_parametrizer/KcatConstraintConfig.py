import pandas as pd
from typing import Optional

class KcatConstraintConfigTable:
    """
    Container for user-defined kcat constraints in enzyme parametrization.

    Stores kcat constraints in a MultiIndex DataFrame with
    ('enzyme_id', 'reaction_id', 'direction') as the index. This makes lookups
    and slicing straightforward while enforcing a structured format.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing user-specified kcat constraints. Must include
        the following columns:
        - 'enzyme_id' : str
        - 'reaction_id' : str
        - 'direction' : {'f', 'b'}
        - 'min_kcat' : float
        - 'max_kcat' : float

    Attributes
    ----------
    df : pandas.DataFrame
        A validated DataFrame with a MultiIndex ('enzyme_id', 'reaction_id', 'direction')
        and columns ['min_kcat', 'max_kcat'] with kcats in 1/s.
    df_model_constraints : pandas.DataFrame
        A validated DataFrame with a MultiIndex ('enzyme_id', 'reaction_id', 'direction')
        and columns ['min_kcat', 'max_kcat'] with kcat in 1/h*1e-6 (compatible with model coefficients).

    Methods
    -------
    get(enzyme_id, reaction_id, direction)
        Return min/max kcat constraints for a given enzyme–reaction–direction.
    get_in_model_constraints(enzyme_id, reaction_id, direction)
        Return min/max kcat constraints for a given enzyme–reaction–direction in model compatible coefficients 1/h*1e-6
    add(enzyme_id, reaction_id, direction, min_kcat, max_kcat)
        Add min/max kcat constraint for a given enzyme–reaction–direction.
    has_constraint(enzyme_id, reaction_id, direction)
        Check if a constraint exists for the given enzyme–reaction–direction.

    Examples
    --------
    >>> import pandas as pd
    >>> data = {
    ...     "enzyme_id": ["E1", "E2"],
    ...     "reaction_id": ["R1", "R2"],
    ...     "direction": ["f", "b"],
    ...     "min_kcat": [0, 10],
    ...     "max_kcat": [100, 200],
    ... }
    >>> df = pd.DataFrame(data)
    >>> config = KcatConstraintConfigTable(df)
    >>> config.get("E1", "R1", "f")
    {'min_kcat': 0, 'max_kcat': 100}
    """

    REQUIRED_COLUMNS = ["enzyme_id", "reaction_id", "direction", "min_kcat", "max_kcat"]
    DIFFUSION_LIMIT = 1e6
    MIN_KCAT = 1e-6

    def __init__(self, df: Optional[pd.DataFrame] = None):
        if df is None:
            df = pd.DataFrame(columns=self.REQUIRED_COLUMNS)

        self._validate_input_df(df=df)
        self.df = df.set_index(["enzyme_id", "reaction_id", "direction"]).sort_index()

    @property
    def df_model_constraints(self) -> pd.DataFrame:
        """Returns the dataframe with the kcat values in model units: [1/h *1e-6] instead of [1/s]
        """
        df = self.df.copy()
        for col in ['min_kcat', 'max_kcat']:
            df[col] = df[col]*3600*1e-6
        return df

    def _validate_input_df(self, df: pd.DataFrame):
        """Check if dataframe has all the columns required to constrain the kcat values"""
        missing = set(self.REQUIRED_COLUMNS) - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Validate numeric ranges
        if (df["min_kcat"] < 0).any():
            raise ValueError("min_kcat must be >= 0")
        if (df["max_kcat"] <= df["min_kcat"]).any():
            raise ValueError("max_kcat must be strictly greater than min_kcat")
        #ToDo check if direction is correctly setup


    def get(self, enzyme_id: str, reaction_id: str, direction: str) -> dict:
        """Return the constraint as a dict for the given enzyme–reaction–direction."""
        try:
            row = self.df.loc[(enzyme_id, reaction_id, direction)]
            return {"min_kcat": row["min_kcat"], "max_kcat": row["max_kcat"]}
        except:
            raise KeyError(
                f"No kcat constraint for enzyme={enzyme_id}, reaction={reaction_id}, direction={direction}"
            )

    def get_in_model_constraints(self, enzyme_id: str, reaction_id: str, direction: str) -> dict:
        """Return the constraint as a dict for the given enzyme–reaction–direction with the
        kcat values in model units: [1/h *1e-6] instead of [1/s]."""
        try:
            row = self.df.loc[(enzyme_id, reaction_id, direction)]
            return {"min_kcat": row["min_kcat"].iloc[0]*3600*1e-6, "max_kcat": row["max_kcat"].iloc[0]*3600*1e-6}
        except:
            return {"min_kcat": self.MIN_KCAT*3600*1e-6, "max_kcat": self.DIFFUSION_LIMIT*3600*1e-6}

    def add(self, enzyme_id:str, reaction_id:str, direction:str,
            min_kcat: Optional[float] = MIN_KCAT, max_kcat: Optional[float] = DIFFUSION_LIMIT) -> dict:
        """Add a new kcat constraint to the KcatConfigTable"""
        new_row = pd.DataFrame(
            [(enzyme_id, reaction_id, direction,min_kcat, max_kcat)],
            columns=["enzyme_id", "reaction_id", "direction", "min_kcat", "max_kcat"]
        ).set_index(["enzyme_id", "reaction_id", "direction"])
        self.df = pd.concat([self.df, new_row])

    def has_constraint(self, enzyme_id: str, reaction_id: str, direction: str) -> bool:
        """Check if a constraint exists for the given enzyme–reaction–direction."""
        return (enzyme_id, reaction_id, direction) in self.df.index

    def ensure_kcats_in_pam_info_file_are_within_bounds(self,
                                                        pam_info_file: str
                                                        ) -> None:
        enzyme_db = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')
        merged_enzyme_db = pd.merge(enzyme_db, self.df,
                                    how='left',
                                    left_on=['rxn_id', 'enzyme_id', 'direction'],
                                    right_on=['reaction_id', 'enzyme_id', 'direction']
                                    )
        merged_enzyme_db['kcat_values'] = merged_enzyme_db['kcat_values'].clip(
            lower=merged_enzyme_db['min_kcat'], upper=merged_enzyme_db['max_kcat']
        )
        updated_enzyme_db = merged_enzyme_db[enzyme_db.columns]
        with pd.ExcelWriter(pam_info_file, mode='a', if_sheet_exists='replace', engine='openpyxl') as writer:
            updated_enzyme_db.to_excel(writer, sheet_name='ActiveEnzymes', index=False)
