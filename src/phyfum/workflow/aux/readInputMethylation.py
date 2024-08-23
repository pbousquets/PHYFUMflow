import pandas as pd
from typing import List


class methylationSample:
    """
    Class for reading samples and casting them into a single comma-separated element
    """

    def __init__(self, methylationArray: pd.Series, precision: int):
        self.methylationArray = methylationArray.round(precision)  # Round floats, we don't need long decimals

    def toString(self) -> str:
        return self.methylationArray.astype("string").str.cat(sep=",")  # Join all sites into a single comma-separated string


class readMethylation:
    """
    Class representing the whole dataset
    """

    def __init__(self, inputStr: str, precision: int, stripRownames: bool, samplesheet: str, age_col: str, age_diagnosis_col: str, sample_col: str, sample_type_col:str, luca_mode: str):
        index_col = 0 if stripRownames else None
        self.input = pd.read_csv(inputStr, index_col=index_col)
        self.precision = precision
        self.luca_mode = luca_mode
        self.samples: List[str] = []
        self.sample_names: List[str] = self.input.columns.tolist()
        self.checkValues()
        self.parseSampleSheet(samplesheet, age_col, age_diagnosis_col, sample_col, sample_type_col)

    def parseSamples(self) -> None:
        for column in self.sample_names:
            sample: methylationSample = methylationSample(self.input[column], self.precision)
            self.samples.append(sample.toString())

    def checkValues(self) -> None:
        assert len(self.input.select_dtypes(include="object").columns) == 0, "\nStrings found in the dataframe. Did you to forget using the --stripRownames flag?\n"
        assert all(self.input.max() <= 1) and all(self.input.min()) >= 0, "\nThe data frame contains values not in (0,1) range. \nDid you to forget using the --stripRownames flag?\n"

    def parseSampleSheet(self, samplesheet: str, age_col: str, age_diagnosis_col: str, sample_col: str, sample_type_col: str) -> None:
        self.samplesheet = pd.read_csv(samplesheet)
        assert sample_col in self.samplesheet.columns, f"Column {sample_col} to estimate samples not found in sample sheet"
        assert age_col in self.samplesheet.columns, f"Column {age_col} to estimate ages not found in sample sheet"
        assert sample_type_col in self.samplesheet.columns, f"Column {sample_type_col} to estimate sample types not found in sample sheet"

        if age_diagnosis_col:
            assert age_diagnosis_col in self.samplesheet.columns, f"Column {age_diagnosis_col} not found in sample sheet"

        # Filter samples with sample_names
        self.samplesheet = self.samplesheet[self.samplesheet[sample_col].isin(self.sample_names)]

        # Check if there are normals
        self.normals = self.samplesheet[sample_type_col].str.lower().str.count('normal|control').sum() > 0

        self.fixed_luca = True if self.luca_mode == "fixed" or (self.luca_mode.lower() == "auto" and self.normals) else False
        # If present, age_diagnosis is the same across all samples. Check it and save it as an int attribute
        if age_diagnosis_col:
            assert self.samplesheet[age_diagnosis_col].nunique() == 1, "Age at diagnosis is not the same across all samples"
            self.age_diagnosis = int(self.samplesheet[age_diagnosis_col].unique()[0])
        else: # If not present, age_diagnosis is the youngest sample in age_col
            self.age_diagnosis = self.samplesheet[age_col].min()
        
        # DOLS is the age of the oldest sample in age_col, DOFS is the age of the first sample
        self.age_dols = self.samplesheet[age_col].max()
        self.age_dofs = self.samplesheet[age_col].min()

        # Save age as a dict of sample_name: age
        self.ages = self.samplesheet.set_index(sample_col)[age_col].to_dict()
        self.tree_settings = {
            "heightRules": "true" if self.fixed_luca else "false",
            "cenancestorHeight": {
                "value": str(self.age_dols) if self.fixed_luca else str(self.age_dols - self.age_diagnosis),
            },
            "cenancestorBranch": {
                "value": "1",
                "lower": "0",
                "upper": str(self.age_dofs),
            },
            "fixed_luca": self.fixed_luca
        }

        if not self.fixed_luca:
            self.tree_settings["cenancestorHeight"]["lower"] = str(self.age_dols - self.age_diagnosis)
            self.tree_settings["cenancestorHeight"]["upper"] = str(self.age_dols)
        self.tree_settings["luca_height_prior"] = { # TODO: when solved bug #6, nest this within the previous ifelse statement
            "lower": str(self.age_dols - self.age_diagnosis),
            "upper": str(self.age_dols)
            }
