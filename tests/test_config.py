import pytest
from pathlib import Path
import snakemake
import src.phyfum.workflow.aux.checkConfig as checkConfig

# Define the base directory for your project
BASE_DIR = Path(__file__).resolve().parent.parent

def test_defaults():
    defaults = BASE_DIR / "src/phyfum/workflow/config.schema.yml"
    cleanConfig = checkConfig.readDefaults(defaults)

    assert isinstance(cleanConfig, dict)
    assert len(cleanConfig) > 0