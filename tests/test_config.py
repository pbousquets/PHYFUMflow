import pytest
from pathlib import Path
import snakemake
import src.phyfum.workflow.aux.checkConfig as checkConfig


def test_defaults(defaults: Path):
    cleanConfig = checkConfig.readDefaults(defaults)

    assert isinstance(cleanConfig, dict)
    assert len(cleanConfig) > 0


def test_config_cleanup(defaults: Path, example_config: Path):
    DEFAULTS = checkConfig.readDefaults(defaults)
    config = snakemake.common.configfile.load_configfile(example_config)
    finalConfig = checkConfig.makeConfig(config, DEFAULTS)
    assert len(finalConfig) > 3 and finalConfig.get("mle_iterations") == 7500
