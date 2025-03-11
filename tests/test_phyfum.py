import shutil
import subprocess as sp
from pathlib import Path, PurePosixPath
import re
import pytest
import os
from tests import common

BASE_DIR = Path(__file__).resolve().parent

@pytest.fixture(scope="module")
def phyfum_output():
    """Fixture to run phyfum once and use its output for all tests"""
    tmpdir = "/tmp/na"
    workdir = Path(tmpdir) / "workdir"
    
    # Create workdir if it doesn't exist
    os.makedirs(workdir, exist_ok=True)
    
    result_path = str(PurePosixPath(f"{workdir}/test_run"))
    inp = str((BASE_DIR / "data/exampleBeta.csv").absolute())
    patientinfo = str((BASE_DIR / "data/metadata.csv").absolute())
    
    # Run the phyfum command once
    sp.check_output([
        "phyfum", 
        "run", 
        "--input", inp, 
        "--patientinfo", patientinfo, 
        "--output", "na", 
        "--iterations", "1000", 
        "--workdir", result_path, 
        "--stemcells", "2-4-1", 
        "--nchains", "1", 
        "--sampling", "2", 
        "--notemp",
        "--mle-sampling", "2",
        "--mle-steps", "2",
        "--mle-ss"
    ])
    
    yield {
        "workdir": workdir,
        "result_path": result_path
    }
    
    # Cleanup after all tests run
    if workdir.exists():
        shutil.rmtree(workdir)

def test_main_output_files(phyfum_output):
    """Test main output directory files match expected"""
    result_path = phyfum_output["result_path"]
    expected_path = PurePosixPath("tests/expected")
    
    # Define patterns to ignore for timestamp/runtime comparisons
    ignore_patterns = [
        r"\[R-package APE, \w+ \w+\s+\d+ \d+:\d+:\d+ \d+\]",  # For overS.tree
        r"# \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d+)?",        # For overS.log
        r"# Generated \w+ \w+ \d+ \d+:\d+:\d+ \w+ \d+ \[seed=\d+\]",  # For seed*.log
        r"# -beagle_off -seed \d+ .+"  # For seed*.log with variable seed and path
    ]
    
    # Create a file checker that ignores runtime patterns
    file_checker = FilePatternChecker(ignore_patterns)
    
    # Check that all expected files exist and have correct content
    check_directory_contents(result_path, expected_path, file_checker)



class FilePatternChecker:
    """Class to compare files while ignoring specific patterns (like timestamps)"""
    def __init__(self, ignore_patterns):
        self.ignore_patterns = ignore_patterns
    
    def are_files_equal(self, file1, file2):
        # Check if the file is binary
        if common.OutputChecker.is_binary(file1) and common.OutputChecker.is_binary(file2):
            hash1 = common.OutputChecker.hash_binary_file(file1)
            hash2 = common.OutputChecker.hash_binary_file(file2)
            return hash1 == hash2
        
        # For text files, read content and compare with pattern ignoring
        with open(file1, 'r', errors='ignore') as f1, open(file2, 'r', errors='ignore') as f2:
            content1 = f1.read()
            content2 = f2.read()
            
            # Apply pattern ignoring
            for pattern in self.ignore_patterns:
                content1 = re.sub(pattern, "", content1)
                content2 = re.sub(pattern, "", content2)
                
            # Compare normalized content
            return content1.strip() == content2.strip()

def check_directory_contents(result_path, expected_path, file_checker):
    """Check that directories contain the same files with same content"""
    # Get all files in the result and expected directories
    result_files = {f.relative_to(result_path) for f in Path(result_path).rglob('*') 
                   if f.is_file() and not f.name.startswith('.') and '.snakemake' not in f.parts and not f.name.endswith('.pdf') and not f.name.endswith('.ops') and not f.name.endswith('Params.csv') and not f.name.endswith('stdout')}
    expected_files = {f.relative_to(expected_path) for f in Path(expected_path).rglob('*') 
                     if f.is_file() and not f.name.startswith('.') and '.snakemake' not in f.parts and not f.name.endswith('.pdf') and not f.name.endswith('.ops') and not f.name.endswith('Params.csv') and not f.name.endswith('stdout')}

    # Check that all expected files exist in the result
    missing_files = expected_files - result_files
    if missing_files:
        pytest.fail(f"Missing expected files in output: {missing_files}")
    
    # Check file content for all expected files
    for file_path in expected_files:
        result_file = result_path / file_path
        expected_file = expected_path / file_path
        
        assert file_checker.are_files_equal(result_file, expected_file), \
            f"Content differs in file: {file_path}"