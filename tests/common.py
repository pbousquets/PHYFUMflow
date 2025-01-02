"""
Common code for unit testing of rules generated with Snakemake 8.18.0.
"""

from pathlib import Path
import subprocess as sp
import os
import hashlib
import string

class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir
        
    def is_binary(file_path):
        # Read a small chunk of the file to determine if it's binary or text
        with open(file_path, 'rb') as file:
            chunk = file.read(1024)  # Read the first 1KB
        # If the chunk contains non-printable characters, treat it as binary
        return any(byte not in string.printable.encode('ascii') for byte in chunk)

    def hash_file(file_path):
        # Determine if the file is binary or text
        if OutputChecker.is_binary(file_path):
            # Handle binary file: hash the raw binary content
            return OutputChecker.hash_binary_file(file_path)
        else:
            # Handle text file: remove comments and hash
            return OutputChecker.hash_text_file(file_path)

    def hash_binary_file(file_path):
        hash_obj = hashlib.sha256()
        with open(file_path, 'rb') as file:
            while chunk := file.read(8192):
                hash_obj.update(chunk)
        return hash_obj.hexdigest()

    def hash_text_file(file_path):
        hash_obj = hashlib.sha256()
        with open(file_path, 'r') as file:
            # Filter out lines that are comments (starting with '#')
            text = ''.join(line for line in file if not line.strip().startswith('#'))
            hash_obj.update(text.encode())
        return hash_obj.hexdigest()


    def check(self):
        input_files = {f.relative_to(self.data_path) for f in self.data_path.rglob('*') if f.is_file()}
        expected_files = {f.relative_to(self.expected_path) for f in self.expected_path.rglob('*') if f.is_file()}

        unexpected_files = expected_files - input_files
        common_files = result_files & expected_files

        for file in common_files:
            result_file = result_path / file
            expected_file = expected_path / file
            compare_files(result_file, expected_file)

        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_files(self, generated_file, expected_file):
        assert OutputChecker.hash_file(generated_file) == OutputChecker.hash_file(expected_file) 
        
