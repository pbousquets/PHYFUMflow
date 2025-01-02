import subprocess as sp

def validate_beast_command():
    """
    Check that beast is available
    """
    try:
        sp.check_output([
            "beast",
            "-version"
        ])
        print("BEAST is available.")
    except (sp.CalledProcessError, FileNotFoundError) as e:
        raise RuntimeError("""
        
        BEAST is not available. Please, make sure that BEAST is installed and available in PATH. 

        See https://phyfum.gitbook.io/tutorial/installation for more help
        
        """)
