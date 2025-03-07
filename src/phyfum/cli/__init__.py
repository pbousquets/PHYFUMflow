from pathlib import Path
import subprocess
from snk_cli import CLI
import typer
import shutil
import os
import sys 
phyfum = CLI(Path(__file__).parent.parent)

@phyfum.app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def betas2xml(ctx: typer.Context):
    """
    Generate XML files for BEAST from methylation beta values CSV file.
    
    This command runs the methylationBetas2xml.py script with the given arguments.
    """
    # Get the path to the methylationBetas2xml.py script
    script_path = Path(__file__).parent.parent / "workflow" / "aux" / "methylationBetas2xml.py"
    
    if not script_path.exists():
        typer.echo(f"Error: Script not found at {script_path}")
        raise typer.Exit(1)
    
    # If no arguments provided, show help
    if len(ctx.args) == 0:
        arguments = ["--help"]
    else:
        arguments = ctx.args
    
    # Build and run the command
    command = [sys.executable, str(script_path)] + arguments
    _ = subprocess.run(command, capture_output=False, text=True)


@phyfum.app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def beast(ctx: typer.Context):
    """
    Run the inference step (BEAST) directly.
    """
    
    if len(ctx.args) == 0:
        arguments = ["-h"]
    else:
        arguments = ctx.args
    command = ["beast"] + arguments
    _ = subprocess.run(command, capture_output=False, text=True)
    

@phyfum.app.command()
def install_beast(
    ctx: typer.Context,
    install_dir: str = typer.Option(None, help="Directory to install Beast (default: ~/phyfum-beast)"),
    add_to_path: bool = typer.Option(False, "--add-to-path", help="Add Beast bin directory to PATH in .bashrc"),
    force: bool = typer.Option(False, "--force", help="Force installation without confirmation"),

):
    """
    Install PHYFUM inference framework (our modified version of BEAST for fluctuating methylation).
    
    This command clones the repository, compiles it using Ant,
    and optionally adds the bin directory to your PATH by modifying .bashrc.
    
    Remember: This process requires git and ant to be installed on your system.
    """

    if not any(param for param in [add_to_path, force, install_dir]):
        ctx.get_help()
        raise typer.Exit(0)

    repo_url = "git@github.com:pbousquets/PHYFUM.git"
    # Set default install directory if not provided
    if install_dir is None:
        install_dir = str(Path.home() / "phyfum-beast")
    
    # Check if git is installed
    try:
        subprocess.run(["git", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.SubprocessError, FileNotFoundError):
        typer.echo("Error: Git is not installed or not available in PATH")
        raise typer.Exit(1)
    
    # Check if ant is installed
    try:
        subprocess.run(["ant", "-version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.SubprocessError, FileNotFoundError):
        typer.echo("Error: Ant is not installed or not available in PATH")
        typer.echo("Please install Ant before proceeding.")
        raise typer.Exit(1)
    
    # Prepare installation directory
    install_path = Path(install_dir).expanduser().resolve()
    repo_dir = install_path / "PHYFUM"
    
    if repo_dir.exists() and not force:
        overwrite = typer.confirm(f"Directory {repo_dir} already exists. Overwrite?", default=False)
        if not overwrite:
            typer.echo("Installation aborted.")
            raise typer.Exit(1)
        shutil.rmtree(repo_dir)
    
    os.makedirs(install_path, exist_ok=True)
    
    # Clone the repository
    typer.echo(f"Cloning repository from {repo_url} to {install_path}...")
    try:
        subprocess.run(["git", "clone", repo_url, str(repo_dir)], check=True)
    except subprocess.CalledProcessError:
        typer.echo(f"Error: Failed to clone repository from {repo_url}")
        typer.echo("If using SSH, ensure your SSH keys are set up correctly.")
        raise typer.Exit(1)
    
    # Compile using ant
    typer.echo(f"Compiling Beast in {repo_dir}...")
    try:
        os.chdir(repo_dir)
        subprocess.run(["ant", "linux"], check=True)
    except subprocess.CalledProcessError:
        typer.echo("Error: Failed to compile Beast")
        raise typer.Exit(1)
    
    # Determine bin directory path
    bin_dir = repo_dir / "release" / "Linux" / "Phyfumv1.0_RC1" / "bin"
    if not bin_dir.exists():
        typer.echo(f"Warning: Expected bin directory {bin_dir} not found")
        # Try to find bin directory
        bin_candidates = list(repo_dir.glob("release/Linux/*/bin"))
        if bin_candidates:
            bin_dir = bin_candidates[0]
            typer.echo(f"Using found bin directory: {bin_dir}")
        else:
            typer.echo("Could not find bin directory, installation may be incomplete")
            bin_dir = repo_dir  # Fallback
    
    # Add to PATH if requested
    if add_to_path and bin_dir.exists():
        bashrc_path = Path.home() / ".bashrc"
        export_line = f'export PATH="{bin_dir}:$PATH"'
        
        # Check if line already exists in .bashrc
        add_line = True
        if bashrc_path.exists():
            with open(bashrc_path, "r") as f:
                if export_line in f.read():
                    add_line = False
        
        if add_line:
            with open(bashrc_path, "a") as f:
                f.write(f"\n# Added by phyfum install_beast\n{export_line}\n")
            typer.echo(f"Added Beast bin directory to PATH in {bashrc_path}")
            typer.echo("You'll need to restart your terminal or run 'source ~/.bashrc' for this to take effect")
        else:
            typer.echo(f"Beast bin directory is already in PATH")
    
    typer.echo("\n✅ Beast MCMC Flipflop installation completed successfully!")
    typer.echo(f"Binary location: {bin_dir}")
    
    if not add_to_path:
        typer.echo("\nTo use Beast without adding to PATH, run:")
        typer.echo(f'export PATH="{bin_dir}:$PATH"')
    
    typer.echo("\n📋 Uninstallation instructions:")
    typer.echo(f"1. Remove the installation directory: rm -rf {install_path}")
    if add_to_path:
        typer.echo(f"2. Remove the PATH entry from {bashrc_path}:")
        typer.echo(f"   Edit {bashrc_path} and remove the lines:")
        typer.echo(f"   # Added by phyfum install_beast")
        typer.echo(f"   {export_line}")


if __name__ == "__main__":
    phyfum.app()
