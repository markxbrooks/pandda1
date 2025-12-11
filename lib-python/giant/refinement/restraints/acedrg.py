import os

import giant.logs as lg
from giant.dispatcher import Dispatcher

logger = lg.getLogger(__name__)


def generate_restraints(smiles, name="LIG", prefix="ligand", verbose=False):
    """Generate pdb and cif wrappers from smiles string"""

    assert len(name) == 3

    # Common fixes
    smiles = smiles.replace("CL", "Cl")
    smiles = smiles.replace("BR", "Br")

    out_pdb = prefix + ".pdb"
    out_cif = prefix + ".cif"
    out_log = prefix + "-acedrg.log"

    if os.path.exists(out_pdb):
        raise IOError("Output PDB file already exists: {}".format(out_pdb))
    if os.path.exists(out_cif):
        raise IOError("Output CIF file already exists: {}".format(out_cif))

    # Run acedrg
    acedrg = Dispatcher("acedrg")

    acedrg.extend_args(
        [
            "--smi={}".format(smiles),
            "-r",
            name,
            "-o",
            prefix,
        ]
    )

    if verbose is True:
        logger(acedrg.as_string())

    acedrg.run()
    acedrg.write_output(out_log)

    # Check which files actually exist
    pdb_exists = os.path.exists(out_pdb)
    cif_exists = os.path.exists(out_cif)
    
    # Get the directory and base name for searching
    prefix_dir = os.path.dirname(os.path.abspath(prefix)) if os.path.dirname(prefix) else os.getcwd()
    prefix_base = os.path.basename(prefix)
    
    # acedrg might create files with the residue name instead of prefix
    alt_pdb = os.path.join(prefix_dir, f"{name}.pdb")
    alt_cif = os.path.join(prefix_dir, f"{name}.cif")
    
    # Also check in current working directory if acedrg changed it
    cwd_pdb = os.path.join(os.getcwd(), f"{prefix_base}.pdb")
    cwd_cif = os.path.join(os.getcwd(), f"{prefix_base}.cif")
    
    # Use alternative paths if main paths don't exist
    if not pdb_exists:
        if os.path.exists(alt_pdb):
            out_pdb = alt_pdb
            pdb_exists = True
        elif os.path.exists(cwd_pdb):
            out_pdb = cwd_pdb
            pdb_exists = True
            
    if not cif_exists:
        if os.path.exists(alt_cif):
            out_cif = alt_cif
            cif_exists = True
        elif os.path.exists(cwd_cif):
            out_cif = cwd_cif
            cif_exists = True

    if not (pdb_exists and cif_exists):
        stdout = str(acedrg.result.stdout) if acedrg.result.stdout else ""
        stderr = str(acedrg.result.stderr) if acedrg.result.stderr else ""
        returncode = getattr(acedrg.result, "returncode", None)
        
        # List files in the output directory and current directory for debugging
        files_in_prefix_dir = []
        files_in_cwd = []
        if os.path.exists(prefix_dir):
            try:
                files_in_prefix_dir = [f for f in os.listdir(prefix_dir) if f.endswith(('.pdb', '.cif'))]
            except (OSError, PermissionError):
                pass
        try:
            files_in_cwd = [f for f in os.listdir(os.getcwd()) if f.endswith(('.pdb', '.cif'))]
        except (OSError, PermissionError):
            pass
        
        error_msg = "acedrg failed during ligand generation"
        if returncode is not None:
            error_msg += f" (return code: {returncode})"
        error_msg += f"\nExpected PDB: {out_pdb} (exists: {pdb_exists})"
        error_msg += f"\nExpected CIF: {out_cif} (exists: {cif_exists})"
        error_msg += f"\nPrefix directory: {prefix_dir}"
        error_msg += f"\nCurrent working directory: {os.getcwd()}"
        if files_in_prefix_dir:
            error_msg += f"\nFiles in prefix directory: {', '.join(files_in_prefix_dir)}"
        if files_in_cwd:
            error_msg += f"\nFiles in current directory: {', '.join(files_in_cwd)}"
        if stdout:
            error_msg += f"\nSTDOUT:\n{stdout}"
        if stderr:
            error_msg += f"\nSTDERR:\n{stderr}"
        
        logger(error_msg)
        raise Exception(error_msg)

    return out_pdb, out_cif
