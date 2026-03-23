from tempfile import TemporaryDirectory
import os

from .run_external import run_external
from .sam import SAM


def map_reads_to_ref(ref: str, r1: str, r2: str) -> SAM:
    """Map provided reads to ref and return path to resulting SAM

    Args:
        ref (str): path to reference sequence(s)
        r1 (str): path to read 1 file
        r2 (str): path to read 2 file

    Returns:
        SAM: magnum_opus.SAM instance containing read mappings
    """
    
    mapping_args = [
        "minimap2",
        "-ax", "sr",
        "-k", "10",
        "-B", "1",
        "--score-N", "0"
    ]
    mapping_args += [ref, r1, r2]

    map_cmd = " ".join([str(a) for a in mapping_args])
    map_cmd += "| samtools view -hF4"

    with TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "tmp.sam")
        map_cmd += f" > {outfile}"
        run_external(map_cmd, shell=True)
        sam = SAM.from_sam(outfile)
    
    return sam
