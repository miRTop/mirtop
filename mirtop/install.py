"""
Some commands to install common databases like mirbase and some human/mouse annotation
"""

import os.path as op
import subprocess

try:
    import bcbio
except:
    print("Probably this will fail, you need bcbio-nextgen for many installation functions.")
    pass

REMOTES = {
            "requirements": "https://raw.github.com/lpantano/seqcluster/master/requirements.txt",
            "gitrepo": "https://github.com/lpantano/seqcluster.git",
            "gitrepo-bcbio": "https://github.com/chapmanb/bcbio-nextgen.git"
           }


def _get_miraligner():
    try:
        tool = bcbio.pipeline.config_utils.get_program("miraligner", {}, "cmd")
    except ImportError:
        tool = None
        pass
    if not tool:
        url = "https://github.com/lpantano/seqbuster/raw/master/modules/miraligner/miraligner.jar"
        subprocess.check_call(["wget", "-O", "miraligner.jar", "--no-check-certificate", url])
        tool = "java -jar {opts} %s" % op.abspath("miraligner.jar")
    else:
        tool = "%s {opts}" % tool
    return tool
