from openmmxtb.openmmxtb import XtbForce, XtbPointCharge

def _get_forcefield_dir():
    from pkg_resources import resource_filename
    return resource_filename('openmmxtb', 'forcefields')
