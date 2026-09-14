"""The one error BSReadSim raises for a command that cannot be carried out."""


class BSReadSimError(RuntimeError):
    """Invalid settings or inputs, a failed htsim, or an output that cannot be written."""
