"""Dependency-free errors shared by output encoders and publication."""


class OutputError(RuntimeError):
    """An artifact cannot be encoded or published safely."""


__all__ = ["OutputError"]
