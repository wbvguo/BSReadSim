"""Immutable data resources bundled with BSReadSim releases."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from importlib import resources as importlib_resources
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any


_RESOURCE_NAME = re.compile(r"^[a-z0-9][a-z0-9-]*$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_RESOURCE_STRING_FIELDS = (
    "name",
    "kind",
    "format",
    "path",
    "sha256",
    "source",
    "license",
    "description",
)
_RESOURCE_FIELDS = frozenset((*_RESOURCE_STRING_FIELDS, "size_bytes"))


class ResourceError(ValueError):
    """Bundled resource metadata or bytes are unavailable or invalid."""


@dataclass(frozen=True)
class BundledResource:
    """One immutable resource declared by the release registry."""

    name: str
    kind: str
    format: str
    path: str
    size_bytes: int
    sha256: str
    source: str
    license: str
    description: str


def _join(root: Any, relative_path: str) -> Any:
    target = root
    for part in PurePosixPath(relative_path).parts:
        target = target.joinpath(part)
    return target


def _data_root() -> Any:
    """Return installed package data, with a source-checkout fallback."""
    packaged = importlib_resources.files("bsreadsim").joinpath("data")
    try:
        if packaged.joinpath("registry.json").is_file():
            return packaged
    except OSError:
        pass

    source_data = Path(__file__).resolve().parents[2] / "data"
    if (source_data / "registry.json").is_file():
        return source_data
    raise ResourceError("the BSReadSim resource registry is not installed")


def _registry_document() -> dict[str, Any]:
    registry = _data_root().joinpath("registry.json")
    try:
        document = json.loads(registry.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ResourceError("cannot read the BSReadSim resource registry") from error
    if not isinstance(document, dict):
        raise ResourceError("resource registry must be a JSON object")
    if set(document) != {"resources"}:
        raise ResourceError("resource registry fields are invalid")
    return document


def _resource_from_document(document: Any) -> BundledResource:
    if not isinstance(document, dict):
        raise ResourceError("each resource registry entry must be an object")
    if set(document) != _RESOURCE_FIELDS:
        raise ResourceError("resource registry entry fields are invalid")
    values: dict[str, Any] = {}
    for field in _RESOURCE_STRING_FIELDS:
        value = document.get(field)
        if not isinstance(value, str) or not value:
            raise ResourceError("resource field {!r} must be a nonempty string".format(field))
        values[field] = value

    name = values["name"]
    if _RESOURCE_NAME.fullmatch(name) is None:
        raise ResourceError("resource name is invalid: {}".format(name))
    relative = PurePosixPath(values["path"])
    if (
        relative.is_absolute()
        or any(part in ("", ".", "..") for part in relative.parts)
        or relative.as_posix() != values["path"]
    ):
        raise ResourceError("resource path is not canonical: {}".format(values["path"]))
    if _SHA256.fullmatch(values["sha256"]) is None:
        raise ResourceError("resource SHA-256 is invalid: {}".format(name))

    size_bytes = document.get("size_bytes")
    if isinstance(size_bytes, bool) or not isinstance(size_bytes, int) or size_bytes < 0:
        raise ResourceError("resource size_bytes is invalid: {}".format(name))
    return BundledResource(size_bytes=size_bytes, **values)


def list_resources() -> tuple[BundledResource, ...]:
    """List resources declared by the installed release."""
    entries = _registry_document().get("resources")
    if not isinstance(entries, list):
        raise ResourceError("resource registry resources must be an array")
    resources = tuple(_resource_from_document(entry) for entry in entries)
    names = [resource.name for resource in resources]
    if len(names) != len(set(names)):
        raise ResourceError("resource registry contains duplicate names")
    return resources


def get_resource(name: str) -> BundledResource:
    """Return metadata for a named bundled resource."""
    if not isinstance(name, str) or not name:
        raise ResourceError("resource name must be a nonempty string")
    for resource in list_resources():
        if resource.name == name:
            return resource
    raise ResourceError("unknown bundled resource: {}".format(name))


def read_resource(name: str) -> bytes:
    """Read and verify the exact bytes of a named resource."""
    resource = get_resource(name)
    path = _join(_data_root(), resource.path)
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise ResourceError("cannot read bundled resource: {}".format(name)) from error
    if len(payload) != resource.size_bytes:
        raise ResourceError("bundled resource size mismatch: {}".format(name))
    observed = hashlib.sha256(payload).hexdigest()
    if observed != resource.sha256:
        raise ResourceError("bundled resource checksum mismatch: {}".format(name))
    return payload


def copy_resource(name: str, output: str | Path) -> Path:
    """Copy a verified resource to a new user-owned path without overwriting."""
    payload = read_resource(name)
    destination = Path(output).expanduser().resolve(strict=False)
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as output_file:
            output_file.write(payload)
    except FileExistsError as error:
        raise ResourceError("resource output already exists: {}".format(destination)) from error
    except OSError as error:
        try:
            if destination.is_file():
                destination.unlink()
        except OSError:
            pass
        raise ResourceError(
            "cannot write resource output {}: {}".format(destination, error)
        ) from error
    return destination


__all__ = [
    "BundledResource",
    "ResourceError",
    "copy_resource",
    "get_resource",
    "list_resources",
    "read_resource",
]
