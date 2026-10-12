"""Handler-derived output and scratch paths without implicit data/store writes."""
from pathlib import Path
import re

from . import instance


def _component(value):
    if not isinstance(value, str) or value.casefold() == "store" or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError("A component must be a simple directory name.")
    return value


class _ManagedPaths:
    """Initialize legacy path fields lazily and permit only the designated root."""

    def __init__(self, root, component):
        self._initializing = True
        super().__init__(str(root), _component(component))
        self._initializing = False

    def check_folder(self, folder):
        if self._initializing:
            return
        allowed = getattr(self, "_" + self.kind)
        if str(folder) != str(allowed):
            raise PermissionError("Use the designated output/scratch handler; data and store are read-only.")
        configured = Path(
            instance.get_data_folder() if self.kind == "out" else instance.get_compute_folder()
        )
        if not configured.is_dir():
            raise FileNotFoundError(
                f"Configured {self.kind} root {configured} is unavailable; request access "
                "to this destination and its missing parents. No alternate root was selected."
            )
        Path(folder).mkdir(parents=True, exist_ok=True)

    def file_return(self, base_folder, filename):
        allowed = Path(getattr(self, "_" + self.kind))
        name = Path(str(filename))
        if name.is_absolute() or ".." in name.parts or any(part.casefold() == "store" for part in name.parts):
            raise ValueError("Output filenames must stay inside their handler directory.")
        configured = Path(instance.get_data_folder() if self.kind == "out" else instance.get_compute_folder())
        relative = allowed.relative_to(configured)
        expected_root = configured.resolve() / relative
        if allowed.resolve() != expected_root:
            raise ValueError("A handler directory link escapes its configured project/component root.")
        target = allowed / name
        if not target.resolve().is_relative_to(allowed.resolve()):
            raise ValueError("A compatibility link points outside the designated handler directory.")
        self.check_folder(base_folder)
        target.parent.mkdir(parents=True, exist_ok=True)
        return str(target)

    def store_file(self, filename):
        raise PermissionError("Store paths require explicit user permission.")

    def data_file(self, filename):
        raise PermissionError("Shared input data must be opened read-only without creating its directory.")


class ManagedFilePaths(_ManagedPaths, instance.FilePathHandler):
    """Saved files beneath DataDir/<project>/<component>/out."""
    kind = "out"


class ManagedComputePaths(_ManagedPaths, instance.ComputePathHandler):
    """Scratch files beneath ComputeDir/<project>/<component>/tmp."""
    kind = "tmp"


def make_file_paths(root, component):
    return ManagedFilePaths(root, component)


def make_compute_paths(root, component):
    return ManagedComputePaths(root, component)
