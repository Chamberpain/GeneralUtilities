"""Lazy legacy cache paths using FilePathHandler without creating store folders."""


def make_download_paths(root, name):
    from GeneralUtilities.Data.Filepath.instance import FilePathHandler

    class DownloadPaths(FilePathHandler):
        def check_folder(self, folder):
            if folder not in (self._store, self._data):
                super().check_folder(folder)

        def store_file(self, filename):
            raise PermissionError("Use of store paths requires explicit permission")

    return DownloadPaths(root, name)


class LazyDownloadPaths:
    def __init__(self, root, name):
        self.root = root
        self.name = name
        self._handler = None

    def __getattr__(self, name):
        if name.startswith('__'):
            raise AttributeError(name)
        if self._handler is None:
            self._handler = make_download_paths(self.root, self.name)
        return getattr(self._handler, name)
