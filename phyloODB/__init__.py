from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("phyloODB")
except PackageNotFoundError:
    # Source-tree fallback; keep in sync with pyproject.toml.
    __version__ = "0.2.0"


def main(*args, **kwargs):
    from .cli.main import main as cli_main

    return cli_main(*args, **kwargs)


__all__ = ["__version__", "main"]
