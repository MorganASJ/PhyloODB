def main(*args, **kwargs):
    from .cli.main import main as cli_main

    return cli_main(*args, **kwargs)


__all__ = ["main"]
