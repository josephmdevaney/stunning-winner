"""Starts the vendored `interpret_paraphase.py` with a default timeout on its HTTP requests.

    python launch.py <path to interpret_paraphase.py> [the script's own arguments]

PacBio's script looks every variant up on Ensembl VEP and NCBI ClinVar through `requests`, and it passes no timeout.
A connection that goes quiet half way (a network change, the machine sleeping) is then waited on forever: a
CDS48001 run sat on one for 17 hours, using two seconds of CPU, with nothing in the log to say so. This launcher
gives every request that has no timeout of its own a connect and a read timeout. The script's retry logic then
treats a timed-out request like any other failed one, so a bad connection is retried and, if it stays bad, the run
fails instead of hanging.

`socket.setdefaulttimeout` would be the obvious way, and does nothing here: `requests` hands urllib3 an explicit
`timeout=None`, which overrides the socket default. The timeout is therefore set on the HTTP adapter's `send`.

Stdlib only; `requests` is imported if it is installed (the script cannot run without it either).
"""

from __future__ import annotations

import os
import runpy
import sys

# Seconds. Generous for a single-variant lookup, short enough that a dead connection is noticed in minutes.
CONNECT_TIMEOUT = float(os.environ.get("PARAPHASE_REQUEST_CONNECT_TIMEOUT", 15))
READ_TIMEOUT = float(os.environ.get("PARAPHASE_REQUEST_READ_TIMEOUT", 90))


def install_default_timeout() -> bool:
    """Make `requests` use (CONNECT_TIMEOUT, READ_TIMEOUT) for any request that sets no timeout. A timeout the caller
    does set is left alone. Returns False when `requests` is not installed."""
    try:
        from requests.adapters import HTTPAdapter
    except ImportError:
        return False
    send = HTTPAdapter.send

    def send_with_default_timeout(self, request, **kwargs):
        if kwargs.get("timeout") is None:
            kwargs["timeout"] = (CONNECT_TIMEOUT, READ_TIMEOUT)
        return send(self, request, **kwargs)

    HTTPAdapter.send = send_with_default_timeout
    return True


def main(argv: list[str]) -> None:
    if len(argv) < 2:
        raise SystemExit(__doc__.split("\n\n")[0])
    script = argv[1]
    install_default_timeout()
    sys.argv = [script, *argv[2:]]       # the script sees its own path and arguments, as if run directly
    runpy.run_path(script, run_name="__main__")


if __name__ == "__main__":
    main(sys.argv)
