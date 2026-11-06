"""The launcher that gives the vendored annotation script's HTTP requests a timeout (paraphase_annotate/launch.py).
Each case runs the launcher in a subprocess against a local server that accepts connections and never answers --
the failure that hung a real run for 17 hours -- so no network and no real sample data are involved."""

import os
import socket
import subprocess
import sys
import textwrap
import threading
import time

import pytest

pytest.importorskip("requests")

from paraphase_annotate import launch  # noqa: E402

LAUNCHER = launch.__file__


@pytest.fixture
def silent_server():
    """A port that accepts connections and then says nothing."""
    srv = socket.socket()
    srv.bind(("127.0.0.1", 0))
    srv.listen(16)
    held = []

    def accept():
        try:
            while True:
                held.append(srv.accept()[0])
        except OSError:
            pass

    threading.Thread(target=accept, daemon=True).start()
    yield srv.getsockname()[1]
    srv.close()
    for c in held:
        c.close()


def _run(tmp_path, body, **env):
    script = tmp_path / "script.py"
    script.write_text(textwrap.dedent(body))
    started = time.time()
    done = subprocess.run(
        [sys.executable, LAUNCHER, str(script), "--flag", "value"],
        capture_output=True, text=True, timeout=60, cwd=tmp_path, env={**os.environ, **{k: str(v) for k, v in env.items()}},
    )
    return done, time.time() - started


def test_a_request_with_no_timeout_of_its_own_gives_up_instead_of_waiting_forever(tmp_path, silent_server):
    done, seconds = _run(tmp_path, f"""
        import requests
        try:
            requests.get("http://127.0.0.1:{silent_server}/")
        except requests.exceptions.RequestException as e:
            print(type(e).__name__)
    """, PARAPHASE_REQUEST_READ_TIMEOUT=1, PARAPHASE_REQUEST_CONNECT_TIMEOUT=1)
    assert done.stdout.strip() == "ReadTimeout", done.stderr
    assert seconds < 15


def test_the_scripts_own_timeout_is_left_alone(tmp_path, silent_server):
    done, seconds = _run(tmp_path, f"""
        import requests
        try:
            requests.get("http://127.0.0.1:{silent_server}/", timeout=1)
        except requests.exceptions.RequestException as e:
            print(type(e).__name__)
    """, PARAPHASE_REQUEST_READ_TIMEOUT=45)
    assert done.stdout.strip() == "ReadTimeout" and seconds < 20       # not the launcher's 45 s


def test_the_scripts_retry_loop_sees_a_timeout_like_any_other_failed_request(tmp_path, silent_server):
    # the same shape as the vendored script: a Session with urllib3 retries, which treats read timeouts as retryable
    done, seconds = _run(tmp_path, f"""
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
        session = requests.Session()
        session.mount("http://", HTTPAdapter(max_retries=Retry(total=1, backoff_factor=0)))
        try:
            session.request("GET", "http://127.0.0.1:{silent_server}/")
        except requests.exceptions.RequestException as e:
            print(type(e).__name__)
    """, PARAPHASE_REQUEST_READ_TIMEOUT=1, PARAPHASE_REQUEST_CONNECT_TIMEOUT=1)
    assert done.stdout.strip() in ("ConnectionError", "ReadTimeout") and seconds < 20


def test_the_script_runs_as_if_started_directly(tmp_path):
    done, _ = _run(tmp_path, """
        import sys
        print(__name__, sys.argv[0].endswith("script.py"), sys.argv[1:], __file__.endswith("script.py"))
    """)
    assert done.stdout.strip() == "__main__ True ['--flag', 'value'] True", done.stderr


def test_a_script_that_fails_still_fails(tmp_path):
    done, _ = _run(tmp_path, "raise SystemExit(3)")
    assert done.returncode == 3


def test_without_requests_the_launcher_does_nothing_and_says_so(monkeypatch):
    monkeypatch.setitem(sys.modules, "requests.adapters", None)      # makes `from requests.adapters import ...` fail
    assert launch.install_default_timeout() is False


def test_the_default_can_be_set_in_the_environment():
    done = subprocess.run(
        [sys.executable, "-c", f"import sys; sys.path.insert(0, {os.path.dirname(LAUNCHER)!r}); import launch; print(launch.CONNECT_TIMEOUT, launch.READ_TIMEOUT)"],
        capture_output=True, text=True, env={**os.environ, "PARAPHASE_REQUEST_CONNECT_TIMEOUT": "7", "PARAPHASE_REQUEST_READ_TIMEOUT": "33"},
    )
    assert done.stdout.split() == ["7.0", "33.0"]
