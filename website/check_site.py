"""Serve the built site locally and check every page in headless Chrome.

For each notebook folder the page must load without a page error, a console
error or a failed request, and must show the text named by ``browser_text``
in its ``meta.toml``. That text should be something the page only shows once
Python has run in the browser, so a pre-computed preview does not satisfy it.

The site is served under a sub-path, as GitHub Pages serves a project site
(``https://<owner>.github.io/<repository>/``). A link that only works when
the site sits at the root of a server therefore fails here too.

Build the site first, then:

    pixi run -e e2e check-site
"""

from __future__ import annotations

import argparse
import functools
import http.server
import json
import sys
import threading
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import build  # noqa: E402
import check_page  # noqa: E402

# Stands in for the repository name in a GitHub Pages address.
SUB_PATH = "/project/"


class SubPathHandler(http.server.SimpleHTTPRequestHandler):
    """Serve files under ``SUB_PATH`` only, without logging every request."""

    def log_message(self, format, *args):
        pass

    def translate_path(self, path):
        if not path.startswith(SUB_PATH):
            # Outside the sub-path nothing exists, as on GitHub Pages.
            return str(Path(self.directory) / "outside-the-site")
        return super().translate_path("/" + path[len(SUB_PATH) :])


def serve(directory: Path) -> tuple[http.server.ThreadingHTTPServer, str]:
    """Start a server for ``directory`` on a free local port.

    Returns
    -------
    server : http.server.ThreadingHTTPServer
        Call ``shutdown()`` on it when done.
    address : str
        The base URL of the site, ending with a slash.
    """
    handler = functools.partial(SubPathHandler, directory=str(directory))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}{SUB_PATH}"


def main(argv: list[str] | None = None) -> int:
    """Check the site. Returns 0 when every page passes, 1 otherwise."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--site",
        type=Path,
        default=build.ROOT / "_site",
        help="folder holding the built site (default: _site)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=600.0,
        help="seconds to wait for each page's text (default: 600)",
    )
    parser.add_argument(
        "--screenshots",
        type=Path,
        help="folder for a full-page PNG of every page",
    )
    args = parser.parse_args(argv)
    if not (args.site / "index.html").is_file():
        parser.error(f"{args.site} holds no index.html; run `pixi run build` first")
    if args.screenshots:
        args.screenshots.mkdir(parents=True, exist_ok=True)

    server, address = serve(args.site)
    reports = []
    try:
        targets = [("index", address, ["Notebooks and tutorials"])]
        for page in build.load_pages():
            expect = [page["browser_text"]] if "browser_text" in page else []
            targets.append((page["slug"], f"{address}{page['slug']}/", expect))
        for name, url, expect in targets:
            screenshot = (
                str(args.screenshots / f"{name}.png") if args.screenshots else None
            )
            report = check_page.check_page(
                url, expect, timeout_s=args.timeout, screenshot=screenshot
            )
            report["page"] = name
            reports.append(report)
            print(
                f"{'ok  ' if report['ok'] else 'FAIL'} {name}: {report['seconds_to']}",
                flush=True,
            )
    finally:
        server.shutdown()
    failed = [report for report in reports if not report["ok"]]
    if failed:
        print(json.dumps(failed, indent=4))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
