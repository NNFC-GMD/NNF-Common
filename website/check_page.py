"""Open a built page in headless Chrome and report whether it works.

The page passes when every ``--expect`` text appears before the timeout, the
browser reports no page error, console error or failed request, and no
notebook cell fails.

It needs the ``playwright`` package and an installed Google Chrome:

    pixi run -e e2e python website/check_page.py http://127.0.0.1:8000/hello/ --expect emscripten
"""

from __future__ import annotations

import argparse
import json
import time

from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
from playwright.sync_api import sync_playwright

# Console errors that marimo 0.25.1 logs on a page that works: its sidebar
# reports these while the page is being laid out.
HARMLESS_CONSOLE_ERRORS = (
    "[handleFillUpdated] component was expected to be defined",
    "[handleFillUnmount] component was expected to be defined",
)

# How a failed notebook cell shows up. In a read-only page marimo hides the
# cell and everything that depends on it, so the page itself says nothing;
# the error is written to the browser console, as an ordinary log message.
CELL_ERROR_MARKS = (
    '"type":"exception"',
    "Traceback (most recent call last)",
)


def check_page(
    url: str,
    expect: list[str],
    timeout_s: float = 240.0,
    screenshot: str | None = None,
    slider_steps: int = 0,
    settle_s: float = 5.0,
    width: int = 1280,
    height: int = 900,
) -> dict:
    """Load ``url`` and collect what a reader's browser would run into.

    Parameters
    ----------
    url
        Address of the page, served over HTTP.
    expect
        Texts that must become visible, in the order given.
    timeout_s
        Seconds to wait for each expected text.
    screenshot
        Path for a full-page PNG taken at the end, or None for no screenshot.
    slider_steps
        Number of right-arrow presses to send to the first slider, to check
        that the page reacts.
    settle_s
        Seconds to keep watching after the expected texts have appeared.
    width, height
        Size of the browser window in pixels.

    Returns
    -------
    dict
        ``ok``, the seconds until each expected text appeared, the texts that
        never appeared, and the page errors, console errors, failed requests
        and failed cells seen.
    """
    report: dict = {
        "url": url,
        "seconds_to": {},
        "missing": [],
        "page_errors": [],
        "console_errors": [],
        "failed_requests": [],
        "cell_errors": [],
    }
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(channel="chrome", headless=True)
        page = browser.new_page(viewport={"width": width, "height": height})
        page.on("pageerror", lambda error: report["page_errors"].append(str(error)))

        def on_console(message):
            text = message.text
            if any(mark in text for mark in CELL_ERROR_MARKS):
                report["cell_errors"].append(text[:500])
            elif message.type == "error" and not text.startswith(
                HARMLESS_CONSOLE_ERRORS
            ):
                report["console_errors"].append(text[:500])

        page.on("console", on_console)
        page.on(
            "requestfailed",
            lambda request: report["failed_requests"].append(
                f"{request.url} ({request.failure})"
            ),
        )
        page.on(
            "response",
            lambda response: (
                report["failed_requests"].append(f"{response.url} ({response.status})")
                if response.status >= 400
                else None
            ),
        )
        started = time.monotonic()
        page.goto(url, wait_until="domcontentloaded", timeout=timeout_s * 1000)
        for text in expect:
            try:
                # The page's text, not what is on screen: marimo does not lay
                # out cells that are far below the visible part of the page.
                page.wait_for_function(
                    "text => document.body.innerText.includes(text)",
                    arg=text,
                    timeout=timeout_s * 1000,
                )
                report["seconds_to"][text] = round(time.monotonic() - started, 1)
            except PlaywrightTimeoutError:
                report["missing"].append(text)
        if slider_steps:
            slider = page.get_by_role("slider").first
            slider.focus()
            for _ in range(slider_steps):
                page.keyboard.press("ArrowRight")
            page.wait_for_timeout(3000)
        # Give the cells that draw the last outputs time to finish, or to fail.
        page.wait_for_timeout(settle_s * 1000)
        if screenshot:
            page.screenshot(path=screenshot, full_page=True)
        browser.close()
    report["ok"] = not (
        report["missing"]
        or report["page_errors"]
        or report["console_errors"]
        or report["failed_requests"]
        or report["cell_errors"]
    )
    return report


def main(argv: list[str] | None = None) -> int:
    """Check one page. Returns 0 when it passes, 1 otherwise."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("url", help="address of the page, served over HTTP")
    parser.add_argument(
        "--expect",
        action="append",
        default=[],
        metavar="TEXT",
        help="text that must become visible; repeat for several",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=240.0,
        help="seconds to wait for each expected text (default: 240)",
    )
    parser.add_argument("--screenshot", help="path for a full-page PNG")
    parser.add_argument(
        "--slider-steps",
        type=int,
        default=0,
        help="right-arrow presses to send to the first slider (default: 0)",
    )
    args = parser.parse_args(argv)
    report = check_page(
        args.url,
        args.expect,
        timeout_s=args.timeout,
        screenshot=args.screenshot,
        slider_steps=args.slider_steps,
    )
    print(json.dumps(report, indent=4))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
