#!/usr/bin/env python3
"""Read-only post-deploy Shiny session smoke: two distinct browser contexts.

No patient identifier, dose, biological observation, or custom C++ is submitted.
A successful result proves websocket-connected session separation, NOT server
filesystem/log non-retention or any clinical correctness.
"""
import argparse
import sys
from playwright.sync_api import sync_playwright

PUBLIC = "https://tdmhub.shinyapps.io/MIPD_Engine/"


def connected(page, base, lang):
    response = page.goto(
        base.rstrip("/") + "/?lang=" + lang,
        wait_until="domcontentloaded",
        timeout=120000
    )
    if response is None or response.status != 200:
        raise AssertionError("Public Shiny page returned no HTTP 200 response")
    page.wait_for_function(
        """() => Boolean(window.Shiny?.shinyapp?.config?.sessionId &&
                         window.Shiny.shinyapp.$socket?.readyState === 1)""",
        timeout=120000
    )
    session = page.evaluate("() => window.Shiny.shinyapp.config.sessionId")
    if not isinstance(session, str) or not session:
        raise AssertionError("No real Shiny session ID")
    if page.locator(".shiny-output-error:visible").count():
        raise AssertionError("Shiny output error visible on startup")
    return session


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", default=PUBLIC)
    args = p.parse_args()
    if not args.base.startswith("https://"):
        p.error("The public app must use HTTPS")
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            first = browser.new_context(locale="fr-FR")
            second = browser.new_context(locale="en-US")
            a = first.new_page()
            b = second.new_page()
            a.set_default_timeout(120000)
            b.set_default_timeout(120000)
            first_id = connected(a,args.base,"fr")
            second_id = connected(b,args.base,"en")
            if first_id == second_id:
                raise AssertionError("Different browser contexts share a Shiny session ID")
            first.close()
            b.wait_for_function(
                "() => window.Shiny.shinyapp.$socket?.readyState === 1",
                timeout=15000
            )
            print("PUBLIC SHINY BROWSER PASS: FR and EN websocket sessions are distinct;")
            print("closing one session leaves the independent second connection alive.")
            return 0
        finally:
            browser.close()


if __name__ == "__main__":
    sys.exit(main())
