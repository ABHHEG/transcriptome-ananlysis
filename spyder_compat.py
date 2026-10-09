"""Small console interface that lets the Streamlit workflow run from Spyder."""
from pathlib import Path
import os
import re


class _ConsoleProgress:
    def __init__(self):
        self._last = None

    def progress(self, value):
        try:
            percent = int(float(value) * 100)
        except (TypeError, ValueError):
            percent = 0
        if percent != self._last:
            print(f"Progress: {percent}%")
            self._last = percent

    def write(self, value=""):
        print(value)

    def success(self, value=""):
        print(f"SUCCESS: {value}")


class _Column:
    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


class _ConsoleUI:
    """Implements the small subset of Streamlit calls used by app.py."""
    def __init__(self):
        self.session_state = {}
        self.secrets = {
            "NCBI_API_KEY": os.getenv("NCBI_API_KEY", ""),
            "NCBI_EMAIL": os.getenv("NCBI_EMAIL", ""),
        }
        self.output_dir = Path(__file__).resolve().parent / "spyder_outputs"
        self._web_windows = []

    def set_page_config(self, **kwargs):
        pass

    def cache_data(self, func=None, **kwargs):
        # Streamlit's cache is optional for a single local Spyder run.
        def decorate(target):
            return target
        return decorate(func) if func is not None else decorate

    def title(self, value): print(f"\n{value}\n{'=' * 70}")
    def header(self, value): print(f"\n{value}\n{'-' * 70}")
    def subheader(self, value): print(f"\n{value}")
    def caption(self, value): print(value)
    def write(self, *values, **kwargs):
        for value in values:
            if value is not None:
                print(value)

    def info(self, value): print(f"INFO: {value}")
    def warning(self, value): print(f"WARNING: {value}")
    def error(self, value): print(f"ERROR: {value}")
    def success(self, value): print(f"SUCCESS: {value}")
    def code(self, value, **kwargs): print(value)

    def markdown(self, value, **kwargs):
        # Print a readable version of Markdown links in the console.
        print(re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1: \2", str(value)))

    def save_dataframe(self, name, value):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        destination = self.output_dir / f"{name}.csv"
        value.to_csv(destination, index=False)
        print(f"Saved full table: {destination}")

    def dataframe(self, value, **kwargs):
        if hasattr(value, "shape"):
            rows, cols = value.shape
            print(f"DataFrame: {rows} rows × {cols} columns")
            if rows > 40:
                print(value.head(40).to_string(index=False))
                print(f"(Showing first 40 rows; full table has {rows} rows.)")
            else:
                print(value.to_string(index=False))
        else:
            print(value)

    def selectbox(self, label, options, index=0, format_func=None, **kwargs):
        options = list(options)
        if not options:
            raise ValueError(f"No choices available for: {label}")
        print(f"\n{label}")
        for number, option in enumerate(options, start=1):
            shown = format_func(option) if format_func else option
            print(f"  {number}. {shown}")
        default_index = min(max(index, 0), len(options) - 1)
        raw = input(f"Choose 1-{len(options)} [default {default_index + 1}]: ").strip()
        if not raw:
            return options[default_index]
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1]
        for option in options:
            if str(option) == raw:
                return option
        print("Invalid choice; using the default.")
        return options[default_index]

    def slider(self, label, min_value, max_value, value, step=None, **kwargs):
        raw = input(f"{label} [{value}]: ").strip()
        if not raw:
            return value
        try:
            answer = type(value)(raw)
            if min_value <= answer <= max_value:
                return answer
        except (ValueError, TypeError):
            pass
        print(f"Invalid value; using {value}.")
        return value

    def button(self, label, **kwargs):
        answer = input(f"\n{label}? [y/N]: ").strip().lower()
        return answer in {"y", "yes"}

    def columns(self, spec, **kwargs):
        count = spec if isinstance(spec, int) else len(spec)
        return [_Column() for _ in range(count)]

    def metric(self, label, value, delta=None, **kwargs):
        extra = f" ({delta})" if delta is not None else ""
        print(f"{label}: {value}{extra}")

    def progress(self, value=0):
        progress = _ConsoleProgress()
        progress.progress(value)
        return progress

    def empty(self):
        return _ConsoleProgress()


    def open_kegg_pathway(self, url):
        """Open a KEGG pathway in a Qt WebEngine window attached to Spyder's Qt session."""
        try:
            from spyder_analysis import render_colored_kegg_map
            render_colored_kegg_map(url, self.output_dir)
        except Exception as exc:
            print(f"Could not render KEGG image in the Spyder Plots pane: {exc}")
        try:
            from qtpy.QtCore import QUrl
            from qtpy.QtWidgets import QApplication, QMainWindow
            from qtpy.QtWebEngineWidgets import QWebEngineView
            app = QApplication.instance()
            if app is None:
                print("Spyder's Qt application was not found; cannot open the embedded KEGG viewer.")
                print(url)
                return False
            window = QMainWindow()
            window.setWindowTitle("KEGG Pathway Viewer")
            window.resize(1200, 850)
            view = QWebEngineView(window)
            view.setUrl(QUrl(url))
            window.setCentralWidget(view)
            window.show()
            self._web_windows.append(window)
            print("Opened the interactive KEGG map in a Qt viewer from this Spyder session.")
            return True
        except Exception as exc:
            print("Could not start Spyder's embedded KEGG viewer. Check that QtWebEngine is installed in Spyder's Python environment.")
            print(f"Viewer error: {exc}")
            print(f"KEGG map URL: {url}")
            return False
    def download_button(self, label, data, file_name, **kwargs):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        destination = self.output_dir / file_name
        if isinstance(data, str):
            data = data.encode("utf-8")
        destination.write_bytes(data)
        print(f"{label}: saved {destination}")
        return destination

    def stop(self):
        raise SystemExit("Stopped because the input or selected options are invalid.")


st = _ConsoleUI()









