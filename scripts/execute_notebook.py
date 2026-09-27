"""Execute this project's trusted notebook cells sequentially without a server.

Uses one in-process IPython session and records text and rich display outputs.
This is useful in environments where a separate Jupyter kernel cannot open its
communication sockets. Ordinary Jupyter/VS Code execution is also supported.
Only execute notebooks you trust: their cells contain Python code.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys

import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "notebooks/01_factorial_anova.ipynb"


def main() -> None:
    os.chdir(ROOT)
    notebook = nbformat.read(PATH, as_version=4)
    shell = InteractiveShell.instance()
    shell.display_formatter.active_types = ["text/plain", "text/html", "image/png", "image/svg+xml"]
    count = 0
    for cell in notebook.cells:
        if cell.cell_type != "code":
            continue
        count += 1
        cell.outputs = []
        cell.execution_count = count
        with capture_output(stdout=True, stderr=True, display=True) as captured:
            result = shell.run_cell(cell.source, store_history=False)
        result.raise_error()
        if captured.stdout:
            cell.outputs.append(nbformat.v4.new_output("stream", name="stdout", text=captured.stdout))
        if captured.stderr:
            cell.outputs.append(nbformat.v4.new_output("stream", name="stderr", text=captured.stderr))
        for output in captured.outputs:
            cell.outputs.append(nbformat.v4.new_output("display_data", data=output.data, metadata=output.metadata))
    notebook.metadata["execution_method"] = "All code cells executed sequentially in one in-process IPython session."
    notebook.metadata["language_info"]["version"] = sys.version.split()[0]
    nbformat.validate(notebook)
    nbformat.write(notebook, PATH)
    print(f"Executed {count} code cells successfully; notebook outputs saved.")


if __name__ == "__main__":
    main()
