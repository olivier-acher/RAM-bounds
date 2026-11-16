"""Reading reflection curves from text files."""
import re

import numpy as np


def read_curve(path):
    """Read a two-column text file: frequency (GHz) and attenuation (dB).

    Separators may be commas, semicolons, tabs or spaces; lines that do not start with a number
    (headers, comments) are skipped; decimal commas are accepted when the separator is ';' or a tab.
    The second column may be the attenuation Gamma_dB (positive) or the reflection |R| in dB
    (negative): negative values are turned into attenuations.
    Returns (f, gamma_dB) as numpy arrays sorted by frequency.
    """
    rows = []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if not line or not re.match(r"^[-+]?(\d|\.\d)", line):
                continue
            if ";" in line or "\t" in line:
                parts = re.split(r"[;\t]+", line)
                parts = [p.replace(",", ".") for p in parts]
            else:
                parts = re.split(r"[,\s]+", line)
            vals = [float(p) for p in parts if p != ""]
            if len(vals) >= 2:
                rows.append(vals[:2])
    if not rows:
        raise ValueError("no numeric data found in %s" % path)
    data = np.array(rows)
    f, y = data[:, 0], data[:, 1]
    if np.median(y) < 0:
        y = -y
    order = np.argsort(f)
    return f[order], y[order]
