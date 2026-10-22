"""STARsolo's run time in update_logs.py (needs pandas; run with pytest)."""

import os, sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "helper_py_scripts"))
from update_logs import star_run_time


def test_star_run_time():
    start = pd.Series([
        "2025-05-13 19:51:55",  # with a year
        "May/13/19:51:55",      # as update_logs.py writes STAR's log
        "Dec/31/23:00:00",      # over a new year, without a year
        "2025-05-13 10:00:00",  # over a day
        "",                     # missing
        "garbage",              # unparseable
    ])
    finish = pd.Series([
        "2025-05-14 01:34:05",
        "May/13/20:00:00",
        "Jan/01/01:30:00",
        "2025-05-14 15:15:00",
        "2025-05-13 10:00:00",
        "2025-05-13 10:00:00",
    ])
    assert list(star_run_time(start, finish)) == [
        "05:42:10", "00:08:05", "02:30:00", "29:15:00", "", "",
    ]
