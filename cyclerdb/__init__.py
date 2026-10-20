# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""cyclerdb — search, validate and catalogue planetary cycler trajectories.

Note: pykep must be imported before pandas/pyarrow. With the reverse order the
interpreter can abort at exit ("corrupted double-linked list"), a shared
library teardown conflict. Importing this package first takes care of it.
"""
import pykep as _pk  # noqa: F401  (import order matters, see above)

__version__ = "0.1.0"
