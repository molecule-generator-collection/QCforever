"""Resolve the scratch paths shared by GAMESS execution and cleanup."""

import os
import re
import shlex
from pathlib import Path


_ASSIGNMENT = re.compile(
    r"^\s*(?:set\s+(SCR|USERSCR)\s*=\s*"
    r"|setenv\s+(SCR|USERSCR)\s+"
    r"|export\s+(SCR|USERSCR)\s*=\s*)(.*)$"
)
_VARIABLE = re.compile(r"\$(?:\{(\w+)\}|(\w+))")


def resolve_scratch_directories(rungms_path, work_directory=None):
    """Use script assignments, inherited variables, then job-local defaults.

    This reads simple assignments without executing the rungms shell script.
    The returned absolute paths can also be passed to the child environment.
    """
    work_directory = Path(work_directory or Path.cwd()).resolve()
    variables = dict(os.environ)

    def expand_path(value):
        value = _VARIABLE.sub(
            lambda match: variables.get(match[1] or match[2], match[0]), value
        )
        if "$" in value or "`" in value:
            raise ValueError(
                f"Cannot resolve GAMESS scratch path {value!r}; use a simple "
                "path in rungms or export SCR and USERSCR before starting Python."
            )
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = work_directory / path
        return str(path.resolve())

    for variable in ("SCR", "USERSCR"):
        default = work_directory / ".gamess_scratch" / variable.lower()
        variables[variable] = expand_path(variables.get(variable) or str(default))

    if rungms_path is not None:
        with open(rungms_path, errors="replace") as script:
            for line in script:
                match = _ASSIGNMENT.match(line)
                if match is None:
                    continue
                variable = next(name for name in match.groups()[:3] if name)
                tokens = shlex.split(match[4], comments=True)
                if not tokens or tokens == [""]:
                    continue
                if len(tokens) != 1:
                    raise ValueError(
                        f"Cannot resolve GAMESS {variable} assignment: {line.strip()}"
                    )
                variables[variable] = expand_path(tokens[0])

    return tuple(Path(variables[name]) for name in ("SCR", "USERSCR"))
