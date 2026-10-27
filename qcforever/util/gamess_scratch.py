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
_SHELL_CONTROL = re.compile(
    r"^\s*(?:(?:if|else|endif|foreach|while|switch|case)\b|source(?:\s|$)|\.\s)"
)
_REPORTED_PATH = re.compile(
    r"^\s*GAMESS (temporary binary|supplementary output) files "
    r"will be written to[ \t]+([^\r\n]+)"
)


def read_scratch_log(log_path):
    """Read the paths reported by the running launcher, including partial logs."""
    paths = {}
    if log_path is None:
        return paths
    log_path = Path(log_path)
    try:
        with log_path.open(errors="replace") as log:
            for line in log:
                match = _REPORTED_PATH.match(line)
                if match is None:
                    continue
                name = "SCR" if match[1] == "temporary binary" else "USERSCR"
                path = Path(match[2].strip())
                if not path.is_absolute():
                    path = log_path.resolve().parent / path
                paths[name] = path.resolve()
    except OSError:
        pass
    return paths


def resolve_scratch_directories(rungms_path, work_directory=None, log_path=None):
    """Use reported paths, simple script settings, environment, then defaults.

    Runtime shell expressions and conditional assignments cannot be evaluated
    safely before launch. Defer those to the launcher's own execution and log.
    """
    reported = read_scratch_log(log_path)
    if len(reported) == 2:
        return tuple(reported[name] for name in ("SCR", "USERSCR"))
    work_directory = Path(work_directory or Path.cwd()).resolve()
    variables = dict(os.environ)

    def expand_path(value):
        value = _VARIABLE.sub(
            lambda match: variables.get(match[1] or match[2]) or match[0], value
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
    fallback = {name: variables[name] for name in ("SCR", "USERSCR")}

    if rungms_path is not None:
        with open(rungms_path, errors="replace") as script:
            for line in script:
                # Do not guess which branch a shell will take or execute a
                # sourced script merely to discover its settings.
                if _SHELL_CONTROL.match(line):
                    break
                match = _ASSIGNMENT.match(line)
                if match is None:
                    continue
                variable = next(name for name in match.groups()[:3] if name)
                try:
                    tokens = shlex.split(match[4], comments=True)
                    if not tokens or tokens == [""]:
                        continue
                    if len(tokens) != 1:
                        variables[variable] = None
                    else:
                        variables[variable] = expand_path(tokens[0])
                except ValueError:
                    # $$, local variables such as $PID, and command
                    # substitutions are resolved only by the real launcher.
                    variables[variable] = None

    return tuple(
        reported.get(name, Path(variables[name] or fallback[name]))
        for name in ("SCR", "USERSCR")
    )
