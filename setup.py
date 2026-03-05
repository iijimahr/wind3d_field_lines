from __future__ import annotations

import importlib.machinery
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class F2pyBuildExt(build_ext):
    """Build the Fortran extension with f2py and place it in setuptools output."""

    def run(self) -> None:
        try:
            import numpy  # noqa: F401
        except ModuleNotFoundError as exc:
            raise RuntimeError("numpy is required. Please install it first.") from exc

        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        source = (
            Path(__file__).parent
            / "src"
            / "wind3d_field_lines"
            / "fortran"
            / "field_line_integrator.f90"
        )
        module_name = "_bbtobln"

        cmd = [
            sys.executable,
            "-m",
            "numpy.f2py",
            "-c",
            "-m",
            module_name,
            str(source),
        ]
        subprocess.check_call(cmd, cwd=build_temp)

        built_ext = self._find_built_extension(build_temp, module_name)
        if built_ext is None:
            raise RuntimeError("Built f2py extension artifact was not found.")

        target = Path(self.get_ext_fullpath("wind3d_field_lines._bbtobln"))
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built_ext, target)

    @staticmethod
    def _find_built_extension(build_temp: Path, module_name: str) -> Path | None:
        candidates = []
        for suffix in importlib.machinery.EXTENSION_SUFFIXES:
            candidates.extend(build_temp.glob(f"{module_name}*{suffix}"))
        if not candidates:
            return None
        return max(candidates, key=lambda p: p.stat().st_mtime)


setup(
    ext_modules=[Extension("wind3d_field_lines._bbtobln", sources=[])],
    cmdclass={"build_ext": F2pyBuildExt},
)
