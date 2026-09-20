from pathlib import Path
import os
import shutil
import sys


build_root = Path(os.environ["MESON_BUILD_ROOT"])
install_root = Path(
    os.environ.get("MESON_INSTALL_DESTDIR_PREFIX", os.environ["MESON_INSTALL_PREFIX"])
)
destination = install_root / sys.argv[1]
destination.mkdir(parents=True, exist_ok=True)
for name in ("disprs.mod", "dftd3.mod", "dftd4.mod", "multicharge.mod"):
    shutil.copy2(next(build_root.rglob(name)), destination)