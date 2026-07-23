#!/usr/bin/env bash
# One-shot environment setup for this workspace.
#   - Python deps (numpy/scipy/sympy/mpmath) via pip
#   - Singular (needed only for the D2 system regeneration / certificates)
# Safe to re-run. Verifies each dependency at the end.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "== Python dependencies =="
python3 -m pip install --quiet -r "$here/requirements.txt"

echo "== Singular =="
if command -v Singular >/dev/null 2>&1; then
  echo "Singular already installed: $(command -v Singular)"
elif command -v apt-get >/dev/null 2>&1; then
  echo "installing Singular via apt-get (needs sudo/root)..."
  if [ "$(id -u)" = "0" ]; then
    apt-get update -qq && apt-get install -y -qq singular
  else
    sudo apt-get update -qq && sudo apt-get install -y -qq singular
  fi
else
  echo "WARNING: no apt-get; install Singular manually (only needed for D2 regeneration)."
fi

echo "== Verify =="
python3 - <<'PY'
import numpy, scipy, sympy, mpmath
print(f"numpy {numpy.__version__}, scipy {scipy.__version__}, "
      f"sympy {sympy.__version__}, mpmath {mpmath.__version__}")
PY
command -v Singular >/dev/null 2>&1 && Singular --version 2>/dev/null | head -1 || \
  echo "Singular: not available (D2 regeneration will be unavailable)"
echo "Setup complete."
