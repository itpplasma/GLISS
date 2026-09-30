"""Write a VMEC input for the TERPSICHORE QAS3 benchmark family.

usage: make_input.py TEMPLATE OUTPUT PRESSURE_FACTOR CURRENT_FACTOR

TEMPLATE is input.QAS3_fx_cur75_bench from the TERPSICHORE repository. The
factors scale the pressure profile (PRES_SCALE) and the net toroidal current
CURTOR = -75 kA; the text wout is requested for vmecv92terps.
"""

import re
import sys

template, output, pressure, current = sys.argv[1:5]
text = open(template).read()
text = re.sub(r"(?im)^\s*(LWOUTTXT|PRES_SCALE)\s*=.*\n", "", text)
text = re.sub(
    r"(?im)^(\s*CURTOR\s*=\s*)\S+",
    lambda match: f"{match[1]}{-7.5e4 * float(current):.4e}",
    text,
)
text = re.sub(
    r"(?i)&INDATA\s*\n",
    f"&INDATA\n  LWOUTTXT = T\n  PRES_SCALE = {float(pressure)}\n",
    text,
    count=1,
)
# The template closes the namelist with "/" followed by a stray "&end".
text = re.sub(r"(?im)^\s*&end\s*$\n?", "", text)
open(output, "w").write(text)
