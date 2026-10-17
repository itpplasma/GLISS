"""Independent analytical finite-stencil truncation, without GLISS outputs."""
import json
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 80
profiles = [-mp.mpf(9) / 5, 4 - 4 * mp.sqrt(2), -mp.mpf(1)]
steps = [mp.mpf(1) / 16, mp.mpf(1) / 32, mp.mpf(1) / 64]


def coefficients(r, c):
    z = mp.sqrt(1 + c * (r * r - 1))
    btheta, k, m = r, -mp.mpf(2), mp.mpf(2)
    fieldline = m * btheta / r + k * z
    wave2 = m * m / r**2 + k * k
    cross_field = m * z / r + k * btheta
    a = fieldline * fieldline / wave2
    b = (k * k * z * z - (m / r)**2 * btheta * btheta) / (r * wave2)
    cc = (fieldline * fieldline + (z * z - btheta * btheta) / r**2
          - 2 * btheta / r - cross_field * cross_field / (r * r * wave2))
    return a, b, cc


rows = []
for c in profiles:
    alpha, beta = c * c / 2, -(2 + c)
    delta = 1 + 4 * beta / alpha
    estimates = []
    for h in steps:
        ap, bp, cp = coefficients(1 + h, c)
        am, bm, cm = coefficients(1 - h, c)
        ae = ((1 + h) * ap + (1 - h) * am) / (2 * h * h)
        be = ((1 + h) * cp + (1 - h) * cm) / 2
        be -= ((1 + h) * bp - (1 - h) * bm) / (2 * h)
        estimates.append((ae, be))
    ar = (4 * estimates[2][0] - estimates[1][0]) / 3
    br = (4 * estimates[2][1] - estimates[1][1]) / 3
    dr = 1 + 4 * br / ar
    rows.append({
        "c": str(c), "alpha": str(alpha), "beta": str(beta), "delta": str(delta),
        "central_estimates": [[str(a), str(b)] for a, b in estimates],
        "richardson": [str(ar), str(br), str(dr)],
        "richardson_errors": [str(abs(ar - alpha) / alpha), str(abs(br - beta)),
                              str(abs(dr - delta))],
        "central_error_ratios": [
            str(abs((estimates[0][j] - [alpha, beta][j])
                    / (estimates[1][j] - [alpha, beta][j]))) for j in range(2)],
    })

output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_suffix(".json")
output.write_text(json.dumps({"mpmath_dps": 80, "rows": rows}, indent=2) + "\n")
