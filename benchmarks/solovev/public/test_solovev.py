"""Independent analytic checks; run with python test_solovev.py."""
from decimal import Decimal, localcontext
import unittest

import numpy as np

import solovev


def series_reference(psi_n, q0):
    """80-digit hypergeometric oracle, independent of contour quadrature.

    q/q0 = 2F1(3/4, 5/4; 1; 4 a² psi_n / R0²).
    Phi = pi F e a² psi_n/R0 * 2F1(3/4, 5/4; 2; 4 a² psi_n/R0²).
    Constants are the GPEC sol_run parameters, not measured generator values.
    """
    with localcontext() as context:
        context.prec = 80
        d = Decimal
        pi = d("3.141592653589793238462643383279502884197169399375105820"
               "974944592307816406286208998628")
        a, elongation, major, field = d("0.33"), d("1.6"), d(1), d(1)
        psi = d(str(psi_n))
        argument = 4 * a**2 * psi / major**2

        def hypergeometric(c):
            total = term = d(1)
            for n in range(1, 350):
                term *= ((d(n) - d("0.25")) * (d(n) + d("0.25"))
                         * argument / (d(n) * (d(n) - 1 + c)))
                total += term
            return total

        q = d(str(q0)) * hypergeometric(d(1))
        flux = pi * field * elongation * a**2 * psi / major * hypergeometric(d(2))
        return float(q), float(flux)


class SolovevTests(unittest.TestCase):
    def test_axis_and_analytic_surfaces(self):
        for q0 in (1.035, 1.5, 2.0):
            self.assertEqual(solovev.safety_factor_and_flux(0.0, q0), (q0, 0.0))
            previous_flux = 0.0
            for psi in (1.0e-8, 0.01, 0.4, 1.0):
                with self.subTest(q0=q0, psi_n=psi):
                    q, flux = solovev.safety_factor_and_flux(psi, q0)
                    reference_q, reference_flux = series_reference(psi, q0)
                    self.assertLess(abs(q - reference_q), 1.0e-13)
                    self.assertLess(abs(flux - reference_flux), 1.0e-14)
                    self.assertGreater(flux, previous_flux)
                    previous_flux = flux

    def test_angular_quadrature_convergence(self):
        reference = np.array(series_reference(1.0, 1.5))
        coarse = np.abs(np.array(solovev.safety_factor_and_flux(1.0, 1.5, 16))
                        - reference)
        fine = np.abs(np.array(solovev.safety_factor_and_flux(1.0, 1.5, 64))
                      - reference)
        self.assertTrue(np.all(coarse > 1.0e-8))
        self.assertTrue(np.all(fine < 1.0e-14))
        self.assertTrue(np.all(fine < 1.0e-4 * coarse))

    def test_profile_scaling(self):
        first = solovev.profiles(1.5)
        second = solovev.profiles(3.0)
        self.assertTrue(np.all(np.diff(first[0]) > 0.0))
        self.assertTrue(np.all(np.diff(first[1]) < 0.0))
        np.testing.assert_allclose(first[0], second[0], rtol=0.0, atol=1.0e-14)
        np.testing.assert_allclose(first[1], 2.0 * second[1], rtol=1.0e-14)
        np.testing.assert_allclose(first[2], 4.0 * second[2], rtol=1.0e-14)
        self.assertEqual(first[2][-1], 0.0)
        self.assertLess(abs(first[4] - series_reference(1.0, 1.5)[1]), 1.0e-14)


if __name__ == "__main__":
    unittest.main()
