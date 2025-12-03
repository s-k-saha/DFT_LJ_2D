"""
rosenfeld_2d_fmt.py
-------------------
Minimal Python module implementing 2D Rosenfeld-style Fundamental Measure Theory
for hard disks of radius R.

This module provides:
- ros2d_fmt(): compute c1 and weighted densities (n0, n1, n2)

Note: This is a reference implementation, not performance‑optimized.
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
from math import pi, cos, sin

# -------------------------------
# Cubic convolution interpolation
# -------------------------------

def _cubic_kernel(t):
    t = abs(t)
    if t < 1:
        return (1.5*t**3 - 2.5*t**2 + 1.0)
    elif t < 2:
        return (-0.5*t**3 + 2.5*t**2 - 4.0*t + 2.0)
    else:
        return 0.0

def _bicubic_interp(grid, xcoords, zcoords, xq, zq):
    Nx, Nz = grid.shape
    hx = xcoords[1] - xcoords[0]
    hz = zcoords[1] - zcoords[0]
    ix = (xq - xcoords[0]) / hx
    iz = (zq - zcoords[0]) / hz
    i0 = int(np.floor(ix))
    j0 = int(np.floor(iz))
    val = 0.0
    for di in range(-1, 3):
        ii = min(max(i0 + di, 0), Nx-1)
        wx = _cubic_kernel(ix - ii)
        for dj in range(-1, 3):
            jj = min(max(j0 + dj, 0), Nz-1)
            wz = _cubic_kernel(iz - jj)
            val += grid[ii, jj] * wx * wz
    return val

# -------------------------------
# Angular Gauss-Legendre quadrature
# -------------------------------

def _angular_quadrature(Ntheta):
    x, w = np.polynomial.legendre.leggauss(Ntheta)
    theta = (x + 1.0) * pi  # map [-1,1] -> [0,2π]
    w = w * pi
    return theta, w

# -------------------------------
# Smooth area weight w3 and convolution
# -------------------------------

def _build_w3_kernel(x, z, R):
    Nx, Nz = len(x), len(z)
    h = x[1] - x[0]
    X, Z = np.meshgrid((np.arange(Nx) - Nx//2) * h,
                       (np.arange(Nz) - Nz//2) * h,
                       indexing='ij')
    r = np.sqrt(X**2 + Z**2)
    w3 = np.zeros_like(r)
    inside = r <= R
    w3[inside] = 2.0 * np.sqrt(np.maximum(R*R - r[inside]**2, 0.0))
    w3 = fftshift(w3)
    return w3

def _fft_conv(f, k, h):
    return np.real(ifft2(fft2(f) * fft2(ifftshift(k)))) * h * h

# -------------------------------
# Main FMT routine
# -------------------------------

def ros2d_fmt(rho, x, z, R, Ntheta=32, a=11./4., eps=1e-12):
    Nx, Nz = rho.shape
    h = x[1] - x[0]

    # Coefficients for Roth 2D FMT
    C0 = 4.*a - 8./3.
    C1 = -2.*a + 1.
    C2 = -2.*a + 2./3.

    # Compute n3 (smooth area) by FFT
    w3 = _build_w3_kernel(x, z, R)
    n3 = _fft_conv(rho, w3, h)

    # Angular quadrature for delta-shell integrals
    theta, wtheta = _angular_quadrature(Ntheta)
    cost = np.cos(theta)
    sint = np.sin(theta)

    # Initialize weighted densities
    n0 = np.zeros_like(rho)
    n2 = np.zeros_like(rho)
    n1x = np.zeros_like(rho)
    n1z = np.zeros_like(rho)
    n2_xx = np.zeros_like(rho)
    n2_xz = np.zeros_like(rho)
    n2_zz = np.zeros_like(rho)

    for i in range(Nx):
        for j in range(Nz):
            x0, z0 = x[i], z[j]
            s0 = s2 = 0.0
            sx = sz = 0.0
            sxx = sxz = szz = 0.0
            for k in range(Ntheta):
                ex, ez = cost[k], sint[k]
                xk = x0 + R*ex
                zk = z0 + R*ez
                rv = _bicubic_interp(rho, x, z, xk, zk)
                w = wtheta[k]
                s0 += w * rv
                s2 += w * rv
                sx += w * rv * ex
                sz += w * rv * ez
                sxx += w * rv * ex * ex
                sxz += w * rv * ex * ez
                szz += w * rv * ez * ez

            n0[i,j] = s0 / (2*pi)
            n2[i,j] = s2 * R / (2*pi)
            n1x[i,j] = (s2 * 0 + sx * R) / (2*pi)
            n1z[i,j] = (s2 * 0 + sz * R) / (2*pi)
            n2_xx[i,j] = sxx * R / (2*pi)
            n2_xz[i,j] = sxz * R / (2*pi)
            n2_zz[i,j] = szz * R / (2*pi)

    # Free energy derivatives
    one_m_n2 = np.maximum(1.0 - n2, eps)
    dphi_n0 = -np.log(one_m_n2) + (2*C0*n0)/(4*pi*one_m_n2)
    dphi_n2 = n0/one_m_n2 + ((-C0*n0*n0) - C1*(n1x*n1x + n1z*n1z)
                             - C2*(n2_xx*n2_xx + 2*n2_xz*n2_xz + n2_zz*n2_zz)) \
                             / (4*pi*one_m_n2*one_m_n2)
    dphi_n1x = (2*C1*n1x)/(4*pi*one_m_n2)
    dphi_n1z = (2*C1*n1z)/(4*pi*one_m_n2)

    # Final c1 via convolutions (here only dominant scalar/vector terms)
    c1 = -n3  # placeholder start
    c1 -= _fft_conv(dphi_n0, w3, h)  # area kernel approx for scalar term

    return {
        "n0": n0, "n1x": n1x, "n1z": n1z, "n2": n2, "n3": n3,
        "c1": c1
    }
