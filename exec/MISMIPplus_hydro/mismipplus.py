"""
MISMIP+ geometry module for BISICLES Python IBC.

Implements the standard MISMIP+ domain with a 24 km half-width fjord channel,
as described in Asay-Davis et al. (2016).

Reference:
    Asay-Davis, X.S. et al. (2016), Experimental design for three interrelated
    marine ice sheet and ocean model intercomparison projects,
    Geosci. Model Dev., 9, 2471-2497, doi:10.5194/gmd-9-2471-2016.

Domain:
    x in [0, 800 km]  –  x = 0 is the upstream ice divide
    y in [0,  80 km]  –  y = 40 km is the along-flow channel centreline

Usage (in BISICLES input file):
    geometry.problem_type       = Python
    PythonIBC.module            = mismipplus
    PythonIBC.thicknessFunction = thickness_init
    PythonIBC.topographyFunction = topography
    PythonIBC.bc_lo = 1 0
    PythonIBC.bc_hi = 0 0

    geometry.beta_type          = Python
    PythonBasalFriction.module  = mismipplus
    PythonBasalFriction.function = friction

    basalFlux.floating.module   = mismipplus
    basalFlux.floating.function = melt_spinup   # or melt_control / melt_ocean
"""

import math as m

# ---------------------------------------------------------------------------
# Bedrock polynomial coefficients (Asay-Davis et al. 2016 Table 1)
# ---------------------------------------------------------------------------
_B0 = -150.0    # m
_B2 = -728.8    # m
_B4 =  343.91   # m
_B6 =  -50.57   # m

# ---------------------------------------------------------------------------
# Fjord-channel cross-sectional parameters (24 km half-width)
# ---------------------------------------------------------------------------
_fc = 4.0e3    # wall steepness scale [m]
_dc = 500.0    # wall elevation amplitude [m]
_wc = 24.0e3   # channel half-width [m]

# ---------------------------------------------------------------------------
# Physical constants – must match values in the BISICLES input file
# ---------------------------------------------------------------------------
RHO_ICE   = 918.0    # kg m-3
RHO_OCEAN = 1028.0   # kg m-3

# ---------------------------------------------------------------------------
# Domain geometry parameters
# ---------------------------------------------------------------------------
X_CALVE  = 640.0e3   # m – calving front held fixed by large melt beyond here
MELT_CF  = 1.0e4     # m yr-1 – large melt rate beyond calving front


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _Bx(x):
    """Along-flow bed-elevation polynomial (centreline profile)."""
    xx  = x / 300.0e3
    xx2 = xx * xx
    xx4 = xx2 * xx2
    xx6 = xx4 * xx2
    return _B0 + _B2*xx2 + _B4*xx4 + _B6*xx6


def _By(yc):
    """
    Across-flow valley-wall elevation perturbation.

    yc : float
        y measured from the channel centreline (= y_domain - 40 km) [m].
    """
    return (  _dc / (1.0 + m.exp(-2.0*(yc - _wc) / _fc))
            + _dc / (1.0 + m.exp( 2.0*(yc + _wc) / _fc)))


def _cfsource(x):
    """
    Large melt rate applied beyond the calving front to pin it at X_CALVE.
    Returns a negative mass-balance rate (melting) in m yr-1.
    """
    return -MELT_CF if x > X_CALVE else 0.0


# ---------------------------------------------------------------------------
# Thickness initial condition
# ---------------------------------------------------------------------------

def thickness_init(x, y, *etc):
    """
    Uniform 100 m ice slab for x < X_CALVE, zero beyond.

    This thin initial condition quickly evolves to the grounded-ice steady
    state under accumulation.  For a faster spin-up, replace with a warmer
    (thicker) initial guess derived from shallow-ice theory.
    """
    return 100.0 if x < X_CALVE else 0.0


# ---------------------------------------------------------------------------
# Bedrock topography
# ---------------------------------------------------------------------------

def topography(x, y, *etc):
    """
    MISMIP+ bed topography: polynomial along-flow profile with a
    24 km half-width Gaussian valley-wall perturbation.

    Returns bed elevation in metres relative to sea level.
    The bed is floored at -720 m (MISMIP+ specification).
    """
    yc  = y - 40.0e3            # y relative to centreline
    bed = _Bx(x) + _By(yc)
    return max(bed, -720.0)     # MISMIP+ bedrock lower limit


# ---------------------------------------------------------------------------
# Basal friction coefficient
# ---------------------------------------------------------------------------

def friction(x, y, t, thck, topg, *etc):
    """
    Spatially uniform basal friction coefficient C = 10^4 Pa (m/yr)^{-1/3}.

    Appropriate for use with the Tsai (2015) pressure-limited sliding law
    and a Weertman power-law exponent m = 1/3.  Under the pressure-limited
    law the effective coefficient is C * (N / N_hydrostatic)^{2/3}, so the
    prescribed C here is an upper bound reached when N equals the full
    overburden (fully drained bed).
    """
    return 1.0e4


# ---------------------------------------------------------------------------
# Basal melt / mass balance functions (floating shelf only)
# ---------------------------------------------------------------------------

def melt_spinup(x, y, *etc):
    """
    Zero ocean melt with calving-front sink beyond X_CALVE.

    Use this for the BISICLES spin-up prior to coupling with SUHMO.
    The large melt beyond the calving front keeps the ice front fixed at
    x = 640 km during the transient spin-up.
    """
    return _cfsource(x)


def melt_control(x, y, t, thck, topg, *etc):
    """
    MISMIP+ 'melt0' control scenario: no ocean melt perturbation.

    Identical to melt_spinup.  Use this for a coupled run that continues
    from the spin-up geometry without any imposed ocean forcing.
    """
    return _cfsource(x)


def melt_ocean(x, y, t, thck, topg, *etc):
    """
    MISMIP+ 'melt4' parameterised ocean melt (Asay-Davis 2016, Eq. 9).

    Active for 10 yr < t < 110 yr.  The melt is proportional to the
    depth of the ice-shelf base below -100 m and tapers to zero at the
    grounding line via a tanh function over a 75 m scale.

    Equivalent to Omega = 0.2 m yr-1 m-1.
    """
    z0    = -100.0
    zb    = -thck * RHO_ICE / RHO_OCEAN    # ice-shelf base elevation
    Omega = -0.2
    mr    = 0.0
    if zb < z0 and 10.0 <= t < 110.0:
        wct = zb - topg
        mr  = -Omega * (zb - z0) * m.tanh(wct / 75.0)
    return mr + _cfsource(x)
