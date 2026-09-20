Damping functions
=================

Current two-body damping support in disprs:

.. list-table::
   :header-rows: 1
   :widths: 60 20 20

   * - Scheme
     - D3
     - D4
   * - Rational (BJ)
     - Yes
     - Yes
   * - Zero
     - Yes
     - No
   * - Modified zero
     - Yes
     - No
   * - Modified rational (BJ)
     - Yes
     - No
   * - Optimized power
     - Yes
     - No
   * - C6-only (CSO)
     - Yes
     - No
   * - Z
     - Yes
     - No

These entries describe implemented damping families, not the availability of fitted parameters for every density functional.
ATM three-body damping and real-space cutoff switching are separate from the functions below.

For atoms :math:`A` and :math:`B` separated by :math:`R=R_{AB}>0`, the two-body pair energy is

.. math::

   E_{AB}^{(2)} = -\sum_{n\in\{6,8\}} s_n \frac{C_n^{AB}}{R^n} f_n^{AB}(R).

Here :math:`s_n` are fitted scale factors.
Define :math:`\rho_{AB}=\sqrt{C_8^{AB}/C_6^{AB}}`, while :math:`R_{\mathrm{vdW}}^{AB}` denotes the tabulated D3 pair reference radius.
Distances and radius offsets are in bohr.
C6-only damping replaces the pair-energy expression as given below.


Rational damping function
-------------------------

.. math::

   R_0^{AB} &= a_1\rho_{AB}+a_2, \\
   f_{n,\mathrm{BJ}}^{AB}(R)
      &= \frac{R^n}{R^n+(R_0^{AB})^n}, \qquad n\in\{6,8\}.

The modified rational (modified BJ) scheme uses the same function with a different fitted parameter set.


Zero damping function
---------------------

.. math::

   f_{n,\mathrm{zero}}^{AB}(R)
      &= \left[1+6\left(\frac{s_{r,n}R_{\mathrm{vdW}}^{AB}}{R}
         \right)^{\alpha_n}\right]^{-1}, \\
   \alpha_6 &= \alpha, \qquad \alpha_8=\alpha+2.

The dimensionless radius multipliers :math:`s_{r,6}` and :math:`s_{r,8}` correspond to ``rs6`` and ``rs8``.
They are distinct from the energy scales :math:`s_6` and :math:`s_8`.


Modified zero damping function
------------------------------

Modified zero damping shifts the dimensionless distance before applying the power law:

.. math::

   f_{n,\mathrm{modified\ zero}}^{AB}(R)
      = \left[1+6\left(
         \frac{R}{s_{r,n}R_{\mathrm{vdW}}^{AB}}
         +\beta R_{\mathrm{vdW}}^{AB}
         \right)^{-\alpha_n}\right]^{-1}.

The exponents are again :math:`\alpha_6=\alpha` and :math:`\alpha_8=\alpha+2`.
Here :math:`\beta` has units of inverse bohr, :math:`\beta=0` recovers zero damping.


Optimized power damping function
--------------------------------

.. math::

   R_0^{AB} &= a_1\rho_{AB}+a_2, \\
   f_{n,\mathrm{OP}}^{AB}(R)
      &= \frac{R^{n+\beta}}{R^{n+\beta}+(R_0^{AB})^{n+\beta}},
         \qquad n\in\{6,8\}.

Here :math:`\beta` is a dimensionless additional exponent, unrelated to the inverse-length shift in modified zero damping.
Setting :math:`\beta=0` recovers rational damping.


C6-only damping function
------------------------

The CSO scheme has no explicit :math:`C_8/R^8` energy term.
It combines rational sixth-power damping with a distance-dependent scale:

.. math::

   R_0^{AB} &= a_3\rho_{AB}+a_4, \\
   S_{AB}(R) &= s_6+
      \frac{a_1}{1+\exp\left[(R-a_2\rho_{AB})/\ell\right]},
      \qquad \ell=1\,\mathrm{bohr}, \\
   f_{6,\mathrm{CSO}}^{AB}(R)
      &= \frac{R^6}{R^6+(R_0^{AB})^6}, \\
   E_{AB}^{(2)} &= -S_{AB}(R)\frac{C_6^{AB}}{R^6}
      f_{6,\mathrm{CSO}}^{AB}(R)
      = -\frac{S_{AB}(R)C_6^{AB}}{R^6+(R_0^{AB})^6}.

The sigmoid width :math:`\ell` is fixed, not a fitted parameter.
The coefficient :math:`C_8^{AB}` enters only through :math:`\rho_{AB}`.


Z damping function
------------------

.. math::

   \lambda_{AB} &= \frac{a_1}{Z_A+Z_B+2}, \\
   f_{n,\mathrm{Z}}^{AB}(R)
      &= \frac{R^n}{R^n+\lambda_{AB}C_n^{AB}},
         \qquad n\in\{6,8\}, \\
   E_{AB}^{(2)} &= -\sum_{n\in\{6,8\}}
      \frac{s_n C_n^{AB}}{R^n+\lambda_{AB}C_n^{AB}}.

Here :math:`a_1` is supplied in atomic units (inverse energy).
Z damping requires explicit parameters; it has no named parameter-table entry.