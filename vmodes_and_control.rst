V-Modes and the OFC Control Loop
=================================

This chapter provides a mathematical description of the Optical Feedback
Control (OFC) pipeline implemented in
`ts_ofc <https://github.com/lsst-ts/ts_ofc>`_ and
`ts_mtaos <https://github.com/lsst-ts/ts_mtaos>`_.  It covers the forward
optical model, the state estimation inverse problem, the definition and
physical interpretation of **v-modes**, and how PID gains are applied in
both DOF space and v-mode space.


The forward optical model
-------------------------

The relationship between the telescope's degrees of freedom (DOFs) and the
resulting wavefront error across the focal plane is captured by the
**sensitivity matrix**.  For each sensor position on the focal plane, the
sensitivity matrix encodes how a unit change in each DOF affects each
Zernike coefficient of the wavefront.

If we stack the Zernike coefficients from all sensors into a single vector,
the forward model is:

.. math::

   \mathbf{y} = \mathbf{A}\,\mathbf{x}

where:

- :math:`\mathbf{y} \in \mathbb{R}^{m}` is the wavefront error vector
  (Zernike coefficients across all sensors),
  with :math:`m = n_\text{sensors} \times n_\text{Zernikes}`.
  For four corner wavefront sensors and Zernikes Z4--Z22 (19 coefficients),
  :math:`m = 4 \times 19 = 76`.

- :math:`\mathbf{A} \in \mathbb{R}^{m \times n}` is the sensitivity matrix,
  evaluated at the sensor field positions and the current camera rotation
  angle.

- :math:`\mathbf{x} \in \mathbb{R}^{n}` is the DOF state vector, with
  :math:`n = n_\text{used DOFs}`.  A typical configuration uses
  :math:`n = 22` DOFs: 5 M2 hexapod, 5 camera hexapod, 7 M1M3 bending
  modes, and 5 M2 bending modes.

The sensitivity matrix :math:`\mathbf{A}` is constructed from a double
Zernike expansion of the optical model and is evaluated at each sensor's
field angle.  It depends on the camera rotation angle because the Zernike
basis rotates with the instrument.


DOF normalization
-----------------

The physical DOFs span different scales: hexapod translations are in
micrometers, hexapod rotations in degrees, and bending mode coefficients are
dimensionless with varying magnitudes.  To give each DOF comparable
influence during the inversion, we introduce a diagonal normalization
matrix:

.. math::

   \mathbf{N} = \operatorname{diag}\!\bigl(w_1, w_2, \ldots, w_n\bigr)

where the weights :math:`w_i` are read from the
``normalization_weights`` configuration file (for details see
`Megias Homar, G. et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024ApJ...974..108M>`_, Eqs. 9--11).

The **normalized sensitivity matrix** is:

.. math::

   \tilde{\mathbf{A}} = \mathbf{A}\,\mathbf{N}

so that the forward model in normalized DOF coordinates
:math:`\tilde{\mathbf{x}} = \mathbf{N}^{-1}\mathbf{x}` becomes:

.. math::

   \mathbf{y} = \tilde{\mathbf{A}}\,\tilde{\mathbf{x}}


State estimation: the inverse problem
--------------------------------------

Given a measured wavefront error :math:`\mathbf{y}_\text{meas}`, we want
to estimate the DOF state :math:`\hat{\mathbf{x}}`.  This is an
overdetermined linear inverse problem (:math:`m > n`), which can be solved
via regularized pseudo-inversion.  There are two main approaches:
**SVD truncation** and **noise-weighted least squares**.  The OFC combines
both for optimal state estimation.


.. _Problem setup:

Problem setup
^^^^^^^^^^^^^

The residual wavefront error used as input is:

.. math::

   \mathbf{y} = \mathbf{y}_\text{meas} - \mathbf{y}_\text{intrinsic}
                - \mathbf{y}_2

These subtractions operate in **measurement space** (Zernike coefficients at
each sensor, :math:`\mathbf{y} \in \mathbb{R}^m`):

- :math:`\mathbf{y}_\text{intrinsic}` contains the design (intrinsic) Zernike
  coefficients for the current filter and sensor positions.  These represent
  the wavefront error that the optical design produces even when the optics are
  in perfect alignment (:math:`\mathbf{x} = \mathbf{0}`).  When
  ``subtract_intrinsics = True`` (the default for the main closed-loop path),
  these are subtracted so the residual :math:`\mathbf{y}` reflects only the
  misalignment component.

- :math:`\mathbf{y}_2` is a **static per-sensor correction** in measurement
  space.  In the deployed configuration (``lsst_y2.yaml``), only the Z4
  (focus) entries for the corner wavefront sensors are non-zero (e.g.,
  R00: -0.139 um, R44: -0.145 um).  These values were determined from on-sky
  data as the focus offsets at which the best image quality is achieved.

The purpose of these subtractions is to ensure that when the telescope is at
the desired alignment state (:math:`\mathbf{x} = \mathbf{0}` in DOF space),
the residual :math:`\mathbf{y} \approx \mathbf{0}`.  The state estimator
then inverts :math:`\mathbf{y}` to produce the DOF estimate :math:`\hat{\mathbf{x}}`.
The PID controller operates in **DOF space** (:math:`\mathbf{x} \in \mathbb{R}^n`)
with a setpoint of zero (see `PID control`_ below).  So the full chain is:

1. Subtract known contributions in y-space:
   :math:`\mathbf{y} = \mathbf{y}_\text{meas} - \mathbf{y}_\text{intrinsic} - \mathbf{y}_2`
2. Estimate the DOF state: :math:`\hat{\mathbf{x}} = f(\mathbf{y})`
3. PID drives :math:`\hat{\mathbf{x}} \to \mathbf{0}` (the setpoint in x-space)

The y-space subtractions and the x-space setpoint work together: the
subtractions define what "zero misalignment" looks like in measurement space,
and the PID targets zero in DOF space.  If the y2 values were wrong (e.g.,
the focus offsets didn't correspond to best image quality), the system would
converge to a suboptimal state even though the PID error would be zero.

We seek to solve:

.. math::

   \mathbf{y} = \tilde{\mathbf{A}}\,\tilde{\mathbf{x}} + \mathbf{n}

where :math:`\tilde{\mathbf{A}}` is the normalized sensitivity matrix,
:math:`\tilde{\mathbf{x}}` is the normalized DOF state, and
:math:`\mathbf{n}` is measurement noise with covariance
:math:`\mathbf{C}_n = \mathbb{E}[\mathbf{n}\mathbf{n}^T]`.

Here :math:`\mathbb{E}[\cdot]` denotes the **expectation operator** (statistical
mean), and :math:`\mathbf{n}\mathbf{n}^T` is the outer product of the noise
vector with itself.  The covariance matrix :math:`\mathbf{C}_n` therefore
contains:

- **Diagonal entries** :math:`(\mathbf{C}_n)_{ii} = \mathbb{E}[n_i^2] = \sigma_i^2`:
  the variance of the :math:`i`-th Zernike measurement.
- **Off-diagonal entries** :math:`(\mathbf{C}_n)_{ij} = \mathbb{E}[n_i n_j]`:
  the covariance between Zernikes :math:`i` and :math:`j`.

In practice, :math:`\mathbf{C}_n` is estimated from stability data by computing
the **sample covariance** of the residual Zernikes across :math:`N` exposures:

.. math::

   \hat{\mathbf{C}}_n = \frac{1}{N-1}
   \sum_{k=1}^{N} (\mathbf{y}_k - \bar{\mathbf{y}})
                     (\mathbf{y}_k - \bar{\mathbf{y}})^T

where :math:`\mathbf{y}_k` are the measured Zernike vectors and
:math:`\bar{\mathbf{y}}` is their mean.  This matrix tells the state estimator
which Zernikes and corners are noisy (large diagonal entries) and which are
correlated (large off-diagonal entries), enabling appropriate down-weighting
during state estimation.  See `SITCOMTN-177 <https://sitcomtn-177.lsst.io>`_
for details on how this covariance is estimated from on-sky data.


Approach 1: SVD with truncation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The standard least-squares solution uses the pseudo-inverse.  Start by
computing the singular value decomposition (SVD):

.. math::

   \tilde{\mathbf{A}} = \mathbf{U}\,\operatorname{diag}(\boldsymbol{\Sigma})\,\mathbf{V}_h

In the **full SVD** of an :math:`m \times n` matrix with :math:`m > n`:

- :math:`\mathbf{U} \in \mathbb{R}^{m \times m}` is an orthogonal matrix (left singular vectors),
- :math:`\boldsymbol{\Sigma} \in \mathbb{R}^{n}` is a vector of :math:`n` singular values
  :math:`(\sigma_1, \ldots, \sigma_n)`, and :math:`\operatorname{diag}(\boldsymbol{\Sigma})`
  is the :math:`m \times n` rectangular matrix with these values on the main diagonal,
- :math:`\mathbf{V}_h \in \mathbb{R}^{n \times n}` is an orthogonal matrix (right singular vectors).

The singular values are ordered :math:`\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n \geq 0`.

In practice, the **thin (economy) SVD** is used, which discards the extra
columns of :math:`\mathbf{U}` corresponding to zero rows of the diagonal matrix:

- :math:`\mathbf{U} \in \mathbb{R}^{m \times n}` — only the first :math:`n` columns,
- :math:`\boldsymbol{\Sigma} \in \mathbb{R}^{n}` — vector of :math:`n` singular values
  :math:`(\sigma_1, \sigma_2, \ldots, \sigma_n)`,
- :math:`\mathbf{V}_h \in \mathbb{R}^{n \times n}` — unchanged.

In this case, :math:`\operatorname{diag}(\boldsymbol{\Sigma})` is an
:math:`n \times n` square diagonal matrix.

The pseudo-inverse is:

.. math::

   \tilde{\mathbf{A}}^{+} = \mathbf{V}_h^T\,\operatorname{diag}(1/\boldsymbol{\Sigma})\,\mathbf{U}^T

where :math:`\operatorname{diag}(1/\boldsymbol{\Sigma}) = \operatorname{diag}(1/\sigma_1, \ldots, 1/\sigma_n)`.

**Truncation regularization.**  Small singular values amplify noise because
:math:`1/\sigma_i \to \infty` as :math:`\sigma_i \to 0`.  To regularize,
truncate at index :math:`k` by keeping only the first :math:`k` singular values:

.. math::

   \tilde{\mathbf{A}}_k^{+}
   = \mathbf{V}_{h,k}^T\,\operatorname{diag}(1/\boldsymbol{\Sigma}_k)\,\mathbf{U}_k^T

where :math:`\mathbf{U}_k` contains the first :math:`k` columns of :math:`\mathbf{U}`,
:math:`\mathbf{V}_{h,k}` contains the first :math:`k` rows of :math:`\mathbf{V}_h`,
and :math:`\boldsymbol{\Sigma}_k` contains the first :math:`k` elements of
:math:`\boldsymbol{\Sigma}`.

The truncated least-squares solution is:

.. math::

   \tilde{\mathbf{x}}_\text{trunc} = \tilde{\mathbf{A}}_k^{+}\,\mathbf{y}

This discards the :math:`n - k` least-constrained DOF combinations,
retaining only directions in DOF space that produce measurable wavefront
changes above the noise floor.

**Concrete example.**  Consider the typical OFC configuration with
:math:`m = 76` measurements (4 corners × 19 Zernikes) and :math:`n = 22` DOFs.
The thin SVD of :math:`\tilde{\mathbf{A}} \in \mathbb{R}^{76 \times 22}` produces:

- :math:`\mathbf{U} \in \mathbb{R}^{76 \times 22}` — 22 left singular vectors (columns)
- :math:`\boldsymbol{\Sigma} \in \mathbb{R}^{22}` — vector of 22 singular values :math:`(\sigma_1, \ldots, \sigma_{22})`
- :math:`\mathbf{V}_h \in \mathbb{R}^{22 \times 22}` — 22 right singular vectors (the v-modes)

With ``truncation_index = 12``, we keep only the first :math:`k = 12` modes:

- :math:`\mathbf{U}_{12} \in \mathbb{R}^{76 \times 12}` — first 12 columns of U
- :math:`\boldsymbol{\Sigma}_{12} \in \mathbb{R}^{12}` — first 12 singular values :math:`(\sigma_1, \ldots, \sigma_{12})`
- :math:`\mathbf{V}_{h,12} \in \mathbb{R}^{12 \times 22}` — first 12 rows of Vₕ

The truncated pseudo-inverse is then:

.. math::

   \tilde{\mathbf{A}}_{12}^{+}
   = \underbrace{\mathbf{V}_{h,12}^T}_{22 \times 12}
     \underbrace{\operatorname{diag}(1/\boldsymbol{\Sigma}_{12})}_{12 \times 12}
     \underbrace{\mathbf{U}_{12}^T}_{12 \times 76}
   \;\in\; \mathbb{R}^{22 \times 76}

The final state estimate is:

.. math::

   \underbrace{\tilde{\mathbf{x}}}_{22 \times 1}
   = \underbrace{\tilde{\mathbf{A}}_{12}^{+}}_{22 \times 76}
     \cdot \underbrace{\mathbf{y}}_{76 \times 1}

The output :math:`\tilde{\mathbf{x}}` has 22 components (one per DOF), but only
12 independent directions are actually constrained by the data.  Internally,
the matrix multiplication proceeds right-to-left:

1. :math:`\mathbf{U}_{12}^T \cdot \mathbf{y}` projects the 76 measurements onto
   12 left singular vectors, producing a 12-element vector.
2. :math:`\operatorname{diag}(1/\boldsymbol{\Sigma}_{12})` scales each component
   by the inverse of its singular value.
3. :math:`\mathbf{V}_{h,12}^T` expands this 12-element vector back to 22 DOFs
   using the first 12 v-modes.

The result lies entirely in the subspace spanned by v-modes 0--11.  The
coefficients of v-modes 12--21 are implicitly zero, preventing noise
amplification from small singular values :math:`\sigma_{13}` through
:math:`\sigma_{22}`.


Approach 2: Noise-weighted pseudo-inverse
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the measurement noise covariance :math:`\mathbf{C}_n` is known and
non-uniform, we want to minimize the weighted least-squares objective:

.. math::

   \min_{\tilde{\mathbf{x}}} \;
   \bigl(\mathbf{y} - \tilde{\mathbf{A}}\,\tilde{\mathbf{x}}\bigr)^T
   \mathbf{C}_n^{-1}
   \bigl(\mathbf{y} - \tilde{\mathbf{A}}\,\tilde{\mathbf{x}}\bigr)

The solution is:

.. math::

   \tilde{\mathbf{x}}_\text{wls}
   = \bigl(\tilde{\mathbf{A}}^T \mathbf{C}_n^{-1} \tilde{\mathbf{A}}\bigr)^{-1}
     \tilde{\mathbf{A}}^T \mathbf{C}_n^{-1} \mathbf{y}

Equivalently, **whiten** the problem by defining:

.. math::

   \bar{\mathbf{y}} = \mathbf{C}_n^{-1/2}\,\mathbf{y}, \qquad
   \bar{\mathbf{A}} = \mathbf{C}_n^{-1/2}\,\tilde{\mathbf{A}}

Then solve the standard (unweighted) least-squares:

.. math::

   \tilde{\mathbf{x}}_\text{wls} = \bar{\mathbf{A}}^{+}\,\bar{\mathbf{y}}

The noise covariance down-weights Zernike modes and sensors with high
measurement uncertainty.  See `SITCOMTN-177 <https://sitcomtn-177.lsst.io>`_
for details on how :math:`\mathbf{C}_n` is estimated from on-sky stability data.


Combined approach (as implemented in ts_ofc)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``StateEstimator`` in ts_ofc combines both approaches:

**Step 1 — Whiten by noise covariance:**

.. math::

   \bar{\mathbf{y}} = \mathbf{C}_n^{-1/2}\,\mathbf{y}, \qquad
   \bar{\mathbf{A}} = \mathbf{C}_n^{-1/2}\,\tilde{\mathbf{A}}

**Step 2 — Truncated pseudo-inverse of whitened system:**

.. math::

   \tilde{\mathbf{x}} = \bar{\mathbf{A}}_k^{+}\,\bar{\mathbf{y}}

where :math:`\bar{\mathbf{A}}_k^{+}` is the truncated pseudo-inverse of the
whitened sensitivity matrix (SVD with :math:`k` retained singular values).

**Step 3 — De-normalize to physical DOFs:**

.. math::

   \hat{\mathbf{x}} = \mathbf{N}\,\tilde{\mathbf{x}}

In compact form:

.. math::

   \hat{\mathbf{x}}
   = \mathbf{N} \cdot
     \bigl(\mathbf{C}_n^{-1/2}\,\tilde{\mathbf{A}}\bigr)_k^{+} \cdot
     \mathbf{C}_n^{-1/2}\,\mathbf{y}

For implementation details, see Equation (10) in
`arXiv:2406.04656 <https://arxiv.org/abs/2406.04656>`_ and
`Martínez-Galarza et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024ApJ...974..108M>`_.


Definition of v-modes
---------------------

The **v-modes** are defined as the right singular vectors of the normalized
sensitivity matrix.  From the thin SVD:

.. math::

   \tilde{\mathbf{A}} = \mathbf{U}\,\operatorname{diag}(\boldsymbol{\Sigma})\,\mathbf{V}_h

the rows of :math:`\mathbf{V}_h` form an orthonormal basis for the
normalized DOF space.  Row :math:`i` of :math:`\mathbf{V}_h` is v-mode
:math:`i`.

**Each v-mode is a linear combination of all (normalized) DOFs.**
Concretely, v-mode :math:`i` is the vector
:math:`\mathbf{v}_i = (\mathbf{V}_h)_{i,:} \in \mathbb{R}^{n}`, which
specifies a coordinated adjustment of hexapod positions and bending modes.

The associated singular value :math:`\sigma_i` quantifies the **optical
leverage** of v-mode :math:`i`: the amount of wavefront change (in Zernike
amplitude) produced per unit of that DOF combination.  Because the singular
values are sorted in decreasing order, v-mode 0 is the most optically
important mode and v-mode :math:`n{-}1` is the least.

**The number of v-modes equals the number of used DOFs.**  With
:math:`n = 22` used DOFs, the SVD of the :math:`(76 \times 22)` matrix
produces :math:`\mathbf{V}_h \in \mathbb{R}^{22 \times 22}`, yielding
exactly 22 v-modes.


Coordinate transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^

**DOFs to v-modes.**  Given a DOF vector :math:`\mathbf{x}`, the v-mode
representation is:

.. math::

   \mathbf{v} = \mathbf{V}_h\,\mathbf{N}^{-1}\,\mathbf{x}

**V-modes to DOFs.**  Given a v-mode vector :math:`\mathbf{v}`, the DOF
representation is:

.. math::

   \mathbf{x} = \mathbf{N}\,\mathbf{V}_h^{T}\,\mathbf{v}

These transformations are exact and lossless because
:math:`\mathbf{V}_h^{T}\mathbf{V}_h = \mathbf{I}`.

In the code, these correspond to
``StateEstimator.get_vmodes_from_dofs()`` and
``StateEstimator.get_dofs_from_vmodes()`` in
`state_estimator.py <https://github.com/lsst-ts/ts_ofc/blob/develop/python/lsst/ts/ofc/state_estimator.py>`_.


Physical interpretation of v-modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although each v-mode is a linear combination of *all* DOFs, each v-mode
produces a characteristic wavefront pattern described by the corresponding
left singular vector :math:`\mathbf{u}_i = \mathbf{U}_{:,i}`.  This vector
contains the Zernike coefficients (across all sensors) that result from
commanding v-mode :math:`i`.

If :math:`\mathbf{u}_i` has most of its energy in the Zernike coefficients
corresponding to a particular aberration — say, focus (Z4) — then v-mode
:math:`i` is effectively "the focus mode": the coordinated DOF adjustment
that most efficiently changes the telescope's focus.  Similarly, if
:math:`\mathbf{u}_j` is dominated by spherical aberration (Z11), then
v-mode :math:`j` is "the spherical mode."

This is why operators associate specific v-modes with specific aberrations:

- V-mode :math:`i` is named after the **dominant Zernike** in
  :math:`\mathbf{u}_i`.
- The association is approximate — each v-mode typically produces a
  dominant Zernike plus smaller contributions to other Zernikes.
- The SVD discovers these "optical eigenmodes" automatically from the
  sensitivity matrix; they are not prescribed by the user.


Truncation and v-mode importance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The singular values :math:`\sigma_i` determine which v-modes are
well-constrained by the wavefront measurements:

- **V-modes** :math:`0` **through** :math:`k{-}1` (large :math:`\sigma_i`):
  These DOF combinations produce large, easily measured wavefront changes.
  The state estimator can reliably reconstruct their amplitudes.

- **V-modes** :math:`k` **through** :math:`n{-}1` (small :math:`\sigma_i`):
  These combinations produce tiny wavefront changes that are buried in
  measurement noise.  The truncation at index :math:`k` discards these
  during state estimation.

With ``truncation_index = 12`` and :math:`n = 22` DOFs, the first 12
v-modes are retained and the remaining 10 are discarded during state
estimation.  This does not reduce the dimension of the correction — the PID
controller still operates on all 22 v-modes (or DOFs) — but the state
estimate for v-modes 12--21 is effectively zero because their singular
values were zeroed in the truncated sensitivity matrix.


.. _PID control:

PID control
-----------

The PID controller computes a correction :math:`\mathbf{u}_k` at each
iteration :math:`k` based on the estimated state.  The correction is
applied with **negative feedback**: the telescope moves in the opposite
direction of the estimated error.

The overall flow in ``OFC.calculate_corrections()`` is:

1. Estimate the optical state :math:`\hat{\mathbf{x}}_k` via the state estimator.
2. Pass :math:`\hat{\mathbf{x}}_k` to the PID controller, which computes an
   internal control signal.
3. **Negate** the control signal: ``uk = -controller.control_step(...)``
4. Aggregate the correction into the telescope state.

The negation in step 3 ensures that if the estimated state indicates a positive
defocus, the correction drives it back toward zero.


The setpoint: what the PID targets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The PID controller drives the estimated state toward a **setpoint vector**
:math:`\mathbf{s} \in \mathbb{R}^{50}`, which specifies the desired value for
each of the 50 possible degrees of freedom.  In the default configuration
(``pid_controller.yaml``), all setpoint entries are zero:

.. math::

   \mathbf{s} = \mathbf{0}

This means the controller targets a state where all DOF deviations are zero —
i.e., the telescope optics are in their nominal alignment.  A non-zero setpoint
would allow the controller to target a deliberate offset (e.g., intentional
defocus for engineering purposes).

**The effective target of the closed loop is not simply "zero physical DOFs."**
The PID setpoint lives in DOF space (:math:`\mathbf{x}`), but what "zero" means
depends on the measurement preprocessing (see `Problem setup`_ above).  Before
state estimation, the intrinsic Zernikes and y2 correction are subtracted in
measurement space (:math:`\mathbf{y}`).  The y2 correction encodes per-sensor
focus offsets determined from on-sky data.

The system converges when :math:`\hat{\mathbf{x}} = \mathbf{0}`, which
corresponds to the measurement residual :math:`\mathbf{y} \approx \mathbf{0}`,
which in turn means :math:`\mathbf{y}_\text{meas} \approx \mathbf{y}_\text{intrinsic} + \mathbf{y}_2`.
So the telescope is driven to the state where its measured wavefront matches the
intrinsic design prediction plus the empirical y2 offsets — i.e., the state that
produces the best image quality as determined from on-sky calibration.

Only the entries corresponding to the active DOFs (``dof_idx``) are used.


PID formula (as implemented)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At each iteration :math:`k`, the controller receives the estimated state
:math:`\hat{\mathbf{x}}_k` (either in DOF or v-mode space) and computes:

**Error signal:**

.. math::

   \mathbf{e}_k = \mathbf{s}[\text{dof\_idx}] - \hat{\mathbf{x}}_k

**Integral with anti-windup:**

.. math::

   \mathbf{I}_k &= \mathbf{I}_{k-1} + \mathbf{e}_k \\
   \mathbf{I}_k &\leftarrow \operatorname{clip}\bigl(
       \mathbf{I}_k,\;
       -\mathbf{I}_\text{max}[\text{dof\_idx}],\;
       +\mathbf{I}_\text{max}[\text{dof\_idx}]
   \bigr)

The ``max_integral`` vector (default: 1.0 for all DOFs) prevents the integral
from growing without bound during sustained errors (anti-windup protection).

**Derivative with low-pass filter:**

.. math::

   \mathbf{d}_k^\text{raw} &= \mathbf{e}_k - \mathbf{e}_{k-1} \\
   \mathbf{d}_k &= \alpha\,\mathbf{d}_k^\text{raw}
                  + (1 - \alpha)\,\mathbf{d}_{k-1}

where :math:`\alpha` is the ``derivative_filter_coeff``.  When
:math:`\alpha = 1` (the default), no filtering is applied and the derivative
is simply the difference between consecutive errors.  When
:math:`\alpha < 1`, the derivative is smoothed with an exponential moving
average to reduce sensitivity to noise.

**Control output:**

.. math::

   \mathbf{u}_k^\text{pid} = \mathbf{K}_p\,\mathbf{e}_k
                            + \mathbf{K}_i\,\mathbf{I}_k
                            + \mathbf{K}_d\,\mathbf{d}_k

where:

- :math:`\mathbf{K}_p = \operatorname{diag}\!\bigl(k_{p,1}, \ldots, k_{p,n}\bigr)`
  is the proportional gain (default: 0.7 for all DOFs),
- :math:`\mathbf{K}_i = \operatorname{diag}\!\bigl(k_{i,1}, \ldots, k_{i,n}\bigr)`
  is the integral gain (default: 0.0 — integral action disabled),
- :math:`\mathbf{K}_d = \operatorname{diag}\!\bigl(k_{d,1}, \ldots, k_{d,n}\bigr)`
  is the derivative gain (default: 0.1).

The gains can be a scalar (applied uniformly to all DOFs) or an array (per-DOF
tuning).  They are indexed by ``dof_idx`` to select only the active DOFs.

The final correction applied to the telescope is:

.. math::

   \mathbf{u}_k = -\mathbf{u}_k^\text{pid}

This negation, applied in ``OFC.calculate_corrections()``, converts the
PID output (which points from state toward setpoint) into a physical correction
(which drives the optics back toward the setpoint).


Default configuration
^^^^^^^^^^^^^^^^^^^^^

The default PID parameters from ``pid_controller.yaml``:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Description
   * - ``kp``
     - 0.7
     - Proportional gain (scalar or per-DOF array)
   * - ``ki``
     - 0.0
     - Integral gain (disabled by default)
   * - ``kd``
     - 0.1
     - Derivative gain
   * - ``derivative_filter_coeff``
     - 1.0
     - Low-pass filter coefficient (1.0 = no filtering)
   * - ``setpoint``
     - all zeros (50 elements)
     - Target state for each DOF
   * - ``max_integral``
     - 1.0 (all DOFs)
     - Anti-windup clipping limit

**Physical meaning of gains:** :math:`k_p = 0.7` means "correct 70% of the
estimated error per iteration."  :math:`k_d = 0.1` provides damping that
reduces overshoot by opposing rapid changes in the error.  With
:math:`k_i = 0`, there is no integral action by default — steady-state errors
are not automatically eliminated across iterations.


Controller state and history
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The controller maintains internal state between iterations:

- **aggregated_state** (:math:`\mathbf{x}_\text{agg}`): the cumulative
  correction applied to the telescope (initialized from ``dof_state0``).
- **integral** (:math:`\mathbf{I}`): accumulated integral of the error.
- **previous_error** (:math:`\mathbf{e}_{k-1}`): error from the previous step.
- **filtered_derivative** (:math:`\mathbf{d}_{k-1}`): filtered derivative from
  the previous step.

The ``reset_history()`` method:

1. Saves the current aggregated state as the new reference:
   :math:`\mathbf{x}_0 \leftarrow \mathbf{x}_\text{agg}`.
2. Zeros the integral, previous error, and filtered derivative.

This is called when the closed loop is restarted with a new ``comp_dof_idx``
configuration, effectively starting the PID fresh from the current telescope
state.


PID in DOF space
^^^^^^^^^^^^^^^^

When ``control_vmodes = False``, the state estimator returns the optical state
directly in DOF space:

.. math::

   \hat{\mathbf{x}}_k = \mathbf{N} \cdot
     \bigl(\mathbf{C}_n^{-1/2}\,\tilde{\mathbf{A}}\bigr)_k^{+} \cdot
     \mathbf{C}_n^{-1/2}\,\mathbf{y}

The PID operates on this 22-dimensional (or however many DOFs are active)
state vector.  Each gain entry corresponds directly to a physical DOF.

**Physical meaning of DOF gains:** each diagonal entry of
:math:`\mathbf{K}_p` controls the aggressiveness of correction for one
physical DOF.  For example, :math:`k_{p,0} = 0.7` means "correct 70% of
the estimated M2 hexapod dZ error per iteration."


PID in v-mode space
^^^^^^^^^^^^^^^^^^^

When ``control_vmodes = True``, the state estimator returns the optical
state directly in v-mode space.  Internally, the state estimator solves
the least-squares problem and returns the v-mode coefficients (zero-padded
to the full DOF dimension):

.. math::

   \hat{\mathbf{v}}_k = \text{pad}_{n}\bigl(\mathbf{v}_\text{state}\bigr)

where :math:`\mathbf{v}_\text{state}` has dimension equal to the number of
retained modes (truncation index), padded with zeros for the discarded modes.

The PID then operates on this v-mode vector using the **same** ``calculate_pid_step()``
method — the controller does not know whether it is operating on DOFs or v-modes.
The error, integral, and derivative are all computed in v-mode space:

.. math::

   \mathbf{e}^{(v)}_k &= \mathbf{s}[\text{dof\_idx}] - \hat{\mathbf{v}}_k \\
   \mathbf{u}^{(v)}_k &= \mathbf{K}_p[\text{dof\_idx}]\,\mathbf{e}^{(v)}_k
                        + \mathbf{K}_i[\text{dof\_idx}]\,\mathbf{I}^{(v)}_k
                        + \mathbf{K}_d[\text{dof\_idx}]\,\mathbf{d}^{(v)}_k

After the PID produces :math:`\mathbf{u}^{(v)}_k`, the OFC converts it back
to physical DOF corrections:

.. math::

   \mathbf{u}_k = -\mathbf{N}\,\mathbf{V}_{h,k}^{T}\,\mathbf{u}^{(v)}_k[\text{:truncate\_index}]

This uses ``get_dofs_from_vmodes()``, which takes the first ``truncate_index``
elements of the v-mode correction, multiplies by :math:`\mathbf{V}_{h,k}^T` to
project to normalized DOF space, then de-normalizes with :math:`\mathbf{N}`.

**Important implementation detail:** because the PID uses the same gain
arrays indexed by ``dof_idx`` regardless of whether it operates on DOFs or
v-modes, the gain at position :math:`i` applies to v-mode :math:`i` (ordered
by singular value).  This means v-mode 0 (highest optical leverage) gets the
gain at ``dof_idx[0]``, v-mode 1 gets the gain at ``dof_idx[1]``, and so on.

**Physical meaning of v-mode gains:** each entry of
:math:`\mathbf{K}_p` controls how aggressively a specific optical
mode is corrected.  For example, assigning a higher gain to v-mode 0
(focus) means the system corrects focus errors more aggressively than
higher-order aberrations.

**The advantage of v-mode control** is the ability to assign different gains
based on optical importance:

- **High gain** on well-constrained v-modes (large :math:`\sigma_i`,
  reliable state estimates) for fast convergence.
- **Low gain** on poorly-constrained v-modes (small :math:`\sigma_i`,
  noisy estimates) to avoid injecting noise into the corrections.


Equivalent DOF-space gain
^^^^^^^^^^^^^^^^^^^^^^^^^

Composing the three steps, the effective DOF-space gain when controlling in
v-mode space is:

.. math::

   \mathbf{K}_\text{eff}
   = \mathbf{N}\,\mathbf{V}_h^{T}\,\mathbf{K}^{(v)}_p\,
     \mathbf{V}_h\,\mathbf{N}^{-1}

When the v-mode gain is a **uniform scalar** (:math:`\mathbf{K}^{(v)}_p = k_p \mathbf{I}`):

.. math::

   \mathbf{K}_\text{eff}
   = k_p\,\mathbf{N}\,\mathbf{V}_h^{T}\,\mathbf{V}_h\,\mathbf{N}^{-1}
   = k_p\,\mathbf{I}

The v-mode transformation cancels entirely; v-mode control with uniform
gains is identical to DOF-space control with the same scalar gain.

When the v-mode gain has **different values per mode**, the matrix
:math:`\mathbf{V}_h^{T}\,\mathbf{K}^{(v)}_p\,\mathbf{V}_h` is a full
(non-diagonal) matrix, creating **cross-coupling between DOFs**.  This is
the desired behavior: correcting an optical eigenmode (e.g., focus) requires
coordinated adjustments to multiple physical DOFs.

For example, with the two-level gain structure
:math:`k_p^{(v)} = 0.18` for v-modes :math:`0`--:math:`10` and
:math:`k_p^{(v)} = 0.045` for v-modes :math:`11`--:math:`21`, the gain
decomposes as:

.. math::

   \mathbf{K}^{(v)}_p
   = 0.045\,\mathbf{I}
   + 0.135\,\mathbf{P}_{11}

where :math:`\mathbf{P}_{11}` is the projection onto the subspace of the
first 11 v-modes.  This applies a uniform base gain plus a boost on the
most optically important modes.


Correction aggregation
----------------------

After each PID step, the DOF correction :math:`\mathbf{u}_k` is
accumulated into the telescope state:

.. math::

   \mathbf{x}_{k+1} = \mathbf{x}_k + \mathbf{u}_k

This aggregated state :math:`\mathbf{x}_{k+1}` is then decomposed into
physical corrections for the individual components:

- **M2 Hexapod:** DOFs 0--4 (dZ, dX, dY, rX, rY), transformed via a
  rotation matrix to hexapod command coordinates.
- **Camera Hexapod:** DOFs 5--9, same transformation.
- **M1M3 Bending Modes:** DOFs 10--29, converted from bending mode
  coefficients to actuator forces via the bending-mode-to-force matrix.
- **M2 Bending Modes:** DOFs 30--49, same conversion.

These corrections are then issued to the respective subsystem CSCs by
``ts_mtaos``.


Summary
-------

The following table summarizes the key properties of DOF-space and v-mode
control:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Property
     - DOF space
     - V-mode space
   * - State vector
     - Physical DOFs (hexapod positions, bending modes)
     - Optical eigenmodes ordered by singular value
   * - Gain meaning
     - Correction fraction per physical DOF
     - Correction fraction per optical mode
   * - Cross-coupling
     - None (diagonal gain matrix)
     - Coordinated multi-DOF corrections
   * - Advantage
     - Simple, physically intuitive
     - Gain tuning by optical importance
   * - Number of modes
     - :math:`n` = number of used DOFs
     - :math:`n` (same; from SVD of :math:`m \times n` matrix)
   * - Truncation
     - Applied during state estimation only
     - Same; PID still operates on all :math:`n` modes
