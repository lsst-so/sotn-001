V-Modes and the OFC Control Loop
=================================

This chapter provides a mathematical description of the Optical Feedback
Control (OFC) pipeline implemented in
`ts_ofc <https://github.com/lsst-ts/ts_ofc>`_, and of the way the
`ts_mtaos <https://github.com/lsst-ts/ts_mtaos>`_ closed loop wraps it.
It covers the forward optical model, the preprocessing of the wavefront
measurements and the state estimation inverse problem.  It then defines
the **v-modes**, explains their physical meaning, and shows how PID gains
are applied in both DOF space and v-mode space.  It closes with the
decisions that MTAOS makes around each OFC call.


Scope and conventions
---------------------

**What is covered.**  The core of this note is the mathematics implemented
in ``ts_ofc``: how the wavefront error vector :math:`\mathbf{y}` is turned
into an estimate of the telescope state :math:`\hat{\mathbf{x}}`, and how a
PID controller converts that estimate into a correction.  A dedicated
section, `Where OFC sits in the closed loop`_, describes only the parts of
``ts_mtaos`` that modify the control law or its inputs (gain overrides,
image skipping, the large-defocus bypass, stress limits and the pointing
correction).  The CSC state machine, event publishing and pipeline
management are documented in the
`MTAOS documentation <https://ts-mtaos.lsst.io>`_ and are not repeated
here.

**What is not covered.**  The measurement side of the loop is out of scope:
how `ts_wep <https://github.com/lsst-ts/ts_wep>`_ estimates Zernike
coefficients from donut images, and how those results are visualized with
``donut_viz``.  Here the wavefront error is treated as an input.  At the
OFC boundary it is:

- one vector of Zernike coefficients per corner wavefront sensor, in
  micrometers, for Noll indices Z4 up to Z28 (``znmin = 4``,
  ``znmax = 28``), read from the ``zernikes`` butler table produced by WEP;
- the average of the intra- and extra-focal estimates of each corner;
- expressed in the camera frame (not de-rotated);
- with sensors whose estimate failed (``NaN``) dropped before inversion.

In the deployed configuration, MTAOS reads the ``deviation_columns`` of
that table, that is, the measured Zernikes **minus the intrinsic (design)
Zernikes**.  It therefore sets ``subtract_intrinsics = False`` so that OFC
does not subtract them a second time.  This affects the subtraction
described in `Measurement preprocessing`_.

**Which configuration the numbers come from.**  ``ts_ofc`` ships a set of
policy files, but the values actually used on the summit live in
`ts_config_mttcs <https://github.com/lsst-ts/ts_config_mttcs>`_:

- ``MTAOS/ofc/`` holds the OFC data (sensitivity matrix, intrinsic
  Zernikes, normalization weights, ``y2`` offsets, noise covariance, and
  the controller files).  ``OFCData`` loads ``configurations/init.yaml``
  by default, and MTAOS points ``OFCData`` to this directory.
- ``MTAOS/v13/_init.yaml`` (at the time of writing) holds the CSC
  configuration: which DOFs are used, thresholds, stress limits, and the
  closed-loop parameters.

Unless stated otherwise, every default quoted in this note comes from
those two files.  The copies packaged with ``ts_ofc`` differ, for example
in the gains, and should not be used as a reference.

**Notation.**  Vectors are bold lower case, matrices bold upper case.
:math:`k` is the closed-loop iteration index.  :math:`r` is the truncation
index (the number of v-modes retained).  :math:`\mathcal{D}` is the set of
active DOF indices (``dof_idx``); a subscript :math:`\mathcal{D}` on a
50-element vector or on a gain selects those entries.


The forward optical model
-------------------------

The **sensitivity matrix** relates the telescope's degrees of freedom
(DOFs) to the wavefront error they produce across the focal plane.  For
each sensor position, it encodes how a unit change in each DOF affects
each Zernike coefficient of the wavefront.

If we stack the Zernike coefficients from all sensors into a single vector,
the forward model is:

.. math::

   \mathbf{y} = \mathbf{A}\,\mathbf{x}

where:

- :math:`\mathbf{y} \in \mathbb{R}^{m}` is the wavefront error vector
  (Zernike coefficients across all sensors, in micrometers),
  with :math:`m = n_\mathrm{sensors} \times n_\mathrm{Zernikes}`.
  For four corner wavefront sensors and Zernikes Z4--Z22 (19 coefficients),
  :math:`m = 4 \times 19 = 76`.

- :math:`\mathbf{A} \in \mathbb{R}^{m \times n}` is the sensitivity matrix,
  evaluated at the sensor field positions and the camera rotation angle.

- :math:`\mathbf{x} \in \mathbb{R}^{n}` is the DOF state vector, with
  :math:`n = |\mathcal{D}|` the number of used DOFs.

The full DOF vector has 50 entries: M2 hexapod (dZ, dX, dY in
micrometers; rX, rY in arcseconds), camera hexapod (same), 20 M1M3 bending
modes and 20 M2 bending modes (mode amplitudes in micrometers).  The
deployed configuration (``used_dofs`` in ``MTAOS/v13/_init.yaml``) uses
:math:`n = 20`: the 10 hexapod DOFs, M1M3 bending modes 1--5 (DOF indices
10--14) and M2 bending modes 1--5 (DOF indices 30--34).  Earlier
configurations used 22 DOFs (7 M1M3 modes).  All dimensions in this note
assume the deployed 20.

.. _From the double Zernike table to A:

From the double Zernike table to the matrix A
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sensitivity matrix is **not stored** as the :math:`m \times n` matrix
:math:`\mathbf{A}` of the forward model.  It is stored as a
**double Zernike** coefficient array of shape :math:`31 \times 29 \times 50`
(``sensitivity_matrix/lsst_sensitivity_dz_31_29_50.yaml``):

- the first axis runs over 31 **field** Zernike terms (Noll indices 0--30),
  describing how an aberration varies across the field of view;
- the second axis runs over 29 **pupil** Zernike terms (Noll indices
  0--28), the aberration itself;
- the third axis runs over the 50 DOFs.

Index 0 of each axis is unused by convention, so the pupil terms are
Z1--Z28, which is why ``znmax = 28``.  For each DOF the
:math:`31 \times 29` slice is a double Zernike: a function of field
position :math:`(u, v)` that returns the vector of pupil Zernike
coefficients produced by a unit change of that DOF at that field
position.  This single array describes the optical response over the
whole field of view, independently of where the sensors are.

.. note::

   The header comment of the YAML file labels the first two axes the other
   way round.  The generator script ``generate_sensitivity_matrix.py`` and
   the galsim ``DoubleZernike`` convention used to evaluate the table both
   put the field index first, which is the convention used here.

The code turns this 3-D table into a 2-D matrix in **two different
ways**, for two different purposes.

**Path 1: evaluate at the sensors (the A that is inverted).**  This is
done for every exposure by ``SensitivityMatrix.evaluate`` and
``StateEstimator.get_sensitivity_matrix``:

1. Evaluate the double Zernike of each DOF at the field position of each
   sensor (``sample_points``, rotated by minus the rotator angle, see
   `Measurement preprocessing`_).  The 31-term field axis disappears and
   is replaced by one pupil Zernike vector per sensor:
   :math:`(4, 29, 50)`.
2. Keep pupil terms Z4--Z28, then select ``zn_selected``:
   :math:`(4, 19, 50)`.  The deployed controller file selects Z4--Z22; in
   the closed loop, MTAOS overrides this with the indices that WEP
   actually produced.
3. Flatten the sensor and Zernike axes into a single row index:
   :math:`(4 \times 19, 50) = (76, 50)`.  Rows are sensor-major, so rows
   0--18 are Z4--Z22 of the first sensor, rows 19--37 those of the second
   sensor, and so on.  The wavefront vector :math:`\mathbf{y}` and the
   noise covariance :math:`\mathbf{C}_n` are ordered the same way.
4. Select the columns of the active DOFs :math:`\mathcal{D}`:
   :math:`(76, 20)`.

The result is the matrix :math:`\mathbf{A} \in \mathbb{R}^{m \times n}`
of the forward model.  This is the matrix that is normalized, projected
onto the v-modes, whitened and inverted in
`State estimation: the inverse problem`_.  Because it is built from the
list of sensors that reported, it shrinks automatically when a corner is
missing, for example to :math:`57 \times 20`.

**Path 2: flatten the table (the matrix that defines the v-modes).**  This
is done once per configuration by ``StateEstimator._update_from_ofc_data``
and involves no evaluation at all.  The table is reshaped so that every
(field term, pupil term) coefficient pair becomes a row:

.. math::

   (31, 29, 50) \;\rightarrow\; (31 \times 29,\, 50) = (899, 50)
   \;\rightarrow\; (899, 20)

where the last step selects the active DOFs.  The field axis is kept as
coefficients rather than collapsed onto sensor positions, so this matrix
describes the response over the whole field.  Its SVD defines the v-modes
(see `Definition of v-modes`_), which is why they do not depend on the
rotator angle, on how many sensors reported, or on the Zernike subset.
Rows corresponding to index 0 of either axis, or to the piston, tip and
tilt pupil terms, are zero; they are harmless in the SVD.

In short, the per-exposure :math:`76 \times 20` matrix answers "what do
*these sensors* see?", and the :math:`899 \times 20` matrix answers "what
does the *whole field* see?".  The first is inverted.  The second only
chooses the basis in which the inversion is done.


DOF normalization
-----------------

The physical DOFs span very different scales.  Hexapod translations are
in micrometers and hexapod rotations in arcseconds.  Bending mode
amplitudes are also in micrometers, but their optical leverage varies
widely.  To give each DOF comparable influence in the inversion, we
introduce a diagonal normalization matrix:

.. math::

   \mathbf{N} = \operatorname{diag}\!\bigl(w_1, w_2, \ldots, w_n\bigr)

where the weights :math:`w_i` are read from the ``normalization_weights``
file named in the controller configuration.  The deployed ``init.yaml``
uses ``range0.5_fwhm-0.15.yaml``.  Each weight combines a range-based
factor and an FWHM-based factor, :math:`w_i = r_i^{\alpha} f_i^{\beta}`,
with the exponents encoded in the file name.

This normalization follows
`Megias Homar et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024ApJ...974..108M>`_
(ApJ 974, 108; also available as
`arXiv:2406.04656 <https://arxiv.org/abs/2406.04656>`_), Eqs. 9 and 11.
The paper calls the matrix :math:`\boldsymbol{\Omega}` and defines each
diagonal entry as the operational range of the DOF times its effect on
the FWHM.

The normalized sensitivity matrix and normalized state are:

.. math::

   \tilde{\mathbf{A}} = \mathbf{A}\,\mathbf{N}, \qquad
   \tilde{\mathbf{x}} = \mathbf{N}^{-1}\,\mathbf{x}

so that :math:`\mathbf{y} = \mathbf{A}\,\mathbf{x}
= \tilde{\mathbf{A}}\,\tilde{\mathbf{x}}`.  All estimation is done in the
normalized space, and the result is de-normalized with
:math:`\mathbf{x} = \mathbf{N}\,\tilde{\mathbf{x}}` at the end.


.. _Measurement preprocessing:

Measurement preprocessing
-------------------------

Before the inversion, the state estimator (``StateEstimator.dof_state``)
brings the measurement and the model into the same frame and removes
known contributions from the measurement.

Camera rotation
^^^^^^^^^^^^^^^

WEP reports the Zernikes in the camera frame, while the double Zernike
model is defined in the telescope (field) frame.  Let
:math:`\theta = \theta_\mathrm{rot} + \theta_\mathrm{off}`, where
:math:`\theta_\mathrm{rot}` is the camera rotator angle at the time of the
exposure (passed by MTAOS as ``rotation_angle``) and
:math:`\theta_\mathrm{off}` is the ``rotation_offset`` of the controller
file.  The offset is 0 degrees in the deployed ``init.yaml``; the packaged
``pid_controller.yaml`` uses 90 degrees.  Two things are rotated:

1. **The measurement.**  Each sensor's Zernike vector is rotated by
   :math:`\theta` as a pupil Zernike series (``galsim.zernike.Zernike.rotate``).
   Rotating a Zernike series mixes the :math:`\pm m` pairs
   (e.g. Z5/Z6, Z7/Z8) and leaves the rotationally symmetric terms
   (Z4, Z11, Z22) unchanged.

2. **The model.**  The sensor field angles are rotated by :math:`-\theta`
   before the double Zernike is evaluated, so that :math:`\mathbf{A}`
   describes where the sensors actually were on the sky.  The intrinsic
   Zernikes are evaluated at the same rotated positions.

Subtraction of known terms
^^^^^^^^^^^^^^^^^^^^^^^^^^

The residual wavefront error used as input to the inversion is:

.. math::

   \mathbf{y} = \mathbf{y}_\mathrm{meas} - \mathbf{y}_\mathrm{intrinsic}
                - \mathbf{y}_2

These subtractions operate in **measurement space** (Zernike coefficients at
each sensor, :math:`\mathbf{y} \in \mathbb{R}^m`):

- :math:`\mathbf{y}_\mathrm{intrinsic}` contains the design (intrinsic)
  Zernike coefficients for the current filter and sensor positions.  They
  represent the wavefront error that the optical design produces even when
  the optics are perfectly aligned (:math:`\mathbf{x} = \mathbf{0}`).
  They are subtracted only when ``subtract_intrinsics = True``.  **In the
  deployed closed loop this flag is False**, because the
  ``deviation_columns`` read from the WEP table already have the
  intrinsics removed (see `Scope and conventions`_).  Setting the flag to
  ``True`` in that case would subtract them twice.

- :math:`\mathbf{y}_2` is a **static per-sensor correction** in measurement
  space, read from the instrument's ``y2`` file (``lsst_y2.yaml``).  In
  the deployed file, the corner sensors carry a focus (Z4) offset of
  -0.11 µm (R00, R04, R40) or -0.19 µm (R44).  All four corners also carry
  -0.07 µm of spherical aberration (Z11) and +0.13 µm of Z14.  The science
  sensors are zero.  These values were determined from on-sky data as the
  residuals at which the best image quality is achieved.

These subtractions ensure that the residual :math:`\mathbf{y}` is close
to zero when the telescope is at the desired alignment state
(:math:`\mathbf{x} = \mathbf{0}` in DOF space).  The state estimator then
inverts :math:`\mathbf{y}` to produce the DOF estimate
:math:`\hat{\mathbf{x}}`, and the PID controller drives that estimate to a
setpoint of zero (see `PID control`_).  The full chain is therefore:

1. Rotate the measurement and the model into a common frame.
2. Subtract known contributions in y-space:
   :math:`\mathbf{y} = \mathbf{y}_\mathrm{meas} - \mathbf{y}_\mathrm{intrinsic} - \mathbf{y}_2`
   (or without the intrinsic term when WEP already removed it).
3. Estimate the DOF state: :math:`\hat{\mathbf{x}} = f(\mathbf{y})`.
4. Let the PID drive :math:`\hat{\mathbf{x}} \to \mathbf{0}`, the setpoint in x-space.

The y-space subtractions and the x-space setpoint work together.  The
subtractions define what "zero misalignment" looks like in measurement
space, and the PID targets zero in DOF space.  If the :math:`\mathbf{y}_2`
values were wrong, for example if the focus offsets did not correspond to
the best image quality, the system would converge to a suboptimal state
even though the PID error would be zero.


State estimation: the inverse problem
--------------------------------------

Given the residual wavefront error :math:`\mathbf{y}`, the goal is to
estimate the DOF state :math:`\hat{\mathbf{x}}`.  Including measurement
noise, the forward model in normalized coordinates reads:

.. math::

   \mathbf{y} = \tilde{\mathbf{A}}\,\tilde{\mathbf{x}} + \mathbf{n}

where :math:`\tilde{\mathbf{A}}` is the normalized sensitivity matrix,
:math:`\tilde{\mathbf{x}}` the normalized DOF state and :math:`\mathbf{n}`
measurement noise with covariance
:math:`\mathbf{C}_n = \mathbb{E}[\mathbf{n}\mathbf{n}^T]`.  This is an
overdetermined linear inverse problem (:math:`m > n`), solved by
regularized least squares.  Two regularization ideas are combined:
**truncation**, which restricts the solution to the well-constrained
directions of DOF space, and **noise weighting**, which down-weights noisy
measurements.  Both are introduced below in their textbook form, followed
by what the code actually does.

Noise covariance
^^^^^^^^^^^^^^^^

In the definition of :math:`\mathbf{C}_n` above,
:math:`\mathbb{E}[\cdot]` denotes the **expectation operator** (the
statistical mean), and :math:`\mathbf{n}\mathbf{n}^T` is the outer
product of the noise vector with itself.  The covariance matrix :math:`\mathbf{C}_n` therefore
contains:

- **Diagonal entries** :math:`(\mathbf{C}_n)_{ii} = \mathbb{E}[n_i^2] = \sigma_i^2`:
  the variance of the :math:`i`-th Zernike measurement.
- **Off-diagonal entries** :math:`(\mathbf{C}_n)_{ij} = \mathbb{E}[n_i n_j]`:
  the covariance between Zernikes :math:`i` and :math:`j`.

In practice, :math:`\mathbf{C}_n` can be estimated from stability data by
computing the **sample covariance** of the residual Zernikes across
:math:`N` exposures:

.. math::

   \hat{\mathbf{C}}_n = \frac{1}{N-1}
   \sum_{k=1}^{N} (\mathbf{y}_k - \bar{\mathbf{y}})
                     (\mathbf{y}_k - \bar{\mathbf{y}})^T

where :math:`\mathbf{y}_k` are the measured Zernike vectors and
:math:`\bar{\mathbf{y}}` is their mean.  See
`SITCOMTN-177 <https://sitcomtn-177.lsst.io>`_ for how this covariance is
estimated from on-sky data.

.. note::

   The deployed configuration does **not** name a
   ``noise_covariance_filename``, and the only file shipped in
   ``ts_config_mttcs`` (``default_noise_covariance.yaml``) is an identity
   matrix.  The noise-weighting path is therefore implemented but has no
   effect on the summit today.  Running the loop with a non-identity
   covariance has not been validated yet.  It will be validated with the
   OFC/PID simulator, currently under development, before it is enabled.

Approach 1: SVD with truncation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The standard least-squares solution uses the pseudo-inverse.  Start by
computing the singular value decomposition (SVD):

.. math::

   \tilde{\mathbf{A}} = \mathbf{U}\,\operatorname{diag}(\boldsymbol{\Sigma})\,\mathbf{V}_h

In the **thin (economy) SVD** of an :math:`m \times n` matrix with :math:`m > n`:

- :math:`\mathbf{U} \in \mathbb{R}^{m \times n}` holds the left singular
  vectors as columns,
- :math:`\boldsymbol{\Sigma} \in \mathbb{R}^{n}` is the vector of singular
  values :math:`(\sigma_1, \ldots, \sigma_n)`, ordered
  :math:`\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n \geq 0`,
- :math:`\mathbf{V}_h \in \mathbb{R}^{n \times n}` is an orthogonal matrix
  whose rows are the right singular vectors.

The pseudo-inverse is:

.. math::

   \tilde{\mathbf{A}}^{+} = \mathbf{V}_h^T\,\operatorname{diag}(1/\boldsymbol{\Sigma})\,\mathbf{U}^T

**Truncation regularization.**  Small singular values amplify noise because
:math:`1/\sigma_i \to \infty` as :math:`\sigma_i \to 0`.  To regularize,
keep only the :math:`r` largest singular values:

.. math::

   \tilde{\mathbf{A}}_r^{+}
   = \mathbf{V}_{h,r}^T\,\operatorname{diag}(1/\boldsymbol{\Sigma}_r)\,\mathbf{U}_r^T

where :math:`\mathbf{U}_r` contains the first :math:`r` columns of
:math:`\mathbf{U}`, :math:`\mathbf{V}_{h,r}` the first :math:`r` rows of
:math:`\mathbf{V}_h`, and :math:`\boldsymbol{\Sigma}_r` the first :math:`r`
singular values.  The truncated solution
:math:`\tilde{\mathbf{x}}_\mathrm{trunc} = \tilde{\mathbf{A}}_r^{+}\,\mathbf{y}`
discards the :math:`n - r` least-constrained DOF combinations and lies
entirely in the subspace spanned by the first :math:`r` right singular
vectors.

Approach 2: Noise-weighted least squares
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the measurement noise covariance :math:`\mathbf{C}_n` is known and
non-uniform, we minimize the weighted least-squares objective:

.. math::

   \min_{\tilde{\mathbf{x}}} \;
   \bigl(\mathbf{y} - \tilde{\mathbf{A}}\,\tilde{\mathbf{x}}\bigr)^T
   \mathbf{C}_n^{-1}
   \bigl(\mathbf{y} - \tilde{\mathbf{A}}\,\tilde{\mathbf{x}}\bigr)

Equivalently, **whiten** the problem by defining:

.. math::

   \bar{\mathbf{y}} = \mathbf{C}_n^{-1/2}\,\mathbf{y}, \qquad
   \bar{\mathbf{A}} = \mathbf{C}_n^{-1/2}\,\tilde{\mathbf{A}}

and solve the ordinary least-squares problem
:math:`\tilde{\mathbf{x}}_\mathrm{wls} = \bar{\mathbf{A}}^{+}\,\bar{\mathbf{y}}`.
Whitening down-weights the Zernike terms and sensors with the largest
measurement uncertainty.

What ts_ofc actually computes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The implementation combines the two ideas, but **not** by taking the
truncated SVD of the per-exposure matrix.  Instead, the truncation
subspace is fixed once from the full-field optical model, and the
least-squares fit is done inside that subspace using the per-exposure
matrix.  The steps in ``StateEstimator`` are as follows.

**Step 0 -- Define the v-modes (once per configuration).**  Take the
normalized double Zernike sensitivity array, reshaped to one row per
(field term, pupil term) coefficient pair (Path 2 in
`From the double Zernike table to A`_), and restricted to the active
DOFs: :math:`\tilde{\mathbf{A}}_\mathrm{full} \in \mathbb{R}^{899 \times n}`.
Its thin SVD,

.. math::

   \tilde{\mathbf{A}}_\mathrm{full}
   = \mathbf{U}\,\operatorname{diag}(\boldsymbol{\Sigma})\,\mathbf{V}_h ,

defines the v-modes as the rows of :math:`\mathbf{V}_h`.  This happens in
``StateEstimator._update_from_ofc_data`` and does not depend on which
sensors reported, on the rotation angle, or on ``zn_selected``.

**Step 1 -- Project the per-exposure matrix onto the first r v-modes.**
With :math:`\tilde{\mathbf{A}} \in \mathbb{R}^{m \times n}` the normalized
sensitivity matrix evaluated at the (rotated) sensor positions (Path 1 in
`From the double Zernike table to A`_):

.. math::

   \tilde{\mathbf{A}}_v = \tilde{\mathbf{A}}\,\mathbf{V}_{h,r}^T
   \;\in\; \mathbb{R}^{m \times r}

Column :math:`i` of :math:`\tilde{\mathbf{A}}_v` is the wavefront pattern
seen by the sensors when v-mode :math:`i` is commanded with unit amplitude.

**Step 2 -- Whiten by the noise covariance.**

.. math::

   \bar{\mathbf{y}} = \mathbf{C}_n^{-1/2}\,\mathbf{y}, \qquad
   \bar{\mathbf{A}}_v = \mathbf{C}_n^{-1/2}\,\tilde{\mathbf{A}}_v

The inverse square root is computed with
``scipy.linalg.fractional_matrix_power``.  In the deployed configuration
:math:`\mathbf{C}_n` is the identity, so this step changes nothing.

**Step 3 -- Ordinary least squares in v-mode space.**

.. math::

   \hat{\mathbf{v}} = \arg\min_{\mathbf{v} \in \mathbb{R}^{r}}
   \bigl\lVert \bar{\mathbf{A}}_v\,\mathbf{v} - \bar{\mathbf{y}} \bigr\rVert^2
   = \bar{\mathbf{A}}_v^{+}\,\bar{\mathbf{y}}

This is solved with ``numpy.linalg.lstsq``, without any further
truncation.  The result has :math:`r` components, one per retained
v-mode.

**Step 4 -- Return the state in the requested basis.**

- DOF basis (``control_vmodes = False``, the deployed setting):

  .. math::

     \hat{\mathbf{x}} = \mathbf{N}\,\mathbf{V}_{h,r}^T\,\hat{\mathbf{v}}
     \;\in\; \mathbb{R}^{n}

  The :math:`r` v-mode amplitudes are expanded into the :math:`n` DOFs and
  then de-normalized.  The estimate lies entirely in the subspace spanned
  by the first :math:`r` v-modes.

- V-mode basis (``control_vmodes = True``): :math:`\hat{\mathbf{v}}` is
  zero-padded from :math:`r` to :math:`n` elements and returned as is.

In compact form, for the DOF basis:

.. math::

   \hat{\mathbf{x}}
   = \mathbf{N}\,\mathbf{V}_{h,r}^T\,
     \bigl(\mathbf{C}_n^{-1/2}\,\tilde{\mathbf{A}}\,\mathbf{V}_{h,r}^T\bigr)^{+}\,
     \mathbf{C}_n^{-1/2}\,\mathbf{y}

**Relation to the textbook forms.**  If the v-modes had been computed from
the per-exposure :math:`\tilde{\mathbf{A}}` itself and
:math:`\mathbf{C}_n = \mathbf{I}`, the expression above would reduce
exactly to the truncated pseudo-inverse of Approach 1.  Because the
v-modes come from the full-field model instead, the implementation is
best read as *fitting the measurements with the best linear combination of
the first r full-field optical eigenmodes*.  A practical advantage is that
the meaning of "v-mode :math:`i`" does not change from one exposure to the
next, which is what makes per-mode gains well defined (see
`PID in v-mode space`_).  Note also that the projection is applied before
the whitening.  With a non-identity :math:`\mathbf{C}_n`, the result is
therefore not the same as whitening first and truncating afterwards.

For background, see
`Megias Homar et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024ApJ...974..108M>`_.
Eq. 2 gives the SVD pseudo-inverse, Eqs. 4--6 the noise amplification
that motivates truncation, Eqs. 7--8 the truncation criterion and the
truncated effective response, and Eqs. 9--10 the rescaled (normalized)
problem and its de-normalized estimate.

Concrete example
^^^^^^^^^^^^^^^^

Consider the deployed configuration, with :math:`m = 76` measurements
(4 corners × 19 Zernikes), :math:`n = 20` DOFs and ``truncation_index = 12``
(:math:`r = 12`):

- :math:`\mathbf{V}_h \in \mathbb{R}^{20 \times 20}` from the SVD of the
  :math:`899 \times 20` full-field matrix; :math:`\mathbf{V}_{h,12} \in \mathbb{R}^{12 \times 20}`.
- :math:`\tilde{\mathbf{A}} \in \mathbb{R}^{76 \times 20}`,
  :math:`\tilde{\mathbf{A}}_v = \tilde{\mathbf{A}}\,\mathbf{V}_{h,12}^T \in \mathbb{R}^{76 \times 12}`.
- :math:`\hat{\mathbf{v}} = \bar{\mathbf{A}}_v^{+}\,\bar{\mathbf{y}} \in \mathbb{R}^{12}`.
- :math:`\hat{\mathbf{x}} = \mathbf{N}\,\mathbf{V}_{h,12}^T\,\hat{\mathbf{v}} \in \mathbb{R}^{20}`.

With dimensions written out:

.. math::

   \underbrace{\hat{\mathbf{x}}}_{20 \times 1}
   = \underbrace{\mathbf{N}}_{20 \times 20}\,
     \underbrace{\mathbf{V}_{h,12}^T}_{20 \times 12}\,
     \underbrace{\bar{\mathbf{A}}_v^{+}}_{12 \times 76}\,
     \underbrace{\bar{\mathbf{y}}}_{76 \times 1}

The output has 20 components, one per DOF, but the data constrain only
12 independent directions.  The amplitudes of v-modes 12--19 are
identically zero in the estimate.


.. _Definition of v-modes:

Definition of v-modes
---------------------

The **v-modes** are the right singular vectors of the normalized
**full-field** sensitivity matrix
:math:`\tilde{\mathbf{A}}_\mathrm{full} \in \mathbb{R}^{899 \times n}`
(Step 0 above).  From its thin SVD,

.. math::

   \tilde{\mathbf{A}}_\mathrm{full}
   = \mathbf{U}\,\operatorname{diag}(\boldsymbol{\Sigma})\,\mathbf{V}_h ,

the rows of :math:`\mathbf{V}_h \in \mathbb{R}^{n \times n}` form an
orthonormal basis for the normalized DOF space.  Row :math:`i` of
:math:`\mathbf{V}_h` is v-mode :math:`i`.  These are the "characteristic
modes" :math:`\tilde{\mathbf{v}}_i` of
`Megias Homar et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024ApJ...974..108M>`_
(Eq. 2 for the unscaled matrix, Eq. 10 for the rescaled one).

**Each v-mode is a linear combination of all (normalized) DOFs.**
Concretely, v-mode :math:`i` is the vector
:math:`\mathbf{v}_i = (\mathbf{V}_h)_{i,:} \in \mathbb{R}^{n}`, which
specifies a coordinated adjustment of hexapod positions and bending modes.

The associated singular value :math:`\sigma_i` quantifies the **optical
leverage** of v-mode :math:`i` over the whole field: the amount of
wavefront change (in double Zernike amplitude) produced per unit of that
DOF combination.  Because the singular values are sorted in decreasing
order, v-mode 0 is the most optically important mode and v-mode
:math:`n{-}1` is the least.

**The number of v-modes equals the number of used DOFs.**  With
:math:`n = 20` used DOFs, :math:`\mathbf{V}_h \in \mathbb{R}^{20 \times 20}`
and there are exactly 20 v-modes.

**V-modes are a property of the optical model, not of the exposure.**
They are computed from the whole-field double Zernike array and are
therefore independent of the number of sensors that reported, of the
rotator angle, and of the Zernike subset in use.  They change only when
the sensitivity matrix file, the normalization weights, or the set of used
DOFs changes.


Coordinate transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^

The code provides two conversions, both restricted to the :math:`r`
retained modes.

**DOFs to v-modes** (``StateEstimator.get_vmodes_from_dofs``).  Given a
DOF vector :math:`\mathbf{x}`:

.. math::

   \mathbf{v} = \mathbf{V}_{h,r}\,\mathbf{N}^{-1}\,\mathbf{x}_{\mathcal{D}}
   \;\in\; \mathbb{R}^{r}

**V-modes to DOFs** (``StateEstimator.get_dofs_from_vmodes``).  Given a
v-mode vector :math:`\mathbf{v}` (of length :math:`r` or longer; only the
first :math:`r` entries are used):

.. math::

   \mathbf{x} = \mathbf{N}\,\mathbf{V}_{h,r}^{T}\,\mathbf{v}_{1:r}
   \;\in\; \mathbb{R}^{n}

Because :math:`\mathbf{V}_{h,r}\mathbf{V}_{h,r}^{T} = \mathbf{I}_r`, going
from v-modes to DOFs and back is lossless.  Going from DOFs to v-modes and
back is lossless only for DOF vectors that already lie in the span of the
first :math:`r` v-modes; any component along the discarded modes is
removed.  With :math:`r = n` both directions are exact.

Both methods live in
`state_estimator.py <https://github.com/lsst-ts/ts_ofc/blob/develop/python/lsst/ts/ofc/state_estimator.py>`_.


Physical interpretation of v-modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although each v-mode is a linear combination of *all* DOFs, it produces a
characteristic wavefront pattern, described by the corresponding left
singular vector :math:`\mathbf{u}_i = \mathbf{U}_{:,i}`.  Since the SVD is
taken over the full-field array, :math:`\mathbf{u}_i` is a set of double
Zernike coefficients.  It describes both which pupil aberration the mode
produces and how that aberration varies across the field.

If :math:`\mathbf{u}_i` has most of its energy in the coefficients of a
particular pupil Zernike, say focus (Z4), then v-mode :math:`i` is
effectively "the focus mode": the coordinated DOF adjustment that most
efficiently changes the telescope's focus.  Similarly, if
:math:`\mathbf{u}_j` is dominated by spherical aberration (Z11), then
v-mode :math:`j` is "the spherical mode."

This is why operators associate specific v-modes with specific aberrations:

- V-mode :math:`i` is named after the **dominant pupil Zernike** in
  :math:`\mathbf{u}_i`.
- The association is approximate.  Each v-mode typically produces a
  dominant Zernike plus smaller contributions to others, and the name
  ignores the field dependence.
- The SVD discovers these "optical eigenmodes" automatically from the
  sensitivity matrix; they are not prescribed by the user.


Truncation and v-mode importance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The singular values :math:`\sigma_i` determine which v-modes are
well-constrained by the wavefront measurements:

- **V-modes** :math:`0` **through** :math:`r{-}1` (large :math:`\sigma_i`):
  these DOF combinations produce large, easily measured wavefront changes.
  The state estimator can reliably reconstruct their amplitudes.

- **V-modes** :math:`r` **through** :math:`n{-}1` (small :math:`\sigma_i`):
  these combinations produce tiny wavefront changes that are buried in
  measurement noise.  They are excluded from the fit.

With ``truncation_index = 12`` and :math:`n = 20` DOFs, the first 12
v-modes are retained and the remaining 8 are never estimated.  The state
vector handed to the controller still has :math:`n` entries in either
basis.  However, its components along the discarded modes are exactly
zero, so the controller never commands them.


.. _PID control:

PID control
-----------

The PID controller computes a correction :math:`\mathbf{u}_k` at each
iteration :math:`k` based on the estimated state.  The correction is
applied with **negative feedback**: the telescope moves in the opposite
direction of the estimated error.

The overall flow in ``OFC.calculate_corrections()`` is:

1. Estimate the optical state :math:`\hat{\mathbf{x}}_k` (or
   :math:`\hat{\mathbf{v}}_k`) via the state estimator.
2. Pass it to the PID controller, which computes an internal control
   signal :math:`\mathbf{u}_k^{\mathrm{pid}}`.
3. **Negate** the control signal: ``uk = -controller.control_step(...)``.
4. If controlling in v-mode space, convert :math:`\mathbf{u}_k` back to
   DOFs.
5. Aggregate the correction into the telescope state.

The negation in step 3 ensures that if the estimated state indicates a
positive defocus, the correction drives it back toward zero.

.. note::

   ``ts_ofc`` also provides an **Optimal Integral Controller**
   (``OICController``).  Instead of a proportional step on the estimated
   state, it chooses the correction that minimizes a quadratic cost
   combining the predicted image-quality-weighted wavefront error with a
   penalty on actuator motion, referenced to the initial state
   (``xref`` = ``x00``, ``x0`` or ``0``); see Eq. 12 of
   `Megias Homar et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024ApJ...974..108M>`_.
   The ``name`` key of the controller file selects the controller, and the
   deployed ``init.yaml`` uses ``name: PID``.  The OIC is out of scope for
   this note.


The setpoint: what the PID targets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The PID controller drives the estimated state toward a **setpoint vector**
:math:`\mathbf{s} \in \mathbb{R}^{50}`, which specifies the desired value for
each of the 50 possible degrees of freedom.  In the deployed configuration
all setpoint entries are zero:

.. math::

   \mathbf{s} = \mathbf{0}

Only the entries :math:`\mathbf{s}_{\mathcal{D}}` that correspond to the
active DOFs are used.  A zero setpoint means the controller targets a
state where all DOF deviations are zero, that is, the optics are in their
nominal alignment.  A non-zero setpoint would let the controller target a
deliberate offset, for example an intentional defocus for engineering
purposes.

**The effective target of the closed loop is not simply "zero physical DOFs."**
The PID setpoint lives in DOF space, but what "zero" means depends on the
measurement preprocessing (see `Measurement preprocessing`_).  The system
converges when :math:`\hat{\mathbf{x}} = \mathbf{0}`.  This corresponds
to a residual :math:`\mathbf{y} \approx \mathbf{0}`, and therefore to
:math:`\mathbf{y}_\mathrm{meas} \approx \mathbf{y}_\mathrm{intrinsic} + \mathbf{y}_2`.
The telescope is thus driven to the state where its measured wavefront
matches the intrinsic design prediction plus the empirical
:math:`\mathbf{y}_2` offsets.  That is the state of best image quality,
as determined from on-sky calibration.


PID formula (as implemented)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At each iteration :math:`k`, the controller
(``BaseController.calculate_pid_step``) receives the estimated state
:math:`\hat{\mathbf{x}}_k` and computes the terms below.  The formulas are
the same in DOF space and in v-mode space.  Note that
`Megias Homar et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024ApJ...974..108M>`_
write the PID in incremental (velocity) form in their Eq. 13.  The
implementation uses the positional form given here, which acts on the
error with respect to the setpoint.

**Error signal:**

.. math::

   \mathbf{e}_k = \mathbf{s}_{\mathcal{D}} - \hat{\mathbf{x}}_k

**Integral term.**  Two variants exist, selected by
``use_leaky_integrator``.  The deployed configuration uses the leaky one.

*Leaky integrator* (``use_leaky_integrator = true``).  This is a
finite-history sum over the last :math:`n_I` errors (``n_iterms``), in
which older errors are down-weighted geometrically by :math:`\lambda`
(``i_factor``):

.. math::

   \mathbf{I}_k^{\mathrm{raw}}
   = \sum_{j=0}^{\min(k,\,n_I)-1} \lambda^{j}\,\mathbf{e}_{k-j}

The newest error has weight 1, the previous one :math:`\lambda`, the one
before that :math:`\lambda^2`, and so on.  With the deployed values
:math:`n_I = 10` and :math:`\lambda = 0.5`, a constant error gives a sum
bounded by
:math:`\lVert\mathbf{e}\rVert/(1-\lambda) = 2\,\lVert\mathbf{e}\rVert`.
The integral therefore cannot wind up indefinitely, even before
clipping.

*Plain integrator* (``use_leaky_integrator = false``).  The errors are
accumulated without decay:

.. math::

   \mathbf{I}_k^{\mathrm{raw}} = \mathbf{I}_{k-1} + \mathbf{e}_k

*Anti-windup clipping.*  Both variants are then clipped:

.. math::

   \mathbf{I}_k = \operatorname{clip}\!\bigl(
       \mathbf{I}_k^{\mathrm{raw}},\;
       -\mathbf{I}_{\max,\mathcal{D}},\;
       +\mathbf{I}_{\max,\mathcal{D}}
   \bigr)

where :math:`\mathbf{I}_{\max}` is the ``max_integral`` vector from the
controller file.  The deployed ``init.yaml`` sets per-DOF values (see
`Deployed configuration`_); the code defaults to ones when the key is
absent.

**Derivative with low-pass filter:**

.. math::

   \mathbf{d}_k^{\mathrm{raw}} = \mathbf{e}_k - \mathbf{e}_{k-1}

.. math::

   \mathbf{d}_k = \alpha\,\mathbf{d}_k^{\mathrm{raw}}
                  + (1 - \alpha)\,\mathbf{d}_{k-1}

where :math:`\alpha` is the ``derivative_filter_coeff``.  When
:math:`\alpha = 1` (the deployed value), no filtering is applied and the
derivative is simply the difference between consecutive errors.  When
:math:`\alpha < 1`, the derivative is smoothed with an exponential moving
average to reduce sensitivity to noise.  After a history reset,
:math:`\mathbf{e}_{k-1}` and :math:`\mathbf{d}_{k-1}` are zero.

**Control output:**

.. math::

   \mathbf{u}_k^{\mathrm{pid}} = \mathbf{K}_{p,\mathcal{D}}\,\mathbf{e}_k
                               + \mathbf{K}_{i,\mathcal{D}}\,\mathbf{I}_k
                               + \mathbf{K}_{d,\mathcal{D}}\,\mathbf{d}_k

where :math:`\mathbf{K}_p`, :math:`\mathbf{K}_i` and :math:`\mathbf{K}_d`
are diagonal gain matrices built from ``kp``, ``ki`` and ``kd``.  Each can
be a scalar, applied uniformly to all DOFs, or a 50-element array for
per-DOF tuning.  The subscript :math:`\mathcal{D}` selects the active
entries.

The final correction applied to the telescope is:

.. math::

   \mathbf{u}_k = -\mathbf{u}_k^{\mathrm{pid}}

This negation, applied in ``OFC.calculate_corrections()``, turns the PID
output, which points from the state toward the setpoint, into a physical
correction that drives the optics back toward the setpoint.


.. _Deployed configuration:

Deployed configuration
^^^^^^^^^^^^^^^^^^^^^^

The controller parameters from ``MTAOS/ofc/configurations/init.yaml`` in
``ts_config_mttcs``:

.. list-table::
   :header-rows: 1
   :widths: 28 24 48

   * - Parameter
     - Deployed value
     - Description
   * - ``name``
     - ``PID``
     - Controller type (``PID`` or ``OIC``)
   * - ``kp``
     - 0.30 (all DOFs)
     - Proportional gain (scalar or per-DOF array)
   * - ``ki``
     - 0.05 (all DOFs)
     - Integral gain
   * - ``kd``
     - 0.0
     - Derivative gain (disabled)
   * - ``derivative_filter_coeff``
     - 1.0
     - Low-pass filter coefficient (1.0 = no filtering)
   * - ``use_leaky_integrator``
     - true
     - Finite-history integral with geometric decay
   * - ``n_iterms``
     - 10
     - Number of error terms kept in the leaky integral
   * - ``i_factor``
     - 0.5
     - Decay factor :math:`\lambda` per iteration of age
   * - ``setpoint``
     - all zeros (50 elements)
     - Target state for each DOF
   * - ``max_integral``
     - hexapods: 500, 100, 100, 0.01, 0.01; bending modes 1--10: 0.5; modes 11--20: 0
     - Anti-windup clipping limit per DOF (µm or arcsec)
   * - ``zn_selected``
     - Z4--Z22
     - Zernikes used (overridden by WEP output in the closed loop)
   * - ``normalization_weights_filename``
     - ``range0.5_fwhm-0.15.yaml``
     - DOF normalization weights
   * - ``truncation_index``
     - 12
     - Number of v-modes retained (:math:`r`)
   * - ``rotation_offset``
     - 0.0 deg
     - Added to the rotator angle before de-rotation
   * - (``noise_covariance_filename``)
     - absent
     - Identity covariance is used

**Physical meaning of the gains.**  :math:`k_p = 0.30` means "correct 30%
of the estimated error per iteration".  With :math:`k_i = 0.05` and the
leaky integral bounded by twice the error, the integral term adds at most
another 10% of a persistent error.  For a constant error, the correction
per iteration therefore approaches 40% of that error.  With
:math:`k_d = 0` there is no derivative action.  MTAOS can also override
``kp`` for an individual image (user gain, elevation scaling) and all
three gains after a filter change; see `Where OFC sits in the closed loop`_.

The packaged ``pid_controller.yaml`` in ``ts_ofc`` (``kp = 0.7``,
``ki = 0``, ``kd = 0.1``) is **not** what runs on the summit.


Controller state and history
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The controller maintains internal state between iterations:

- **aggregated_state** (:math:`\mathbf{x}_\mathrm{agg}`): the cumulative
  correction applied to the telescope, initialized from ``dof_state0``.
- **previous_error**: a queue of the last :math:`n_I` errors, used by both
  the leaky integral and the derivative.
- **integral** (:math:`\mathbf{I}`): the clipped integral; the unclipped
  value is kept alongside for diagnostics.
- **filtered_derivative** (:math:`\mathbf{d}_{k-1}`): filtered derivative
  from the previous step.

The ``reset_history()`` method:

1. Saves the current aggregated state as the new reference:
   :math:`\mathbf{x}_0 \leftarrow \mathbf{x}_\mathrm{agg}`.
2. Empties the error queue and zeros the integral and derivative terms.

It is called when the closed loop is started or stopped with a new
``comp_dof_idx`` configuration, so that the PID restarts from the current
telescope state.


PID in DOF space
^^^^^^^^^^^^^^^^

When ``control_vmodes = False`` (the deployed setting), the state estimator
returns the optical state in DOF space:

.. math::

   \hat{\mathbf{x}}_k = \mathbf{N}\,\mathbf{V}_{h,r}^T\,\hat{\mathbf{v}}_k

The PID operates on this :math:`n`-dimensional state vector.  Each gain
entry corresponds directly to a physical DOF.

**Physical meaning of DOF gains.**  Each diagonal entry of
:math:`\mathbf{K}_p` sets how aggressively one physical DOF is corrected.
For example, :math:`k_{p,0} = 0.30` means "correct 30% of the estimated
M2 hexapod dZ error per iteration".


.. _PID in v-mode space:

PID in v-mode space
^^^^^^^^^^^^^^^^^^^

When ``control_vmodes = True``, the state estimator returns the optical
state in v-mode space, zero-padded to :math:`n` elements:

.. math::

   \hat{\mathbf{v}}_k = \bigl(\hat{v}_{k,1}, \ldots, \hat{v}_{k,r},
                              0, \ldots, 0\bigr)

The PID then operates on this vector with the **same**
``calculate_pid_step()`` method.  The controller does not know whether it
is operating on DOFs or v-modes.  The error, integral and derivative are all
computed in v-mode space:

.. math::

   \mathbf{e}^{(v)}_k = \mathbf{s}_{\mathcal{D}} - \hat{\mathbf{v}}_k

.. math::

   \mathbf{u}^{(v)}_k = \mathbf{K}_{p,\mathcal{D}}\,\mathbf{e}^{(v)}_k
                      + \mathbf{K}_{i,\mathcal{D}}\,\mathbf{I}^{(v)}_k
                      + \mathbf{K}_{d,\mathcal{D}}\,\mathbf{d}^{(v)}_k

After the PID produces :math:`\mathbf{u}^{(v)}_k`, OFC converts it back
to physical DOF corrections with ``get_dofs_from_vmodes()``:

.. math::

   \mathbf{u}_k = -\mathbf{N}\,\mathbf{V}_{h,r}^{T}\,\mathbf{u}^{(v)}_{k,1:r}

Only the first :math:`r` entries of the v-mode correction are used.

**Important implementation detail.**  The PID uses the same gain arrays,
indexed by :math:`\mathcal{D}`, regardless of the basis.  The gain at
position :math:`i` of the active list therefore applies to v-mode
:math:`i`, with v-modes ordered by singular value.  V-mode 0, which has
the highest optical leverage, gets the gain of the first active DOF,
v-mode 1 the gain of the second, and so on.  The setpoint is interpreted
positionally in the same way.  With an all-zero setpoint this is harmless,
but a non-zero setpoint only makes sense in DOF space.

**Physical meaning of v-mode gains.**  Each entry of :math:`\mathbf{K}_p`
sets how aggressively one optical mode is corrected.  For example, if
v-mode 0 is the focus mode, giving it a higher gain makes the system
correct focus errors more aggressively than higher-order aberrations.

**The advantage of v-mode control** is the ability to assign different
gains based on optical importance:

- **High gain** on well-constrained v-modes, which have large
  :math:`\sigma_i` and reliable estimates, for fast convergence.
- **Low gain** on poorly constrained v-modes, which have small
  :math:`\sigma_i` and noisy estimates, to avoid injecting noise into the
  corrections.


Equivalent DOF-space gain
^^^^^^^^^^^^^^^^^^^^^^^^^

Composing the three steps (DOF to v-mode, gain, v-mode to DOF), the
effective DOF-space proportional gain when controlling in v-mode space is:

.. math::

   \mathbf{K}_\mathrm{eff}
   = \mathbf{N}\,\mathbf{V}_{h,r}^{T}\,\mathbf{K}^{(v)}_{p}\,
     \mathbf{V}_{h,r}\,\mathbf{N}^{-1}

When the v-mode gain is a **uniform scalar**
(:math:`\mathbf{K}^{(v)}_p = k_p \mathbf{I}_r`):

.. math::

   \mathbf{K}_\mathrm{eff}
   = k_p\,\mathbf{N}\,\mathbf{P}_r\,\mathbf{N}^{-1}, \qquad
   \mathbf{P}_r = \mathbf{V}_{h,r}^{T}\,\mathbf{V}_{h,r}

where :math:`\mathbf{P}_r` is the orthogonal projector onto the retained
subspace.  Since the DOF-space estimate :math:`\hat{\mathbf{x}}_k` already
lies in that subspace,
:math:`\mathbf{K}_\mathrm{eff}\hat{\mathbf{x}}_k = k_p\hat{\mathbf{x}}_k`.
V-mode control with a uniform gain is therefore identical to DOF-space
control with the same scalar gain.

When the v-mode gain has **different values per mode**, the matrix
:math:`\mathbf{V}_{h,r}^{T}\,\mathbf{K}^{(v)}_p\,\mathbf{V}_{h,r}` is a
full (non-diagonal) matrix that **couples the DOFs**.  This is the
intended behavior, since correcting an optical eigenmode such as focus
requires coordinated adjustments of several physical DOFs.

For example, a hypothetical two-level gain structure with
:math:`k_p^{(v)} = 0.30` for v-modes 0--5 and :math:`k_p^{(v)} = 0.15` for
v-modes 6--11 decomposes as:

.. math::

   \mathbf{K}^{(v)}_p
   = 0.15\,\mathbf{I}_r
   + 0.15\,\mathbf{P}_{6}

where :math:`\mathbf{P}_{6}` is the projection onto the first 6 v-modes.
In words, a uniform base gain is applied to all retained modes, with an
extra boost on the most optically important ones.


Correction aggregation
----------------------

After each PID step, the DOF correction :math:`\mathbf{u}_k` is added to
the telescope state (``BaseController.aggregate_state``):

.. math::

   \mathbf{x}_{k+1} = \mathbf{x}_k + \mathbf{u}_k

Only the active entries :math:`\mathcal{D}` of the 50-element aggregated
state are updated.  The per-visit correction :math:`\mathbf{u}_k` is also
stored on its own (``lv_dof``) so that it can be undone if a subsystem
rejects it.

The aggregated state is then decomposed into physical corrections for the
individual components (``OFC.get_correction``):

- **M2 Hexapod:** DOFs 0--4 (dZ, dX, dY, rX, rY), mapped to hexapod
  command coordinates through the (pseudo-)inverse of the hexapod rotation
  matrix in ``comp_dof_idx``.
- **Camera Hexapod:** DOFs 5--9, same transformation.
- **M1M3 Bending Modes:** DOFs 10--29, converted from bending mode
  amplitudes to actuator forces with the bending-mode-to-force influence
  matrix (``BendModeToForce``).
- **M2 Bending Modes:** DOFs 30--49, same conversion.

``ts_mtaos`` issues the resulting ``Correction`` objects to the subsystem
CSCs, after the checks described in the next section.


.. _Where OFC sits in the closed loop:

Where OFC sits in the closed loop
---------------------------------

This section summarizes the parts of the MTAOS closed loop
(``MTAOS.run_closed_loop``) that change what OFC computes or what is
finally applied.  The full operational description, with diagrams of the
main flow and of ``_execute_ofc``, is in the Closed Loop Operations
section of the `MTAOS documentation <https://ts-mtaos.lsst.io>`_.  The
configuration values below come from ``MTAOS/v13/_init.yaml`` in
``ts_config_mttcs``.

Per-image sequence
^^^^^^^^^^^^^^^^^^

For every image ingested by the Observatory Operations Data Service
(OODS), the loop:

1. Checks that the image is one the loop is following, and that it is not
   stale.
2. Reads the filter and the mean elevation of the exposure, and the mean
   rotator angle recorded during the exposure.
3. Skips the image if the telescope has moved since it was taken: more
   than ``rotation_delta_limit`` (9 deg) in rotator angle or more than
   ``elevation_delta_limit_max`` (9 deg) in elevation relative to the
   current position.
4. Runs or polls the wavefront estimation (WEP) and collects the Zernikes.
5. Computes the gain for this image (see `Gain overrides`_).  A zero gain
   skips the OFC step entirely.
6. Calls ``_execute_ofc``, which temporarily overrides the controller
   gains, runs ``Model.calculate_corrections`` (large-defocus check,
   then ``OFC.calculate_corrections``), publishes the DOF and correction
   events, and restores the gains.
7. Waits for the shutter to close, then applies the corrections (see
   `Stress limits, force thresholds and pointing`_).

The OFC data parameters passed with ``startClosedLoop``, for instance
``truncation_index`` or ``comp_dof_idx``, are applied once for the whole
session and reset the controller history.  For each image, only
``filter_name`` and ``rotation_angle`` are passed through to OFC.

Gain overrides
^^^^^^^^^^^^^^

Three mechanisms can change the gains that the PID uses for a given
image.  MTAOS applies all of them on top of the controller file and
restores the original gains after the call.

- **Elevation scaling.**  ``Model.get_correction_gain`` compares the
  elevation of this image with that of the previous *processed* image.
  If the change is at or above ``elevation_delta_limit_max`` the gain is
  zero and the image is skipped; at or below ``elevation_delta_limit_min``
  the gain is the controller's ``kp``; in between ``kp`` is scaled
  linearly to zero.  In the deployed configuration both limits are
  9 degrees, so this is a binary skip rather than a ramp.  The resulting
  gain replaces ``kp`` only; ``ki`` and ``kd`` are untouched.
- **User gain.**  The ``runOFC`` command carries a ``userGain``.  When it
  is non-zero, it replaces ``kp`` for that call in the same way.
- **Filter-change override.**  ``closed_loop_filter_change_gain`` can
  replace ``kp``, ``ki`` and ``kd`` for ``n_iter`` successful OFC
  iterations after a filter change.  The deployed default is
  ``n_iter = 0``, which disables it.  An optional override file enables it
  with ``n_iter = 4`` and ``kp = 0.75``.

Large-defocus bypass
^^^^^^^^^^^^^^^^^^^^

Before calling OFC, ``Model.calculate_corrections`` takes the median of
the defocus :math:`dz` inferred from the donut radii on the corner
sensors.  If :math:`|dz|` exceeds ``dz_threshold_min`` (300 µm):

- with ``raise_on_large_defocus = True`` (deployed), the computation
  raises an error, no correction is computed, and the failure counts
  towards ``max_ofc_consecutive_failures`` (3);
- with ``raise_on_large_defocus = False``, the value is clipped to
  ``dz_threshold_max`` (1500 µm).  Half of it is applied directly to the
  M2 hexapod dZ and half to the camera hexapod dZ through ``offset_dof``,
  **bypassing the state estimator and the PID**.  The offset is still
  added to the aggregated state.

Stress limits, force thresholds and pointing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the corrections are applied (``MTAOS.handle_corrections``):

- **Mirror stress limits.**  The bending-mode stresses of the aggregated
  state are combined as a root sum of squares and multiplied by
  ``stress_scale_factor`` (1.25).  If the result exceeds
  ``m1m3_stress_limit`` (137900 Pa) or ``m2_stress_limit`` (344737 Pa),
  the bending modes of that mirror are reduced.  The ``scale`` approach
  (deployed) reduces them proportionally; the ``truncate`` approach zeroes
  the highest modes one by one.  The aggregated state in the controller
  is overwritten with the reduced values, and the corrections are
  recomputed.  This is the only place where the applied state can differ
  from what the PID produced.
- **Force deltas.**  M1M3 and M2 force commands are skipped when the
  change from the currently applied forces is below
  ``m1m3_min_forces_to_apply`` or ``m2_min_forces_to_apply``,
  respectively.  Hexapod commands are always issued.
- **Pointing correction.**  When ``enable_pointing_correction`` is true
  (deployed), the aggregated DOF vector is multiplied by the
  :math:`50 \times 2` ``pointing_correction_matrix``.  The hexapod rows
  give a pointing origin offset in millimeters, and the remaining rows
  give a collimation offset in arcseconds after the plate-scale factor.
  Both are sent to the pointing component so that the AOS corrections do
  not move the field.
- **Undo on failure.**  If any subsystem rejects its correction, the
  successful ones are undone and the last-visit correction is subtracted
  from the aggregated state.


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
     - :math:`n` = number of used DOFs (20)
     - :math:`n` (same; from the SVD of the full-field model)
   * - Truncation
     - Estimate confined to the first :math:`r` v-modes
     - Same; only the first :math:`r` v-mode corrections are converted
   * - Deployed
     - Yes (``control_vmodes = False``)
     - Available (``test_control_vmodes.yaml``)
