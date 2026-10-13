#########################################
OFC Control Loop Mathematical Description
#########################################

.. abstract::

   This technical note provides a mathematical description of the
   Optical Feedback Control (OFC) system used in the Vera C. Rubin
   Observatory's Active Optics System (AOS), as implemented in ``ts_ofc``,
   together with the decisions that the ``ts_mtaos`` closed loop makes
   around each OFC call.  The measurement side (wavefront estimation in
   ``ts_wep``) is out of scope; the wavefront error is treated as an input.
   All quoted defaults come from the summit configuration in
   ``ts_config_mttcs``.

   The document begins with the **forward optical model**, explaining how
   the double Zernike sensitivity matrix relates telescope degrees of
   freedom (DOFs) to wavefront Zernike coefficients and how DOF
   normalization balances their contributions.  It then describes how the
   measurements are de-rotated and corrected for known terms before
   inversion.

   The core of the note addresses the **state estimation inverse problem**:
   recovering DOF states from noisy wavefront measurements.  The textbook
   regularizations (SVD truncation and noise-weighted least squares) are
   introduced and then contrasted with what ``ts_ofc`` actually computes:
   a least-squares fit confined to the leading **v-modes** of the
   full-field optical model.  A concrete numerical example gives the
   matrix dimensions for the deployed 76-measurement, 20-DOF,
   12-mode configuration.

   A dedicated section defines v-modes, explains their physical
   interpretation as coordinated DOF adjustments that produce characteristic
   aberration patterns, and describes the coordinate transformations
   between DOF space and v-mode space.

   The note then covers **PID control** as implemented, including the
   leaky integrator, anti-windup clipping and the derivative filter, in
   both DOF space and v-mode space, and shows how v-mode gains enable
   tuning by optical importance.  A final section describes where OFC sits
   in the MTAOS closed loop: gain overrides, the large-defocus bypass,
   mirror stress limits and the pointing correction.

   This document serves as a reference for developers working on AOS
   software and operators seeking to understand the control loop
   mathematics.

.. include:: vmodes_and_control.rst
