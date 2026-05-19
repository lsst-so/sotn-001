#########################################
OFC Control Loop Mathematical Description
#########################################

.. abstract::

   This technical note provides a comprehensive mathematical description of the
   Optical Feedback Control (OFC) system used in the Vera C. Rubin Observatory's
   Active Optics System (AOS), as implemented in ``ts_ofc`` and ``ts_mtaos``.

   The document begins with the **forward optical model**, explaining how the
   sensitivity matrix relates telescope degrees of freedom (DOFs) to wavefront
   Zernike coefficients, and how DOF normalization ensures balanced contributions
   during state estimation.

   The core of the note addresses the **state estimation inverse problem**:
   recovering DOF states from noisy wavefront measurements. Two complementary
   regularization approaches are presented — SVD truncation and noise-weighted
   least squares — along with a concrete numerical example showing matrix
   dimensions for the typical 76-measurement, 22-DOF configuration.

   A dedicated section defines **v-modes** (optical eigenmodes from the SVD),
   explains their physical interpretation as coordinated DOF adjustments that
   produce characteristic aberration patterns, and describes coordinate
   transformations between DOF space and v-mode space.

   Finally, the note covers **PID control** in both DOF space and v-mode space,
   showing how v-mode gains enable tuning by optical importance, and how
   corrections are aggregated and dispatched to the telescope subsystems.

   This document serves as a reference for developers working on AOS software
   and operators seeking to understand the control loop mathematics.

.. include:: vmodes_and_control.rst
