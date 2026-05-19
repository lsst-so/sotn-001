.. image:: https://img.shields.io/badge/sotn--001-lsst.io-brightgreen.svg
   :target: https://sotn-001.lsst.io
.. image:: https://github.com/lsst-so/sotn-001/workflows/CI/badge.svg
   :target: https://github.com/lsst-so/sotn-001/actions/

#########################################
OFC Control Loop Mathematical Description
#########################################

SOTN-001
========

This technical note provides a mathematical description of the Optical Feedback Control (OFC) pipeline implemented in ``ts_ofc`` and ``ts_mtaos``. It covers the forward optical model, the state estimation inverse problem, the definition and physical interpretation of v-modes, and how PID gains are applied in both DOF space and v-mode space.

**Links:**

- Publication URL: https://sotn-001.lsst.io
- Alternative editions: https://sotn-001.lsst.io/v
- GitHub repository: https://github.com/lsst-so/sotn-001
- Build system: https://github.com/lsst-so/sotn-001/actions/


Build this technical note
=========================

You can clone this repository and build the technote locally if your system has Python 3.11 or later:

.. code-block:: bash

   git clone https://github.com/lsst-so/sotn-001
   cd sotn-001
   make init
   make html

Repeat the ``make html`` command to rebuild the technote after making changes.
If you need to delete any intermediate files for a clean build, run ``make clean``.

The built technote is located at ``_build/html/index.html``.

Publishing changes to the web
=============================

This technote is published to https://sotn-001.lsst.io whenever you push changes to the ``main`` branch on GitHub.
When you push changes to a another branch, a preview of the technote is published to https://sotn-001.lsst.io/v.

Editing this technical note
===========================

The main content of this technote is in ``index.rst`` (a reStructuredText file).
Metadata and configuration is in the ``technote.toml`` file.
For guidance on creating content and information about specifying metadata and configuration, see the Documenteer documentation: https://documenteer.lsst.io/technotes.
