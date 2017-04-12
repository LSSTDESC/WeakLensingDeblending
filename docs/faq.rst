Frequently Asked Questions
##########################

.. contents::
    :local:
    :depth: 2

How do I...
===========

…get started with this package?
-------------------------------

Start with the :doc:`quickstart tutorial<quickstart>` to get an overview of this package, then select one of the :doc:`installation methods<install>`.

…simulate my own source catalogs?
---------------------------------

If your source objects are galaxy-like, you can format you catalog with the columns :ref:`described here<catalog-format>`.  Otherwise, :ref:`create an issue<faq-feedback>` describing your catalog and what you are trying to do.

…simulate other instruments?
----------------------------

If the instrument you would like to simulate is similar to one of the defaults (LSST, DES, CFHT), you can simply modify a few parameters.  See the :doc:`quickstart tutorial<quickstart>` for examples.  You can also define a new default survey configuration by modifying the file ``descwl/survey.py``.  Instruments are described with a simple model that may not be sufficient for your needs: in this case, please :ref:`create an issue<faq-feedback>` describing what you are trying to do.

…add my own pixel-level analysis?
---------------------------------

In case the :ref:`existing output catalog<output>` does not already calculate the quantities you need, you can either write your own post-processor using the `skeleton code provided <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending/blob/master/skeleton.py>`_, or else modify the ``descwl/analysis.py`` file that generates the output catalog.  In either case, feel free to :ref:`create an issue<faq-feedback>` to let us know what you are trying to do and get feedback.

.. _faq-feedback:

…ask questions or report bugs?
------------------------------

The code is hosted on `github <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending>`_.  Please use the `issue tracker <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending/issues>`_ to let us know about any issues you have with installing or running this code, or to request new features.

…contribute to the code?
------------------------

This software is open source and your contributions are welcome! General information for package developers is :doc:`here<developer>`. If you would like to add a new feature, please start by :ref:`creating a new issue<faq-feedback>` to describe it.
