wind3d_field_lines
==================

This package provides a Fortran-backed magnetic field-line tracer for wind3d data.

Status
======

This package is under active development for research and learning purposes.

Installation
============

To install the package in editable mode:

.. code-block:: bash

   pip install -e ".[dev,docs]"

Examples
========

Basic import:

.. doctest::

   >>> import wind3d_field_lines as wfl
   >>> hasattr(wfl, "trace_field_lines")
   True

API Reference
=============

.. toctree::
   :maxdepth: 2

   api
   demo_arcade
