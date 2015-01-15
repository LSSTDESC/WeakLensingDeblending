Source Catalogs
===============

Galaxy Catalog Format
---------------------

The input catalog is a text file (or an :ref:`equivalent FITS table <catalog-to-fits>`) consisting of an initial header line and followed by one line of numbers per catalog object, for example a catalog with two objects might look like this::

	galtileid ra dec redshift fluxnorm_bulge fluxnorm_disk fluxnorm_agn a_b a_d b_b b_d pa_bulge pa_disk u_ab g_ab r_ab i_ab z_ab y_ab
	2200871446 0.418319702147 -0.000148399994941 0.496377289295 0.0 1.4144730572e-17 0.0 0.0 0.278649687767 0.0 0.221303001046 0.0 307.344329834 25.9418621063 25.129743576 23.9588813782 23.3607368469 23.0723800659 22.9095973969
	2205921112 0.420028448104 -0.00100259995088 1.89508104324 0.0 1.91501907101e-18 0.0 0.0 0.358063697815 0.0 0.313674807549 0.0 137.791702271 25.848903656 25.867565155 25.9179477692 25.9851398468 25.8779563904 25.7642536163

The header line specifies a list of names, separated by white space, associated with each catalog parameter. Each catalog entry is represented by a a list of numbers, separated by white space, giving the corresponding parameter values.

The parameter names below are required by the :ref:`prog-simulate` program, but other parameters may also be present in the file. Parameters can be listed in any order as long as the order of values matches the header line. The names in the table below are a bit idiosyncratic (e.g., mixing under_scores with CamelCase and using different schemes to denote bulge vs. disk) but are chosen to match the names used in the `LSST DM galaxy catalog schema <https://confluence.lsstcorp.org/display/SIM/Database+Schema>`_.

==================== ==========================================================================
Name                 Description
==================== ==========================================================================
id                   Unique integer identifier for this object (stored procedure `galtileid`)
ra                   Object centroid right ascension (degrees)
dec                  Object centroid declination (degrees)
redshift             Object cosmological redshift
fluxnorm_bulge       Multiplicative scaling factor to apply to the bulge SED
fluxnorm_disk        Multiplicative scaling factor to apply to the disk SED
fluxnorm_agn         Multiplicative scaling factor to apply to the AGN SED
a_b                  Semi-major axis of bulge 50% isophote (arcseconds)
b_b                  Semi-minor axis of bulge 50% isophote (arcseconds)
a_d                  Semi-major axis of disk 50% isophote (arcseconds)
b_d                  Semi-minor axis of disk 50% isophote (arcseconds)
pa_bulge             Position angle of bulge (degrees)
pa_disk              Position angle of disk (degrees)
u_ab                 Apparent AB magnitude in the LSST u-band, including extinction effects 
g_ab                 Apparent AB magnitude in the LSST g-band, including extinction effects 
r_ab                 Apparent AB magnitude in the LSST r-band, including extinction effects 
i_ab                 Apparent AB magnitude in the LSST i-band, including extinction effects 
z_ab                 Apparent AB magnitude in the LSST z-band, including extinction effects 
y_ab                 Apparent AB magnitude in the LSST y-band, including extinction effects 
==================== ==========================================================================

The catalog file is read using a :py:class:`astropy.io.ascii.Basic` reader (created with default options) so can be embellished with comments and blank lines for readability, as supported by that class.

.. _catalog-create:

Create Galaxy Catalog From LSST Catalog Database
------------------------------------------------

This section documents the process for creating an input catalog derived from the `LSST Data Management <http://dm.lsst.org>`_ (DM) simulation database, described in `these SPIE proceedings <http://dx.doi.org/10.1117/12.2054953>`_ and being `documented here <https://confluence.lsstcorp.org/display/SIM/Catalog+Simulations+Documentation>`_. This is not something you will normally need (or want) to do since the resulting catalog is already :doc:`available to download</products>`. However, if you want to know exactly how the catalog was created or want to create your own, read on.

.. warning::
	In case you are using the DM galaxy database directly, you should **avoid using the DiskHalfLightRadius and BulgeHalfLightRadius columns** since these do not have the usual definition where `hlr == sqrt(a*b)` and instead satisfy `hlr == a`. To avoid this confusion, use the `a,b` parameters directly and calculate `hlr == sqrt(a*b)`.

The `dbquery.py` program automates the process of connecting to the database, extracting the relevant data, and writing a catalog file.  The program uses the `pymsql <http://pymssql.org/en/stable/>`_ python interface to Microsoft SQL Server (which LSST DM uses), so you will need to install that and its dependencies (`FreeTDS <http://www.freetds.org>`_ and `Cython <http://cython.org>`_) in order to run `dbquery.py`.

The `OneDegSq.dat` catalog file was created using::

	dbquery.py -o OneDegSq.dat --dec-min -0.5 --dec-max +0.5 --ra-min -0.5 --ra-max +0.5 --verbose

with FreeTDS v0.91, Cython v0.21, pymsql v2.1.1 under OS-X 10.10.1.  The program takes about 2 minutes to run and saves 858502 rows to the 200 Mbyte output text file. The set of rows returned by a query is invariant, but they are returned in an arbitrary order so the output files obtained by running this program twice can be different (but only by a permutation of lines). Note that the "1 sq.deg." here is defined by RA and DEC intervals of 1 degree centered on (RA,DEC) = (0,0), which is not exactly 1 sq.deg. of 2D imaging because of projection effects.

The query uses a `custom stored procedure <https://listserv.lsstcorp.org/mailman/private/lsst-imsim/2013-July/42.html>`_ that tiles the full sky by repeating the same 4deg x 4deg patch of about 14M galaxies (the catalog actually contains a 4.5deg x 4.5deg patch, but only the smaller patch is tiled). The stored routine synthesizes the returned `ra`, `dec`, and `galtileid` values, where::

	galtileid = TTTTGGGGGGGG

is a 64-bit integer that combines the tile number `TTTT` (0-4049) with the unique galaxy ID `GGGGGGGG` (0-17,428,283).  The ID that `dbquery` writes to its output file is `galtileid`. For the `OneDegSq.dat` example above, all galaxies are unique and located in tileid = 22 or 4027.

Note that the LSST database uses standard port assignments for its Microsoft SQL Server. However, since these ports are frequently targets of network attacks, many organizations block outbound packets to these ports from internal IP addresses. In addition, the UW hosts of the database only allow inbound packets to their SQL server from IP address ranges included in their whitelist. Therefore, if you are unable to connect, the most likely reason is packet filtering on one or both ends of your connection. One option is to use a host that has already been configured for access to the LSST catalog database (on the NCSA cluster, for example).

You might find the following interactive python snippet useful for debugging connection problems::

	import _mssql
	conn = _mssql.connect(server='fatboy.npl.washington.edu', user='LSST-2', password='L$$TUser', database='LSST', port=1433)
	print conn.tds_version
	conn.execute_query("Select name from sysobjects where type like 'u'")
	for row in conn: print row['name']
	conn.execute_query("select COLUMN_NAME from INFORMATION_SCHEMA.COLUMNS where TABLE_NAME = 'galaxy'")
	for row in conn: print row[0]
	conn.execute_scalar("select count(*) from galaxy")
	conn.close()

The second line will fail with a connection error after about 30 seconds if your packets are being filtered on either end::

	MSSQLDatabaseException: (20009, 'DB-Lib error message 20009, severity 9:\nUnable to connect: Adaptive Server is unavailable or does not exist\nNet-Lib error during Operation now in progress (36)\n')

.. _catalog-to-fits:

Convert ASCII Catalog to FITS Table
-----------------------------------

The catalog file can also be read as a FITS file containing a single table. A catalog in text format can be converted to a FITS file using, for example (this will not overwrite an existing output file)::

	import astropy.table
	catalog = astropy.table.Table.read('OneDegSq.dat',format='ascii.basic')
	catalog.write('OneDegSq.fits')

The resulting FITS file will be somewhat smaller (by about 30%) than the text original and significantly faster for the :ref:`prog-simulate` program to read, but may be less convenient for other programs to read or generate.
