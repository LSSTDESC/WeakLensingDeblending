Source Catalogs
===============

Galaxy Catalog Format
---------------------

The input catalog is a text file consisting of an initial header line and followed by one line of numbers per catalog object, for example a catalog with two objects might look like this::

	id ra dec redshift fluxnorm_bulge fluxnorm_disk fluxnorm_agn a_b a_d b_b b_d pa_bulge pa_disk u_ab g_ab r_ab i_ab z_ab y_ab
	11962286 0.593881845474 0.000286099995719 2.02162647247 0.0 5.9757610912e-19 0.0 0.0 0.132357493043 0.0 0.130808100104 0.0 144.235595703 26.9993686676 26.9311523438 26.9816989899 27.0538825989 27.0064048767 26.9156532288
	12663653 0.594969153404 0.00124180002604 0.459954589605 0.0 3.85022202105e-19 0.0 0.0 0.46328830719 0.0 0.430979907513 0.0 325.863311768 29.0297031403 28.4987335205 27.4711303711 27.1866264343 27.030462265 26.9395713806

The header line specifies a list of names, separated by white space, associated with each catalog parameter. Each catalog entry is represented by a a list of numbers, separated by white space, giving the corresponding parameter values.

The parameter names below are required by the :ref:`prog-simulate` program, but other parameters may also be present in the file. Parameters can be listed in any order as long as the order of values matches the header line. The names in the table below are a bit idiosyncratic (e.g., mixing under_scores with CamelCase and using different schemes to denote bulge vs. disk) but are chosen to match the names used in the `LSST DM galaxy catalog schema <https://confluence.lsstcorp.org/display/SIM/Database+Schema>`_.

==================== ===========
Name                 Description
==================== ===========
id                   Unique integer identifier for this object (the primary key for objects in the LSST db)
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
==================== ===========

The catalog file is read using a :py:class:`astropy.io.ascii.Basic` reader (created with default options) so can be embellished with comments and blank lines for readability, as supported by that class.

Convert ASCII Catalog to FITS Table
-----------------------------------

The catalog file can also be read as a FITS file containing a single table. A catalog in text format can be converted to a FITS file using, for example::

	import astropy.table
	catalog = astropy.table.Table.read('OneDegSq.dat',format='ascii.basic')
	catalog.write('OneDegSq.fits')

The resulting FITS file will be somewhat smaller (by about 30%) than the text original and significantly faster for the :ref:`prog-simulate` program to read, but may be less convenient for other programs to read or generate.

.. _catalog-create:

Create Galaxy Catalog From LSST Catalog Database
------------------------------------------------

This section documents the process for creating a galaxy catalog from the LSST Data Management database. This is not something you will normally need (or want) to do since suitable catalogs are already provided. However, if you want to know exactly how the provided catalogs were created or do need to create your own, read on.

The `dbquery.py` program automates the process of connecting to the database, extracting the relevant data, and writing a catalog file.  The program uses the `pymsql <http://pymssql.org/en/stable/>`_ python interface to Microsoft SQL Server (which LSST DM uses), so you will need to install that and its dependencies (`FreeTDS <http://www.freetds.org>`_ and `Cython <http://cython.org>`_) in order to run `dbquery.py`.

The `OneDegSq.dat` catalog file was created using::

	dbquery.py -o OneDegSq.dat --dec-min -0.5 --dec-max +0.5 --ra-min 0.0 --ra-max 1.0 --verbose

with FreeTDS v0.91, Cython v0.21, pymsql v2.1.1 under OS-X 10.10.1.  The program takes about ? seconds to run. The set of 

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
