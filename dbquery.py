#!/usr/bin/env python
"""Query the LSST DM simulation galaxy catalog.

Documentation for this program is available at
http://weaklensingdeblending.readthedocs.org/en/latest/programs.html#dbquery
"""

import argparse
import math

import _mssql

def main():

    # Parser command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action = "store_true",
        help = "provide more verbose output")
    parser.add_argument("-o", "--output", default = "gcat.dat", metavar = 'FILE',
        help = "name of output catalog file to write")
    parser.add_argument("--dec-min", type = float, default = -0.5, metavar = 'DEC',
        help = "minimum DEC value to fetch (deg)")
    parser.add_argument("--dec-max", type = float, default = +0.5, metavar = 'DEC',
        help = "maximum DEC value to fetch (deg)")
    parser.add_argument("--ra-min", type = float, default = 0.0, metavar = 'RA',
        help = "minimum RA value to fetch (deg)")
    parser.add_argument("--ra-max", type = float, default = 1.0, metavar = 'RA',
        help = "maximum RA value to fetch (deg)")
    parser.add_argument("--null-sub", type = float, default = -1, metavar = 'SUB',
        help = "numeric value to substitute for any SQL NULLs")
    args = parser.parse_args()

    # Try to open the output file
    try:
        # The filename to write
        f = open(args.output,'w')
    except IOError,e:
        print 'Cannot open output %r for writing' % args.output
        print str(e)
        return -2

    # The ra,dec window to retrieve
    window = { 'RAmin':args.ra_min, 'RAmax':args.ra_max, 'DECmin':args.dec_min, 'DECmax':args.dec_max }
    if args.ra_min >= args.ra_max or args.dec_min >= args.dec_max:
        print 'Invalid RA-DEC window %r' % window
        return -2

    def addColumns(patterns,types):
        text = ''
        for p in patterns:
            text += ',' + ','.join([p%t for t in types])
        return text

    # Specify the header columns we need
    columns = 'id,ra,dec,redshift'

    # Add bulge and disk specific columns
    columns += addColumns(('fluxnorm_%s',),('bulge','disk','agn'))
    columns += addColumns(('a_%s','b_%s'),('b','d'))
    columns += addColumns(('pa_%s',),('bulge','disk'))

    # Add filter-specific columns
    columns += addColumns(('%s_ab',),"ugrizy")

    # calculate the radius in arcmins of a circle enclosing our search box
    radius = 60*0.5*math.sqrt((args.ra_max-args.ra_min)**2 + (args.dec_max-args.dec_min)**2)
    # use the stored procedure described at
    # http://listserv.lsstcorp.org/mailman/private/lsst-imsim/2013-July/42.html
    query = "GalaxySearchSpecColsConstraint2013 @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f, @ColumnNames = '%s', @WhereClause = ''" % (
        0.5*(args.ra_min+args.ra_max),0.5*(args.dec_min+args.dec_max),radius,columns)

    if args.verbose:
        print 'using query: "%s"' % query

    clist = columns.split(',')
    print >>f, ' '.join(clist)
    conn = None
    nulls = { }
    clipCount = 0
    try:
        conn = _mssql.connect(
            server='fatboy.npl.washington.edu', port=1433,
            user='LSST-2', password='L$$TUser',
            database='LSST')
        conn.execute_query(query)
        nrows = 0
        for row in conn:
            # Filter out any SQL NULLs
            for col in clist:
                if row[col] is None:
                    if col not in nulls:
                        nulls[col] = 0
                    nulls[col] += 1
                    row[col] = args.null_sub
            # Skip any objects outside the requested bounding box
            ra = row['ra']
            dec = row['dec']
            if ra < args.ra_min or ra > args.ra_max or dec < args.dec_min or dec > args.dec_max:
                clipCount += 1
                continue
            # Dump this row to our output file
            print >>f, ' '.join([str(row[col]) for col in clist])
            nrows += 1
        if args.verbose:
            print 'Dumped',nrows,'rows to',args.output
            if nulls:
                print 'Replaced NULLs with',args.null_sub,'for:'
                for col in nulls:
                    print '%10d %s' % (nulls[col],col)
            print '%d rows with (ra,dec) outside window were clipped' % clipCount
    except _mssql.MssqlDatabaseException,e:
        print 'Database Exception'
        raise
    finally:
        if conn: conn.close()

    f.close()

if __name__ == "__main__":
    main()
