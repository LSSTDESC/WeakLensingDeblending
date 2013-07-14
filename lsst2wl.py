#!/usr/bin/env python
#######################################################################################################
## Program to query the official LSST simulation galaxy catalog and write a simple text catalog
## file suitable for weak lensing studies.
##
## Uses the pymsql package for talking to the remote Microsoft SQL server, which in turns requires
## that Cython and FreeTDS are installed.
##
## Usage example: lsst2wl.py --maxrows 100 --output galaxies.dat
##
## Created 07-Nov-2012 by David Kirkby <dkirkby@uci.edu>
#######################################################################################################
import _mssql
import sys
import argparse

def main():

    # Parser command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action = "store_true",
        help = "provide more verbose output")
    parser.add_argument("-o", "--output", default = "gcat.dat",
        help = "name of output catalog file to write")
    parser.add_argument("--dec-min", type = float, default = -1,
        help = "minimum DEC value to fetch (deg)")
    parser.add_argument("--dec-max", type = float, default = 0,
        help = "maximum DEC value to fetch (deg)")
    parser.add_argument("--ra-min", type = float, default = 0,
        help = "minimum RA value to fetch (deg)")
    parser.add_argument("--ra-max", type = float, default = 1,
        help = "maximum RA value to fetch (deg)")
    parser.add_argument("--maxrows", type = int, default = 100,
        help = "maximum number of rows to fetch from the server")
    parser.add_argument("--null-sub", type = float, default = -1,
        help = "numeric value to substitute for any SQL NULLs")
    args = parser.parse_args()

    # Check for a valid maxrows (total rows in the DB is about 17M)
    if args.maxrows <= 0:
        print 'Invalid maxrows %r <= 0'
        sys.exit(-2)

    # Try to open the output file
    try:
        # The filename to write
        f = open(args.output,'w')
    except IOError,e:
        print 'Cannot open output %r for writing' % args.output
        print str(e)
        sys.exit(-2)

    # The ra,dec window to retrieve
    window = { 'RAmin':args.ra_min, 'RAmax':args.ra_max, 'DECmin':args.dec_min, 'DECmax':args.dec_max }
    if args.ra_min >= args.ra_max or args.dec_min >= args.dec_max:
        print 'Invalid RA-DEC window %r' % window
        sys.exit(-2)

    def addColumns(patterns,types):
        text = ''
        for p in patterns:
            text += ',' + ','.join([p%t for t in types])
        return text

    # Specify the header columns we need
    columns = 'id,ra,dec,redshift,absmag_r_total'

    # Add bulge and disk specific columns
    columns += addColumns(('%sHalfLightRadius',),('Bulge','Disk'))
    columns += addColumns(('pa_%s','magnorm_%s','sedid_%s','sedname_%s'),('bulge','disk'))
    columns += addColumns(('a_%s','b_%s'),('b','d'))

    # Add filter-specific columns
    columns += addColumns(('%s_ab',),"ugrizy")

    # SQL filter to use
    filter = 'WHERE ra BETWEEN %(RAmin)s AND %(RAmax)s AND dec BETWEEN %(DECmin)s and %(DECmax)s' % window

    clist = columns.split(',')
    print >>f, ' '.join(clist)
    conn = None
    nulls = { }
    try:
        conn = _mssql.connect(
            server='fatboy.npl.washington.edu', port=1433,
            user='LSST-2', password='L$$TUser',
            database='LSST')
        conn.execute_query("SELECT TOP %d %s FROM galaxy %s" % (args.maxrows,columns,filter))
        nrows = 0
        for row in conn:
            # Filter out any SQL NULLs
            for col in clist:
                if row[col] is None:
                    if col not in nulls:
                        nulls[col] = 0
                    nulls[col] += 1
                    row[col] = args.null_sub
            # Dump this row to our output file
            print >>f, ' '.join([str(row[col]) for col in clist])
            nrows += 1
        if args.verbose:
            print 'Dumped',nrows,'rows to',args.output
            if nulls:
                print 'Replaced NULLs with',args.null_sub,'for:'
                for col in nulls:
                    print '%10d %s' % (nulls[col],col)
    except _mssql.MssqlDatabaseException,e:
        print 'Database Exception'
        raise
    finally:
        if conn: conn.close()

    f.close()

if __name__ == "__main__":
    main()
