#!/usr/bin/env python
#######################################################################################################
## Program to query the official LSST simulation galaxy catalog and write a simple text catalog
## file suitable for weak lensing studies.
##
## Uses the pymsql package for talking to the remote Microsoft SQL server, which in turns requires
## that Cython and FreeTDS are installed.
##
## Usage: lsst2wl.py <maxrows> <fname>
##
## Created 07-Nov-2012 by David Kirkby <dkirkby@uci.edu>
#######################################################################################################
import _mssql
import sys
import argparse

def main():

    if len(sys.argv) != 3:
        print 'usage: %s <maxrows> <filename>' % sys.argv[0]
        sys.exit(-1)
        
    try:
        # The maximum number of rows to retrieve (total number in the catalog is ~17M)
        maxrows = int(sys.argv[1])
        if maxrows < 0:
            maxrows = 100
    except ValueError:
        print 'Expected an integer for maxrows'
        sys.exit(-1)

    try:
        # The filename to write
        fname = sys.argv[2]
        f = open(fname,'w')
    except IOError,e:
        print 'Cannot write to filename "%s"' % fname
        print str(e)
        sys.exit(-2)

    # The numeric value to substitute for any SQL NULLs
    nullSub = -1

    # The ra,dec window to retrieve
    window = { 'RAmin':0, 'RAmax':1, 'DECmin':-1, 'DECmax':0 }

    def addColumns(patterns,types):
        text = ''
        for p in patterns:
            text += ',' + ','.join([p%t for t in types])
        return text

    # Specify the header columns we need
    columns = 'id,ra,dec,redshift,r_ab,absmag_r_total'

    # Add bulge and disk specific columns
    columns += addColumns(('%sHalfLightRadius',),('Bulge','Disk'))
    columns += addColumns(('pa_%s','magnorm_%s','sedid_%s','sedname_%s'),('bulge','disk'))
    columns += addColumns(('a_%s','b_%s'),('b','d'))

    # Add filter-specific columns
    columns += addColumns(('%s_ab',),"ugrizy")

    # Filename to write
    fname = 'gcat.dat'

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
        conn.execute_query("SELECT TOP %d %s FROM galaxy %s" % (maxrows,columns,filter))
        nrows = 0
        for row in conn:
            # Filter out any SQL NULLs
            for col in clist:
                if row[col] is None:
                    if col not in nulls:
                        nulls[col] = 0
                    nulls[col] += 1
                    row[col] = nullSub
            # Dump this row to our output file
            print >>f, ' '.join([str(row[col]) for col in clist])
            nrows += 1
        print 'Dumped',nrows,'rows to',fname
        if nulls:
            print 'Replaced NULLs with',nullSub,'for:'
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
