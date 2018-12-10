#!/usr/bin/env python
# $Id: 1d.py  hjm $
# by harry mangalam, uc irvine. 2008

"""
This Program is a wrapper for AJ Shaka's & Vladimir Mandelshtam's
fd_rr1d.f 1D NMR analysis program originally writ in Fortran,
trivially converted to a subroutine and thence converted, along with
its subsidiary Fortran subroutines via f2py (part of the numpy pkg).
try:
    "./1d.py --help"     for how to use it
    "./1d.py --help1d"   for the information included wth the original Fortran code

Explanation of how it was wrapped in Python and the features added, see:
<http:moo.nac.uci.edu/~hjm/fd_rrt1d/index.html>
"""

from __future__ import division # true floating point div is nice
import sys
from types import *
from pydoc import pipepager
from math import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from GUI_1D import *
import sys, os, commands, md5, re, fnmatch, getopt, errno, psycopg2
# this is the magic resulting from f2py ...?
from fd_rrt1d import *
from configobj import ConfigObj# for the configuration module
import time # for the hi-res timer..
import socket
from IPython.Shell import IPShellEmbed # for superior debugging support, iPython is unbeatable
dbg = IPShellEmbed(['']) # and insert 'dbg()' where you want to break.


def usage(code, msg=''):
    #read in the help file
    try:
        help_fp = file("1d_help.txt", "r")
        help_txt = help_fp.read()
    except:
        print "Can't find the help file - should be called '1d_help.txt' - Did you rename it?"
        sys.exit(code)
    pipepager(help_txt, '/usr/bin/less -NS') # pipe help text into 'less -NS'
    sys.exit(code)


# this should definitely be externalized into a file and read in as needed.
def usage1d(code):
    try:
        help_fp = file("1d_orig_help.txt", "r")
        help_txt = help_fp.read()
    except:
        print "Can't find the help file - should be called '1d_orig_help.txt' - Did you rename it?"
        sys.exit(code)
    pipepager(help_txt, '/usr/bin/less -NS') # pipe help text into 'less -NS'
    sys.exit(code)

def ppdict(dct):
    max_wid = 60
    nbr_spcs = 2
    max_nme_len = 20
    dotchar = '.'
    for key, val in dct.iteritems():
        # calc len of (left) key
        key_nme_len = len(key)
        if key_nme_len < max_nme_len:
            rstr = key
        else:
            rstr = key[max_nme_len]
            key_nme_len = max_nme_len
        # calc len of (right) val - tricky due to dif var types.
        #assume only 3 types: IntType, FloatType, StringType (default)
        if type(val) == IntType: val = "%d" % (val)
        elif  type(val) == FloatType: val = "%10.5e" % (val)
        val_nme_len = len(val)
        if val_nme_len > max_nme_len:
            val_nme_len = max_nme_len
        rstr += " " * nbr_spcs  # add left spaces
        # add dots
        nbr_dots = max_wid - (key_nme_len + val_nme_len + (nbr_spcs*2))
        rstr += dotchar * nbr_dots
        rstr += " " * nbr_spcs  # add right spaces
        rstr += val  # add (right) name
        print rstr  # and print it


# the following class is a dummy which should be removed
# before it's released into the wild.  It does nothing.
class Form(QDialog):

    def __init__(self, parent=None):
        super(Form, self).__init__(parent)
        self.browser = QTextBrowser()
        self.lineedit = QLineEdit("Type an expression and press Enter")
        self.lineedit.selectAll()
        layout = QVBoxLayout()
        layout.addWidget(self.browser)
        layout.addWidget(self.lineedit)
        self.setLayout(layout)
        self.lineedit.setFocus()
        self.connect(self.lineedit, SIGNAL("returnPressed()"),
                     self.updateUi)
        self.setWindowTitle("Calculate")


    def updateUi(self):
        try:
            text = unicode(self.lineedit.text())
            self.browser.append("%s = <b>%s</b>" % (text, eval(text)))
        except:
            self.browser.append(
                    "<font color=red>%s is invalid!</font>" % text)

# the class multiply inherits from the library prototype and the specific interface
# class defined in the designer ui -> py
class Form1D(QDialog,Ui_Dialog):

    # to pop it up, it only needs to __init__ itself, declare its parent (itself) as
    # it's a top-level dialog, and then call the designer -> pyuic4-generated setupUi()
    # to make it do anything useful, I have to write all the glue code to pass
    # the params, connect signals & slots, do error-checking, etc. but his pops it up
    # don't forget to erase the no-longer needed class and defs when finished.

    def __init__(self, parent=None):
        super(Form1D, self).__init__(parent)
        self.setupUi(self)

def gui():
    """
    to pop up the designer-built form, it only needs to declare an instance of the
    QtApplication, ditto the form itself, and show it
    To make it do anything useful, still have to write all the glue code to pass
    the params, connect signals & slots, do error-checking, etc. but this pops it up
    don't forget to erase the no-longer needed class and defs when finished.
    """

    app = QApplication(sys.argv)
    form = Form1D()
    form.show()
    app.exec_()



# Typedef and initialize ALL vars as parts of a dict
# the named vars are separate from the params for fd_rrt1d
help      = 0 # no help required
help1d    = 0 # no extended help required
paramfile =""  # no default - if null, then have to read the params
DEBUG     = 0  # 1 = debug mode.
USEDB     = 1;
# note: no ','s tokens in the connection string below; whitespace only!
# following dsn is for postgresql running on the same machine - your user and pass will be different of course.
dsn = "dbname='shaka1d' port='5432' user='hjm' host='localhost' password='gulahula'"
# following dsn is for postgresql running on a remote  machine
#dsn = "dbname='shaka1d' port='5433' user='hjm' host='128.200.34.95' password='flibberty'"
# Following are the system defaults for the fd_rrt1d parameters.
# if any values are not defined by the commandline or file, these are the last resort.

clcfg={}  # the dict for the params from the cmd line
fcfg={} # the dict for the params from the specified file
cfg={}    # the dict that is finally passed to fd_rrt1d

# the dicts should contain ONLY those params that are specific to the fortan code
# other params or variables should be stored in other variables.
cfg['signal']    = "2p-no-noise.txt"          # signal file
cfg['theta']     = 1.5            # zero, order phase correction
cfg['nmr_method']    = "FDM"           # FDM, RRT or DFT
cfg['Nsig']      = 40960             # NONdefault
cfg['wmino']     = -9000   # window size min value
cfg['wmaxo']     =  4500     # window size max value
cfg['par']       = "FDM_par.out"  # linelist output file
cfg['threshhold']= 0.0001  # threshhold for par above
cfg['ReSp']      = "ReSp_spectra.data"   # spectra output files
cfg['ImSp']      = "ImSp_spectra.data"   # spectra output files
cfg['AbsSp']     = "AbsSp_spectra.data"  # spectra output files
cfg['rho']       = 2.0       # basis density
cfg['Nb0']       = 100    # basis size
cfg['Nbc']       = -20     # coarse basis
cfg['Nsp']       = 20000              # plotting pnts
cfg['Gamm']      = 5e-02              # smoothing (5d-2)
cfg['cheat']     = 1.0                # used in FDM, see discriptions.
cfg['cheatmore'] = "F"                # used in FDM, see discriptions.
cfg['ros']       = 1e-08              # RRT regularization parameter (1d-8)

try:
    opts, args = getopt.getopt(sys.argv[1:], 'hD', ['help', 'debug', 'help1d', 'gui', 'nodb', 'paramfile=', 'signal=', 'theta=', 'method=', 'Nsig=', 'wmino=', 'wmaxo=', 'par=', 'threshhold=', 'ReSp=', 'ImSp=', 'AbsSp=', 'rho=', 'Nb0=', 'Nbc=', 'Nsp=', 'Gamm=', 'cheat=', 'cheatmore=', 'ros='])
except getopt.GetoptError:
    # print help information and exit:
    print "There was an error specifying an option flag.  Here's the correct usage:"
    usage(1)

#dbg()

# set up the options required
for opt, arg in opts:
    if opt in ('-h', '--help'):   usage(1)
    elif opt in ('-D', '--debug'):      DEBUG = 1
    elif opt in ('--help1d'):     usage1d(1)
    elif opt in ('--gui'):        gui()
    elif opt in ('--nodb'):       USEDB = 0
    elif opt in ('--paramfile'):  paramfile = arg
    elif opt in ('--signal'):     clcfg['signal'] = arg         # file name
    elif opt in ('--theta'):      clcfg['theta'] = float(arg)   # float
    elif opt in ('--method'):     clcfg['nmr_method'] = arg         # FDM, RRT or DFT
    elif opt in ('--Nsig'):       clcfg['Nsig'] = int(round(float(arg)))      #int
    elif opt in ('--wmino'):      clcfg['wmino'] = int(round(float(arg)))   #int
    elif opt in ('--wmaxo'):      clcfg['wmaxo'] = int(round(float(arg)))   #int
    elif opt in ('--par'):        clcfg['par'] =  arg           # linelist output file
    elif opt in ('--threshhold'): clcfg['threshhold'] = float(arg)
    elif opt in ('--ReSp'):       clcfg['ReSp'] = arg           # file name
    elif opt in ('--ImSp'):       clcfg['ImSp'] = arg           # file name
    elif opt in ('--AbsSp'):      clcfg['AbsSp'] = arg          # file name
    elif opt in ('--rho'):        clcfg['rho'] = float(arg)   #
    elif opt in ('--Nb0'):        clcfg['Nb0'] = int(round(float(arg)))   #
    elif opt in ('--Nbc'):        clcfg['Nbc'] = int(round(float(arg)))   #
    elif opt in ('--Nsp'):        clcfg['Nsp'] = int(round(float(arg)))   #
    elif opt in ('--Gamm'):       clcfg['Gamm'] = float(arg)   #
    elif opt in ('--cheat'):
        clcfg['cheat'] = float(arg)   # float
        if clcfg['cheat'] != 1 or clcfg['cheat'] != 0:
            print >> sys.stderr, "cheat must be '1' or '0'"
            sys.exit(1)
    elif opt in ('--cheatmore'):
        clcfg['cheatmore'] = arg   # T or F
        if clcfg['cheatmore'] != 'T' or clcfg['cheatmore'] != 'F':
            print >> sys.stderr, "cheatmore must be 'T' or 'F'"
            sys.exit(1)
    elif opt in ('--ros'):        clcfg['ros'] = float(arg)   #

#dbg()
# the precedence of parameter input is cmdline > file > default, so if any values are entered from
# the cmdline, they supercede the ones from any other input.
# if there's a param file specified as well (--paramfile="/path/to/parameter/file")
# then the cmdline params supercede the ones from the file.
# if any params are missing, then we finally go to the default dict dcfg{} to get them


if paramfile != "": # If there's a param file named, try to get params from it.
    fcfg = ConfigObj(file(paramfile))

    # have to coerce everything from the param file that is not a
    # string to the correct type
    # following members used to iterate over to coerce into int or float
    int_params = ("Nsig", "Nsp", "wmino", "wmaxo", "Nb0", "Nbc")
    float_params = ("cheat", "theta", "threshhold", "rho", "Gamm", "ros")

    for ip in int_params: fcfg[ip]=int(fcfg[ip])
    for fp in float_params: fcfg[fp]=float(fcfg[fp])

    if DEBUG:
        print "\nParameters from the config FILE are:"
        for k, v in fcfg.iteritems():
            print k, '\t', v

# we take most values from the file, but if there were also values from the cl, the cl values
# take precendence, so overwrite the fgfg values with them. (populate the fcfg dict)
for k, v in cfg.iteritems(): # nb: iterating thru the k,v from the /defaults/ dict
    # load cfg with parameters from the file.
    # print "default: ", k, '\t',  v

    #print "Updating default keys from file [", paramfile, "]"
    if k in fcfg: # test if there are any params from named file
        if fcfg[k] != "":
            print "default Key ", k, " -> file value", fcfg[k]
            cfg[k] = fcfg[k] # if so, they take prec over defaults stored in cfg{}

    #print "Updating default keys from commandline parameters"
    if k in clcfg: # test if there are any params from cmdline
        if clcfg[k] != "":
            print "default Key ", k, " -> cmdline value ", clcfg[k]
            cfg[k] = clcfg[k] # if so, they take prec over the rest


# and finally print the params that ARE going to be passed to fd_rrt1d()
print "\nFINAL Parameters passed to fd_rrt1d():"
#for k, v in cfg.iteritems():
#    print k, '=', v
ppdict(cfg)

# and now punt all these to fd_rrt1d
if DEBUG:
    tt = raw_input("hit Enter to call fd_rrt1d() with the params above..")

time_start = time.time()  #Start the timer
# and run the fortran code
fd_rrt1d(cfg['signal'], cfg['theta'], cfg['nmr_method'], cfg['Nsig'], cfg['wmino'], cfg['wmaxo'], cfg['par'], cfg['threshhold'], cfg['ReSp'], cfg['ImSp'], cfg['AbsSp'], cfg['rho'], cfg['Nb0'], cfg['Nbc'], cfg['Nsp'], cfg['Gamm'], cfg['cheat'], cfg['cheatmore'], cfg['ros'])

time_end = time.time()
runtime = time_end - time_start

if DEBUG: print "Wall clock time taken(s) = ", runtime

# now write all the info we need to the database
# using the dsn "DB_name,DB_Server,port#,DB_user,DB_user_password"
# but 1st, make sure that all the params are available:
date = time.asctime()
#username = os.getenv("USER")
opsys,hostname,j,jj,CPU = os.uname()
ipnbr = socket.gethostbyname(hostname)
# get the OS-specific hardware info
print "Following Warning messages are from the Hardware profiler:"
if opsys == "Linux":    f= os.popen("lshw -short")
elif opsys == "Darwin": f= os.popen("system_profiler -detailLevel mini SPHardwareDataType")
# read the entire text blob into 'sysinfo'.  This should be trimmed to include only
# lines that have the following terms: memory, processor (1st line), display, eth0
sysinfo = f.read()


# all the rest is contained in the cfg{} hash

if DEBUG:
    print "date:\t", date
    print "host:\t", hostname
    print "ipnbr:\t", ipnbr
    print "OS:\t", opsys
    print "sysinfo:\t", sysinfo
    print "runtime:\t", runtime
    for k,v in cfg.iteritems():
        print k, ":\t", v

response = raw_input("Do you want to send this data back to the author? [y,N]")
if response.lower() == "y": USEDB=1
else: USEDB=0

if USEDB:
    try:
        conn = psycopg2.connect(dsn)
        print "Success: DB connection made!"
        #compose the insert string
        # apologies - this is fairly awful to debug.  There's probably a better way.
        insert_str = "INSERT INTO all_data (date, hostname, ipnbr, os, sysinfo, runtime, signal, theta, nmr_methd, Nsig, wmino, wmaxo, par, threshhold, ReSp, ImSp, AbsSp, rho, Nb0, Nbc, Nsp, Gamm, cheat, cheatmore, ros) VALUES ('%s',  '%s', '%s', '%s', '%s', %f, '%s', %f, '%s', %d, %d, %d, '%s', %f, '%s', '%s', '%s', %f, %d, %d, %d, %f, %f, '%c', %f) RETURNING all_data.id;" % (date, hostname, ipnbr, opsys, sysinfo, runtime, cfg['signal'], cfg['theta'], cfg['nmr_method'], cfg['Nsig'], cfg['wmino'], cfg['wmaxo'], cfg['par'], cfg['threshhold'], cfg['ReSp'], cfg['ImSp'], cfg['AbsSp'], cfg['rho'], cfg['Nb0'], cfg['Nbc'], cfg['Nsp'], cfg['Gamm'], cfg['cheat'], cfg['cheatmore'], cfg['ros'])

        if DEBUG:
            print "\n\nSending this information to author database: ", insert_str
        all_data_id = 0 # initialize var


        try:
            curs = conn.cursor()
            curs.execute(insert_str)
            conn.commit()
            all_data_id = curs.fetchone()[0]
            print "inserted data @ index = ",  all_data_id
        except: # insertion fails for some reason so report it
            print "INSERT FAILURE at curs.execute"
    except:
        print "ERROR: DB connection NOT made!"
        sys.exit(1)
#    curs = conn.cursor()


else:
    if DEBUG: print "...skipping DB connection..."

print "Normal exit.  Thanks for using fd_rr1d."
sys.exit(1)
