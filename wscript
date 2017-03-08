# -*- mode: python; -*- 
# 

import os
import os.path as op
import sys
import Configure 
import Options
import commands

APPNAME = 'bao'
VERSION = '0.1.0'
top  = '.'
out  = 'build'
description = "BAO tools"
requirements = [('cfitsio', '3.280', True)]
debug="True"
optimize=3
omp="True"

@Configure.conftest
def check_packages(conf, pkg_list):
    """
    - Check whether the packages passed in argument 
      can be found by pkg-config. 
    - Parse the cflags and libs and store the information 
      in the uselib_store.
    """
    
    load_pkg_config(conf)

    # check for the packages passed in arguments 
    for pkg in pkg_list:
        try:
            pkg_name, pkg_version, pkg_mandatory = pkg
            conf.check_cfg(args='--cflags --libs', 
                           package = pkg_name + '-' + pkg_version, 
                           mandatory=pkg_mandatory, 
                           uselib_store=pkg_name.upper())
        except Configure.ConfigurationError:
            conf.fatal('unable to locate %s-%s (mandatory)' % (pkg_name, pkg_version))

def load_pkg_config(conf):
    """
    Load pkg-config (used to read the package metadata)
    """
    
    # first, we need pkg-config 
    if not conf.env['PKG_CONFIG'] or \
       not conf.env['PKG_CONFIG_PATH']:
        try:
            conf.find_program('pkg-config', var='PKG_CONFIG')
            pkgcpath = op.join( conf.env.PREFIX, 'lib', 'pkgconfig')
            if os.environ.has_key('PKG_CONFIG_PATH'):
                os.environ['PKG_CONFIG_PATH'] = pkgcpath + ":" + os.environ['PKG_CONFIG_PATH']
            else:
                os.environ['PKG_CONFIG_PATH'] = pkgcpath
        except Configure.ConfigurationError:
            conf.fatal('pkg-config not found')
    
    
def install_headers(bld, headers):
    """
    Just hide the install header commands 
    """
    name    = APPNAME
    version = VERSION
    install_dir = op.join('$PREFIX', 'include', '%s-%s' % (name, version))
    
    bld.install_files(install_dir, headers)






def options(opt):
    opt.load('compiler_cc')
    opt.load('compiler_cxx')
    opt.load('compiler_fc')
    
    
    opt.add_option('--with-cfitsio', 
                   action='store', 
                   help='Path to the cfitsio root dir')

    opt.add_option('--debug', 
                   action='store', 
                   default=debug, 
                   dest='debug', 
                   help='True or False, turn on/off the -g / -DNDEBUG option, (True by default)')

    opt.add_option('--omp', 
                   action='store', 
                   default=omp,
                   dest='omp', 
                   help='turn on open mp')
    
    opt.add_option('--optimize', 
                   action='store', 
                   default=optimize, 
                   dest='optimize', 
                   type=int,
                   help='0, 1, 2 or 3, default is %d, i.e. -O%d'%(optimize,optimize))

def configure(conf):
    # conf.load('frogs')    # c compiler 
    conf.load( 'compiler_c' )
    conf.env['CCFLAGS'] = ['-fPIC', '-DPIC']
    
    # c++ compiler 
    conf.load( 'compiler_cxx' )
    conf.env['CXXFLAGS'] = ['-fPIC', '-DPIC','-Wuninitialized','-Wunused-value','-Wunused-variable']
    conf.env['PKG_INCDIR'] = op.join('include', '%s-%s' % (APPNAME,VERSION))
    
    # fortran compilers 
    # apparently, we need to check for a c-compiler before.
    conf.load( 'compiler_fc')
    conf.check_fortran()
    conf.check_fortran_verbose_flag()
    conf.check_fortran_clib()
    conf.check_fortran_dummy_main()
    conf.check_fortran_mangling()
    conf.env['FCFLAGS'] = ['-fPIC']
    
    
    if conf.options.debug == "True" :
        print "with debug symbols"
        conf.env['CXXFLAGS'].append('-g')
    
        
    if conf.options.omp  == "True" :
        print "with openmp"
        conf.env['CXXFLAGS'].append('-fopenmp')
        conf.env['CXXFLAGS'].append('-DUSEOMP')
        conf.env['LINKFLAGS'].append('-fopenmp')
        conf.env['LINKFLAGS'].append('-lgomp')
        
    
    # Optimizer ? 
    if conf.options.optimize:
        print "optimizing"
        conf.env['CXXFLAGS'].append('-O%d'%conf.options.optimize)
        conf.env['LINKFLAGS'].append('-O%d'%conf.options.optimize)

    
    
    conf.check_cc( lib='lapack', msg='checking for lapack' )
    conf.check_packages(requirements)
    conf.env['CXXFLAGS'].append('-DUSELAPACK')
    conf.env['CXXFLAGS'].append('-DUSECFITSIO')
    conf.write_config_header('config.h')
        
def build(bld):
    bld.add_subdirs( ['src','apps'] )    
    
