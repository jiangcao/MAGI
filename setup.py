import os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class f2py_Extension(Extension):

    def __init__(self, name, sourcedirs):
        Extension.__init__(self, name, sources=[],f2py_options=['--quiet'])
        self.sourcedirs = [os.path.abspath(sourcedir) for sourcedir in sourcedirs]
        self.dirs = sourcedirs

class f2py_Build(build_ext):

    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        # compile
        F90FLAGS='"-march=native -ffree-line-length-none -fbounds-check -fbacktrace -ffast-math -fopenmp -fexternal-blas"'
        MKLROOT="/usr/pack/intel_compiler-2020-af/x64/mkl/"
        F2PYFLAGS="-c --fcompiler=gnu95 -lgomp -llapack -lblas --f90flags=%s" % (F90FLAGS)
        outdir='bin/'        
        
        for ind,to_compile in enumerate(ext.sourcedirs):
            module_loc = os.path.split(ext.dirs[ind])[0]
            module_name = os.path.split(to_compile)[1].split('.')[0]
            os.system('rm %s/*.mod %s/*.o ' % (module_loc,module_loc))
            os.system('cd %s;f2py %s %s/include/mkl_dfti.f90 -m %s %s skip: trimul_c mul_c : ' % (module_loc,F2PYFLAGS,MKLROOT,module_name,to_compile))
            os.system('mv %s/%s*.so %s' % (module_loc,module_name,outdir))

setup(
    name="magi",
    ext_modules=[f2py_Extension('fortran_external',['magi/core/interface/util.f95',
                                                    'magi/core/interface/wannier.f95',
                                                    'magi/core/interface/negf.f95'])],
    cmdclass=dict(build_ext=f2py_Build),
)