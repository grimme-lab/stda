# This file is part of stda.
#
# Copyright (C) 2019 Sebastian Ehlert
#
# stda is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# stda is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb4stda.  If not, see <https://www.gnu.org/licenses/>.

# modified to include options for Windows with Intel Fortran only (Shoubhik Maiti, April 2022)

project('stda', 'fortran',
        version: '1.6.1',
        license: 'LGPL3',
        meson_version: '>=0.51')

if get_option('interface') == '64' and get_option('la_backend') != 'mkl'
  error('64 bit integer interface only supported for MKL backend')
endif

fc = meson.get_compiler('fortran')

if fc.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  add_project_arguments('-fbacktrace', language: 'fortran')
  if get_option('interface') == '64'
    add_project_arguments('-fdefault-integer-8', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments('-axAVX2',    language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  if get_option('interface') == '64'
    add_project_arguments('-i8', language: 'fortran')
  endif
  if get_option('static')
    add_project_link_arguments('-static', language: 'fortran')
  endif
elif fc.get_id() == 'intel-cl' # for Windows with Intel Fortran
  add_project_arguments('-QaxCORE-AVX2',    language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  if get_option('interface') == '64'
    add_project_arguments('-integer-size:64', language: 'fortran')
  endif
  if get_option('static')
    add_project_link_arguments('-MT', language: 'fortran')
  endif
endif

dependencies = []

la_backend = get_option('la_backend')
if la_backend == 'mkl'
  if host_machine.system() != 'windows'
    libmkl = [fc.find_library('pthread')]
    libmkl += fc.find_library('m')
    libmkl += fc.find_library('dl')
  endif
  if (fc.get_id() == 'intel') or (fc.get_id() == 'intel-cl')
    if get_option('interface') == '64'
      libmkl += fc.find_library('mkl_intel_ilp64')
    else
      libmkl += fc.find_library('mkl_intel_lp64')
    endif
    libmkl += fc.find_library('mkl_intel_thread')
  else
    if get_option('interface') == '64'
      libmkl += fc.find_library('mkl_gf_ilp64')
    else
      libmkl += fc.find_library('mkl_gf_lp64')
    endif
    libmkl += fc.find_library('mkl_gnu_thread')
  endif
  libmkl += fc.find_library('mkl_core')
  libmkl += fc.find_library('iomp5')
  dependencies += libmkl
elif la_backend == 'openblas'
  dependencies += fc.find_library('openblas', required : true)
  dependencies += fc.find_library('lapack', required : true)
elif la_backend == 'custom'
  foreach lib: get_option('custom_libraries')
    dependencies += fc.find_library(lib)
  endforeach
else
  dependencies += fc.find_library('blas', required : true)
  dependencies += fc.find_library('lapack', required : true)
endif

if get_option('openmp')
  if fc.get_id() == 'intel'
    add_project_arguments('-qopenmp', language : 'fortran')
    add_project_link_arguments('-qopenmp', language : 'fortran')
  elif fc.get_id() == 'intel-cl'
    add_project_arguments('-Qopenmp', language : 'fortran')
    add_project_link_arguments('-Qopenmp', language : 'fortran')
  else
    add_project_arguments('-fopenmp', language : 'fortran')
    add_project_link_arguments('-fopenmp', language : 'fortran')
  endif
endif

stda_srcs = [
    'apbtrafo.f',
    'block.f',
    'header.f',
    'intpack.f90',
    'intslvm.f',
    'io.f',
    'linal.f',
    'linear_response.f',
    'main.f',
    'molden.f',
    'normalize.f',
    'onetri.f',
    'pckao.f',
    'print_nto.f',
    'printvec.f',
    'prmat.f',
    'readbasa.f',
    'readbasmold.f',
    'readl.f',
    'readxtb.f',
    'sfstda.f',
    'sosor.f',
    'srpapack.f',
    'stdacommon.f90',
    'stda.f',
    'stda-rw.f',
    'stringmod.f90',
    'sutda.f',
    'velo.f'
]

g_spec_srcs = [
    'g_spec/g_spec.f'
]

g2molden_srcs = [
    'g2molden/main.f',
    'g2molden/stringmod.f90'
]

stda_exe = executable(meson.project_name(), stda_srcs,
                      dependencies: dependencies,
                      install: true)

g2molden_exe = executable('g2molden', g2molden_srcs,
                          dependencies: dependencies,
                          install: true)


g_spec_exe = executable('g_spec', g_spec_srcs,
                        dependencies: dependencies,
                        install: true)


