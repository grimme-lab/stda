# This file is part of stda.
#
# Copyright (C) 2019 Sebastian Ehlert
#
# xtb4stda is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb4stda is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb4stda.  If not, see <https://www.gnu.org/licenses/>.

project('stda', 'fortran',
        version: '1.6.1',
        license: 'LGPL3',
        meson_version: '>=0.51')

fc = meson.get_compiler('fortran')

if fc.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  add_project_arguments('-fbacktrace', language: 'fortran')
elif fc.get_id() == 'intel'
  add_project_arguments('-axAVX2',    language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  if get_option('static')
    add_project_link_arguments('-static', language: 'fortran')
  endif
endif

dependencies = []

la_backend = get_option('la_backend')
if la_backend == 'mkl'
  libmkl = [fc.find_library('pthread')]
  libmkl += fc.find_library('m')
  libmkl += fc.find_library('dl')
  if fc.get_id() == 'intel'
    libmkl += [fc.find_library('mkl_intel_lp64')]
    libmkl += fc.find_library('mkl_intel_thread')
  else
    libmkl += [fc.find_library('mkl_gf_lp64')]
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


