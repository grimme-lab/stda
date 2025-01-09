# This file is part of std2.
#
# Copyright (C) 2019 Sebastian Ehlert
# Copyright (C) 2024 Marc de Wergifosse

# Modified by P. Beaujean

project('std2', 'fortran',
  version: '1.6.1',
  meson_version: '>=0.51',
  default_options : ['warning_level=0', 'fortran_std=legacy'],
)

if get_option('interface') == '64' and (get_option('la_backend') == 'netlib' or get_option('la_backend') == 'openblas')
  error('64 bit integer interface not supported by OpenBLAS/netlib backends, use -Dinterface=32')
endif

# set compiler options
fc = meson.get_compiler('fortran')

build_args = [
  '-DPROJECT_NAME="' + meson.project_name() + '"',
  '-DPROJECT_VERSION="' + meson.project_version() + '"',
]

if fc.has_argument('-march=native')
    build_args += ['-march=native', '-fno-math-errno']
    if fc.get_id() == 'intel'
        build_args += ['-unroll-aggressive', '-ipo']
    else
        build_args += ['-funroll-loops', '-ftree-vectorize']
    endif
    message('added -march=native and loop unroling flags')
endif

# configure project_dep
project_dep = []
libcint_options = ['with_fortran=true', 'with_cint2_interface=true']

la_backend = get_option('la_backend')
if la_backend == 'mkl'
  # build a pkg-config compatible string, list all possibilities with `pkg-config --list-all | grep "mkl"`
  mkl_kind = 'mkl'
  if get_option('static')
    mkl_kind += '-static'
    build_args += '-static'
  else
    mkl_kind += '-dynamic'
  endif

  if get_option('interface') == '64'
    mkl_kind += '-ilp64'
    # specific requirements for ilp64 libraries, see https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-windows/2023-0/using-the-ilp64-interface-vs-lp64-interface.html
    build_args += '-DMKL_ILP64'
    libcint_options += 'i8=true'
    if fc.get_id() == 'intel' or fc.get_id() == 'intel-llvm'
      build_args += '-i8'
    else
      build_args += '-fdefault-integer-8'
    endif
  else
    mkl_kind += '-lp64'
  endif

  if get_option('openmp')
    mkl_kind += '-iomp'
  else
    mkl_kind += '-seq'
  endif

  message('MKL kind: ' +  mkl_kind)

  # backups for system that do not provide pkg-config
  mkl_libraries = {
    'mkl-dynamic-ilp64-iomp': ['mkl_core', 'mkl_intel_ilp64', 'mkl_gnu_thread', 'iomp5', 'pthread', 'm', 'dl'],
    'mkl-dynamic-ilp64-seq': ['mkl_core', 'mkl_intel_ilp64', 'mkl_sequential', 'm', 'dl'],
    'mkl-dynamic-lp64-iomp': ['mkl_core', 'mkl_intel_lp64', 'mkl_gnu_thread', 'iomp5', 'pthread', 'm', 'dl'],
    'mkl-dynamic-lp64-seq': ['mkl_core', 'mkl_intel_lp64', 'mkl_sequential', 'm', 'dl']
    # ... some of them are missing at the moment
  }

  mkl_dep = dependency(mkl_kind, required: false)
  if mkl_dep.found() # pkg-config
    project_dep += mkl_dep
  else # back up to finding libraries one per one
    foreach lib: mkl_libraries[mkl_kind]
      project_dep += fc.find_library(lib)
    endforeach
  endif
elif la_backend == 'openblas'
  project_dep += fc.find_library('openblas', required : true)
  project_dep += fc.find_library('lapack', required : true)
elif la_backend == 'netlib'
  project_dep += fc.find_library('blas', required : true)
  project_dep += fc.find_library('lapack', required : true)
elif la_backend == 'custom'
  foreach lib: get_option('custom_libraries')
    project_dep += fc.find_library(lib)
  endforeach
endif

# add openMP if requested but MKL (and thus intel-openMP) was not used
if la_backend != 'mkl' and get_option('openmp')
  dep_openmp = dependency('openmp', required: false)
  if dep_openmp.found()
    project_dep += dep_openmp
  else
    if fc.get_id() == 'intel' or fc.get_id() == 'intel-llvm'
      build_args += '-qopenmp'
    else
      build_args += '-fopenmp'
    endif
  endif
endif

libcint_dep = fc.find_library('libcint', required: false)
if not libcint_dep.found()
  libcint_proj = subproject('libcint', default_options: libcint_options)
  libcint_dep = libcint_proj.get_variable('libcint_dep')
endif
project_dep += libcint_dep

message('Build args are: ' + ', '.join(build_args))

# Sources
std2_srcs = [
    '2PA.f90',
    'apbtrafo.f',
    'block.f',
    'full.f',
    'header.f',
    'intpack.f90',
    'intslvm.f',
    'io.f',
    'libcint.f',
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
    'stda.f',
    'stdacommon.f90',
    'stda-rw.f',
    'stda-rw_dual.f',
    'stringmod.f90',
    'sutda.f',
    'velo.f',
    'xstd.f90',
]

g_spec_srcs = [
    'g_spec/g_spec.f'
]

g2molden_srcs = [
    'g2molden/main.f',
    'g2molden/stringmod.f90'
]

# Executables
std2_exe = executable(
  meson.project_name(),
  std2_srcs,
  dependencies: project_dep,
  fortran_args : build_args,
  link_language : 'fortran',
  install: true
)

g2molden_exe = executable(
  'g2molden',
  g2molden_srcs,
  dependencies: project_dep,
  fortran_args : build_args,
  link_language : 'fortran',
  install: true
)

g_spec_exe = executable(
  'g_spec', g_spec_srcs,
  dependencies: project_dep,
  fortran_args : build_args,
  link_language : 'fortran',
  install: true
)

# add test
test('valid input', std2_exe, args: ['-f', '../tests/water_sto3g.molden', '-sty', '3', '-ax', '1', '-e', '20'], suite: 'app')
