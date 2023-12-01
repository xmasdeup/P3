project(
    'Práctica 3 de PAV - detección de pitch', 'cpp',
    default_options: 'cpp_std=c++11',
    version: 'v2.0'
)

inc = include_directories(['src/include', 'src/pav/docopt_cpp','src/include/ffft',])

# Librería PAV
subdir('src/pav')

# Programas del proyecto
subdir('src/get_pitch')

# Scripts de la práctica
subdir('scripts')

# Documentación del proyecto
subdir('src/doxyfile')
