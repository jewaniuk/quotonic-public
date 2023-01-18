import os

# Check if Cython is installed
try:
    from Cython.Build import cythonize

# If Cython is not installed, print message and do nothing
except ImportError:
    print("Cython is not correctly installed.")

    # This function is necessary, otherwise, Poetry will fail
    def build(setup_kwargs):
        pass


# If Cython is installed, then compile the extensions
else:
    from distutils.command.build_ext import build_ext

    from numpy import get_include
    from setuptools import Extension

    # This function will be executed in the setup.py created from the pyproject.toml
    def build(setup_kwargs):
        # Generate the extensions for each Cython module
        extensions = [Extension("quotonic.fock", ["quotonic/fock.pyx"]), Extension("quotonic.aa", ["quotonic/aa.pyx"])]

        # gcc arguments hack: enable optimizations
        os.environ["CFLAGS"] = "-O3"

        # Build
        setup_kwargs.update(
            {
                "ext_modules": cythonize(extensions, language_level=3, compiler_directives={"linetrace": True}),
                "include_dirs": [get_include()],
                "cmdclass": {"build_ext": build_ext},
            }
        )
