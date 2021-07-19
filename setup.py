from distutils.core import setup
from Cython.Build import cythonize

exec(open("chemfold/version.py").read())

setup(
    version     = __version__,
    name        = "chemfold",
    description = "Train-validation-test splits for molecular data, including federated ML setups",
    long_description = "Implementations of sphere exclusion clustering, scaffold trees, locality sensitive hashing",
    packages    = ["chemfold"],
    ext_modules = cythonize("./chemfold/sphere.pyx"),
    setup_requires = ['numpy>=1.0', 'scipy', 'cython>=0.24.1'],
    install_requires = ["numpy>=1.0", "scipy", "cython>=0.24.1", "statsmodels"],
)
