from distutils.core import setup  
from Cython.Build import cythonize  

setup(name='vcf_simplify',
    ext_modules= cythonize(['metadata_parser/*.py',
    'records_parser/*.py',
    'records_parser/simplifyvcf/*.py',
    'records_parser/buildvcf/*.py']),
)

# python3 setup.py build_ext --inplace