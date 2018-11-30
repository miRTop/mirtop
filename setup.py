"""small RNA-seq annotation"""

import os
from setuptools import setup, find_packages

version = '0.4.14a'
url = 'http://github.com/mirtop/mirtop'


def readme():
    with open('README.md') as f:
        return f.read()


def write_version_py():
    version_py = os.path.join(os.path.dirname(__file__), 'mirtop',
                              'version.py')
    with open(version_py, "w") as out_handle:
        out_handle.write("\n".join(['__version__ = "%s"' % version,
                                    '__url__ = "%s"' % url]))


write_version_py()

setup(name='mirtop',
      version=version,
      description='Small RNA-seq annotation',
      long_description=readme(),
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        "Programming Language :: Python :: 3",
        'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      keywords='RNA-seq miRNA isomiRs annotation',
      url=url,
      author='Lorena Pantano',
      author_email='lorena.pantano@gmail.com',
      license='MIT',
      packages=find_packages(),
      test_suite='nose',
      entry_points={
          'console_scripts': ['mirtop=mirtop.command_line:main'],
      },
      include_package_data=True,
      zip_safe=False)
