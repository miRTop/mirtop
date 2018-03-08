"""small RNA-seq annotation"""

from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

# with open("reqs.txt", "r") as f:
#     install_requires = [x.strip() for x in f.readlines() if not x.startswith("#")]


setup(name='mirtop',
      version='0.1.11a',
      description='Small RNA-seq annotation',
      long_description=readme(),
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      keywords='RNA-seq miRNA isomiRs annotation',
      url='http://github.com/mirtop/mirtop',
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
