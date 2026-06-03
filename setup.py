from setuptools import setup

setup(name='ldsc',
      version='1.0',
      description='LD Score Regression (LDSC)',
      url='http://github.com/bulik/ldsc',
      author='Brendan Bulik-Sullivan and Hilary Finucane',
      author_email='',
      license='GPLv3',
      packages=['ldscore'],
      scripts=['ldsc.py', 'munge_sumstats.py'],
      python_requires='>=3.8',
      install_requires=[
            'bitarray>=2.6,<3',
            'nose>=1.3,<2',
            'pybedtools>=0.10,<1',
            'scipy>=1.10,<2',
            'numpy>=1.24,<2',
            'pandas>=1.5,<2'
      ]
)
