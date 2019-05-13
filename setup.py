from setuptools import setup, find_packages
import glob

setup(
        name='autovps',
        version="0.0.1",
        packages=find_packages(),
        url='https://github.com/rhancockn/autovps.git',
        license='MIT',
        author='Roeland Hancock',
        author_email='rhancock@gmail.com',
        description='Preprocessing utilities for CMRR SVS data',
        classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
            'Topic :: Scientific/Engineering :: Physics',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: MIT License',

            # Specify the Python versions you support here.
            # Indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',

        ],
        install_requires=['numpy', 'pydicom', 'nibabel'],
        test_requires=['pytest'],
        python_requires= '>=3',
        scripts=glob.glob('bin/*.py'),
        entry_points={
          'console_scripts': [
              'process_cmrr.py = autovps.bin.process_cmrr:main',
              'tarquin_to_dict.py = autovps.bin.tarquin_to_dict:main'
          ]
      },
)
