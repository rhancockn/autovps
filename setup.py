from setuptools import setup, find_packages

setup(
        name='autovps',
        version="0.0.1",
        packages=find_packages(),
        url='https://github.com/rhancockn/autovps.git',
        license='MIT',
        author='Roeland Hancock',
        author_email='rhancock@gmail.com',
        description='',
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

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',

        ],
        install_requires=['numpy', 'pydicom', 'nibabel'],
        test_requires=['pytest']
)