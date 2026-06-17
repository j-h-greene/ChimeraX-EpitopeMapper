from setuptools import setup, find_packages

setup(
    name='chimerax-epitope_mapper',
    version='1.0',
    description='A ChimeraX bundle for epitope mapping',
    author='Sophia Greene',
    author_email='sophia.greene@case.edu',
    url='https://github.com/j-h-greene/ChimeraX-EpitopeMapper',
    packages=find_packages(),
    package_data={
        'chimerax.epitope_mapper': ['bundle_info.xml'],
    },
    entry_points={
        'chimerax.bundle': ['epitope_mapper = chimerax.epitope_mapper:bundle_api']
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy',
        'scipy',
        'chimerax-core',
        'chimerax-map'
    ],
)