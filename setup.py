from setuptools import setup

setup(
    name='hcpasl', # is hyphen ok here and below or should they all be underscores?
    entry_points={
        'console_scripts': ['hcp_asl = hcpasl.scripts.run_pipeline:main']
    }
    python_requires=['>=3.4'] # for pathlib.Path
    setup_requires=[
        'numpy',
        'fslpy',
        'pyfab'
    ]
)