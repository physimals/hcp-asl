from setuptools import setup

setup(
    name='hcp-asl', # is hyphen ok here and below or should they all be underscores?
    entry_points={
        'console_scripts': ['hcp_asl = hcp-asl.scripts.pipeline:main']
    }
    python_requires=['>=3.4'] # for pathlib.Path
    setup_requires=[
        'numpy',
        'fslpy',
        'pyfab'
    ]
)