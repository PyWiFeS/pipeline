from setuptools import find_packages, setup

setup(
    name="pywifes",
    version="1.0.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires="==3.10.*",
    install_requires=[
        "astropy",
        "matplotlib",
        "numpy<2.0",
        "pyjson5",
        "pandas",
        "photutils==1.8.0",
        "scipy==1.9.1",
        "setuptools",
        "wheel",
    ],
    description="A Python package for optical data reduction pipeline.",
)
