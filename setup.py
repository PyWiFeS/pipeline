from setuptools import find_packages, setup

setup(
    name="pywifes",
    version="0.0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "wheel",
        "setuptools",
        "astropy",
        "scipy",
        "numpy>=1.23",
        "matplotlib",
    ],
    description="A Python package for optical data reduction pipeline.",
    author="Timothy Davies",
    author_email="tim.davies@uwa.edu.au",
)
