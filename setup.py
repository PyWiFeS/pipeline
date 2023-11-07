from setuptools import find_packages, setup

setup(
    name="pywifes",
    version="0.7.4",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "wheel",
        "setuptools",
        "astropy",
        "scipy",
        "numpy",
        "matplotlib",
    ],
    description="A Python package for optical data reduction pipeline.",
)
