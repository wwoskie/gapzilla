from setuptools import setup, find_packages

setup(
    name="gapzilla",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pyyaml",
        "tqdm",
        "viennarna",
    ],
    entry_points={
        "console_scripts": [
            "gapzilla=gapzilla.main:main",
        ],
    },
)
