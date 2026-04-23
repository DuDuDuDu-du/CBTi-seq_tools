from setuptools import setup, find_packages

setup(
    name="CBTi-seq_tools",
    version="1.0.0",
    description="CBTi-seq complete analysis pipeline",
    author="Wenyi Zhang",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "pysam",
        "biopython",
        "pandas",
        "pyyaml"
    ],
    entry_points={
        "console_scripts": [
            "CBTi-seq_tools=cbtiseq_tools.cli:main"
        ]
    }
)
