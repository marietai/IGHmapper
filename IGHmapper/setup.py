from setuptools import setup, find_packages

setup(
    name="IGHmapper",
    version="0.1",
    packages=find_packages(include=["IGHmapper", "IGHmapper.*"]),
    install_requires=[
        "pandas",
    ],
    package_data={
        "IGHmapper.igh_converter": ["igh_to_hg19.txt", "igh_to_hg38.txt"]
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "igh_convert=IGHmapper.igh_converter.conversion:convert_to_custom_reference",
            "igh_validate=IGHmapper.igh_converter.validation:validate_mapping",
            "igh_identify=IGHmapper.igh_converter.mapping:identify_build",
            "igh_extract=IGHmapper.igh_converter.mapping:extract_region",
        ],
    },
    description="A Python package for converting and validating SNPs using custom IGH reference files.",
    author="marietai",
    url="https://github.com/marietai/IGH",
    license="MIT",
    python_requires=">=3.6",
)
