# IGHmapper

A Python package for converting and validating SNPs using custom IGH reference files.

---

## Installation

Clone the repository and install the package locally:
```bash
git clone https://github.com/marietai/IGH.git
cd IGHmapper
pip install -e .



## Usage

### Importing Modules
```python
from igh_converter import conversion, validation, mapping
```

### Example Usage
```python
from igh_converter import mapping

# Identify the build (GRCh37 or GRCh38)
from IGHmapper.igh_converter import identify_build

bim_file = r"path_to_your_bim_file.bim"
build = identify_build(bim_file)
print(f"Identified Build: {build}")


# Extract region chr14:105000000-107500000
from IGHmapper.igh_converter import extract_region

bim_file = r"path_to_your_bim_file"
output_prefix = r"path_to_your_output_folder\IGH_extracted"

extract_region(
    bim_file=bim_file,
    bed_file=f"{bim_file}.bed",
    fam_file=f"{bim_file}.fam",
    output_prefix=output_prefix
)

# Convert to Custom Reference
from IGHmapper.igh_converter import convert_to_custom_reference

bim_file = r"path_to_your_bim_file.bim"
output_bim_file = r"path_to_your_output_folder\converted.bim"

convert_to_custom_reference(
    bim_file=bim_file,
    build_type="GRCh38",  # Specify either "GRCh37" or "GRCh38"
    output_bim_file=output_bim_file
)


# Generate Converted Files (.bim, .bed, .fam)
from IGHmapper.igh_converter import generate_binary_files

prefix = r"path_to_your_bim_prefix"
output_prefix = r"path_to_your_output_folder\final_output"

generate_binary_files(
    prefix=prefix,
    converted_bim=output_bim_file,
    output_prefix=output_prefix,
    plink_path="plink2"  # Use "plink" or "plink2" depending on your installation
)


# Validate Mapping
from IGHmapper.igh_converter import validate_mapping

validate_mapping(
    custom_bim_file=r"path_to_your_converted_bim.bim",
    build_type="GRCh38",  # Specify either "GRCh37" or "GRCh38"
    output_success=r"path_to_your_output_folder\successfully_mapped_variants.txt",
    output_failure=r"path_to_your_output_folder\unsuccessfully_mapped_variants.txt"
)


```

---

## Modules

### `conversion.py`
- Handles the conversion of `.bim` files to the custom reference using mapping files.

### `mapping.py`
- Identifies genome build (GRCh37 or GRCh38).
- Extracts relevant loci from input files.

### `validation.py`
- Validates successful and unsuccessful mappings.
- Generates logs of successfully and unsuccessfully mapped variants.

### `utils.py`
- Utility functions to support the other modules (e.g., file handling, logging).

---

## API Reference
### `conversion.py`
```python
def convert_to_custom_reference(bim_file, map_file, output_bim_file):
```
- Converts `.bim` files to custom reference.
- Inputs: `.bim` file, mapping file, output path.
- Outputs: Converted `.bim` file.

### `mapping.py`
```python
def identify_build(bim_file):
```
- Identifies whether the input is `GRCh37` or `GRCh38`.
- Returns: `GRCh37`, `GRCh38`, or `Unknown`.

```python
def extract_region(bim_file, bed_file, fam_file, output_prefix):
```
- Extracts region `chr14:105000000-107500000`.
- Saves `.bim`, `.bed`, and `.fam` files.

### `validation.py`
```python
def validate_mapping(custom_bim_file, map_file, output_success, output_failure):
```
- Validates mapping process.
- Outputs: Successfully and unsuccessfully mapped variants.

---
## Reference File Handling

The package uses embedded reference files for both conversion and validation:

- **GRCh37 Reference File:** `igh_to_hg19.txt`
- **GRCh38 Reference File:** `igh_to_hg38.txt`

These files are automatically accessed from the `igh_converter` directory within the package.

### No User-Provided Reference Files Needed
Simply specify the build type (`GRCh37` or `GRCh38`) when running conversion or validation, and the package will automatically use the correct reference file.


## Contributing

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -m 'Add new feature'`).
4. Push to the branch (`git push origin feature-branch`).
5. Create a Pull Request.

---

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.


