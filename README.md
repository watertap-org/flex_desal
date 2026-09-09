This is a repository for WaterTAP development and analysis for the NAWI 3.25 project.

## Installation

### 1. Create the conda environment

Open command prompt or terminal and from the location of this repo execute:

```bash
conda env create -f flex_desal_analysis_setup.yml
conda activate flex_desal_analysis
```

### 2. Install WaterTAP (development version)

This project requires the development version of WaterTAP installed in editable mode:

```bash
# Clone watertap repository (if not already done)
git clone https://github.com/watertap-org/watertap.git
cd watertap

# Install in editable mode without dependencies to avoid version conflicts
pip install --no-deps -e .
cd ..
```

### 3. Install idaes-pse 2.10

Force install the required version:

```bash
pip install --upgrade --force-reinstall idaes-pse==2.10.0
```

### 4. Install this package

Return to the flex_desal directory and install in editable mode:

```bash
pip install --no-deps -e .
```

### Notes

- This project requires Python 3.11
- WaterTAP 1.5.dev0+ is required for compatibility
- idaes-pse 2.10.0 is required (WaterTAP may try to downgrade this)
- After installation, you may see dependency warnings from pip - these can typically be ignored
