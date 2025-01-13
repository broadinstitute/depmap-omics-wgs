# PureCN postprocess

Collect outputs from PureCN as a single JSON file and call chromosomal instability and whole genome doubling.

## Installation

1. Install the required system dependencies:
    - [pyenv](https://github.com/pyenv/pyenv)
    - [Poetry](https://python-poetry.org/)

2. Install the required Python version (developed with 3.12.3, but other 3.12+ versions should work):
   ```shell
   pyenv install "$(cat .python-version)"
   ```

3. Confirm that `python` maps to the correct version:
   ```
   python --version
   ```

4. Set the Poetry interpreter and install the Python dependencies:
   ```shell
   poetry env use "$(pyenv which python)"
   poetry install
   ```

 ## Usage

```shell
poetry run python -m purecn_postprocess \
    --solution="purecn_solution.csv" \
    --loh="purecn_loh.csv" \
    --out="out.json"
```
