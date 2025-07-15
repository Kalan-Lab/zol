# PEP 8 Compliance Guide for ZOL Project

## Overview
This guide provides comprehensive standards for maintaining PEP 8 compliance across the ZOL codebase. It covers both automated tools and manual practices to ensure consistent, readable, and maintainable Python code.

## Table of Contents
1. [Automated Tools](#automated-tools)
2. [Code Style Standards](#code-style-standards)
3. [File Organization](#file-organization)
4. [Naming Conventions](#naming-conventions)
5. [Documentation Standards](#documentation-standards)
6. [Common Issues and Fixes](#common-issues-and-fixes)
7. [Pre-commit Hooks](#pre-commit-hooks)
8. [Continuous Integration](#continuous-integration)

## Automated Tools

### Required Tools
```bash
# Install required tools
pip install black flake8 mypy isort

# Or add to requirements.txt
black>=23.0.0
flake8>=6.0.0
mypy>=1.0.0
isort>=5.12.0
```

### Configuration Files

#### `.flake8` Configuration
```ini
[flake8]
max-line-length = 88
extend-ignore = 
    E203,  # whitespace before ':'
    E501,  # line too long (handled by black)
    W503,  # line break before binary operator
    W504   # line break after binary operator
exclude = 
    .git,
    __pycache__,
    .venv,
    venv,
    .mypy_cache,
    build,
    dist,
    *.egg-info
per-file-ignores =
    __init__.py:F401
```

#### `pyproject.toml` for Black
```toml
[tool.black]
line-length = 88
target-version = ['py38', 'py39', 'py310', 'py311']
include = '\.pyi?$'
extend-exclude = '''
/(
  # directories
  \.eggs
  | \.git
  | \.hg
  | \.mypy_cache
  | \.tox
  | \.venv
  | build
  | dist
)/
'''
```

#### `pyproject.toml` for isort
```toml
[tool.isort]
profile = "black"
multi_line_output = 3
line_length = 88
known_first_party = ["zol"]
known_third_party = ["numpy", "pandas", "scipy", "Bio", "ete3"]
sections = ["FUTURE", "STDLIB", "THIRDPARTY", "FIRSTPARTY", "LOCALFOLDER"]
```

## Code Style Standards

### 1. Indentation
- **Use 4 spaces per indentation level**
- **Never use tabs**
- **Maximum line length: 88 characters** (Black default)

```python
# ✅ Correct
def long_function_name(
    var_one,
    var_two,
    var_three,
):
    print(var_one)

# ❌ Incorrect
def long_function_name(var_one, var_two, var_three):
    print(var_one)
```

### 2. Imports
- **Group imports in this order:**
  1. Standard library imports
  2. Third-party imports
  3. Local application imports
- **Use absolute imports when possible**
- **Use `isort` to automatically organize imports**

```python
# ✅ Correct
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from Bio import SeqIO

from zol import util
from zol.orthologs import findOrthologs
```

### 3. String Formatting
- **Prefer f-strings over other formatting methods**
- **Use `.format()` for complex formatting**
- **Avoid `%` formatting**

```python
# ✅ Preferred
name = "World"
print(f"Hello, {name}!")

# ✅ For complex cases
print("Value: {:.2f}".format(3.14159))

# ❌ Avoid
print("Hello, %s!" % name)
```

### 4. Variable and Function Naming
- **Use `snake_case` for variables and functions**
- **Use `UPPER_CASE` for constants**
- **Use `PascalCase` for classes**

```python
# ✅ Correct
def calculate_gene_cluster_score():
    MAX_SCORE = 100
    gene_cluster_name = "cluster_001"
    
class GeneClusterAnalyzer:
    pass
```

### 5. Type Hints
- **Use type hints for function parameters and return values**
- **Import types from `typing` module**

```python
# ✅ Correct
from typing import Dict, List, Optional, Union

def process_genbank_file(
    file_path: str,
    output_dir: Optional[str] = None,
) -> Dict[str, List[str]]:
    """Process GenBank file and return results."""
    pass
```

## File Organization

### 1. File Structure
```
src/zol/
├── __init__.py
├── data_dictionary.py
├── fai.py
├── util.py
├── zol.py
└── orthologs/
    ├── __init__.py
    ├── findOrthologs.py
    └── ...
```

### 2. Module Organization
Each Python file should follow this structure:
1. **Module docstring**
2. **Imports** (organized by isort)
3. **Constants**
4. **Classes**
5. **Functions**
6. **Main execution block** (if applicable)

```python
"""
Module for processing GenBank files.

This module provides utilities for parsing and analyzing GenBank format files
used in the ZOL pipeline.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from Bio import SeqIO

# Constants
DEFAULT_THREADS = 1
MAX_MEMORY_GB = 16

# Classes
class GenBankProcessor:
    """Process GenBank files for analysis."""
    
    def __init__(self, input_dir: str):
        self.input_dir = Path(input_dir)

# Functions
def parse_genbank_file(file_path: str) -> Dict[str, str]:
    """Parse a GenBank file and return sequence data."""
    pass

# Main execution
if __name__ == "__main__":
    main()
```

## Naming Conventions

### 1. Variables and Functions
```python
# ✅ Good names
gene_cluster_count = 0
calculate_ortholog_groups()
process_diamond_results()

# ❌ Poor names
gc_count = 0
calc_og()
process_dr()
```

### 2. Constants
```python
# ✅ Constants in UPPER_CASE
DEFAULT_THREADS = 1
MAX_MEMORY_GB = 16
DIAMOND_SENSITIVITY = "very-sensitive"
```

### 3. Classes
```python
# ✅ PascalCase for classes
class GeneClusterAnalyzer:
    pass

class OrthologGroupFinder:
    pass
```

## Documentation Standards

### 1. Docstrings
Use Google-style docstrings for all public functions and classes:

```python
def process_genbank_files(
    input_directory: str,
    output_directory: str,
    threads: int = 1,
) -> Dict[str, List[str]]:
    """Process GenBank files in the input directory.
    
    Args:
        input_directory: Path to directory containing GenBank files.
        output_directory: Path to directory for output files.
        threads: Number of threads to use for processing.
        
    Returns:
        Dictionary mapping file names to processed results.
        
    Raises:
        FileNotFoundError: If input directory doesn't exist.
        ValueError: If threads is less than 1.
    """
    pass
```

### 2. Inline Comments
- **Use comments sparingly**
- **Explain "why" not "what"**
- **Keep comments up to date**

```python
# ✅ Good comment
# Skip processing if no valid files found
if not valid_files:
    return

# ❌ Poor comment
# Set x to 5
x = 5
```

## Common Issues and Fixes

### 1. Indentation Issues
```python
# ❌ Mixed tabs and spaces
def function():
	# This line uses tabs
    # This line uses spaces
    pass

# ✅ Consistent spaces
def function():
    # All lines use 4 spaces
    pass
```

### 2. Line Length Issues
```python
# ❌ Too long
result = very_long_function_name_with_many_parameters(param1, param2, param3, param4, param5, param6, param7, param8)

# ✅ Properly formatted
result = very_long_function_name_with_many_parameters(
    param1,
    param2,
    param3,
    param4,
    param5,
    param6,
    param7,
    param8,
)
```

### 3. Import Issues
```python
# ❌ Wildcard imports
from zol import *

# ✅ Specific imports
from zol import util, fai
from zol.orthologs import findOrthologs
```

### 4. Path Concatenation
```python
# ❌ Old Path object usage (if removing pathlib support)
file_path = Path("dir") / "file.txt"

# ✅ String concatenation
file_path = "dir" + "file.txt"

# ✅ os.path.join (alternative)
import os
file_path = os.path.join("dir", "file.txt")
```

## Pre-commit Hooks

### 1. Install pre-commit
```bash
pip install pre-commit
pre-commit install
```

### 2. `.pre-commit-config.yaml`
```yaml
repos:
  - repo: https://github.com/psf/black
    rev: 23.3.0
    hooks:
      - id: black
        language_version: python3

  - repo: https://github.com/pycqa/isort
    rev: 5.12.0
    hooks:
      - id: isort
        args: ["--profile", "black"]

  - repo: https://github.com/pycqa/flake8
    rev: 6.0.0
    hooks:
      - id: flake8
        args: [--max-line-length=88]

  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.3.0
    hooks:
      - id: mypy
        additional_dependencies: [types-all]
```

## Continuous Integration

### 1. GitHub Actions Workflow
Create `.github/workflows/pep8-check.yml`:

```yaml
name: PEP 8 Compliance Check

on: [push, pull_request]

jobs:
  pep8-check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.9'
          
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install black flake8 mypy isort
          
      - name: Check code formatting with Black
        run: black --check --diff .
        
      - name: Check import sorting with isort
        run: isort --check-only --diff .
        
      - name: Lint with flake8
        run: flake8 .
        
      - name: Type check with mypy
        run: mypy src/
```

## Manual Code Review Checklist

Before committing code, ensure:

- [ ] Code follows PEP 8 style guidelines
- [ ] All imports are properly organized
- [ ] Type hints are used where appropriate
- [ ] Docstrings are present for public functions
- [ ] Variable names are descriptive and follow conventions
- [ ] No hardcoded paths or magic numbers
- [ ] Error handling is appropriate
- [ ] No unused imports or variables
- [ ] Line length is under 88 characters
- [ ] Proper indentation (4 spaces, no tabs)

## Quick Commands

### Format Code
```bash
# Format all Python files
black .

# Sort imports
isort .

# Check formatting without changes
black --check --diff .
isort --check-only --diff .
```

### Lint Code
```bash
# Run flake8
flake8 .

# Run mypy
mypy src/

# Run all checks
black --check --diff . && isort --check-only --diff . && flake8 . && mypy src/
```

### Fix Common Issues
```bash
# Auto-fix formatting
black .
isort .

# Remove unused imports
autoflake --in-place --remove-all-unused-imports --recursive .
```

## Best Practices Summary

1. **Use automated tools** (Black, isort, flake8, mypy)
2. **Follow naming conventions** consistently
3. **Write clear docstrings** for all public functions
4. **Use type hints** to improve code clarity
5. **Keep functions small** and focused
6. **Use meaningful variable names**
7. **Handle errors appropriately**
8. **Test your code** thoroughly
9. **Review code** before committing
10. **Keep dependencies updated**

## Resources

- [PEP 8 Style Guide](https://www.python.org/dev/peps/pep-0008/)
- [Black Documentation](https://black.readthedocs.io/)
- [isort Documentation](https://pycqa.github.io/isort/)
- [flake8 Documentation](https://flake8.pycqa.org/)
- [mypy Documentation](https://mypy.readthedocs.io/)

---

**Remember**: The goal is readable, maintainable code. When in doubt, prioritize clarity over brevity. 