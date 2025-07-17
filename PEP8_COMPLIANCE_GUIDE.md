# PEP 8 Compliance Guide for ZOL Project

## Overview
This guide provides comprehensive standards for maintaining PEP 8 compliance across the ZOL codebase. It covers manual practices to ensure consistent, readable, and maintainable Python code.

## Table of Contents
1. [Code Style Standards](#code-style-standards)
2. [File Organization](#file-organization)
3. [Naming Conventions](#naming-conventions)
4. [Documentation Standards](#documentation-standards)
5. [Common Issues and Fixes](#common-issues-and-fixes)
6. [Manual Code Review Checklist](#manual-code-review-checklist)
7. [Best Practices Summary](#best-practices-summary)

## Code Style Standards

### 1. Indentation
- **Use 4 spaces per indentation level**
- **Never use tabs**
- **Maximum line length: 88 characters**

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
2. **Imports**
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

## Best Practices Summary

1. **Follow naming conventions** consistently
2. **Write clear docstrings** for all public functions
3. **Use type hints** to improve code clarity
4. **Keep functions small** and focused
5. **Use meaningful variable names**
6. **Handle errors appropriately**
7. **Test your code** thoroughly
8. **Review code** before committing
9. **Keep dependencies updated**

## Resources

- [PEP 8 Style Guide](https://www.python.org/dev/peps/pep-0008/)

---

**Remember**: The goal is readable, maintainable code. When in doubt, prioritize clarity over brevity. 