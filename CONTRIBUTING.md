# Contributing to CRC Transcriptomics Analysis

Thank you for your interest in contributing! This document provides guidelines for contributing to this project.

## Ways to Contribute

### 1. Report Issues
- **Bug reports**: Describe the problem, steps to reproduce, expected vs actual behavior
- **Feature requests**: Explain the use case and proposed solution
- **Documentation improvements**: Identify unclear or missing documentation

### 2. Improve Code
- **Code optimization**: Make analysis faster or more memory-efficient
- **Additional analyses**: Add complementary methods (e.g., survival analysis, subtype classification)
- **Visualization enhancements**: Improve or add new plots
- **Testing**: Add unit tests or validation checks

### 3. Enhance Documentation
- **Tutorial notebooks**: Create Jupyter/R Markdown tutorials
- **Use case examples**: Show applications to other datasets
- **Methods explanations**: Expand on statistical methods
- **Troubleshooting guides**: Document common issues and solutions

## Getting Started

### Fork and Clone

```bash
# Fork repository on GitHub
# Clone your fork
git clone https://github.com/YOUR-USERNAME/crc-transcriptomics-analysis.git
cd crc-transcriptomics-analysis

# Add upstream remote
git remote add upstream https://github.com/ORIGINAL-OWNER/crc-transcriptomics-analysis.git
```

### Create Branch

```bash
# Create feature branch
git checkout -b feature/your-feature-name

# Or for bug fixes
git checkout -b fix/bug-description
```

### Set Up Environment

```bash
# Using conda
conda env create -f environment.yml
conda activate crc-analysis

# Test installation
R -e "library(limma); library(ReactomePA)"
```

## Development Workflow

### 1. Make Changes

- Write clear, commented code
- Follow existing code style
- Update documentation as needed
- Add examples for new features

### 2. Test Changes

```r
# Run main analysis to ensure no breaks
source("scripts/main_analysis.R")

# Check all outputs generated
list.files("results/", recursive = TRUE)
list.files("figures/", recursive = TRUE)
```

### 3. Commit Changes

```bash
# Stage changes
git add file1.R file2.md

# Commit with descriptive message
git commit -m "Add feature: description of what was added"
```

**Commit Message Format**:
```
<type>: <short description>

<longer description if needed>

<reference to issue if applicable>
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

Examples:
```
feat: Add survival analysis module

Implemented Cox proportional hazards analysis for pathway scores
Includes visualization of Kaplan-Meier curves

Closes #15
```

```
fix: Correct GSEA permutation seed

Fixed bug where GSEA results were not reproducible
Added set.seed() before gsePathway call
```

### 4. Push and Pull Request

```bash
# Push to your fork
git push origin feature/your-feature-name

# Create pull request on GitHub
# Provide clear description of changes
```

## Code Style Guidelines

### R Code

**General Principles**:
- Use descriptive variable names
- Comment non-obvious operations
- Keep functions focused and modular
- Avoid hard-coded values

**Naming Conventions**:
```r
# Variables and functions: snake_case
sample_info <- data.frame(...)
calculate_enrichment <- function(genes) { ... }

# Constants: UPPER_SNAKE_CASE
MAX_PVALUE <- 0.05
MIN_FOLD_CHANGE <- 1
```

**Function Structure**:
```r
#' Brief description
#'
#' Longer description if needed
#'
#' @param x Description of parameter
#' @return Description of return value
#' @examples
#' example_function(data)
example_function <- function(x) {
  # Validate input
  stopifnot(is.numeric(x))
  
  # Main logic
  result <- process(x)
  
  # Return
  return(result)
}
```

**Code Organization**:
```r
## ---- Section Header ----

# Single line comment for brief explanation

# Multi-line comment
# for longer explanations
# that need multiple lines
```

### Documentation

**README Updates**:
- Keep table of contents current
- Update installation instructions if dependencies change
- Add new results to Results section
- Update version numbers

**Code Comments**:
```r
# Good: Explains WHY
# Use quantile normalization to remove batch effects
expr_norm <- normalizeBetweenArrays(expr)

# Bad: Explains WHAT (obvious from code)
# Normalize expression
expr_norm <- normalizeBetweenArrays(expr)
```

## Testing

### Manual Testing Checklist

Before submitting PR, verify:

- [ ] Code runs without errors
- [ ] All plots generated correctly
- [ ] Results files created with expected columns
- [ ] No new warnings introduced
- [ ] Works with fresh conda environment
- [ ] Documentation updated
- [ ] Examples run successfully

### Test Datasets

When adding features, test with:
1. Original GSE33126 (included)
2. Similar paired microarray dataset
3. Edge cases (small sample size, no significant genes, etc.)

## Pull Request Process

### 1. Pre-submission Checklist

- [ ] Branch is up-to-date with main
- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] Documentation updated
- [ ] Commit messages are clear
- [ ] No merge conflicts

### 2. PR Description Template

```markdown
## Description
Brief summary of changes

## Motivation
Why is this change needed?

## Changes Made
- Change 1
- Change 2
- Change 3

## Testing
How was this tested?

## Checklist
- [ ] Code follows project style
- [ ] Documentation updated
- [ ] Tests added/updated
- [ ] No breaking changes

## Related Issues
Closes #issue_number
```

### 3. Review Process

1. Maintainers review code
2. Automated checks run (if configured)
3. Discussion and revisions
4. Approval and merge

### 4. After Merge

- Delete feature branch
- Update local main branch
- Close related issues

## Reporting Issues

### Bug Reports

Include:
```markdown
**Describe the bug**
Clear description of the problem

**To Reproduce**
Steps to reproduce:
1. Run script '...'
2. With parameters '...'
3. See error

**Expected behavior**
What should happen

**Actual behavior**
What actually happened

**Environment**
- R version:
- Package versions:
- OS:

**Additional context**
Any other relevant information
```

### Feature Requests

Include:
```markdown
**Problem**
What problem would this solve?

**Proposed Solution**
How might this be implemented?

**Alternatives Considered**
Other approaches you've thought about

**Additional Context**
Any other relevant information
```

## Code Review Guidelines

### For Reviewers

**Be Constructive**:
- ✅ "Consider using dplyr::filter() here for clarity"
- ❌ "This code is bad"

**Focus on**:
- Correctness
- Clarity
- Performance (if relevant)
- Maintainability
- Documentation

**Ask Questions**:
- "Could you explain the rationale for this approach?"
- "Have you considered edge case X?"

### For Contributors

- Respond to all comments
- Don't take criticism personally
- Ask for clarification if needed
- Be open to suggestions

## Getting Help

- **Questions**: Open a GitHub Discussion
- **Bug reports**: Create an Issue
- **Real-time chat**: (Add Slack/Discord if available)
- **Email**: (Add maintainer email)

## Recognition

Contributors are recognized in:
- GitHub contributors page
- CHANGELOG.md
- Publication acknowledgments (for significant contributions)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to this project! Your help makes this resource more valuable for the research community.
