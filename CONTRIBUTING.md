# Contributing to OmniSciKit

Thank you for your interest in contributing to OmniSciKit! This document provides guidelines for contributing to the project.

## Code of Conduct

This project and everyone participating in it is governed by our Code of Conduct. By participating, you are expected to uphold this code.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check the existing issues to see if the problem has already been reported. When you are creating a bug report, please include as many details as possible:

- **Use a clear and descriptive title**
- **Describe the exact steps to reproduce the problem**
- **Provide specific examples to demonstrate the steps**
- **Describe the behavior you observed and what behavior you expected**
- **Include code samples and output**

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion, please include:

- **Use a clear and descriptive title**
- **Provide a step-by-step description of the suggested enhancement**
- **Provide specific examples to demonstrate the enhancement**
- **Explain why this enhancement would be useful**

### Pull Requests

1. Fork the repository
2. Create a new branch from `main`
3. Make your changes
4. Add or update tests as necessary
5. Update documentation
6. Submit a pull request

## Development Setup

### Prerequisites

- R (>= 4.1.0)
- Rtools (Windows) or Xcode (macOS)
- devtools package

### Setting Up Development Environment

```r
# Install devtools
install.packages("devtools")

# Clone the repository
git clone https://github.com/yourusername/OmniSciKit.git
cd OmniSciKit

# Load the package in development mode
devtools::load_all()

# Run tests
devtools::test()

# Check the package
devtools::check()
```

## Style Guidelines

### R Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable and function names
- Add comments for complex logic
- Keep functions focused and modular

### Documentation

- Document all exported functions with roxygen2
- Include examples in documentation
- Write clear, concise descriptions
- Use `@param` tags for all parameters

### Testing

- Write tests for all new functions
- Aim for high test coverage
- Use descriptive test names
- Test edge cases and error conditions

## Commit Messages

- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit the first line to 72 characters or less
- Reference issues and pull requests liberally after the first line

## Release Process

1. Update version number in `DESCRIPTION`
2. Update `NEWS.md` with changes
3. Run `devtools::check()` and ensure no errors
4. Create a new release on GitHub
5. Submit to CRAN if applicable

## Questions?

Feel free to open an issue for any questions or concerns.

Thank you for contributing to OmniSciKit!
