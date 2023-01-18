# quotonic-public

[![Release](https://img.shields.io/github/v/release/jewaniuk/quotonic-public)](https://img.shields.io/github/v/release/jewaniuk/quotonic-public)
[![Build status](https://img.shields.io/github/actions/workflow/status/jewaniuk/quotonic-public/main.yml?branch=main)](https://github.com/jewaniuk/quotonic-public/actions/workflows/main.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/jewaniuk/quotonic-public/branch/main/graph/badge.svg)](https://codecov.io/gh/jewaniuk/quotonic-public)
[![Commit activity](https://img.shields.io/github/commit-activity/m/jewaniuk/quotonic-public)](https://img.shields.io/github/commit-activity/m/jewaniuk/quotonic-public)
[![License](https://img.shields.io/github/license/jewaniuk/quotonic-public)](https://img.shields.io/github/license/jewaniuk/quotonic-public)

A platform for performing efficient simulations in the quantum photonic domain.

- **Github repository**: <https://github.com/jewaniuk/quotonic-public/>
- **Documentation** <https://jewaniuk.github.io/quotonic-public/>

## Getting started with your project

First, create a repository on GitHub with the same name as this project, and then run the following commands:

``` bash
git init -b main
git add .
git commit -m "init commit"
git remote add origin git@github.com:jewaniuk/quotonic-public.git
git push -u origin main
```

Finally, install the environment and the pre-commit hooks with 

```bash
make install
```

You are now ready to start development on your project! The CI/CD
pipeline will be triggered when you open a pull request, merge to main,
or when you create a new release.

To finalize the set-up for publishing to PyPi or Artifactory, see
[here](https://fpgmaas.github.io/cookiecutter-poetry/features/publishing/#set-up-for-pypi).
For activating the automatic documentation with MkDocs, see
[here](https://fpgmaas.github.io/cookiecutter-poetry/features/mkdocs/#enabling-the-documentation-on-github).
To enable the code coverage reports, see [here](https://fpgmaas.github.io/cookiecutter-poetry/features/codecov/).

## Releasing a new version



---

Repository initiated with [fpgmaas/cookiecutter-poetry](https://github.com/fpgmaas/cookiecutter-poetry).