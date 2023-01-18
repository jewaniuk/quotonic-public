# -*- coding: utf-8 -*-
from setuptools import setup

from build import build

packages = ["quotonic"]

package_data = {"": ["*"]}

install_requires = ["cython>=0.29.33,<0.30.0", "numpy>=1.24.1,<2.0.0"]

setup_kwargs = {
    "name": "quotonic",
    "version": "0.0.1",
    "description": "A platform for performing efficient simulations in the quantum photonic domain.",
    "long_description": (
        "# quotonic-public\n\n[![Release](https://img.shields.io/github/v/release/jewaniuk/quotonic-public)](https://img.shields.io/github/v/release/jewaniuk/quotonic-public)\n[![Build"
        " status](https://img.shields.io/github/actions/workflow/status/jewaniuk/quotonic-public/main.yml?branch=main)](https://github.com/jewaniuk/quotonic-public/actions/workflows/main.yml?query=branch%3Amain)\n[![Commit"
        " activity](https://img.shields.io/github/commit-activity/m/jewaniuk/quotonic-public)](https://img.shields.io/github/commit-activity/m/jewaniuk/quotonic-public)\n[![License](https://img.shields.io/github/license/jewaniuk/quotonic-public)](https://img.shields.io/github/license/jewaniuk/quotonic-public)\n\nA"
        " platform for performing efficient simulations in the quantum photonic domain.\n\n- **Github repository**:"
        " <https://github.com/jewaniuk/quotonic-public/>\n- **Documentation**"
        " <https://jewaniuk.github.io/quotonic-public/>\n\n---\n\nRepository initiated with"
        " [fpgmaas/cookiecutter-poetry](https://github.com/fpgmaas/cookiecutter-poetry)."
    ),
    "author": "Jacob Ewaniuk",
    "author_email": "fjacob.ewaniuk@queensu.ca",
    "maintainer": "None",
    "maintainer_email": "None",
    "url": "https://github.com/jewaniuk/quotonic-public",
    "packages": packages,
    "package_data": package_data,
    "install_requires": install_requires,
    "python_requires": ">=3.8,<4.0",
}

build(setup_kwargs)

setup(**setup_kwargs)
