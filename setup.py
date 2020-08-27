import sys

import setuptools
try:
    from setuptools_rust import RustExtension, Binding
except ImportError:
    import subprocess

    errno = subprocess.call([sys.executable, "-m", "pip", "install", "setuptools-rust"])
    if errno:
        print("Please install setuptools-rust package")
        raise SystemExit(errno)
    else:
        from setuptools_rust import RustExtension, Binding

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vpolo",
    version="0.2.1",
    author="Avi Srivastava",
    author_email="asrivastava@cs.stonybrook.edu",
    description="Support packages for Alevin tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/k3yavi/vpolo",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    zip_safe=False,
    rust_extensions=[RustExtension("sce.sce", binding=Binding.RustCPython)],
    install_requires=[
       'pandas',
       'scipy',
    ]
)
