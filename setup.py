import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vpolo",
    version="0.1.6",
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
    install_requires=[
       'pandas',
       'scipy'
    ]
)
