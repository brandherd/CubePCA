import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cubePCA",
    version="0.1",
    author="Bernd Husemann",
    author_email="husemann@mpia.de",
    description="A small package for PCA sky subtraction in IFU datacubes",
    url="https://github.com/brandherd/CubePCA",
    packages = setuptools.find_packages(),
    scripts=['bin/subtractPCAsky.py','bin/createPCAsky.py','bin/applyPCAsky.py','bin/createMask.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)
