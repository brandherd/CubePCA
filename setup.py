import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cubePCA",
    version="0.1",
    author="Bernd Husemann",
    author_email="husemann@mpia.de",
    description="A small example package",
    url="https://github.com/pypa/sampleproject",
    packages = setuptools.find_packages(),
    scripts=['bin/subtractPCAsky'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)
