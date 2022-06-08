from setuptools import setup

setup(
    name="diploid_contam",
    packages=[
        "diploid_contam",
    ],
    install_requires=["numpy", "scipy", "pysam", "pydantic", "more_itertools"],
    python_requires=">3.6.1",
)
