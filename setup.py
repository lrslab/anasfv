from setuptools import setup
from anasfv import __version__

with open('readme.md', "r") as f:
    LONG_DESC = f.read()

setup(
    name='anasfv',
    version=__version__,
    packages=["anasfv"],
    url='https://github.com/lrslab/anasfv',
    license='MIT',
    author='like',
    author_email='taxoroger@gmail.com',
    description='The anasfv project focuses on analyzing nanopore-sequenced data of PCR-amplified African Swine Fever Virus (ASFV).',
    long_description=LONG_DESC,
    long_description_content_type='text/markdown',
    install_requires = ["biopython>=1.8.1",
                        "pandas==2.0.3",
                        "pyyaml>=5.3.1",
                        "pycirclize>=1.6.0"],
    include_package_data=True,
    scripts = ['anasfv/download_asfv_genome.py',
               'anasfv/find_near_ref.py',
               'anasfv/completeness.py',
               'anasfv/recombination_test.py',
               'anasfv/recombination_plot.py',
               'anasfv/get_cds_alignments.py',
               'anasfv/mapping_assembly.py',
               'anasfv/polish_asfv.py',
               'anasfv/make_tree.py',
                ]
)
