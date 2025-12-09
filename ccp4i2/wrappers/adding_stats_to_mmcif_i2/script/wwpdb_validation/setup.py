from setuptools import setup, find_packages

setup(
    name='wwPDB validation',
    version='0.1',
    url='https://github.com/berrisfordjohn/wwpdb_validation',
    author='John Berrisford',
    test_suite='tests',
    zip_safe=True,
    packages=['wwpdb_validation'],
    # dependency_links=['https://github.com/project-gemmi/gemmi/tarball/master#egg=gemmi'],
    install_requires=['requests',
                      'onedep_api>=0.15'
                      ],
)
