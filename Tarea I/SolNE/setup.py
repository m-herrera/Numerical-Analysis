from setuptools import setup, find_packages

setup(
    name='SolNE',
    description='This is my package',
    packages = find_packages(include=['SolNE', 'SolNE.*']),
    version='1.0',
    url='https://github.com/m-herrera/Numerical-Analysis',
    author='Kenneth Hernandez, Marco Herrera, Jasson Rodriguez',
    install_requires=[
        'matplotlib',
        'sympy'
    ],
    zip_safe=False
)
