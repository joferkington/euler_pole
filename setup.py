from setuptools import setup

setup(
    name = 'euler_pole',
    version = '0.1',
    description = "Euler poles for plate motion calculations",
    author = 'Joe Kington',
    author_email = 'joferkington@gmail.com',
    license = 'LICENSE',
    url = 'https://github.com/joferkington/euler_pole',
    packages = ['euler_pole'],
    install_requires = [
        'numpy >= 1.1',
        ]
)
