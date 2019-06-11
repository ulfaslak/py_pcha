from setuptools import setup

with open("README.md") as fp:
	long_description = fp.read()

setup(
    name='py_pcha',
    version='0.1.3',
    description='Python implemenation of PCHA algorithm for Archetypal Analysis',
    url='https://github.com/ulfaslak/py_pcha',
    # download_url = 'https://github.com/ulfaslak/py_pcha/tarball/0.1',
    author='Ulf Aslak',
    author_email='ulfjensen@gmail.com',
    license='MIT',
    packages=['py_pcha'],
    zip_safe=False,
    keywords = ['Archetypal Analysis', 'AA', 'PCHA', 'Clustering', 'Machine Learning', 'Signal Processing'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics'
    ]
)
