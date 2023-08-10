import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='revonto',
    author='Vladimir Smrkolj & Aljoša Škorjanc',
    author_email='todo',
    description='Python library for Gene Ontology Reverse Lookup',
    keywords='gene ontology, python, reverse lookup',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ladismrkolj/revonto',
    project_urls={
        'Documentation': 'https://github.com/ladismrkolj/revonto',
        'Bug Reports':
        'https://github.com/ladismrkolj/revonto/issues',
        'Source Code': 'https://github.com/ladismrkolj/revonto',
        # 'Funding': '',
        # 'Say Thanks!': '',
    },
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    classifiers=[
        # see https://pypi.org/classifiers/
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GPL-3.0 License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    # install_requires=['Pillow'],
    extras_require={
        'dev': ['check-manifest'],
        # 'test': ['coverage'],
    },
    # entry_points={
    #     'console_scripts': [  # This can provide executable scripts
    #         'run=revonto:main',
    # You can execute `run` in bash to run `main()` in src/revonto/__init__.py
    #     ],
    # },
)
