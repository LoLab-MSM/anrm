#!/usr/bin/env python

from distutils.core import setup

def main():
    setup(name='anrm',
          version='1.0',
          description='Extrinsic Apoptosis Reaction Model',
          authors=['Michael Irvin','Carlos Lopez']
          author_emaisl=['michael.w.irvin@vanderbilt.edu','c.lopez@vanderbilt.edu']
          packages=['anrm'],
          requires=['pysb'],
          keywords=['systems', 'biology', 'model', 'rules'],
          classifiers=[
            'Intended Audience :: Science/Research',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            ],
          )

if __name__ == '__main__':
    main()
