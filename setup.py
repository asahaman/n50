from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='n50',
      version='1.0',
      description='N50 calculation from bacterial assemblies',
      long_description=readme(),
      classifiers=[
        'Development Status :: 1 - Planning'
        'License :: OSI Approved :: Academic Free License (AFL)'
        'Programming Language :: Python :: 3.6'
        'Topic :: Text Processing :: General'
        'Operating System :: Unix'
      ],
      url='https://github.com/asahaman/n50.git',
      author='Arnab Saha Mandal',
      author_email='arnab.sahamandal@umanitoba.ca',
      license='MIT',
      packages=['n50'],
      install_requires=[
          'biopython',
      ],
      entry_points = {
        'console_scripts': ['contig_dnld=n50.fetch_contigs:main',
                            'contig_n50=n50.calculating_n50_assemblies:main'
                          ],
      },
      zip_safe=False)