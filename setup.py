from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='n50',
      version='0.9',
      description='N50 calculation from bacterial assemblies',
      long_description=readme(),
      url='https://github.com/asahaman/n50.git',
      author='Arnab Saha Mandal',
      author_email='arnab.sahamandal@umanitoba.ca',
      license='MIT',
      packages=['n50'],
      install_requires=[
          'Bio','markdown'
      ],
      zip_safe=False)