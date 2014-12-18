from setuptools import setup

setup(name='compmod',
      version='0.1',
      description="Compmod project",
      long_description="",
      author='Ludovic Charleux, Moustapha Issack, Laurent Bizet',
      author_email='ludovic.charleux@univ-savoie.fr',
      license='GPL v2',
      packages=['compmod'],
      zip_safe=False,
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib'
          ],
      )
