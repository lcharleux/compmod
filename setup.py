from setuptools import setup

setup(name='compmod',
      version='0.1',
      description="Oh please please please, let me compile sphinx",
      long_description="",
      author='Ludovic Charleux',
      author_email='ludovic.charleux@univ-savoie.fr',
      license='GPL',
      packages=['compmod'],
      zip_safe=False,
      install_requires=[
          'numpy',
          'scipy'
          # 'Sphinx',
          # ^^^ Not sure if this is needed on readthedocs.org
          # 'something else?',
          ],
      )
