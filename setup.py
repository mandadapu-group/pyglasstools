from setuptools import setup, find_packages

setup(name='pyglasstools',
      version='0.1',
      description='A set of tools to analyze supercooled liquid and glass trajectories',
      url='http://github.com/muhammadhasyim/pyglasstools',
      author='Muhammad R. Hasyim',
      author_email='muhammad_hasyim@berkeley.edu',
      license='MIT',
      packages=find_packages(),
      package_data={"":["*.so"]},
      install_requires=['numpy','scipy'],
      zip_safe=False)
