from setuptools import setup, find_packages

setup(name='pyglasstools',
      version='0.1',
      description='A set of tools to analyze trajectories of atomistic supercooled liquid and glass formers',
      url='http://github.com/muhammadhasyim/pyglasstools',
      author='Muhammad R. Hasyim',
      author_email='muhammad_hasyim@berkeley.edu',
      license='MIT',
      packages=find_packages(),
      package_data={"":["*.so"]},
      install_requires=['numpy','scipy','gsd','tqdm'],
      zip_safe=False)
