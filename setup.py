from setuptools import setup

setup(name='pyglasstools',
      version='0.1',
      description='A set of tools to analyze supercooled liquid and glass trajectories',
      #url='http://github.com/storborg/funniest',
      author='Muhammad R. Hasyim',
      author_email='muhammad_hasyim@berkeley.edu',
      license='MIT',
      packages=['pyglasstools'],
      install_requires=['numba','scipy','numba','tqdm','pyfftw'],
      zip_safe=False)
