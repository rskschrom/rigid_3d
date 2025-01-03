from setuptools import setup

setup(name='rigidpy', version='0.0.1',
      description='Rigid body dynamics for particles',
      url='https://github.com/snowwxradar/rigid_3d',
      author='Robert Schrom',
      include_package_data=True,
      packages=['rigidpy','rigidpy.pybind11_lib'],
      package_dir={'rigidpy':'src',
                   'rigidpy.pybind11_lib':'src/pybind11_lib'})
