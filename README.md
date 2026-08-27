# Oropharyngeal Cancer Agent-Based Model
This agent-based model (ABM) models oropharyngeal cancer using patient spatial transcritomics data on CompuCell3D.

# Installing CompuCell3D
To run this model, please first install CompuCell3D. 
To install CompuCell3D on DRAC, from a login node:
```bash
cd $SCRATCH

module load python/3.11 vtk/9.4 swig
git clone --depth 1 --branch 4.9.0 https://github.com/CompuCell3D/CompuCell3D.git CompuCell3D-source
cd CompuCell3D-source/CompuCell3D/

cmake --fresh -B _build -DCMAKE_INSTALL_PREFIX=$HOME/CompuCell3D
srun --account=def-nicoleli_cpu -c 16 --mem-per-cpu=3500M -- cmake --build _build/ --parallel 16
cmake --install _build

virtualenv --clear ~/CompuCell3D-ENV && source ~/CompuCell3D-ENV/bin/activate
pip install -U pip
wget https://files.pythonhosted.org/packages/39/8d/e4b927fc3fe6b7f07fb38662702e8ebcc0f5a460676ac647ec2b5256e6b8/libroadrunner-2.9.2-cp311-cp311-manylinux_2_28_x86_64.whl -O ~/libroadrunner-2.9.2-cp311-cp311-linux_x86_64.whl
setrpaths --path ~/libroadrunner-2.9.2-cp311-cp311-linux_x86_64.whl --any_interpreter

pip install lxml deprecated ~/libroadrunner-2.9.2-cp311-cp311-linux_x86_64.whl 
pip freeze --local > ~/compucell3d.reqs
```
# Running the Model
To run the model on DRAC, include the following in your slurm script:
```bash
module load python/3.11 vtk/9.4 cuda/12.6

virtualenv $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r ~/compucell3d.reqs

export PYTHONPATH=$PYTHONPATH:$HOME/CompuCell3D/lib/site-packages

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

python -X faulthandler -m cc3d.run_script -i output_directory/OP_Cancer_2D/OP_Cancer_2D.cc3d -f 10 -o output_directory>
```
