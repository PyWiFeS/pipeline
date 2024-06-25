.. _installation:

Installation
============

To install PyWiFeS, follow these steps:

1. Clone the repository (automation branch):
   
   .. code-block:: bash
   
      git clone -b automation https://github.com/PyWiFeS/pipeline.git
   
2. Navigate to the project directory and install dependencies:
   
   .. code-block:: bash
   
      pip install .
   
3. Set the `PYWIFES_DIR` environment variable to your reference data directory:
   
   .. code-block:: bash
   
      export PYWIFES_DIR=/Users/.../pipeline/reference_data

4. Set up an alias for the main reduction routine:
   
   .. code-block:: bash
   
      alias pywifes-reduce='/Users/.../pipeline/reduction_scripts/reduce_data.py'
