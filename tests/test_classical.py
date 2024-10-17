import os
import subprocess
import glob
from pathlib import Path
import pytest


class TestClassical:
    cwd = Path.cwd()
    reduce_script = '/code/dev/reduction_scripts/reduce_data.py'
    ref_data = os.environ['PYWIFES_DIR']

    @staticmethod
    def get_input():
        return ['classical_20230922']

    # a helper function to make a symlink
    # using Path.symlink_to
    def mksymlink(self,fpath,ifile):
        """
        Make a symlink in a specific folder.

        Parameters
        ----------
        fpath: pathlib.Path
            Folder to place the symlink
        ifile: pathlib.Path
            File to make the symlink to
        """
        symlink = fpath / ifile.name
        # remove if it already exists
        if Path(symlink).exists():
            os.remove(str(Path(symlink)))
        Path(symlink).symlink_to(str(ifile))

    @pytest.fixture
    def run_classical(self):
        def _run_classical(request,tmp_path):
            ifolder = self.cwd / request.param
            # get a list of fits files in the input folder    
            ifiles = glob.glob(str(ifolder) + "/*.fits")
            # make symlinks to the fits files in the tmp folder 
            for ifile in ifiles:
                self.mksymlink(tmp_path,Path(ifile))
            # go into the tmp folder and run the pipeline
            os.chdir(str(tmp_path))

            cmd = self.reduce_script + " ."
            #cmd = """echo 'HELLO BRENT'"""
#           with subprocess.Popen(args=cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,cwd='.',env={'PYWIFES_DIR' : self.ref_data},bufsize=1,universal_newlines=True) as p:
#               for line in p.stdout:
#                   print (line.strip())
            # check that the output products exist
            # TODO
            return tmp_path
        yield _run_classical
        os.chdir(str(self.cwd))


    @pytest.fixture(params=get_input())
    def check_output(self,request,tmp_path,run_classical):
        # check that the folder returned by the run_classical fixture 
        # contains the datacube results we require   
        ofolder = run_classical(request,tmp_path)
        ofiles = [
            'OBK-124128-WiFeS-Blue-UT20230922T090650-1.cube.fits',
            'OBK-525056-WiFeS-Blue-UT20230922T121452-0.cube.fits',
            'OBK-124128-WiFeS-Red--UT20230922T090650-1.cube.fits',
            'OBK-525056-WiFeS-Red--UT20230922T121452-0.cube.fits'
        ]
        for o in ofiles:
            data_product = Path(ofolder) / 'data_products' / o
            assert Path(data_product).exists()

    def test_run_classical(self,check_output):
        assert 1==0
        pass

