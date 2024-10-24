import os
import subprocess
import glob
from pathlib import Path
import pytest


class TestReduction:
    cwd = Path.cwd()
    reduce_script = '/code/dev/reduction_scripts/reduce_data.py'
    ref_data = os.environ['PYWIFES_DIR']

    @staticmethod
    def get_datasets():
        return [
                'classical_20230922',
                'nod-and-shuffle_20230927',
                'raw_data_half_frame',
                'wifes_b3000_560_r3000_1x2',
                'wifes_b7000_615_i7000_1x2',
                'wifes_stellar_b3000_560_r3000_1x2',
                'wifes_stellar_u7000_480_r7000_1x2',
                'wifes_subns_b3000_560_r3000_1x2',
               ]

    def expected_files(self,dataset):
        # returns a list of cubes to check whether they exist...
        cubes = {
            'classical_20230922' : 
            [
                'OBK-124128-WiFeS-Blue-UT20230922T090650-1.cube.fits',
                'OBK-525056-WiFeS-Blue-UT20230922T121452-0.cube.fits',
                'OBK-124128-WiFeS-Red--UT20230922T090650-1.cube.fits',
                'OBK-525056-WiFeS-Red--UT20230922T121452-0.cube.fits'
            ],
            'nod-and-shuffle_20230927' : 
            [
                'OBK-124576-WiFeS-Blue-UT20230926T085813-5.cube.fits',
                'OBK-124576-WiFeS-Red--UT20230926T085813-5.cube.fits',
                'OBK-527648-WiFeS-Blue-UT20230927T152726-5.cube.fits',
                'OBK-527648-WiFeS-Red--UT20230927T152726-5.cube.fits',
                'OBK-527936-WiFeS-Blue-UT20230927T132601-9.cube.fits',
                'OBK-527936-WiFeS-Red--UT20230927T132601-9.cube.fits',
                'OBK-527968-WiFeS-Blue-UT20230927T141355-9.cube.fits',
                'OBK-527968-WiFeS-Red--UT20230927T141355-9.cube.fits',
                'OBK-528000-WiFeS-Blue-UT20230927T161205-5.cube.fits',
                'OBK-528000-WiFeS-Red--UT20230927T161205-5.cube.fits',
                'OBK-528032-WiFeS-Blue-UT20230927T165551-6.cube.fits',
                'OBK-528032-WiFeS-Red--UT20230927T165551-6.cube.fits',
            ],
            'raw_data_half_frame' : 
            [
                'OBK-727168-WiFeS-Blue-UT20240525T191600-2.cube.fits',
                'cut_OBK-124224-WiFeS-Blue-UT20240523T081350-9.cube.fits',
                'OBK-727168-WiFeS-Red--UT20240525T191600-2.cube.fits',
                'cut_OBK-124224-WiFeS-Red--UT20240523T081350-9.cube.fits',
            ],
            'wifes_b3000_560_r3000_1x2' : 
            [
                'OBK-124288-WiFeS-Blue-UT20240624T140040-4.cube.fits',
                'OBK-124288-WiFeS-Red--UT20240624T140040-4.cube.fits',
                'OBK-124352-WiFeS-Blue-UT20240529T083012-7.cube.fits',
                'OBK-124352-WiFeS-Red--UT20240529T083012-7.cube.fits',
                'OBK-238848b-WiFeS-Blue-UT20240529T072137-8b.cube.fits',
                'OBK-238848b-WiFeS-Red--UT20240529T072137-8b.cube.fits',
                'OBK-636224-WiFeS-Blue-UT20240609T193227-6.cube.fits',
                'OBK-636224-WiFeS-Red--UT20240609T193227-6.cube.fits',
            ],
            'wifes_b7000_615_i7000_1x2' :
            [
                'OBK-124416-WiFeS-Blue-UT20240704T082822-2.cube.fits',
                'OBK-124672-WiFeS-Blue-UT20240710T131935-4.cube.fits',
                'OBK-124416-WiFeS-Red--UT20240704T082822-2.cube.fits',
                'OBK-124672-WiFeS-Red--UT20240710T131935-4.cube.fits',
            ],
            'wifes_stellar_b3000_560_r3000_1x2' : 
            [
                'cut_OBK-124224-WiFeS-Blue-UT20230512T081947-9.cube.fits',
                'cut_OBK-124224-WiFeS-Red--UT20230512T081947-9.cube.fits',
                'cut_OBK-124288-WiFeS-Blue-UT20230529T161409-9.cube.fits',
                'cut_OBK-124288-WiFeS-Red--UT20230529T161409-9.cube.fits',
                'cut_OBK-999999-WiFeS-Blue-UT20230514T073157-0.cube.fits',
                'cut_OBK-999999-WiFeS-Red--UT20230514T073157-0.cube.fits',
                'OBK-383520-WiFeS-Blue-UT20230517T124646-3.cube.fits',
                'OBK-383520-WiFeS-Red--UT20230517T124646-3.cube.fits',
            ],
            'wifes_stellar_u7000_480_r7000_1x2' : 
            [
                'cut_OBK-518144-WiFeS-Blue-UT20240801T084638-2.cube.fits',
                'cut_OBK-518144-WiFeS-Red--UT20240801T084638-2.cube.fits',
                'cut_OBK-518176-WiFeS-Blue-UT20240807T083248-3.cube.fits',
                'cut_OBK-518176-WiFeS-Red--UT20240807T083248-3.cube.fits',
            ],
            'wifes_subns_b3000_560_r3000_1x2' : 
            [
                'OBK-124288-WiFeS-Blue-UT20240412T184018-5.cube.fits',
                'OBK-124288-WiFeS-Red--UT20240412T184018-5.cube.fits',
                'OBK-712864-WiFeS-Blue-UT20240410T113528-0.cube.fits',
                'OBK-712864-WiFeS-Red--UT20240410T113528-0.cube.fits',
            ]
        }

        if dataset in cubes:
            return cubes[dataset]
        else:
           # could also raise an exception here...
           return []

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
    def run_reduction(self):
        def _run_reduction(request,tmp_path):
            ifolder = self.cwd / request.param
            # get a list of fits files in the input folder    
            ifiles = glob.glob(str(ifolder) + "/*.fits*")

            # make a subfolder for the dataset
            sub_folder = tmp_path / request.param
            Path(sub_folder).mkdir(parents=True,exist_ok=True)

            # make symlinks to the fits files in the tmp folder 
            for ifile in ifiles:
                self.mksymlink(sub_folder,Path(ifile))
            # go into the tmp folder and run the pipeline
            os.chdir(str(sub_folder))

            cmd = self.reduce_script + " ." + " --reduce-both"
            with subprocess.Popen(args=cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,cwd='.',env={'PYWIFES_DIR' : self.ref_data},bufsize=1,universal_newlines=True) as p:
                for line in p.stdout:
                    print (line.strip())
            # check that the output products exist
            # TODO
            return sub_folder
        yield _run_reduction
        os.chdir(str(self.cwd))


    @pytest.fixture(params=get_datasets())
    def check_output(self,request,tmp_path,run_reduction):
        # check that the folder returned by the run_reduction fixture 
        # contains the datacube results we require   
        ofolder = run_reduction(request,tmp_path)
        ofiles = self.expected_files(request.param)
        for o in ofiles:
            data_product = Path(ofolder) / 'data_products' / o
            assert Path(data_product).exists()

    def test_run_reductions(self,check_output):
        pass

