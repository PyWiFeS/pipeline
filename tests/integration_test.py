import pytest
import subprocess
import os
import shutil
import logging


# def setup_logging(log_file):
#     logging.basicConfig(
#         level=logging.INFO,
#         format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
#         handlers=[
#             logging.FileHandler(log_file),
#             logging.StreamHandler()
#         ]
#     )
def setup_logging(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    file_handler = logging.FileHandler(log_file)
    stream_handler = logging.StreamHandler()

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    def close_loggers():
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

    return close_loggers

@pytest.fixture(scope="module")
def setup_test_environment():
    # Get the directory of the current file
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Define the data processing flavours to be tested: Classical, Nod-and-Shuffle, and Half-Frame
    data_flavours = ['class_B3000_R3000', 'ns_B3000_R3000', 'hf_B3000_R3000']
    data_flavours = ['hf_B3000_R3000']

    # Initialize the list to store test case details
    test_cases = []

    # Loop through each data flavour and prepare its corresponding paths
    for flavour in data_flavours:
        test_case = {
            'test_case':flavour,
            'data_dir': os.path.join(current_dir, "data", flavour),
            'config_file': os.path.join(current_dir, "data", f"{flavour}.json"),
        }
        test_cases.append(test_case)
    yield test_cases

def test_pipeline_integration(setup_test_environment, caplog):
        test_cases = setup_test_environment
        current_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(current_dir, "data_products")
        script_path = os.path.join(current_dir, "../reduction_scripts/reduce_data.py")

        # Ensure the script path is correct
        assert os.path.exists(script_path), f"{script_path} does not exist"



        # Setup logging
        log_file = os.path.join(current_dir, "test_log.log")
        setup_logging(log_file)
        logger = logging.getLogger(__name__)


        for test_case in test_cases:
            data_dir = test_case['data_dir']

            logger.info(f"Running test case with data directory: {data_dir}")


            # Run the reduce_data.py script with subprocess
            result = subprocess.run([
                "python3", script_path, data_dir,
            ], capture_output=True, text=True)
            
            # Print stdout and stderr for debugging
            print("Standard Output:", result.stdout)
            print("Standard Error:", result.stderr)
            
            # # Check the script output
            # assert "Pipeline started" in result.stdout, f"Pipeline output: {result.stdout}\nPipeline error: {result.stderr}"
            # assert "Raw data directory" in result.stdout, f"Pipeline output: {result.stdout}\nPipeline error: {result.stderr}"
            # assert "Red parameters" in result.stdout, f"Pipeline output: {result.stdout}\nPipeline error: {result.stderr}"
            # assert "Blue parameters" in result.stdout, f"Pipeline output: {result.stdout}\nPipeline error: {result.stderr}"
            # assert "Using master calibration files from" in result.stdout, f"Pipeline output: {result.stdout}\nPipeline error: {result.stderr}"
            # assert "Pipeline finished" in result.stdout, f"Pipeline output: {result.stdout}\nPipeline error: {result.stderr}"
            
            # Verify log file content
            log_file_path = os.path.join(output_dir, "pywifes_logger.log")
            assert os.path.exists(log_file_path), f"Log file not found: {log_file_path}"
            with open(log_file_path, "r") as log_file:
                log_content = log_file.read()
                assert "Starting PyWiFeS data reduction pipeline." in log_content
                assert "All done" in log_content
            
            # Cleanup output directory
            shutil.rmtree(output_dir)

            logger.info(f"Completed test case with data directory: {data_dir}")


if __name__ == "__main__":
    pytest.main()
