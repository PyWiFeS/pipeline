To run the tests:
- cd tests
- git clone git@dev.aao.org.au:datacentral/wifes/sample/classical_20230922.git 
- git clone git@dev.aao.org.au:datacentral/wifes/sample/nod-and-shuffle_20230927.git
- git clone git@dev.aao.org.au:datacentral/wifes/sample/raw_data_half_frame.git
- pytest -vv .

TODO: Add other test datasets to the tests and make the pytests more modular to accommodate them.
