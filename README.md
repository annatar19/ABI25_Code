The RNTuple to ndjson converter used for validation can be found here:
  https://github.com/annatar19/ABI25_to_json

The ORC to ndjson converter used for validation can be found here:
  https://github.com/annatar19/ABI25_orc_to_json

The dataset used, DecayTree, can be found here:
  https://amsuni-my.sharepoint.com/:u:/g/personal/stijn_jongbloed_student_uva_nl/EUJdBj8LeN1Kn9i-QugsN5MBJA7UDAuXpEJ0m36VMqHQRQ?e=E5m4Wl
It should be accessible to persons within the UvA. 
  
The CMakeLists assume a folder named root_src to be present in this directory, containing a locally build version of root: https://github.com/root-project/root 

In the case root has been installed through the systems package manager, or have been build elsewhere, the CMakeLists need to be updated. As for ORC, the CMakeLists assume an installation installed through the system package manager, this might require updating on other systems.
