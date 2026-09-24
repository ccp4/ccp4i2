"""The load plan knows the moorhen task's list-shaped inputs."""

from ccp4i2.cootbridge import api_client

PARAMS = """<?xml version='1.0'?>
<ccp4i2><ccp4i2_body id="moorhen"><inputData>
  <XYZIN_LIST><CPdbDataFile><baseName>model.pdb</baseName><dbFileId>a</dbFileId></CPdbDataFile></XYZIN_LIST>
  <MAPIN_LIST><CMapDataFile><baseName>mask.map</baseName><subType>4</subType></CMapDataFile></MAPIN_LIST>
  <DICT_LIST>
    <CDictDataFile><baseName>LIG.cif</baseName></CDictDataFile>
    <CDictDataFile><baseName>DRG.cif</baseName></CDictDataFile>
  </DICT_LIST>
  <FPHIIN_LIST><CMapCoeffsDataFile><baseName>fphi.mtz</baseName><subType>1</subType></CMapCoeffsDataFile></FPHIIN_LIST>
</inputData></ccp4i2_body></ccp4i2>
"""


def test_dictionaries_come_first_and_real_space_maps_are_a_kind():
    entries = api_client.parse_input_params(PARAMS)
    kinds = [(e["param"], e["kind"], e["base_name"]) for e in entries]
    assert kinds == [
        ("DICT_LIST", "dictionary", "LIG.cif"),
        ("DICT_LIST", "dictionary", "DRG.cif"),
        ("XYZIN_LIST", "coordinates", "model.pdb"),
        ("FPHIIN_LIST", "map_2fofc", "fphi.mtz"),
        ("MAPIN_LIST", "map", "mask.map"),
    ]
    assert entries[-1]["sub_type"] == "4"
