from lxml import etree

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def addXMLelement(containerXML, elementname, elementtext):
    #print 'addElement', elementname, type(elementtext), elementtext 
    e2 = etree.Element(elementname)
    e2.text = elementtext
    containerXML.append(e2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class ReflectionDataTypes():
    # Information about data types for I2 converion from mmcif or MTZ to MTZ

    # Dictionary of data types: for each type
    #   List of [mmcif_code, MTZ column name, column type, datasetID]
    REFLECTION_DATA = {
        'I+- anomalous':
        ["pdbx_I_plus Iplus K 1",
         "pdbx_I_plus_sigma SIGIplus M 1",
         "pdbx_I_minus Iminus K 1",
         "pdbx_I_minus_sigma SIGIminus M 1"],
        'Imean':
        ["intensity_meas I J 1",
         "intensity_sigma SIGI Q 1"],
        'F+- anomalous':
        ["pdbx_F_plus Fplus G 1",
         "pdbx_F_plus_sigma SIGFplus L 1",
         "pdbx_F_minus Fminus G 1",
         "pdbx_F_minus_sigma SIGFminus L 1"],
        'Fmean':
        ["F_meas_au F F 1",
         "F_meas_sigma_au SIGF Q 1"],
        'FreeR status':["status FREER  s 0"],
        'FreeR flag': ["pdbx_r_free_flag FREER I 1"]}

    # From CCP4XtalData.py class CObsDataFile
    #    CONTENT_FLAG_IPAIR = 1
    #    CONTENT_FLAG_FPAIR = 2
    #    CONTENT_FLAG_IMEAN = 3
    #    CONTENT_FLAG_FMEAN = 4
    #    CONTENT_ANNOTATION = ['Anomalous Is', 'Anomalous SFs',
    #                          'Mean Is' ,'Mean SFs']
    CONTENT_FLAGS = {
        'I+- anomalous': 1,
        'F+- anomalous': 2,
        'Imean': 3,
        'Fmean': 4}
    # For mtz to I2 mtz, we want things indexed by content type
    CONTENT_TYPES = {
        1: 'I+- anomalous',
        2: 'F+- anomalous',
        3: 'Imean',
        4: 'Fmean'}

    TYPE_LABELS = {
        'I+- anomalous': "Anomalous pairs of intensities",
        'F+- anomalous': "Anomalous pairs of amplitudes",
        'Imean': "Mean intensities",
        'Fmean': "Mean amplitudes"}

    # Order of preference for data types
    DATA_PRIORITY = ['I+- anomalous', 'Imean', 'F+- anomalous', 'Fmean']

    # FreeR types
    FREER_TYPES = ['FreeR flag', 'FreeR status']

