"""One meaning for File.directory. ccp4i2_django_dbapi once defined its own
PATH_FLAG constants as 0/1 against static_data's 1/2, so its relpath put a
job's own files in CCP4_IMPORTED_FILES and gave imported files no path."""


def test_dbapi_uses_the_static_data_path_flags():
    from ccp4i2.db import ccp4i2_django_dbapi as dbapi
    from ccp4i2.db import ccp4i2_static_data as static
    assert dbapi.PATH_FLAG_JOB_DIR == static.PATH_FLAG_JOB_DIR == 1
    assert dbapi.PATH_FLAG_IMPORT_DIR == static.PATH_FLAG_IMPORT_DIR == 2
