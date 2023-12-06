#if defined (_WIN32)
#include <windows.h>
#include <math.h>
#endif
#include <string>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <Python.h>

int main(int argc, char* argv[]){

#if PY_MAJOR_VERSION > 2
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }
    Py_SetProgramName(program);  /* optional but recommended */
    wchar_t **wargs = new wchar_t*[argc];
    for(int i=0;i<argc;i++){
        wargs[i] = Py_DecodeLocale(argv[i], NULL);
    }
    Py_Initialize();
    int err = Py_Main(argc,wargs);
    if (Py_FinalizeEx() < 0) {
        exit(120);
    }
    PyMem_RawFree(program);
    for(int i=0;i<argc;i++){
        PyMem_RawFree(wargs[i]);
    }

    return err;
#else
    Py_Initialize();
    int err = Py_Main(argc,argv);
    Py_Finalize();  
    return err;
#endif
}
