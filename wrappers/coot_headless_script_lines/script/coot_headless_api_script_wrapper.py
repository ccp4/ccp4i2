import sys
import importlib
import pathlib
import coot_headless_api as chapi

if __name__ == "__main__":
    sys.path.insert(0, pathlib.Path().absolute().__str__())
    '''
    print('Path', sys.path)
    print('ScriptFileRoot', self.scriptFileRoot)
    sys.stdout.flush()
    importedModule = importlib.import_module(self.scriptFileRoot)
    '''
    print(sys.path)
    sys.exit(0)