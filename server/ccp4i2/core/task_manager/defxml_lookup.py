import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path


def main():
    outputDir = Path(__file__).parent
    outputPath = outputDir / "defxml_lookup.json"
    result = {}
    for path in outputDir.parent.parent.glob("**/*.def.xml"):
        relPath = os.path.relpath(path, outputDir)
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError as e:
            raise ET.ParseError(f"Error parsing XML file {path}") from e
        pluginName = root.findtext(".//ccp4i2_header/pluginName")
        if pluginName in result:
            raise RuntimeError(
                f"Duplicate def.xml for plugin {pluginName}\n"
                f"  File 1: {result[pluginName]}\n"
                f"  File 2: {relPath}"
            )
        if pluginName:
            result[pluginName] = relPath
        else:
            print(f"Could not determine plugin name for {path}")
    with outputPath.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
