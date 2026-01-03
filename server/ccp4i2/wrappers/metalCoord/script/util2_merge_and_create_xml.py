import json
import glob
import xml.etree.ElementTree as ET
from collections import defaultdict


def main():

    # List all JSON files in the current directory (you can specify a pattern if needed)
    json_files = glob.glob('*[0-9].json')  # This will match all files like 09.json, 11.json, 24.json, etc.

    # Initialize an empty dictionary to store merged data
    merged_data = defaultdict(list)

    # Loop over all the JSON files
    for json_file in json_files:
        with open(json_file, 'r') as file:
            data = json.load(file)
            
            # Function to remove the "cod" key
            def remove_cod(data):
                for key in data:
                    if "cod" in data[key]:
                        del data[key]["cod"]
                return data
            
            # Remove the "cod" key from the data
            data = remove_cod(data)
            
            # Group the data by coordination in the merged_data dictionary
            for key, value in data.items():
                coordination = value["coordination"]
                merged_data[coordination].append({"name": key, "frequency": value["frequency"], "count": value["count"]})

    # Save the merged data to a new JSON file (grouped by coordination)
    with open('merged.json', 'w') as merged_file:
        json.dump(merged_data, merged_file, indent=4)

    print("Merged JSON saved as 'merged.json'")

    # Prepare XML structure
    root = ET.Element("container", id="coordination")

    # Sort the coordination numbers (COORD) from lowest to highest, with leading zeros
    for coordination in sorted(merged_data.keys()):
        content = ET.SubElement(root, "content", id=f"COORD{coordination:02d}")  # Leading zeroes
        class_name = ET.SubElement(content, "className")
        class_name.text = "CString"
        
        qualifiers = ET.SubElement(content, "qualifiers")
        only_enumerators = ET.SubElement(qualifiers, "onlyEnumerators")
        only_enumerators.text = "True"
        
        # Sort the entries by frequency in descending order
        sorted_entries = sorted(merged_data[coordination], key=lambda x: x["frequency"], reverse=True)
        
        enumerators = ET.SubElement(qualifiers, "enumerators")
        enumerators.text = "auto," + ",".join(entry["name"] for entry in sorted_entries)
        
        # Create menuText with frequency percentage
        menu_text = ET.SubElement(qualifiers, "menuText")
        menu_text.text = "auto," + ",".join([f"{entry['name']} ({entry['frequency'] * 100:.1f} %)" for entry in sorted_entries])
        
        # Add the most frequent class to the <default> field
        default = ET.SubElement(qualifiers, "default")
        default.text = "auto"
        # default.text = sorted_entries[0]["name"]  # Most frequent class

    # Create a nicely formatted XML string with indentation
    tree = ET.ElementTree(root)

    # Save the XML file with pretty indentation
    with open("output.xml", "wb") as xml_file:
        tree.write(xml_file, encoding="utf-8", xml_declaration=True)

    # To format nicely with indentation
    import xml.dom.minidom
    with open("output.xml", "r", encoding="utf-8") as xml_file:
        xml_content = xml_file.read()

    # Pretty print the XML content
    pretty_xml = xml.dom.minidom.parseString(xml_content).toprettyxml(indent="  ")

    with open("output.xml", "w", encoding="utf-8") as xml_file:
        xml_file.write(pretty_xml)

    print("XML saved as 'output.xml'")


if __name__ == "__main__":
    main()