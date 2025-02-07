import re
import json

def preprocess_json(file_path, output_path):
    with open(file_path, 'r') as file:
        raw_data = file.read()
    
    # Fix leading zeros in integers
    raw_data = re.sub(r'"\d+":(\d{2,})', lambda m: m.group(0).lstrip('0'), raw_data)
    
    # Convert hexadecimal floats to decimal
    def hex_to_float(match):
        return str(float.fromhex(match.group(0)))

    raw_data = re.sub(r'0x[0-9A-Fa-f]+\.p[-+][0-9]+', hex_to_float, raw_data)

    # Ensure valid JSON formatting
    raw_data = re.sub(r'("delta": .*?)(\s+"data")', r'\1, \2', raw_data)  # Fix missing commas

    try:
        parsed_data = json.loads(raw_data)
        print("JSON validated successfully!")
        with open(output_path, 'w') as out_file:
            json.dump(parsed_data, out_file, indent=2)
            print(f"Fixed JSON saved to {output_path}")
    except json.JSONDecodeError as e:
        print(f"JSON validation failed: {e}")
        print(raw_data)

# Specify input and output file paths
input_path = "all-compare.json"  # Replace with your file
output_path = "fixed-all-compare.json"  # Replace with desired output path

preprocess_json(input_path, output_path)
